"""In-process job manager for the QC and coverage-check services.

Only one API process is assumed.  Every job owns a private directory under
``DATA_DIR/jobs``.  Uploads and partial outputs live there, and are removed
before a stop or failure is reported.  Successful QC jobs keep only the two
filtered gzip files; successful coverage jobs keep only their published
``coverage.json`` report.

Coverage jobs never copy the QC result pair: they read the published gzip
files in place, so they add no large files to the data directory.

Concurrency contract
--------------------
* A stop requested while a job is pending/running atomically moves it to
  ``stopping`` and sets its event; it can never subsequently become
  ``succeeded``.
* A job that already reached ``succeeded``/``failed``/``stopped`` rejects
  stop requests (the API maps that to HTTP 409).
* ``stopping`` is only turned into ``stopped`` after the worker has exited and
  its temporary files have been removed.
* A coverage report is written to ``coverage.json.part`` and atomically
  renamed to ``coverage.json`` under the manager lock, so status queries never
  observe a partial report.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .coverage import CoverageResult, Reference, run_coverage
from .qc import JobStopped, QCError, QCResult, publish_outputs, run_qc

logger = logging.getLogger(__name__)

PENDING = "pending"
RUNNING = "running"
STOPPING = "stopping"
SUCCEEDED = "succeeded"
FAILED = "failed"
STOPPED = "stopped"

KIND_QC = "qc"
KIND_COVERAGE = "coverage"

TERMINAL_STATUSES = frozenset({SUCCEEDED, FAILED, STOPPED})
STOPPABLE_STATUSES = frozenset({PENDING, RUNNING})


class JobNotFound(Exception):
    pass


class JobConflict(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Job:
    job_id: str
    directory: Path
    kind: str = KIND_QC
    r1_filename: str | None = None
    r2_filename: str | None = None
    min_quality: int | None = None
    # Coverage jobs read the published outputs of a successful QC job.
    source_job_id: str | None = None
    source_r1: Path | None = None
    source_r2: Path | None = None
    reference: Reference | None = None
    min_fragment: int | None = None
    max_fragment: int | None = None
    min_depth: int | None = None
    coverage: CoverageResult | None = None
    status: str = PENDING
    created_at: str = field(default_factory=_now)
    started_at: str | None = None
    finished_at: str | None = None
    pairs_processed: int = 0
    total_pairs: int | None = None
    kept_pairs: int | None = None
    discarded_pairs: int | None = None
    error_code: str | None = None
    error_message: str | None = None
    stop_event: threading.Event = field(default_factory=threading.Event)
    thread: threading.Thread | None = None

    # -- paths -----------------------------------------------------------
    @property
    def input_r1(self) -> Path:
        return self.directory / "R1.fastq.gz"

    @property
    def input_r2(self) -> Path:
        return self.directory / "R2.fastq.gz"

    @property
    def part_r1(self) -> Path:
        return self.directory / "R1.kept.fastq.gz.part"

    @property
    def part_r2(self) -> Path:
        return self.directory / "R2.kept.fastq.gz.part"

    @property
    def output_r1(self) -> Path:
        return self.directory / "R1.kept.fastq.gz"

    @property
    def output_r2(self) -> Path:
        return self.directory / "R2.kept.fastq.gz"

    @property
    def part_report(self) -> Path:
        return self.directory / "coverage.json.part"

    @property
    def report_path(self) -> Path:
        return self.directory / "coverage.json"

    @property
    def downloads_ready(self) -> bool:
        return (
            self.kind == KIND_QC
            and self.status == SUCCEEDED
            and self.output_r1.is_file()
            and self.output_r2.is_file()
        )

    @property
    def report_ready(self) -> bool:
        return (
            self.kind == KIND_COVERAGE
            and self.status == SUCCEEDED
            and self.report_path.is_file()
        )

    def _error_dict(self) -> dict[str, str] | None:
        if not self.error_code:
            return None
        return {"code": self.error_code, "message": self.error_message or ""}

    def to_dict(self) -> dict[str, Any]:
        if self.kind == KIND_COVERAGE:
            return self._coverage_dict()
        result: dict[str, Any] = {
            "job_id": self.job_id,
            "kind": KIND_QC,
            "status": self.status,
            "min_quality": self.min_quality,
            "r1_filename": self.r1_filename,
            "r2_filename": self.r2_filename,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "pairs_processed": self.pairs_processed,
            "total_pairs": self.total_pairs,
            "kept_pairs": self.kept_pairs,
            "discarded_pairs": self.discarded_pairs,
            "error": self._error_dict(),
            "downloads": None,
        }
        if self.downloads_ready:
            result["downloads"] = {
                "r1": f"/jobs/{self.job_id}/download/r1",
                "r2": f"/jobs/{self.job_id}/download/r2",
            }
        return result

    def _coverage_dict(self) -> dict[str, Any]:
        coverage = self.coverage
        intervals = (
            [interval.to_dict() for interval in coverage.low_coverage_intervals]
            if coverage is not None
            else None
        )
        reference = (
            {
                "name": self.reference.name,
                "length": len(self.reference.sequence),
            }
            if self.reference is not None
            else None
        )
        result: dict[str, Any] = {
            "job_id": self.job_id,
            "kind": KIND_COVERAGE,
            "status": self.status,
            "qc_job_id": self.source_job_id,
            "source_job_id": self.source_job_id,
            "reference": reference,
            "reference_name": self.reference.name if reference else None,
            "reference_length": reference["length"] if reference else None,
            "min_fragment": self.min_fragment,
            "max_fragment": self.max_fragment,
            "fragment_length_range": (
                [self.min_fragment, self.max_fragment]
                if self.min_fragment is not None
                and self.max_fragment is not None
                else None
            ),
            "min_depth": self.min_depth,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "pairs_processed": self.pairs_processed,
            "total_pairs": coverage.total_pairs if coverage else None,
            "placed_pairs": coverage.placed_pairs if coverage else None,
            "ambiguous_pairs": coverage.ambiguous_pairs if coverage else None,
            "unplaced_pairs": coverage.unplaced_pairs if coverage else None,
            "placed": coverage.placed_pairs if coverage else None,
            "ambiguous": coverage.ambiguous_pairs if coverage else None,
            "unplaced": coverage.unplaced_pairs if coverage else None,
            "unique_pairs": coverage.placed_pairs if coverage else None,
            "multimapping_pairs": (
                coverage.ambiguous_pairs if coverage else None
            ),
            "unmapped_pairs": coverage.unplaced_pairs if coverage else None,
            "pair_counts": (
                {
                    "placed": coverage.placed_pairs,
                    "ambiguous": coverage.ambiguous_pairs,
                    "unplaced": coverage.unplaced_pairs,
                }
                if coverage
                else None
            ),
            "low_coverage_intervals": intervals,
            # Convenience aliases for the same list.
            "low_coverage": intervals,
            "low_coverage_regions": intervals,
            "result": coverage.to_dict() if coverage else None,
            "report": (
                f"/coverage-jobs/{self.job_id}/report"
                if self.report_ready
                else None
            ),
            "error": self._error_dict(),
        }
        return result


class JobManager:
    def __init__(self, data_dir: Path) -> None:
        self.data_dir = Path(data_dir)
        self.jobs_dir = self.data_dir / "jobs"
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        self._remove_stale_directories()

    # -- lifecycle -------------------------------------------------------
    def _remove_stale_directories(self) -> None:
        """Drop directories left behind by a previous process.

        Job state is in memory only, so such directories can no longer be
        downloaded and would otherwise accumulate forever.
        """
        for child in self.jobs_dir.iterdir():
            try:
                if child.is_dir():
                    shutil.rmtree(child, ignore_errors=True)
                else:
                    child.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort
                logger.warning("could not remove stale path %s", child)

    def create_job(
        self,
        r1_filename: str,
        r2_filename: str,
        min_quality: int,
        kind: str = KIND_QC,
    ) -> Job:
        job_id = uuid.uuid4().hex
        directory = self.jobs_dir / job_id
        directory.mkdir(parents=True)
        job = Job(
            job_id=job_id,
            directory=directory,
            kind=kind,
            r1_filename=r1_filename,
            r2_filename=r2_filename,
            min_quality=min_quality,
        )
        with self._lock:
            self._jobs[job_id] = job
        return job

    def create_coverage_job(
        self,
        *,
        source_job_id: str,
        source_r1: Path,
        source_r2: Path,
        reference: Reference,
        min_fragment: int,
        max_fragment: int,
        min_depth: int,
    ) -> Job:
        """Register a coverage job that reads a published QC result pair."""
        job_id = uuid.uuid4().hex
        directory = self.jobs_dir / job_id
        directory.mkdir(parents=True)
        job = Job(
            job_id=job_id,
            directory=directory,
            kind=KIND_COVERAGE,
            r1_filename=source_r1.name,
            r2_filename=source_r2.name,
            source_job_id=source_job_id,
            source_r1=source_r1,
            source_r2=source_r2,
            reference=reference,
            min_fragment=min_fragment,
            max_fragment=max_fragment,
            min_depth=min_depth,
        )
        with self._lock:
            self._jobs[job_id] = job
        return job

    def abort_upload(self, job_id: str) -> None:
        """Discard a job whose upload did not complete."""
        with self._lock:
            job = self._jobs.pop(job_id, None)
        if job is not None:
            shutil.rmtree(job.directory, ignore_errors=True)

    def start(self, job_id: str) -> None:
        job = self.get(job_id)
        thread = threading.Thread(
            target=self._run_job,
            args=(job,),
            name=f"{job.kind}-{job_id[:8]}",
            daemon=True,
        )
        with self._lock:
            if job.status != PENDING:
                raise JobConflict("job has already been started")
            job.thread = thread
        thread.start()

    def shutdown(self) -> None:
        with self._lock:
            jobs = list(self._jobs.values())
        for job in jobs:
            job.stop_event.set()
        for job in jobs:
            thread = job.thread
            if thread is not None and thread.is_alive():
                thread.join(timeout=5)

    # -- queries ---------------------------------------------------------
    def get(self, job_id: str) -> Job:
        with self._lock:
            job = self._jobs.get(job_id)
        if job is None:
            raise JobNotFound(job_id)
        return job

    def list_jobs(self, kind: str | None = None) -> list[Job]:
        with self._lock:
            jobs = list(self._jobs.values())
        if kind is not None:
            jobs = [job for job in jobs if job.kind == kind]
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)

    def snapshot(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            return job.to_dict()

    def request_stop(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            if job.status in STOPPABLE_STATUSES:
                job.status = STOPPING
                job.stop_event.set()
            elif job.status == STOPPING:
                job.stop_event.set()
            elif job.status == STOPPED:
                pass  # idempotent
            else:  # succeeded / failed
                raise JobConflict(
                    f"job is already {job.status}; it can no longer be stopped"
                )
            return job.to_dict()

    def output_paths(self, job_id: str) -> tuple[Path, Path]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            if job.kind != KIND_QC:
                raise JobConflict(
                    "this job has no downloadable QC output; it is a "
                    "coverage job"
                )
            if job.status != SUCCEEDED:
                raise JobConflict(
                    f"job status is {job.status!r}; downloads are only "
                    "available after successful completion"
                )
            if not (job.output_r1.is_file() and job.output_r2.is_file()):
                raise JobConflict("result files are no longer available")
            return job.output_r1, job.output_r2

    def coverage_report_path(self, job_id: str) -> Path:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            if job.kind != KIND_COVERAGE:
                raise JobConflict("this job is not a coverage job")
            if job.status != SUCCEEDED:
                raise JobConflict(
                    f"job status is {job.status!r}; the coverage report is "
                    "only available after successful completion"
                )
            if not job.report_path.is_file():
                raise JobConflict("coverage report is no longer available")
            return job.report_path

    # -- worker ----------------------------------------------------------
    def _run_job(self, job: Job) -> None:
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                # A stop request landed while the worker was still starting.
                # The worker has now effectively exited, so it is responsible
                # for cleaning up and publishing ``stopped``.
                stop_early = True
            elif job.status == PENDING:
                job.status = RUNNING
                job.started_at = _now()
                stop_early = False
            else:  # pragma: no cover - upload aborted / defensive
                return
        if stop_early:
            self._finish_stopped(job)
            return

        try:
            if job.kind == KIND_COVERAGE:
                result = run_coverage(
                    job.source_r1,
                    job.source_r2,
                    job.reference,
                    job.min_fragment,
                    job.max_fragment,
                    job.min_depth,
                    job.stop_event,
                    progress_cb=lambda pairs: self._set_progress(job, pairs),
                )
            else:
                result = run_qc(
                    job.input_r1,
                    job.input_r2,
                    job.part_r1,
                    job.part_r2,
                    job.min_quality,
                    job.stop_event,
                    progress_cb=lambda pairs: self._set_progress(job, pairs),
                )
        except JobStopped:
            self._finish_stopped(job)
        except QCError as exc:
            self._finish_failed(job, exc.code, exc.message)
        except Exception as exc:  # pragma: no cover - unexpected bug
            logger.exception("unexpected error in job %s", job.job_id)
            self._finish_failed(job, "internal_error", str(exc))
        else:
            if job.kind == KIND_COVERAGE:
                self._finish_coverage_completed(job, result)
            else:
                self._finish_completed(job, result)

    def _set_progress(self, job: Job, pairs: int) -> None:
        with self._lock:
            if pairs > job.pairs_processed:
                job.pairs_processed = pairs

    def _finish_completed(self, job: Job, result: QCResult) -> None:
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                # A stop was accepted while the worker was finishing.
                stop_wins = True
            else:
                stop_wins = False
                # Publish only under the same lock used by stop requests and
                # download lookups: once ``succeeded`` is visible the files
                # are already in place, and a later stop is rejected.
                publish_outputs(
                    job.part_r1,
                    job.part_r2,
                    job.output_r1,
                    job.output_r2,
                )
                job.total_pairs = result.total_pairs
                job.kept_pairs = result.kept_pairs
                job.discarded_pairs = result.discarded_pairs
                job.pairs_processed = result.total_pairs
                job.finished_at = _now()
                job.status = SUCCEEDED
        if stop_wins:
            self._finish_stopped(job)
            return
        self._cleanup_inputs(job)

    def _finish_coverage_completed(
        self, job: Job, result: CoverageResult
    ) -> None:
        # A report is only produced after the entire stream was consumed
        # successfully.  It is first written to ``coverage.json.part`` and
        # published by an atomic rename, under the same lock used by stop
        # requests, so a failed/stopped job never exposes a partial report.
        try:
            with open(job.part_report, "w", encoding="utf-8") as handle:
                json.dump(result.to_dict(), handle, ensure_ascii=False)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
        except OSError as exc:
            self._finish_failed(
                job, "internal_error", f"could not write coverage report: {exc}"
            )
            return

        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                stop_wins = True
                publish_error = None
            else:
                stop_wins = False
                publish_error = None
                try:
                    os.replace(job.part_report, job.report_path)
                except OSError as exc:  # pragma: no cover - disk failure
                    publish_error = exc
                else:
                    job.coverage = result
                    job.total_pairs = result.total_pairs
                    job.pairs_processed = result.total_pairs
                    job.finished_at = _now()
                    job.status = SUCCEEDED
        if stop_wins:
            self._finish_stopped(job)
            return
        if publish_error is not None:  # pragma: no cover - disk failure
            self._finish_failed(
                job,
                "internal_error",
                f"could not publish coverage report: {publish_error}",
            )
            return
        # Coverage jobs own no input copies; the coverage report is the result.
        self._cleanup_inputs(job)

    def _finish_failed(self, job: Job, code: str, message: str) -> None:
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                stop_wins = True
            else:
                stop_wins = False
                job.status = FAILED
                job.error_code = code
                job.error_message = message
                job.finished_at = _now()
        if stop_wins:
            self._finish_stopped(job)
            return
        self._cleanup_all(job)

    def _finish_stopped(self, job: Job) -> None:
        # Temporary inputs and partial outputs are removed *before* the job is
        # reported as stopped, as required by the stop contract.
        self._cleanup_all(job)
        with self._lock:
            job.total_pairs = None
            job.kept_pairs = None
            job.discarded_pairs = None
            job.error_code = None
            job.error_message = None
            job.finished_at = _now()
            job.status = STOPPED

    def _cleanup_inputs(self, job: Job) -> None:
        for path in (job.input_r1, job.input_r2, job.part_r1, job.part_r2):
            try:
                path.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort
                logger.warning("could not remove %s", path)

    def _cleanup_all(self, job: Job) -> None:
        shutil.rmtree(job.directory, ignore_errors=True)
