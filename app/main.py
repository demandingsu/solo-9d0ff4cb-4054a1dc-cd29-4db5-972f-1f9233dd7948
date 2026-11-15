"""FastAPI application for paired-end FASTQ quality control and coverage."""

from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Union

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, RedirectResponse

from . import __version__
from .coverage import (
    MAX_REFERENCE_BASES,
    MAX_REFERENCE_UPLOAD_BYTES,
    ReferenceError,
    parse_reference_upload,
)
from .jobs import (
    JobConflict,
    JobManager,
    JobNotFound,
)
from .qc import MAX_PHRED, MIN_PHRED

logger = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path(os.environ.get("QC_DATA_DIR", "./data"))
UPLOAD_CHUNK_BYTES = 1 << 20
MAX_UPLOAD_BYTES = int(
    os.environ.get("QC_MAX_UPLOAD_BYTES", str(20 * 1024**3))
)


def _manager(request: Request) -> JobManager:
    return request.app.state.manager


async def _save_upload(
    upload: UploadFile, destination: Path, label: str
) -> int:
    """Stream one uploaded file to disk without buffering it in memory."""
    total = 0
    try:
        with open(destination, "wb") as handle:
            while True:
                chunk = await upload.read(UPLOAD_CHUNK_BYTES)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_UPLOAD_BYTES:
                    raise HTTPException(
                        status_code=413,
                        detail=(
                            f"{label} exceeds the maximum upload size of "
                            f"{MAX_UPLOAD_BYTES} bytes"
                        ),
                    )
                handle.write(chunk)
    finally:
        await upload.close()
    return total


def _parse_threshold(
    min_quality: str | None,
    quality_threshold: str | None,
    threshold: str | None,
) -> int:
    raw = None
    for candidate in (min_quality, quality_threshold, threshold):
        if candidate is not None and str(candidate).strip() != "":
            raw = str(candidate).strip()
            break
    if raw is None:
        raise HTTPException(
            status_code=422,
            detail=(
                "missing quality threshold: send 'min_quality' as a form "
                f"field (integer {MIN_PHRED}..{MAX_PHRED}, Phred+33)"
            ),
        )
    try:
        value = int(raw, 10)
    except ValueError:
        raise HTTPException(
            status_code=422,
            detail=f"quality threshold {raw!r} is not an integer",
        ) from None
    if not (MIN_PHRED <= value <= MAX_PHRED):
        raise HTTPException(
            status_code=422,
            detail=(
                f"quality threshold must be between {MIN_PHRED} and "
                f"{MAX_PHRED} (Phred+33)"
            ),
        )
    return value


def _first_present(*candidates: str | None) -> str | None:
    for candidate in candidates:
        if candidate is not None and str(candidate).strip() != "":
            return str(candidate).strip()
    return None


def _parse_int_field(
    name: str, *candidates: str | None, required: bool = True
) -> int | None:
    raw = _first_present(*candidates)
    if raw is None:
        if not required:
            return None
        raise HTTPException(
            status_code=422,
            detail=f"missing required field {name!r}",
        )
    try:
        return int(raw, 10)
    except ValueError:
        raise HTTPException(
            status_code=422,
            detail=f"{name} {raw!r} is not an integer",
        ) from None


async def _read_reference_upload(reference) -> bytes:
    """Accept a FASTA upload either as a file part or as a text field."""
    # Starlette may decode multipart parts into its own UploadFile class, so
    # duck-type on ``read`` rather than relying on a specific class identity.
    if reference is not None and hasattr(reference, "read"):
        total = 0
        chunks: list[bytes] = []
        try:
            while True:
                chunk = await reference.read(UPLOAD_CHUNK_BYTES)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_REFERENCE_UPLOAD_BYTES:
                    raise HTTPException(
                        status_code=413,
                        detail=(
                            "reference exceeds the maximum size of "
                            f"{MAX_REFERENCE_UPLOAD_BYTES} bytes"
                        ),
                    )
                chunks.append(chunk)
        finally:
            try:
                await reference.close()
            except Exception:  # pragma: no cover - best effort
                pass
        return b"".join(chunks)
    if isinstance(reference, (bytes, bytearray, memoryview)):
        raw = bytes(reference)
    elif isinstance(reference, str):
        raw = reference.encode("utf-8")
    else:
        raise HTTPException(
            status_code=422,
            detail="missing required field 'reference' (FASTA content)",
        )
    if len(raw) > MAX_REFERENCE_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail=(
                "reference exceeds the maximum size of "
                f"{MAX_REFERENCE_UPLOAD_BYTES} bytes"
            ),
        )
    return raw


def create_app(data_dir: Path | str | None = None) -> FastAPI:
    resolved_data_dir = Path(
        data_dir
        if data_dir is not None
        else os.environ.get("QC_DATA_DIR", str(DEFAULT_DATA_DIR))
    )

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        manager = JobManager(resolved_data_dir)
        app.state.manager = manager
        try:
            yield
        finally:
            manager.shutdown()

    app = FastAPI(
        title="Paired-end FASTQ QC and coverage API",
        version=__version__,
        description=(
            "Streaming quality-control filtering for paired-end FASTQ files, "
            "plus a coverage-verification step for targeted sequencing.\n\n"
            "**QC jobs** (`POST /jobs`): upload `R1`/`R2` `.fastq.gz` files "
            "together with a Phred+33 quality threshold. A pair is kept only "
            "when every base at both ends has quality at least the threshold. "
            "Pair order is preserved, and successful jobs expose both "
            "filtered gzip files for download.\n\n"
            "**Coverage jobs** (`POST /coverage-jobs`): point at a successful "
            "QC job and submit a single-record ACGT reference FASTA (at most "
            "50 kb), an inclusive fragment-length range and a minimum depth. "
            "The published R1/R2 files are streamed again; pairs with exactly "
            "one compliant exact placement on the reference add depth over "
            "their whole fragment (overlapping bases counted once). The "
            "report lists the maximal contiguous intervals whose depth is "
            "below the threshold, in 0-based half-open coordinates, together "
            "with the placed/ambiguous/unplaced pair counts.\n\n"
            "Both kinds of jobs can be polled and stopped. A stopped, failed "
            "or corrupt run never publishes a partial report."
        ),
        lifespan=lifespan,
    )

    @app.get("/", include_in_schema=False)
    async def index() -> RedirectResponse:
        return RedirectResponse(url="/docs")

    @app.get("/health", tags=["meta"])
    async def health() -> dict[str, str]:
        return {"status": "ok"}

    @app.post(
        "/jobs",
        status_code=202,
        tags=["jobs"],
        summary="Upload a read pair and start a QC job",
        responses={
            202: {"description": "Job accepted and queued"},
            422: {"description": "Invalid threshold or missing files"},
        },
    )
    async def create_job(
        request: Request,
        r1: UploadFile = File(..., description="Read 1 FASTQ.gz"),
        r2: UploadFile = File(..., description="Read 2 FASTQ.gz"),
        min_quality: str | None = Form(
            None,
            description=(
                f"Minimum Phred+33 quality ({MIN_PHRED}..{MAX_PHRED}). "
                "Pairs with any base below it are discarded."
            ),
        ),
        quality_threshold: str | None = Form(None, include_in_schema=False),
        threshold: str | None = Form(None, include_in_schema=False),
    ) -> dict:
        threshold_value = _parse_threshold(
            min_quality, quality_threshold, threshold
        )
        manager = _manager(request)
        job = manager.create_job(
            r1_filename=Path(r1.filename or "R1.fastq.gz").name,
            r2_filename=Path(r2.filename or "R2.fastq.gz").name,
            min_quality=threshold_value,
        )
        try:
            size1 = await _save_upload(r1, job.input_r1, "R1")
            size2 = await _save_upload(r2, job.input_r2, "R2")
            if size1 == 0 or size2 == 0:
                raise HTTPException(
                    status_code=400,
                    detail=(
                        "one of the uploaded files is empty: "
                        f"R1={size1} bytes, R2={size2} bytes"
                    ),
                )
            manager.start(job.job_id)
        except HTTPException:
            manager.abort_upload(job.job_id)
            raise
        except Exception:
            manager.abort_upload(job.job_id)
            raise
        return manager.snapshot(job.job_id)

    @app.get("/jobs", tags=["jobs"], summary="List jobs (newest first)")
    async def list_jobs(request: Request) -> dict:
        manager = _manager(request)
        return {"jobs": [job.to_dict() for job in manager.list_jobs()]}

    @app.get(
        "/jobs/{job_id}",
        tags=["jobs"],
        summary="Get job status and counts",
        responses={404: {"description": "Unknown job"}},
    )
    async def get_job(request: Request, job_id: str) -> dict:
        try:
            return _manager(request).snapshot(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None

    @app.post(
        "/jobs/{job_id}/stop",
        tags=["jobs"],
        summary="Request job cancellation",
        responses={
            404: {"description": "Unknown job"},
            409: {"description": "Job already finished"},
        },
    )
    async def stop_job(request: Request, job_id: str) -> dict:
        try:
            return _manager(request).request_stop(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None

    @app.get(
        "/jobs/{job_id}/download/{side}",
        tags=["jobs"],
        summary="Download a filtered gzip file (successful QC jobs only)",
        response_class=FileResponse,
        responses={
            404: {"description": "Unknown job"},
            409: {"description": "Job not successfully completed"},
        },
    )
    async def download(
        request: Request, job_id: str, side: str
    ) -> FileResponse:
        if side not in {"r1", "r2"}:
            raise HTTPException(
                status_code=404, detail="side must be 'r1' or 'r2'"
            )
        try:
            r1_path, r2_path = _manager(request).output_paths(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None
        path, filename = (
            (r1_path, "R1.kept.fastq.gz")
            if side == "r1"
            else (r2_path, "R2.kept.fastq.gz")
        )
        return FileResponse(
            path,
            media_type="application/gzip",
            filename=filename,
        )

    # -- coverage verification ------------------------------------------
    def _submit_coverage(
        request: Request,
        qc_job_id: str | None,
        reference,
        min_fragment: int | None,
        max_fragment: int | None,
        min_depth: int | None,
    ) -> dict:
        if not qc_job_id:
            raise HTTPException(
                status_code=422,
                detail=(
                    "missing successful QC job id: send 'qc_job_id' as a "
                    "form field"
                ),
            )
        if min_fragment is None or max_fragment is None:
            raise HTTPException(
                status_code=422,
                detail=(
                    "missing fragment length range: send 'min_fragment' and "
                    "'max_fragment' (inclusive, in bases)"
                ),
            )
        if min_depth is None:
            raise HTTPException(
                status_code=422,
                detail="missing 'min_depth' (minimum required depth)",
            )
        if min_fragment < 1 or max_fragment < min_fragment:
            raise HTTPException(
                status_code=422,
                detail=(
                    "invalid fragment length range: require "
                    "1 <= min_fragment <= max_fragment"
                ),
            )
        if min_depth < 0:
            raise HTTPException(
                status_code=422, detail="min_depth must be >= 0"
            )

        raw_reference = reference
        try:
            parsed = parse_reference_upload(raw_reference)
        except ReferenceError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from None

        manager = _manager(request)
        try:
            source_r1, source_r2 = manager.output_paths(qc_job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404,
                detail=f"QC job {qc_job_id!r} not found",
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None

        job = manager.create_coverage_job(
            source_job_id=qc_job_id,
            source_r1=source_r1,
            source_r2=source_r2,
            reference=parsed,
            min_fragment=min_fragment,
            max_fragment=max_fragment,
            min_depth=min_depth,
        )
        try:
            manager.start(job.job_id)
        except Exception:
            manager.abort_upload(job.job_id)
            raise
        return manager.snapshot(job.job_id)

    def _coverage_parameters(
        form: dict,
        *,
        body: dict | None = None,
        query: dict | None = None,
    ) -> tuple[int | None, int | None, int | None]:
        body = body or {}
        query = query or {}
        min_fragment = _parse_int_field(
            "min_fragment",
            form.get("min_fragment"),
            form.get("fragment_min"),
            form.get("min_fragment_length"),
            form.get("fragment_length_min"),
            body.get("min_fragment"),
            query.get("min_fragment"),
            required=False,
        )
        max_fragment = _parse_int_field(
            "max_fragment",
            form.get("max_fragment"),
            form.get("fragment_max"),
            form.get("max_fragment_length"),
            form.get("fragment_length_max"),
            body.get("max_fragment"),
            query.get("max_fragment"),
            required=False,
        )
        # Nested convenience range: {"fragment_length_range": [lo, hi]} or
        # {"fragment_length_range": {"min": lo, "max": hi}}.
        nested = body.get("fragment_length_range")
        if isinstance(nested, (list, tuple)) and len(nested) == 2:
            try:
                nested_min, nested_max = (int(value) for value in nested)
            except (TypeError, ValueError):
                raise HTTPException(
                    status_code=422,
                    detail="fragment_length_range must contain two integers",
                ) from None
            if min_fragment is None:
                min_fragment = nested_min
            if max_fragment is None:
                max_fragment = nested_max
        elif isinstance(nested, dict):
            try:
                nested_min = int(nested.get("min"))
                nested_max = int(nested.get("max"))
            except (TypeError, ValueError):
                raise HTTPException(
                    status_code=422,
                    detail=(
                        "fragment_length_range must have integer 'min' and "
                        "'max' members"
                    ),
                ) from None
            if min_fragment is None:
                min_fragment = nested_min
            if max_fragment is None:
                max_fragment = nested_max
        nested_range = body.get("fragment_length")
        if isinstance(nested_range, dict):
            if min_fragment is None and nested_range.get("min") is not None:
                min_fragment = _parse_int_field(
                    "fragment_length.min", str(nested_range["min"])
                )
            if max_fragment is None and nested_range.get("max") is not None:
                max_fragment = _parse_int_field(
                    "fragment_length.max", str(nested_range["max"])
                )

        # Convenience: "fragment_length_range=100-300" or "100..300".
        raw_range = _first_present(
            form.get("fragment_length_range"),
            form.get("fragment_range"),
            query.get("fragment_length_range"),
            query.get("fragment_range"),
        )
        if isinstance(raw_range, str) and (
            min_fragment is None or max_fragment is None
        ):
            cleaned = raw_range.replace("..", "-").replace(",", "-")
            pieces = [piece for piece in cleaned.split("-") if piece.strip()]
            if len(pieces) == 2:
                try:
                    parsed_min, parsed_max = (int(piece, 10) for piece in pieces)
                except ValueError:
                    raise HTTPException(
                        status_code=422,
                        detail=(
                            f"fragment length range {raw_range!r} is not a "
                            "pair of integers"
                        ),
                    ) from None
                if min_fragment is None:
                    min_fragment = parsed_min
                if max_fragment is None:
                    max_fragment = parsed_max
        min_depth = _parse_int_field(
            "min_depth",
            form.get("min_depth"),
            form.get("minimum_depth"),
            form.get("min_coverage"),
            form.get("min_dp"),
            form.get("required_depth"),
            form.get("coverage_threshold"),
            body.get("min_depth"),
            query.get("min_depth"),
            required=False,
        )
        return min_fragment, max_fragment, min_depth

    async def _coverage_submission(
        request: Request,
        *,
        qc_job_id: str | None = None,
        reference: Union[UploadFile, str, None] = None,
        min_fragment: str | None = None,
        max_fragment: str | None = None,
        min_depth: str | None = None,
        qc_job_id_from_path: str | None = None,
    ) -> dict:
        # ``request.form()`` is cached by Starlette, so this also fills in
        # aliases that were not declared as explicit parameters.
        form = await request.form()
        body: dict = {}
        if not form and "application/json" in request.headers.get(
            "content-type", ""
        ):
            # Convenience: allow the same fields in a JSON body, in which
            # case the reference must be the FASTA text itself.
            try:
                parsed_body = await request.json()
            except ValueError:
                parsed_body = None
            if isinstance(parsed_body, dict):
                body = parsed_body

        def field(*names: str):
            for name in names:
                value = body.get(name)
                if value is not None and str(value).strip() != "":
                    return value
                value = form.get(name)
                if value is not None and str(value).strip() != "":
                    return value
                value = request.query_params.get(name)
                if value is not None and str(value).strip() != "":
                    return value
            return None

        job_id_value = _first_present(
            qc_job_id,
            field(
                "qc_job_id",
                "job_id",
                "source_job_id",
                "qc_job",
                "quality_job_id",
                "qcJobId",
                "qc_job_no",
                "qc_job_number",
                "source_job",
                "job",
            ),
            qc_job_id_from_path,
        )

        reference_value = reference
        if reference_value is None:
            for field_name in (
                "reference",
                "fasta",
                "reference_fasta",
                "reference_sequence",
                "reference_file",
                "fasta_file",
                "ref",
                "reference_sequence_fasta",
            ):
                value = body.get(field_name)
                if value is None:
                    value = form.get(field_name)
                if value is None:
                    value = request.query_params.get(field_name)
                if value is not None:
                    reference_value = value
                    break
        if reference_value is None:
            # Fallback: accept any uploaded file part, or any form/query value
            # that looks like FASTA text, regardless of its field name.
            for _, value in form.items():
                if hasattr(value, "read"):
                    reference_value = value
                    break
                if isinstance(value, str) and value.lstrip().startswith(">"):
                    reference_value = value
                    break
        if reference_value is None:
            for value in request.query_params.values():
                if value.lstrip().startswith(">"):
                    reference_value = value
                    break
        if reference_value is None or (
            isinstance(reference_value, str)
            and not reference_value.strip()
        ):
            raise HTTPException(
                status_code=422,
                detail=(
                    "missing reference: send a FASTA upload in the "
                    "'reference' file field (single record, ACGT, "
                    f"<= {MAX_REFERENCE_BASES} bp)"
                ),
            )
        raw_reference = await _read_reference_upload(reference_value)

        range_values = {
            "min_fragment": min_fragment
            or request.query_params.get("min_fragment"),
            "max_fragment": max_fragment
            or request.query_params.get("max_fragment"),
            "fragment_min": field("fragment_min", "fragmentMin"),
            "fragment_max": field("fragment_max", "fragmentMax"),
            "min_fragment_length": field("min_fragment_length", "frag_min"),
            "max_fragment_length": field("max_fragment_length", "frag_max"),
            "fragment_length_min": field(
                "fragment_length_min", "min_frag_len"
            ),
            "fragment_length_max": field(
                "fragment_length_max", "max_frag_len"
            ),
            "fragment_length_range": field(
                "fragment_length_range",
                "fragment_range",
                "fragment_length",
                "insert_size_range",
                "fragment_sizes",
            ),
            "fragment_range": field("fragment_range"),
            "min_depth": min_depth
            or request.query_params.get("min_depth"),
            "minimum_depth": field(
                "minimum_depth",
                "min_coverage",
                "depth",
                "min_dp",
                "required_depth",
                "coverage_threshold",
            ),
            "min_coverage": field("min_coverage"),
            "min_dp": field("min_dp", "required_depth", "coverage_threshold"),
        }
        parsed_min, parsed_max, parsed_depth = _coverage_parameters(
            range_values, body=body, query=dict(request.query_params)
        )
        return _submit_coverage(
            request,
            job_id_value,
            raw_reference,
            parsed_min,
            parsed_max,
            parsed_depth,
        )

    @app.post(
        "/coverage-jobs",
        status_code=202,
        tags=["coverage"],
        summary="Start a coverage verification job for a successful QC job",
        responses={
            202: {"description": "Coverage job accepted and queued"},
            404: {"description": "Unknown QC job"},
            409: {"description": "QC job has not completed successfully"},
            422: {"description": "Invalid reference, range or depth"},
        },
    )
    async def create_coverage_job(
        request: Request,
        qc_job_id: str | None = Form(
            None,
            description="Job id of a QC job that finished successfully",
        ),
        reference: Union[UploadFile, str, None] = File(
            None,
            description=(
                "Reference FASTA: one ACGT sequence, at most "
                f"{MAX_REFERENCE_BASES} bp (gzip is accepted)"
            ),
        ),
        min_fragment: str | None = Form(
            None, description="Minimum fragment length in bases (inclusive)"
        ),
        max_fragment: str | None = Form(
            None, description="Maximum fragment length in bases (inclusive)"
        ),
        min_depth: str | None = Form(
            None, description="Minimum required depth; lower is reported"
        ),
    ) -> dict:
        return await _coverage_submission(
            request,
            qc_job_id=qc_job_id,
            reference=reference,
            min_fragment=min_fragment,
            max_fragment=max_fragment,
            min_depth=min_depth,
        )

    @app.post(
        "/jobs/{job_id}/coverage",
        status_code=202,
        tags=["coverage"],
        summary="Start coverage verification for this successful QC job",
        include_in_schema=True,
        responses={
            202: {"description": "Coverage job accepted and queued"},
            404: {"description": "Unknown QC job"},
            409: {"description": "QC job has not completed successfully"},
            422: {"description": "Invalid reference, range or depth"},
        },
    )
    async def create_coverage_job_for(
        request: Request,
        job_id: str,
        reference: Union[UploadFile, str, None] = File(
            None,
            description=(
                "Reference FASTA: one ACGT sequence, at most "
                f"{MAX_REFERENCE_BASES} bp (gzip is accepted)"
            ),
        ),
        min_fragment: str | None = Form(
            None, description="Minimum fragment length in bases (inclusive)"
        ),
        max_fragment: str | None = Form(
            None, description="Maximum fragment length in bases (inclusive)"
        ),
        min_depth: str | None = Form(
            None, description="Minimum required depth; lower is reported"
        ),
    ) -> dict:
        return await _coverage_submission(
            request,
            reference=reference,
            min_fragment=min_fragment,
            max_fragment=max_fragment,
            min_depth=min_depth,
            qc_job_id_from_path=job_id,
        )

    @app.post(
        "/coverage",
        status_code=202,
        tags=["coverage"],
        summary="Alias of /coverage-jobs",
        include_in_schema=False,
    )
    async def create_coverage_job_alias(request: Request) -> dict:
        return await _coverage_submission(request)

    @app.post(
        "/coverage-check",
        status_code=202,
        tags=["coverage"],
        summary="Alias of /coverage-jobs",
        include_in_schema=False,
    )
    async def create_coverage_check_alias(request: Request) -> dict:
        return await _coverage_submission(request)

    @app.post(
        "/jobs/{job_id}/coverage-check",
        status_code=202,
        tags=["coverage"],
        summary="Alias of /jobs/{job_id}/coverage",
        include_in_schema=False,
    )
    async def create_coverage_check_for(
        request: Request, job_id: str
    ) -> dict:
        return await _coverage_submission(
            request, qc_job_id_from_path=job_id
        )

    @app.get(
        "/coverage-jobs/{job_id}",
        tags=["coverage"],
        summary="Get coverage job status and report",
        responses={404: {"description": "Unknown coverage job"}},
    )
    async def get_coverage_job(request: Request, job_id: str) -> dict:
        try:
            return _manager(request).snapshot(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404,
                detail=f"coverage job {job_id!r} not found",
            ) from None

    @app.post(
        "/coverage-jobs/{job_id}/stop",
        tags=["coverage"],
        summary="Request coverage job cancellation",
        responses={
            404: {"description": "Unknown coverage job"},
            409: {"description": "Job already finished"},
        },
    )
    async def stop_coverage_job(request: Request, job_id: str) -> dict:
        return await stop_job(request, job_id)

    @app.get(
        "/coverage/{job_id}",
        tags=["coverage"],
        summary="Alias of /coverage-jobs/{job_id}",
        include_in_schema=False,
    )
    async def get_coverage_job_alias(request: Request, job_id: str) -> dict:
        return await get_coverage_job(request, job_id)

    @app.get(
        "/coverage-jobs/{job_id}/status",
        tags=["coverage"],
        summary="Alias of /coverage-jobs/{job_id}",
        include_in_schema=False,
    )
    async def get_coverage_status(
        request: Request, job_id: str
    ) -> dict:
        return await get_coverage_job(request, job_id)

    @app.get(
        "/coverage/{job_id}/status",
        tags=["coverage"],
        summary="Alias of /coverage-jobs/{job_id}",
        include_in_schema=False,
    )
    async def get_coverage_status_alias(
        request: Request, job_id: str
    ) -> dict:
        return await get_coverage_job(request, job_id)

    @app.get(
        "/jobs/{job_id}/coverage",
        tags=["coverage"],
        summary="Latest coverage job started from this QC job",
        include_in_schema=True,
        responses={404: {"description": "No coverage job for this QC job"}},
    )
    async def latest_coverage_for_job(
        request: Request, job_id: str
    ) -> dict:
        manager = _manager(request)
        try:
            manager.get(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"QC job {job_id!r} not found"
            ) from None
        for job in manager.list_jobs(kind="coverage"):
            if job.source_job_id == job_id:
                return manager.snapshot(job.job_id)
        raise HTTPException(
            status_code=404,
            detail=f"no coverage job has been started from QC job {job_id!r}",
        )

    @app.get(
        "/coverage-jobs",
        tags=["coverage"],
        summary="List coverage jobs (newest first)",
    )
    async def list_coverage_jobs(request: Request) -> dict:
        manager = _manager(request)
        return {
            "jobs": [
                job.to_dict() for job in manager.list_jobs(kind="coverage")
            ]
        }

    @app.get(
        "/coverage-jobs/{job_id}/report",
        tags=["coverage"],
        summary="Download the published coverage report (successful jobs only)",
        response_class=FileResponse,
        responses={
            404: {"description": "Unknown coverage job"},
            409: {"description": "Report not published"},
        },
    )
    async def coverage_report(
        request: Request, job_id: str
    ) -> FileResponse:
        try:
            path = _manager(request).coverage_report_path(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404,
                detail=f"coverage job {job_id!r} not found",
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None
        return FileResponse(
            path,
            media_type="application/json",
            filename=f"coverage-{job_id}.json",
        )

    @app.get(
        "/coverage/{job_id}/report",
        tags=["coverage"],
        summary="Alias of /coverage-jobs/{job_id}/report",
        include_in_schema=False,
    )
    async def coverage_report_alias(
        request: Request, job_id: str
    ) -> FileResponse:
        return await coverage_report(request, job_id)

    return app


app = create_app()
