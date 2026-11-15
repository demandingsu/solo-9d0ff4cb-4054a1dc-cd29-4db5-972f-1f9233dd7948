"""Tests for the coverage verification engine and API."""

from __future__ import annotations

import gzip
import json
import threading

import pytest
from fastapi.testclient import TestClient

from app import coverage as coverage_module
from app.coverage import (
    Reference,
    ReferenceError,
    parse_reference_upload,
    reverse_complement,
    run_coverage,
)
from app.main import create_app

from test_api import fastq_record, gz_bytes, gz_pair, read_gz, submit, wait_for  # noqa: E402

HIGH = "I" * 8


def q(seq: str) -> str:
    return "I" * len(seq)


def fasta(sequence: str, name: str = "ref", width: int = 70) -> str:
    lines = [f">{name}"]
    for start in range(0, len(sequence), width):
        lines.append(sequence[start : start + width])
    return "\n".join(lines) + "\n"


@pytest.fixture()
def client(tmp_path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as test_client:
        yield test_client


def run_qc_job(client, r1: str, r2: str, threshold: int = 20) -> str:
    response = submit(client, *gz_pair(r1, r2), threshold=threshold)
    assert response.status_code == 202, response.text
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "succeeded", done
    return job_id


def start_coverage(
    client,
    qc_job_id: str,
    reference: str,
    min_fragment: int,
    max_fragment: int,
    min_depth: int,
):
    return client.post(
        "/coverage-jobs",
        files={
            "reference": ("reference.fa", reference.encode("ascii"), "text/plain")
        },
        data={
            "qc_job_id": qc_job_id,
            "min_fragment": str(min_fragment),
            "max_fragment": str(max_fragment),
            "min_depth": str(min_depth),
        },
    )


def coverage(client, qc_job_id: str, reference: str, lo: int, hi: int, depth: int):
    response = start_coverage(client, qc_job_id, reference, lo, hi, depth)
    assert response.status_code == 202, response.text
    payload = response.json()
    done = wait_for(
        client, payload["job_id"], {"succeeded", "failed", "stopped"}, timeout=20
    )
    return done


# ---------------------------------------------------------------------------
# Pair placement rules
# ---------------------------------------------------------------------------


def test_pair_counts_and_intervals(client):
    # Reference: 20 A/C/G/T pattern so reads are unique.
    reference = (
        "ACGTTGCAACGTAGCTAGCT"
        "AGGCTAACGTTCGATCGTAG"
        "CTAGGATCCGATCGATCGTG"
        "ACGTTGCAACGTAGCTAGCT"
    )
    # R1 starts at 10 (forward). R2 must be the reverse complement of the
    # reference segment ending at position 30.
    ref = reference
    r1_seq = ref[10:18]                  # 8 bp forward
    r2_ref = ref[22:30]                  # 8 bp on reference
    r2_read = reverse_complement(r2_ref.encode()).decode()
    name = "p1"
    r1 = fastq_record(f"{name}/1", r1_seq, q(r1_seq))
    r2 = fastq_record(f"{name}/2", r2_read, q(r2_read))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 8, 100, 1)
    assert report["status"] == "succeeded", report
    assert report["total_pairs"] == 1
    assert report["placed_pairs"] == 1
    assert report["ambiguous_pairs"] == 0
    assert report["unplaced_pairs"] == 0
    assert report["low_coverage_intervals"] == [
        {"start": 0, "end": 10},
        {"start": 30, "end": len(ref)},
    ]


def test_overlapping_ends_counted_once(client):
    ref = (
        "ACGTTGCAAGCTTCGATCGA"
        "TTAGGCATCGATCGTAGCTA"
        "GGCCTTAACGTTAGCCATGA"
        "CTAGGTACCATTCGAGCTTA"
        "GCTAACCGTATGGCATTCGA"
    )
    # R1 1..21, R2 read maps 11..31 (overlap 10 bases).
    r1_seq = ref[0:21]
    r2_seq = reverse_complement(ref[10:31].encode()).decode()
    r1 = fastq_record("o/1", r1_seq, q(r1_seq))
    r2 = fastq_record("o/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 30, 40, 1)
    assert report["status"] == "succeeded", report
    assert report["placed_pairs"] == 1
    # Fragment [0, 31) is covered once; 31.. is uncovered.
    assert report["low_coverage_intervals"] == [{"start": 31, "end": 100}]


def test_ambiguous_pair_is_not_counted(client):
    # Reference has the same 8 bp block twice; a pair made of that block
    # therefore has more than one compliant placement.
    block = "ACGTACGT"
    ref = block + "TTTTTTTT" + block + "GGGGGGGG" + block
    r1 = fastq_record("a/1", block, q(block))
    r2 = fastq_record("a/2", reverse_complement(block.encode()).decode(), q(block))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 8, 20, 1)
    assert report["status"] == "succeeded", report
    assert report["ambiguous_pairs"] == 1
    assert report["placed_pairs"] == 0
    assert report["unplaced_pairs"] == 0
    # Nothing is covered when the pair is ambiguous.
    assert report["low_coverage_intervals"] == [
        {"start": 0, "end": len(ref)}
    ]


def test_unplaced_pairs(client):
    ref = "A" * 40
    r1 = fastq_record("u/1", "ACGTACGT", q("ACGTACGT"))
    r2 = fastq_record("u/2", "TTTTTTTT", q("TTTTTTTT"))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 8, 40, 1)
    assert report["status"] == "succeeded", report
    assert report["unplaced_pairs"] == 1
    assert report["placed_pairs"] == 0
    assert report["ambiguous_pairs"] == 0


def test_fragment_range_filters_placements(client):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGTTCGATCGTAG"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("r/1", r1_seq, q(r1_seq))
    r2 = fastq_record("r/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)
    # Fragment is 30 bp: inside range -> placed.
    inside = coverage(client, job, fasta(ref), 30, 30, 1)
    assert inside["placed_pairs"] == 1
    # Outside range -> unplaced.
    outside = coverage(client, job, fasta(ref), 31, 40, 1)
    assert outside["placed_pairs"] == 0
    assert outside["unplaced_pairs"] == 1


def test_min_depth_zero_reports_no_low_intervals(client):
    ref = "A" * 40
    r1 = fastq_record("z/1", "AAAA", "IIII")
    r2 = fastq_record("z/2", "TTTT", "IIII")
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 4, 40, 0)
    assert report["status"] == "succeeded", report
    assert report["low_coverage_intervals"] == []


def test_interval_merging(client):
    # Two placed pairs covering 0..10 and 10..20; depth >= 2 everywhere else
    # is not possible with only two reads, so use min_depth=2 to get one
    # merged interval covering the whole reference.
    ref = "".join(["ACGTACGTTG"] * 10)  # 100 bp, 10 bp repeat unit
    # Use unique blocks instead to keep placement exact.
    ref = (
        "ACGTTGCAAG" "CTAGCTAGCT" "AGGCTAACGT" "TCGATCGTAG"
        "CCTAGGATCC" "GATCGATCGT" "ACGCGTACGT" "AGCTAGCTAG"
        "TTACGGCAAT" "CCGTAACGTT"
    )
    r1a = ref[0:10]
    r2a = reverse_complement(ref[10:20].encode()).decode()
    # Second pair covers 20..30.
    r1b = ref[20:30]
    r2b = reverse_complement(ref[30:40].encode()).decode()
    r1 = fastq_record("m1/1", r1a, q(r1a)) + fastq_record("m2/1", r1b, q(r1b))
    r2 = fastq_record("m1/2", r2a, q(r2a)) + fastq_record("m2/2", r2b, q(r2b))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 20, 20, 2)
    assert report["status"] == "succeeded", report
    assert report["placed_pairs"] == 2
    # 0..40 has depth 1 < 2; 40..100 has depth 0 < 2; contiguous -> merged.
    assert report["low_coverage_intervals"] == [{"start": 0, "end": 100}]


def test_non_overlapping_ends_cover_the_fragment_span(client):
    # R1 at [0, 8), R2 at [22, 30): fragment [0, 30) is counted once as a
    # whole (the 14 bp gap between the ends belongs to the same fragment).
    ref = "ACGTTGCAACGTAGCTAGCTAGGCTAACGTTCGATCGTAG"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("g/1", r1_seq, q(r1_seq))
    r2 = fastq_record("g/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 30, 30, 1)
    assert report["status"] == "succeeded", report
    assert report["placed_pairs"] == 1
    assert report["low_coverage_intervals"] == [
        {"start": 30, "end": len(ref)}
    ]


def test_empty_published_pair_produces_full_reference_interval(client):
    # A QC job that discards every pair publishes empty FASTQ files.
    r1 = fastq_record("e/1", "ACGT", "IIII")
    r2 = fastq_record("e/2", "ACGT", "II##")
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded", done
    job = job_id
    report = coverage(client, job, fasta("A" * 10), 4, 20, 1)
    assert report["status"] == "succeeded", report
    assert report["total_pairs"] == 0
    assert report["low_coverage_intervals"] == [{"start": 0, "end": 10}]


# ---------------------------------------------------------------------------
# API validation
# ---------------------------------------------------------------------------


def test_reference_rules():
    assert parse_reference_upload(b">r\nACGT\n").sequence == b"ACGT"
    assert parse_reference_upload(b">r description\nacgt\n").sequence == b"ACGT"
    assert parse_reference_upload(gzip.compress(b">r\nACGT\n")).sequence == b"ACGT"
    with pytest.raises(ReferenceError):
        parse_reference_upload(b"ACGT\n")
    with pytest.raises(ReferenceError):
        parse_reference_upload(b">r\nACGN\n")
    with pytest.raises(ReferenceError):
        parse_reference_upload(b">r\nACGT\n>r2\nGGGG\n")
    with pytest.raises(ReferenceError):
        parse_reference_upload(b">r\n\n")
    with pytest.raises(ReferenceError):
        parse_reference_upload(b">r\n" + b"A" * 50_001 + b"\n")
    # 50,000 exactly is accepted.
    parsed = parse_reference_upload(b">r\n" + b"A" * 50_000 + b"\n")
    assert len(parsed.sequence) == 50_000


def test_coverage_requires_successful_qc_job(client):
    r1 = fastq_record("a/1", "ACGT", "IIII")
    r2 = fastq_record("a/2", "ACGT", "III")  # seq/qual length mismatch
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job = response.json()["job_id"]
    failed = wait_for(client, job, {"succeeded", "failed"})
    assert failed["status"] == "failed"
    response = start_coverage(client, job, fasta("ACGT"), 4, 20, 1)
    assert response.status_code == 409

    assert (
        start_coverage(client, "missing-job", fasta("ACGT"), 4, 20, 1).status_code
        == 404
    )


def test_coverage_validates_parameters(client):
    r1 = fastq_record("a/1", "AAAA", "IIII")
    r2 = fastq_record("a/2", "TTTT", "IIII")
    job = run_qc_job(client, r1, r2)

    bad_range = start_coverage(client, job, fasta("A" * 20), 10, 5, 1)
    assert bad_range.status_code == 422
    zero_min = start_coverage(client, job, fasta("A" * 20), 0, 5, 1)
    assert zero_min.status_code == 422
    negative_depth = start_coverage(client, job, fasta("A" * 20), 4, 10, -1)
    assert negative_depth.status_code == 422
    bad_reference = start_coverage(client, job, ">r\nACGN\n", 4, 10, 1)
    assert bad_reference.status_code == 422
    two_records = start_coverage(client, job, ">a\nACGT\n>b\nACGT\n", 4, 10, 1)
    assert two_records.status_code == 422
    missing_reference = client.post(
        "/coverage-jobs",
        data={"qc_job_id": job, "min_fragment": "4", "max_fragment": "10", "min_depth": "1"},
    )
    assert missing_reference.status_code == 422
    missing_depth = client.post(
        "/coverage-jobs",
        files={"reference": ("r.fa", b">r\nACGT\n", "text/plain")},
        data={"qc_job_id": job, "min_fragment": "4", "max_fragment": "10"},
    )
    assert missing_depth.status_code == 422


def test_reference_too_long_rejected_by_api(client):
    r1 = fastq_record("a/1", "AAAA", "IIII")
    r2 = fastq_record("a/2", "TTTT", "IIII")
    job = run_qc_job(client, r1, r2)
    response = start_coverage(client, job, fasta("A" * 50_001), 4, 10, 1)
    assert response.status_code == 422


def test_coverage_job_path_route_and_aliases(client):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("a/1", r1_seq, q(r1_seq))
    r2 = fastq_record("a/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)
    response = client.post(
        f"/jobs/{job}/coverage",
        files={"reference": ("r.fa", fasta(ref).encode(), "text/plain")},
        data={"min_fragment": "30", "max_fragment": "30", "min_depth": "1"},
    )
    assert response.status_code == 202, response.text
    cov_id = response.json()["job_id"]
    done = wait_for(client, cov_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "succeeded", done
    assert done["kind"] == "coverage"
    assert done["qc_job_id"] == job
    assert client.get(f"/coverage/{cov_id}").json()["status"] == "succeeded"
    report = client.get(f"/coverage/{cov_id}/report")
    assert report.status_code == 200
    assert report.headers["content-type"] == "application/json"
    assert json.loads(report.content)["placed_pairs"] == 1
    # A QC download on a coverage job is not available.
    assert client.get(f"/jobs/{cov_id}/download/r1").status_code == 409


def test_qc_downloads_still_work_after_coverage(client):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("a/1", r1_seq, q(r1_seq))
    r2 = fastq_record("a/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)
    report = coverage(client, job, fasta(ref), 30, 30, 1)
    assert report["status"] == "succeeded"
    # The source QC job is untouched and still downloadable.
    assert client.get(f"/jobs/{job}/download/r1").status_code == 200
    assert read_gz(client.get(f"/jobs/{job}/download/r1").content) == r1


# ---------------------------------------------------------------------------
# Failure / stop semantics and streaming
# ---------------------------------------------------------------------------


def test_corrupt_source_gzip_fails_without_report(client, monkeypatch):
    """A corrupted tail must fail the coverage job and leave no report."""
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("a/1", r1_seq, q(r1_seq))
    r2 = fastq_record("a/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)

    manager = client.app.state.manager
    source = manager.get(job)
    # Corrupt the published R2 after the QC job completed.
    data = source.output_r2.read_bytes()
    source.output_r2.write_bytes(data[:-6])

    response = start_coverage(client, job, fasta(ref), 30, 30, 1)
    assert response.status_code == 202
    cov_id = response.json()["job_id"]
    done = wait_for(client, cov_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "failed", done
    assert done["error"]["code"] == "gzip_corrupt"
    assert done["low_coverage_intervals"] is None
    assert client.get(f"/coverage-jobs/{cov_id}/report").status_code == 409
    directory = manager.jobs_dir / cov_id
    assert not directory.exists()


def test_coverage_stop_cleans_and_publishes_nothing(client, monkeypatch):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("a/1", r1_seq, q(r1_seq))
    r2 = fastq_record("a/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)

    import app.jobs as jobs_module

    entered = threading.Event()
    release = threading.Event()
    real_run_coverage = jobs_module.run_coverage

    def blocking(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return real_run_coverage(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_coverage", blocking)
    response = start_coverage(client, job, fasta(ref), 30, 30, 1)
    cov_id = response.json()["job_id"]
    assert entered.wait(10)
    stopped = client.post(f"/coverage-jobs/{cov_id}/stop")
    assert stopped.status_code == 200
    assert stopped.json()["status"] == "stopping"
    assert client.get(f"/coverage-jobs/{cov_id}").json()["status"] == "stopping"
    assert client.get(f"/coverage-jobs/{cov_id}/report").status_code == 409
    release.set()
    done = wait_for(client, cov_id, {"stopped", "succeeded", "failed"})
    assert done["status"] == "stopped", done
    assert done["low_coverage_intervals"] is None
    assert done["error"] is None
    directory = client.app.state.manager.jobs_dir / cov_id
    assert not directory.exists()


def test_other_jobs_query_during_coverage(client, monkeypatch):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1a = fastq_record("a/1", r1_seq, q(r1_seq))
    r2a = fastq_record("a/2", r2_seq, q(r2_seq))
    first = run_qc_job(client, r1a, r2a)

    import app.jobs as jobs_module

    entered = threading.Event()
    release = threading.Event()
    real_run_coverage = jobs_module.run_coverage

    def blocking(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return real_run_coverage(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_coverage", blocking)
    response = start_coverage(client, first, fasta(ref), 30, 30, 1)
    cov_id = response.json()["job_id"]
    assert entered.wait(10)
    # While the coverage job is blocked, both QC queries and QC downloads
    # remain responsive.
    assert client.get(f"/jobs/{first}").status_code == 200
    assert client.get(f"/jobs/{first}/download/r1").status_code == 200
    assert client.get("/jobs").status_code == 200
    assert client.get("/health").status_code == 200
    release.set()
    done = wait_for(client, cov_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "succeeded", done


def test_coverage_requires_running_qc_job(client, monkeypatch):
    import app.jobs as jobs_module

    release = threading.Event()
    entered = threading.Event()
    real_run_qc = jobs_module.run_qc

    def blocking_run_qc(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return real_run_qc(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_qc", blocking_run_qc)
    r1 = fastq_record("a/1", "ACGT", "IIII")
    r2 = fastq_record("a/2", "ACGT", "IIII")
    job = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    assert entered.wait(10)
    response = start_coverage(client, job, fasta("ACGT"), 4, 20, 1)
    assert response.status_code == 409
    release.set()
    wait_for(client, job, {"succeeded", "failed"})


def test_computation_error_publishes_nothing(client, monkeypatch):
    ref = "ACGTTGCAAGCTAGCTAGCTAGGCTAACGT"
    r1_seq = ref[0:8]
    r2_seq = reverse_complement(ref[22:30].encode()).decode()
    r1 = fastq_record("a/1", r1_seq, q(r1_seq))
    r2 = fastq_record("a/2", r2_seq, q(r2_seq))
    job = run_qc_job(client, r1, r2)

    import app.jobs as jobs_module

    def boom(*args, **kwargs):
        raise RuntimeError("simulated calculation failure")

    monkeypatch.setattr(jobs_module, "run_coverage", boom)
    response = start_coverage(client, job, fasta(ref), 30, 30, 1)
    assert response.status_code == 202
    cov_id = response.json()["job_id"]
    done = wait_for(client, cov_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "failed", done
    assert done["error"]["code"] == "internal_error"
    assert done["low_coverage_intervals"] is None
    assert client.get(f"/coverage-jobs/{cov_id}/report").status_code == 409
    assert not (client.app.state.manager.jobs_dir / cov_id).exists()


def test_coverage_engine_streams(tmp_path, monkeypatch):
    """The engine reads lines incrementally and never whole files."""
    reference = Reference(name="r", sequence=b"ACGT" * 25)
    r1_text = "".join(fastq_record(f"r{i}/1", "ACGT", "IIII") for i in range(50))
    r2_text = "".join(fastq_record(f"r{i}/2", "ACGT", "IIII") for i in range(50))
    (tmp_path / "R1.gz").write_bytes(gz_bytes(r1_text))
    (tmp_path / "R2.gz").write_bytes(gz_bytes(r2_text))

    calls = {"readline": 0}
    original = coverage_module._GzipLineReader.readline

    def counting(self):
        calls["readline"] += 1
        return original(self)

    monkeypatch.setattr(coverage_module._GzipLineReader, "readline", counting)
    result = run_coverage(
        tmp_path / "R1.gz",
        tmp_path / "R2.gz",
        reference,
        4,
        4,
        1,
        threading.Event(),
    )
    assert result.total_pairs == 50
    # 4 lines per record per side plus one EOF probe per side.
    assert calls["readline"] == 50 * 4 * 2 + 2


def test_docs_list_coverage_endpoints(client):
    schema = client.get("/openapi.json").json()
    paths = schema["paths"]
    assert "/coverage-jobs" in paths
    assert "/coverage-jobs/{job_id}" in paths
    assert "/coverage-jobs/{job_id}/report" in paths
    assert "/jobs/{job_id}/coverage" in paths
    assert "/jobs/{job_id}/download/{side}" in paths
