"""Streaming exact-placement coverage verification for targeted sequencing.

The engine consumes the *published* (quality-filtered) ``R1``/``R2`` gzip
files of a successful QC job, pair by pair.  For each pair:

* R1 is searched in its original (forward) orientation;
* R2 is reverse-complemented and then searched.

A pair contributes to coverage only when exactly one *compliant* combination
of exact placements exists.  A combination ``(p, q)`` is compliant when

* ``p <= q`` (the R1 start is not to the right of the R2 start);
* the fragment ``[p, q + len(R2rc))`` contains both complete reads;
* the fragment length is inside the inclusive ``[min_fragment, max_fragment]``
  range.

No compliant combination means *unplaced*; two or more mean *ambiguous*.
Only uniquely placed pairs increase depth, and each such pair deposits depth
exactly once over its whole fragment ``[R1 start, R2 end)``.  Where the two
ends overlap the bases are therefore counted once, never twice.

The returned low-coverage intervals are maximal contiguous stretches of the
reference whose depth is strictly below ``min_depth``.  They use 0-based,
half-open coordinates (``start`` inclusive, ``end`` exclusive).

Read processing is streaming: at most one pair is held in memory at a time.
Only the reference index and a depth-difference array scale with the
reference, which is capped at 50,000 bases.
"""

from __future__ import annotations

import gzip
import io
import threading
from bisect import bisect_left, bisect_right
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .qc import (
    JobStopped,
    QCError,
    _GzipLineReader,
    _parse_record_name,
    _strip_eol,
)

# A submitted reference FASTA must contain exactly one sequence and that
# sequence may not be longer than 50 kb.
MAX_REFERENCE_BASES = 50_000
# Compressed or plain FASTA upload cap (50 kb of bases plus line wrapping and
# a header fits comfortably, even with one base per line).
MAX_REFERENCE_UPLOAD_BYTES = 8 * 1024 * 1024
MAX_REFERENCE_TEXT_BYTES = 8 * 1024 * 1024

# Index width for exact-placement candidate lookup.
KMER_SIZE = 16
OCCURRENCE_CACHE_ENTRIES = 4_096
# Do not cache highly repetitive patterns: their placement tuples can be large.
OCCURRENCE_CACHE_MAX_HITS = 64
# Above this many k-mer candidates a C-level ``bytes.find`` scan is cheaper
# than testing every candidate position in Python.
FIND_SCAN_THRESHOLD = 32

PROGRESS_EVERY_PAIRS = 256

_ACGT = frozenset(b"ACGT")
_COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


class ReferenceError(ValueError):
    """Raised when an uploaded reference FASTA is not acceptable."""


@dataclass(frozen=True)
class Reference:
    """One validated reference sequence."""

    name: str
    sequence: bytes


def reverse_complement(sequence: bytes) -> bytes:
    """Reverse-complement a read; non-ACGT bytes are left untouched."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _decompress_reference(raw: bytes) -> bytes:
    if raw[:2] != b"\x1f\x8b":
        return raw
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as handle:
            data = handle.read(MAX_REFERENCE_TEXT_BYTES + 1)
    except (EOFError, OSError, gzip.BadGzipFile, ValueError) as exc:
        raise ReferenceError(
            "reference is gzip-compressed but the stream is corrupt"
        ) from exc
    if len(data) > MAX_REFERENCE_TEXT_BYTES:
        raise ReferenceError(
            f"reference text exceeds {MAX_REFERENCE_TEXT_BYTES} bytes"
        )
    return data


def parse_reference_upload(raw: bytes) -> Reference:
    """Parse and validate a single-record, ACGT-only FASTA upload."""
    if not raw:
        raise ReferenceError("reference FASTA is empty")
    if len(raw) > MAX_REFERENCE_UPLOAD_BYTES:
        raise ReferenceError(
            f"reference upload exceeds {MAX_REFERENCE_UPLOAD_BYTES} bytes"
        )

    text = _decompress_reference(raw)
    if len(text) > MAX_REFERENCE_TEXT_BYTES:
        raise ReferenceError(
            f"reference text exceeds {MAX_REFERENCE_TEXT_BYTES} bytes"
        )
    try:
        text.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ReferenceError("reference FASTA must be ASCII text") from exc

    name: str | None = None
    parts: list[bytes] = []
    records: list[tuple[str, bytes]] = []

    for line_number, raw_line in enumerate(text.split(b"\n"), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(b">"):
            if name is not None:
                records.append((name, b"".join(parts)))
            elif parts:
                raise ReferenceError(
                    f"reference FASTA has sequence data before the first "
                    f"header (line {line_number})"
                )
            header = line[1:].strip()
            tokens = header.split()
            name = (
                tokens[0].decode("utf-8", "replace")
                if tokens
                else f"reference_{len(records) + 1}"
            )
            parts = []
            continue
        if name is None:
            raise ReferenceError(
                f"reference FASTA has sequence data before the first "
                f"header (line {line_number})"
            )
        # Whitespace inside a sequence line is tolerated and removed; lower
        # case is normalised to upper case.
        chunk = b"".join(line.split()).upper()
        invalid = set(chunk) - _ACGT
        if invalid:
            shown = ", ".join(
                repr(chr(byte)) for byte in sorted(invalid)[:5]
            )
            raise ReferenceError(
                f"reference sequence contains non-ACGT character(s) "
                f"{shown} (line {line_number})"
            )
        parts.append(chunk)

    if name is None:
        raise ReferenceError("reference FASTA has no header line")
    records.append((name, b"".join(parts)))

    if len(records) != 1:
        raise ReferenceError(
            "reference FASTA must contain exactly one sequence "
            f"(found {len(records)} records)"
        )
    references_name, sequence = records[0]
    if not sequence:
        raise ReferenceError("reference sequence is empty")
    if len(sequence) > MAX_REFERENCE_BASES:
        raise ReferenceError(
            f"reference sequence length {len(sequence)} exceeds the "
            f"maximum of {MAX_REFERENCE_BASES} bases"
        )
    return Reference(name=references_name, sequence=sequence)


class ReferenceIndex:
    """Exact substring lookup over one reference sequence."""

    def __init__(self, sequence: bytes) -> None:
        self.sequence = sequence
        self._k = min(KMER_SIZE, len(sequence))
        index: dict[bytes, list[int]] = {}
        if self._k:
            for position in range(len(sequence) - self._k + 1):
                key = sequence[position : position + self._k]
                index.setdefault(key, []).append(position)
        self._index = index
        self._cache: OrderedDict[bytes, tuple[int, ...]] = OrderedDict()

    def placements(self, pattern: bytes) -> tuple[int, ...]:
        """Return every start position where *pattern* occurs exactly."""
        if not pattern or len(pattern) > len(self.sequence):
            return ()

        cached = self._cache.get(pattern)
        if cached is not None:
            self._cache.move_to_end(pattern)
            return cached

        if len(pattern) < self._k:
            found: list[int] = []
            start = self.sequence.find(pattern)
            while start != -1:
                found.append(start)
                start = self.sequence.find(pattern, start + 1)
            result = tuple(found)
        else:
            candidates = self._index.get(pattern[: self._k])
            if not candidates:
                result = ()
            elif len(candidates) == 1:
                position = candidates[0]
                result = (
                    (position,)
                    if self.sequence.startswith(pattern, position)
                    else ()
                )
            elif len(candidates) > FIND_SCAN_THRESHOLD:
                # Repetitive k-mer: a single C-level scan beats probing each
                # candidate in Python and yields the positions in order.
                found = []
                start = self.sequence.find(pattern)
                while start != -1:
                    found.append(start)
                    start = self.sequence.find(pattern, start + 1)
                result = tuple(found)
            else:
                result = tuple(
                    position
                    for position in candidates
                    if self.sequence.startswith(pattern, position)
                )

        if len(result) <= OCCURRENCE_CACHE_MAX_HITS:
            self._cache[pattern] = result
            if len(self._cache) > OCCURRENCE_CACHE_ENTRIES:
                self._cache.popitem(last=False)
        return result


def _compliant_combinations(
    r1_placements: tuple[int, ...],
    r1_length: int,
    r2_placements: tuple[int, ...],
    r2_length: int,
    min_fragment: int,
    max_fragment: int,
) -> tuple[int, tuple[int, int] | None]:
    """Count compliant ``(R1 start, R2 start)`` combinations (capped at 2).

    Returns ``(count, first_combination)``.  ``count`` is at most 2; callers
    only need to distinguish zero, one and many.
    """
    if not r1_placements or not r2_placements:
        return 0, None

    count = 0
    first: tuple[int, int] | None = None
    for p in r1_placements:
        # Fragment = [p, q + r2_length).
        #   q >= p                         (R1 start not after R2 start)
        #   q + r2_length >= p + r1_length (fragment contains all of R1)
        #   min_fragment <= q + r2_length - p <= max_fragment
        low = max(p, p + r1_length - r2_length, p + min_fragment - r2_length)
        high = p + max_fragment - r2_length
        if high < low:
            continue
        left = bisect_left(r2_placements, low)
        right = bisect_right(r2_placements, high)
        if right <= left:
            continue
        if count == 0:
            first = (p, r2_placements[left])
        count += right - left
        if count > 1:
            return count, first
    return count, first


def _check_published_record(
    lines: list[bytes],
    side: str,
    pair_number: int,
    expected_suffix: bytes,
) -> tuple[bytes, bytes]:
    """Light validation of a record from an already-published FASTQ file.

    The source files were produced by a successful QC job, so quality bytes
    have already been validated.  Re-checking the structural invariants keeps
    a hand-modified or corrupted file from silently producing nonsense.
    Returns ``(identifier, sequence)``.
    """
    header, sequence_line, plus_line, quality_line = lines

    identifier = _parse_record_name(
        header, expected_suffix, side, pair_number
    )
    if not plus_line.startswith(b"+"):
        raise QCError(
            "malformed_record",
            f"{side}: record {pair_number} third line does not start with '+'",
        )
    sequence = _strip_eol(sequence_line)
    quality = _strip_eol(quality_line)
    if len(sequence) != len(quality):
        raise QCError(
            "seq_qual_length_mismatch",
            f"{side}: record {pair_number} sequence length "
            f"({len(sequence)}) differs from quality length ({len(quality)})",
        )
    return identifier, sequence


@dataclass(frozen=True)
class LowCoverageInterval:
    start: int
    end: int

    def to_dict(self) -> dict[str, int]:
        return {"start": self.start, "end": self.end}


@dataclass(frozen=True)
class CoverageResult:
    reference_length: int
    total_pairs: int
    placed_pairs: int
    ambiguous_pairs: int
    unplaced_pairs: int
    low_coverage_intervals: tuple[LowCoverageInterval, ...]

    def to_dict(self) -> dict:
        return {
            "reference_length": self.reference_length,
            "total_pairs": self.total_pairs,
            "placed_pairs": self.placed_pairs,
            "ambiguous_pairs": self.ambiguous_pairs,
            "unplaced_pairs": self.unplaced_pairs,
            "low_coverage_intervals": [
                interval.to_dict() for interval in self.low_coverage_intervals
            ],
        }


ProgressCallback = Callable[[int], None]


def _low_coverage_intervals(
    difference: list[int], min_depth: int
) -> tuple[LowCoverageInterval, ...]:
    """Turn a depth-difference array into maximal low-depth intervals."""
    intervals: list[LowCoverageInterval] = []
    depth = 0
    in_low = False
    low_start = 0
    length = len(difference) - 1  # trailing slot only holds a decrement
    for position in range(length):
        depth += difference[position]
        is_low = depth < min_depth
        if is_low and not in_low:
            low_start = position
            in_low = True
        elif not is_low and in_low:
            intervals.append(LowCoverageInterval(low_start, position))
            in_low = False
    if in_low:
        intervals.append(LowCoverageInterval(low_start, length))
    return tuple(intervals)


def run_coverage(
    r1_path: Path,
    r2_path: Path,
    reference: Reference,
    min_fragment: int,
    max_fragment: int,
    min_depth: int,
    stop_event: threading.Event,
    *,
    progress_cb: ProgressCallback | None = None,
) -> CoverageResult:
    """Stream a published R1/R2 pair and compute its coverage report.

    Nothing is written by this function; the caller publishes the report only
    after a successful return.  A truncated or otherwise corrupt gzip stream,
    a malformed record or a stop request raises instead, so a partial report
    can never escape.
    """
    if min_fragment < 1:
        raise QCError(
            "invalid_fragment_range", "minimum fragment length must be >= 1"
        )
    if max_fragment < min_fragment:
        raise QCError(
            "invalid_fragment_range",
            "maximum fragment length must be >= minimum fragment length",
        )
    if min_depth < 0:
        raise QCError("invalid_min_depth", "minimum depth must be >= 0")

    sequence = reference.sequence
    reference_length = len(sequence)
    if reference_length == 0 or reference_length > MAX_REFERENCE_BASES:
        raise QCError(
            "invalid_reference",
            f"reference length must be between 1 and {MAX_REFERENCE_BASES}",
        )

    index = ReferenceIndex(sequence)
    difference = [0] * (reference_length + 1)

    placed = 0
    ambiguous = 0
    unplaced = 0
    pair_number = 0

    reader1 = reader2 = None
    try:
        reader1 = _GzipLineReader(r1_path, "R1", stop_event)
        reader2 = _GzipLineReader(r2_path, "R2", stop_event)

        while True:
            if stop_event.is_set():
                raise JobStopped()

            header1 = reader1.readline()
            header2 = reader2.readline()
            if not header1 and not header2:
                break
            if not header1 or not header2:
                missing = "R2" if header1 else "R1"
                other = "R1" if missing == "R2" else "R2"
                raise QCError(
                    "pair_count_mismatch",
                    f"{other} has more records than {missing}; "
                    "the two ends contain different numbers of reads",
                )

            pair_number += 1
            lines1 = [
                header1,
                reader1.readline(),
                reader1.readline(),
                reader1.readline(),
            ]
            lines2 = [
                header2,
                reader2.readline(),
                reader2.readline(),
                reader2.readline(),
            ]
            for side, lines in (("R1", lines1), ("R2", lines2)):
                if any(line == b"" for line in lines[1:]):
                    raise QCError(
                        "malformed_record",
                        f"{side}: record {pair_number} has fewer than four "
                        "lines (truncated FASTQ record)",
                    )

            id1, sequence1 = _check_published_record(
                lines1, "R1", pair_number, b"/1"
            )
            id2, raw_sequence2 = _check_published_record(
                lines2, "R2", pair_number, b"/2"
            )
            if id1 != id2:
                raise QCError(
                    "pair_name_mismatch",
                    f"record {pair_number}: R1 read name "
                    f"{id1.decode('utf-8', 'replace')!r} does not match R2 "
                    f"read name {id2.decode('utf-8', 'replace')!r}",
                )

            # The reference is upper case; normalise reads the same way so
            # soft-masked (lower-case) input is matched consistently.
            sequence1 = sequence1.upper()
            sequence2 = reverse_complement(raw_sequence2.upper())

            placements1 = index.placements(sequence1)
            placements2 = index.placements(sequence2)
            combinations, first = _compliant_combinations(
                placements1,
                len(sequence1),
                placements2,
                len(sequence2),
                min_fragment,
                max_fragment,
            )

            if combinations == 0:
                unplaced += 1
            elif combinations > 1:
                ambiguous += 1
            else:
                placed += 1
                assert first is not None
                fragment_start, r2_start = first
                # The fragment is [R1 start, R2 end).  Compliant placement
                # guarantees R2 end >= R1 end, so one difference pair covers
                # the fragment exactly once -- overlap is never double-counted.
                fragment_end = r2_start + len(sequence2)
                difference[fragment_start] += 1
                difference[fragment_end] -= 1

            if progress_cb is not None and (
                pair_number % PROGRESS_EVERY_PAIRS == 0
            ):
                progress_cb(pair_number)

        if stop_event.is_set():
            raise JobStopped()
    finally:
        for reader in (reader1, reader2):
            if reader is not None:
                reader.close()

    return CoverageResult(
        reference_length=reference_length,
        total_pairs=pair_number,
        placed_pairs=placed,
        ambiguous_pairs=ambiguous,
        unplaced_pairs=unplaced,
        low_coverage_intervals=_low_coverage_intervals(difference, min_depth),
    )
