#!/usr/bin/env python3
"""Reference solution for the peak-set choice task. Standard library only.

    python3 reference_solution.py <workspace_dir>

Reads the three BED files in <workspace_dir>/inputs and prints the answer as
one JSON object. It reads nothing else: no answer key, no labels.

Method: sorted boundaries and binary search. For a query interval [qs, qe),
the subject intervals that overlap it are those that start before qe, minus
those that end at or before qs. Every subject interval that ends at or before
qs also starts before qe, so the difference is exactly the overlap count.
That holds when no interval has zero length, which the replay checks.
"""
from __future__ import annotations

import bisect
import gzip
import json
import sys
from pathlib import Path

COUNT_FIELDS = (
    "n_peaks_a",
    "n_peaks_b",
    "n_reference_peaks",
    "n_a_peaks_overlapping_reference",
    "n_b_peaks_overlapping_reference",
    "n_reference_peaks_overlapped_by_a",
    "n_reference_peaks_overlapped_by_b",
)
CHOICES = ("A", "B", "tie", "neither")


def read_bed3(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Chromosome -> list of (start, end), in file order. Reads .bed or .bed.gz."""
    opener = gzip.open if str(path).endswith(".gz") else open
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            by_chrom.setdefault(chrom, []).append((int(start), int(end)))
    return by_chrom


def n_intervals(bed: dict[str, list[tuple[int, int]]]) -> int:
    return sum(len(v) for v in bed.values())


def overlap_counts(query: dict, subject: dict) -> dict[str, list[int]]:
    """For each query interval, how many subject intervals share at least 1 bp with it."""
    out: dict[str, list[int]] = {}
    for chrom, intervals in query.items():
        other = subject.get(chrom, [])
        starts = sorted(s for s, _ in other)
        ends = sorted(e for _, e in other)
        out[chrom] = [
            bisect.bisect_left(starts, q_end) - bisect.bisect_right(ends, q_start)
            for q_start, q_end in intervals
        ]
    return out


def n_overlapping(query: dict, subject: dict) -> int:
    """Unique query intervals that overlap at least one subject interval."""
    return sum(1 for counts in overlap_counts(query, subject).values() for c in counts if c > 0)


def decide(counts: dict) -> str:
    """The lab's rule in whole numbers.

    A set qualifies when strictly more than 75% of its peaks overlap a
    reference peak: 4 * overlapping > 3 * total. Among qualifying sets the
    one that recovers more reference peaks wins.
    """
    ok_a = 4 * counts["n_a_peaks_overlapping_reference"] > 3 * counts["n_peaks_a"]
    ok_b = 4 * counts["n_b_peaks_overlapping_reference"] > 3 * counts["n_peaks_b"]
    if not ok_a and not ok_b:
        return "neither"
    if ok_a != ok_b:
        return "A" if ok_a else "B"
    got_a = counts["n_reference_peaks_overlapped_by_a"]
    got_b = counts["n_reference_peaks_overlapped_by_b"]
    if got_a == got_b:
        return "tie"
    return "A" if got_a > got_b else "B"


def solve(workspace: Path) -> dict:
    inputs = Path(workspace) / "inputs"
    a = read_bed3(inputs / "peaks_A.bed")
    b = read_bed3(inputs / "peaks_B.bed")
    ref = read_bed3(inputs / "reference_peaks.bed")
    answer = {
        "n_peaks_a": n_intervals(a),
        "n_peaks_b": n_intervals(b),
        "n_reference_peaks": n_intervals(ref),
        "n_a_peaks_overlapping_reference": n_overlapping(a, ref),
        "n_b_peaks_overlapping_reference": n_overlapping(b, ref),
        "n_reference_peaks_overlapped_by_a": n_overlapping(ref, a),
        "n_reference_peaks_overlapped_by_b": n_overlapping(ref, b),
    }
    answer["answer"] = decide(answer)
    return answer


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    print(json.dumps(solve(Path(sys.argv[1])), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
