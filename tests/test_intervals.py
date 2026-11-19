"""Unit tests for the interval toolkit. No data needed."""
import numpy as np

from epigenomics_evals.intervals import (
    Intervals,
    clean_peaks,
    count_overlaps,
    intersect_bp,
    keep_standard,
    merge,
    merged_bp,
    overlaps_any,
    subtract,
)


def iv(rows):
    chrom = np.array([r[0] for r in rows], dtype=object)
    start = np.array([r[1] for r in rows], dtype=np.int64)
    end = np.array([r[2] for r in rows], dtype=np.int64)
    return Intervals(chrom, start, end)


def test_overlap_is_at_least_one_bp_half_open():
    q = iv([("chr1", 100, 200)])
    assert count_overlaps(q, iv([("chr1", 199, 300)]))[0] == 1
    assert count_overlaps(q, iv([("chr1", 200, 300)]))[0] == 0  # book-ended
    assert count_overlaps(q, iv([("chr1", 0, 100)]))[0] == 0  # book-ended
    assert count_overlaps(q, iv([("chr1", 0, 101)]))[0] == 1


def test_closed_convention_counts_book_ended_intervals():
    q = iv([("chr1", 100, 200)])
    assert count_overlaps(q, iv([("chr1", 200, 300)]), closed=True)[0] == 1
    assert count_overlaps(q, iv([("chr1", 201, 300)]), closed=True)[0] == 0


def test_counts_are_per_chromosome_and_handle_nesting():
    q = iv([("chr1", 0, 1000), ("chr2", 0, 1000), ("chr3", 0, 10)])
    s = iv([("chr1", 10, 20), ("chr1", 15, 900), ("chr1", 999, 2000), ("chr2", 1000, 2000)])
    assert count_overlaps(q, s).tolist() == [3, 0, 0]
    assert overlaps_any(q, s).tolist() == [True, False, False]


def test_counts_match_brute_force_on_random_intervals():
    rng = np.random.default_rng(0)

    def random_set(n):
        start = rng.integers(0, 5000, n)
        return iv([(f"chr{rng.integers(1, 3)}", int(s), int(s + rng.integers(1, 400))) for s in start])

    q, s = random_set(300), random_set(300)
    brute = [
        sum(1 for j in range(len(s)) if s.chrom[j] == q.chrom[i]
            and s.start[j] < q.end[i] and s.end[j] > q.start[i])
        for i in range(len(q))
    ]
    assert count_overlaps(q, s).tolist() == brute


def test_merged_and_intersect_bp():
    a = iv([("chr1", 0, 100), ("chr1", 50, 150), ("chr2", 0, 10)])
    assert merged_bp(a) == 160
    assert len(merge(a)) == 2
    b = iv([("chr1", 100, 200), ("chr2", 5, 50)])
    assert intersect_bp(a, b) == 50 + 5


def test_filters():
    peaks = iv([("chr1", 0, 100), ("chr1", 500, 600), ("chrM", 0, 100), ("chrUn_gl000220", 0, 50)])
    blacklist = iv([("chr1", 90, 120)])
    assert len(keep_standard(peaks)) == 2
    kept = clean_peaks(peaks, blacklist)
    assert len(kept) == 1 and kept.start[0] == 500
    assert len(subtract(peaks, blacklist)) == 3
