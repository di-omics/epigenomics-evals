"""Small, dependency-light interval toolkit used by the reference solutions.

Coordinates follow BED: 0-based, half-open. Two intervals overlap when they
share at least 1 bp (a.start < b.end and a.end > b.start). That matches
GenomicRanges::countOverlaps on peaks read with ChIPseeker::readPeakFile,
which is what the paper's code uses.

Everything is numpy. No bedtools, no R.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

STANDARD_CHROMS = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")


@dataclass(frozen=True, eq=False)
class Intervals:
    """A set of genomic intervals with optional per-interval columns."""

    chrom: np.ndarray  # object array of chromosome names
    start: np.ndarray  # int64
    end: np.ndarray  # int64
    extra: pd.DataFrame | None = None
    _cache: dict = field(default_factory=dict, repr=False)

    def groups(self) -> dict[str, np.ndarray]:
        """Row indices per chromosome, computed once per set."""
        if "groups" not in self._cache:
            codes, names = pd.factorize(self.chrom)
            order = np.argsort(codes, kind="stable")
            bounds = np.searchsorted(codes[order], np.arange(len(names) + 1))
            self._cache["groups"] = {
                str(name): order[bounds[i]:bounds[i + 1]] for i, name in enumerate(names)
            }
        return self._cache["groups"]

    def __len__(self) -> int:
        return int(self.start.shape[0])

    @property
    def width(self) -> np.ndarray:
        return self.end - self.start

    def take(self, mask: np.ndarray) -> "Intervals":
        extra = None if self.extra is None else self.extra.loc[mask].reset_index(drop=True)
        return Intervals(self.chrom[mask], self.start[mask], self.end[mask], extra)


def read_bed(path: str | Path) -> Intervals:
    """Read a BED-like file (BED3+, narrowPeak, broadPeak, SEACR, gz or plain)."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", quotechar='"', low_memory=False)
    chrom = df[0].astype(str).to_numpy(dtype=object)
    start = df[1].to_numpy(dtype=np.int64)
    end = df[2].to_numpy(dtype=np.int64)
    extra = df.drop(columns=[0, 1, 2]).reset_index(drop=True) if df.shape[1] > 3 else None
    return Intervals(chrom, start, end, extra)


def _sorted_bounds(iv: Intervals) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per chromosome: all starts sorted, all ends sorted (independently)."""
    if "bounds" not in iv._cache:
        iv._cache["bounds"] = {
            chrom: (np.sort(iv.start[idx]), np.sort(iv.end[idx])) for chrom, idx in iv.groups().items()
        }
    return iv._cache["bounds"]


def count_overlaps(query: Intervals, subject: Intervals, closed: bool = False) -> np.ndarray:
    """For each query interval, the number of subject intervals it overlaps.

    A subject interval overlaps the query when its start lies before the query
    end and its end lies after the query start, so the count is
    (starts before the query end) minus (ends at or before the query start).

    closed=True treats coordinates as 1-based closed (book-ended intervals
    count as overlapping). It exists to measure how much that convention
    moves an answer; the default is the BED convention.
    """
    counts = np.zeros(len(query), dtype=np.int64)
    subj = _sorted_bounds(subject)
    for chrom, idx in query.groups().items():
        if chrom not in subj:
            continue
        starts, ends = subj[chrom]
        q_start, q_end = query.start[idx], query.end[idx]
        if closed:
            n_start_before_end = np.searchsorted(starts, q_end, side="right")
            n_end_before_start = np.searchsorted(ends, q_start, side="left")
        else:
            n_start_before_end = np.searchsorted(starts, q_end, side="left")
            n_end_before_start = np.searchsorted(ends, q_start, side="right")
        counts[idx] = n_start_before_end - n_end_before_start
    return counts


def overlaps_any(query: Intervals, subject: Intervals, closed: bool = False) -> np.ndarray:
    return count_overlaps(query, subject, closed=closed) > 0


def drop_blacklisted(iv: Intervals, blacklist: Intervals) -> Intervals:
    return iv.take(~overlaps_any(iv, blacklist))


def keep_standard(iv: Intervals) -> Intervals:
    """Keep chr1-22, chrX, chrY (drops chrM, random, unplaced, haplotype)."""
    return iv.take(pd.Series(iv.chrom).isin(STANDARD_CHROMS).to_numpy())


def clean_peaks(iv: Intervals, blacklist: Intervals) -> Intervals:
    """The paper's peak filter: drop blacklist overlaps, keep standard chromosomes."""
    return keep_standard(drop_blacklisted(iv, blacklist))


def merged_bp(iv: Intervals) -> int:
    """Total base pairs covered by the union of the intervals."""
    total = 0
    order = np.lexsort((iv.start, iv.chrom))
    chrom, start, end = iv.chrom[order], iv.start[order], iv.end[order]
    bounds = np.flatnonzero(np.r_[True, chrom[1:] != chrom[:-1], True])
    for a, b in zip(bounds[:-1], bounds[1:]):
        s, e = start[a:b], end[a:b]
        run_end = np.maximum.accumulate(e)
        new_block = np.r_[True, s[1:] > run_end[:-1]]
        block_id = np.cumsum(new_block) - 1
        block_start = s[new_block]
        block_end = np.zeros(block_start.shape[0], dtype=np.int64)
        np.maximum.at(block_end, block_id, e)
        total += int((block_end - block_start).sum())
    return total


def merge(iv: Intervals) -> Intervals:
    """Union of the intervals as non-overlapping blocks."""
    chroms, starts, ends = [], [], []
    order = np.lexsort((iv.start, iv.chrom))
    chrom, start, end = iv.chrom[order], iv.start[order], iv.end[order]
    bounds = np.flatnonzero(np.r_[True, chrom[1:] != chrom[:-1], True])
    for a, b in zip(bounds[:-1], bounds[1:]):
        s, e = start[a:b], end[a:b]
        run_end = np.maximum.accumulate(e)
        new_block = np.r_[True, s[1:] > run_end[:-1]]
        block_id = np.cumsum(new_block) - 1
        block_start = s[new_block]
        block_end = np.zeros(block_start.shape[0], dtype=np.int64)
        np.maximum.at(block_end, block_id, e)
        chroms.append(np.full(block_start.shape[0], chrom[a], dtype=object))
        starts.append(block_start)
        ends.append(block_end)
    return Intervals(np.concatenate(chroms), np.concatenate(starts), np.concatenate(ends))


def intersect_bp(a: Intervals, b: Intervals) -> int:
    """Base pairs covered by both sets (union of a intersected with union of b)."""
    return merged_bp(a) + merged_bp(b) - merged_bp(
        Intervals(np.r_[a.chrom, b.chrom], np.r_[a.start, b.start], np.r_[a.end, b.end])
    )


def subtract(a: Intervals, b: Intervals) -> Intervals:
    """Intervals of a that overlap nothing in b (whole-interval removal)."""
    return a.take(~overlaps_any(a, b))
