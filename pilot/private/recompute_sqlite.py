#!/usr/bin/env python3
"""Second, independent recomputation of the answer: a range join in SQLite.

    python3 recompute_sqlite.py <workspace_dir>

Shares no code with reference_solution.py. The files are loaded into an
in-memory database and each direction is one declarative query: join the two
tables on "same chromosome, starts before the other ends, ends after the other
starts", then count the distinct peaks on each side. The same query also
returns the number of overlapping pairs, which is what a careless count gives.
The decision is made with exact fractions instead of cross-multiplied integers.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import sys
from fractions import Fraction
from pathlib import Path

BAR = Fraction(3, 4)


def load(con: sqlite3.Connection, table: str, path: Path) -> None:
    con.execute(f"CREATE TABLE {table} (id INTEGER PRIMARY KEY, chrom TEXT, s INTEGER, e INTEGER)")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        rows = [
            (parts[0], int(parts[1]), int(parts[2]))
            for parts in (line.rstrip("\n").split("\t") for line in handle if line.strip())
        ]
    con.executemany(f"INSERT INTO {table} (chrom, s, e) VALUES (?, ?, ?)", rows)
    con.execute(f"CREATE INDEX {table}_pos ON {table} (chrom, s)")


def joined_counts(con: sqlite3.Connection, x: str, r: str) -> tuple[int, int, int]:
    """(distinct x peaks with a partner, distinct r peaks with a partner, overlapping pairs).

    The last condition only narrows the index scan. An r interval that overlaps
    x ends after x starts, so it cannot start more than the widest r interval
    before x starts.
    """
    widest = con.execute(f"SELECT MAX(e - s) FROM {r}").fetchone()[0] or 0
    row = con.execute(
        f"""
        SELECT COUNT(DISTINCT x.id), COUNT(DISTINCT r.id), COUNT(*)
        FROM {x} AS x CROSS JOIN {r} AS r
        WHERE r.chrom = x.chrom
          AND r.s < x.e
          AND r.e > x.s
          AND r.s > x.s - ?
        """,
        (widest,),
    ).fetchone()
    return int(row[0]), int(row[1]), int(row[2])


def recompute(workspace: Path) -> dict:
    inputs = Path(workspace) / "inputs"
    con = sqlite3.connect(":memory:")
    load(con, "a", inputs / "peaks_A.bed")
    load(con, "b", inputs / "peaks_B.bed")
    load(con, "ref", inputs / "reference_peaks.bed")
    n_a, n_b, n_ref = (con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in ("a", "b", "ref"))
    a_hit, ref_by_a, pairs_a = joined_counts(con, "a", "ref")
    b_hit, ref_by_b, pairs_b = joined_counts(con, "b", "ref")
    con.close()

    def qualifies(hit: int, total: int) -> bool:
        return total > 0 and Fraction(hit, total) > BAR

    ok_a, ok_b = qualifies(a_hit, n_a), qualifies(b_hit, n_b)
    if not ok_a and not ok_b:
        choice = "neither"
    elif ok_a and not ok_b:
        choice = "A"
    elif ok_b and not ok_a:
        choice = "B"
    elif Fraction(ref_by_a, n_ref) == Fraction(ref_by_b, n_ref):
        choice = "tie"
    else:
        choice = "A" if Fraction(ref_by_a, n_ref) > Fraction(ref_by_b, n_ref) else "B"
    return {
        "answer": {
            "n_peaks_a": n_a,
            "n_peaks_b": n_b,
            "n_reference_peaks": n_ref,
            "n_a_peaks_overlapping_reference": a_hit,
            "n_b_peaks_overlapping_reference": b_hit,
            "n_reference_peaks_overlapped_by_a": ref_by_a,
            "n_reference_peaks_overlapped_by_b": ref_by_b,
            "answer": choice,
        },
        "overlapping_pairs": {"a": pairs_a, "b": pairs_b},
    }


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    print(json.dumps(recompute(Path(sys.argv[1])), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
