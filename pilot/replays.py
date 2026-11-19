"""Scripted replays: one correct solution and the wrong analyses a grader must reject.

Every replay is a function of the workspace folder only. None of them reads
the answer key or the paper labels; each one computes its own answer from the
3 staged files, the way an agent with that particular misconception would.

kind:
  correct       must pass on every instance
  wrong_counts  reports wrong counts, so it must fail on every instance
  wrong_logic   counts are right, the decision logic is wrong. Whether that
                changes the final choice depends on the instance, which is the
                point of running it on more than one.
  convention    another coordinate convention
  malformed     not a valid answer
"""
from __future__ import annotations

import bisect
import importlib.util
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("reference_solution", HERE / "private" / "reference_solution.py")
ref = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ref)


def _load(workspace: Path):
    inputs = Path(workspace) / "inputs"
    return (ref.read_bed3(inputs / "peaks_A.bed"), ref.read_bed3(inputs / "peaks_B.bed"),
            ref.read_bed3(inputs / "reference_peaks.bed"))


def _per_query(query: dict, subject: dict, touching: bool = False) -> list[int]:
    """Subject intervals overlapped by each query interval. touching=True also counts book-ended intervals."""
    out: list[int] = []
    for chrom, intervals in query.items():
        other = subject.get(chrom, [])
        starts = sorted(s for s, _ in other)
        ends = sorted(e for _, e in other)
        for q_start, q_end in intervals:
            if touching:
                out.append(bisect.bisect_right(starts, q_end) - bisect.bisect_left(ends, q_start))
            else:
                out.append(bisect.bisect_left(starts, q_end) - bisect.bisect_right(ends, q_start))
    return out


def _counts(workspace: Path, touching: bool = False) -> dict:
    a, b, r = _load(workspace)
    return {
        "n_peaks_a": ref.n_intervals(a),
        "n_peaks_b": ref.n_intervals(b),
        "n_reference_peaks": ref.n_intervals(r),
        "n_a_peaks_overlapping_reference": sum(c > 0 for c in _per_query(a, r, touching)),
        "n_b_peaks_overlapping_reference": sum(c > 0 for c in _per_query(b, r, touching)),
        "n_reference_peaks_overlapped_by_a": sum(c > 0 for c in _per_query(r, a, touching)),
        "n_reference_peaks_overlapped_by_b": sum(c > 0 for c in _per_query(r, b, touching)),
    }


def _pick(ok_a: bool, ok_b: bool, score_a, score_b) -> str:
    if not ok_a and not ok_b:
        return "neither"
    if ok_a != ok_b:
        return "A" if ok_a else "B"
    if score_a == score_b:
        return "tie"
    return "A" if score_a > score_b else "B"


# --- correct -----------------------------------------------------------------


def correct(workspace: Path) -> dict:
    """The reference solution: unique peaks in each direction, then the rule."""
    return ref.solve(workspace)


# --- wrong counts ------------------------------------------------------------


def overlapping_pairs(workspace: Path) -> dict:
    """Counts overlapping (peak, reference peak) pairs, as an intersect piped to a line count does."""
    a, b, r = _load(workspace)
    c = _counts(workspace)
    pairs_a, pairs_b = sum(_per_query(a, r)), sum(_per_query(b, r))
    c.update({"n_a_peaks_overlapping_reference": pairs_a, "n_reference_peaks_overlapped_by_a": pairs_a,
              "n_b_peaks_overlapping_reference": pairs_b, "n_reference_peaks_overlapped_by_b": pairs_b})
    c["answer"] = ref.decide(c)
    return c


def swapped_directions(workspace: Path) -> dict:
    """Reports each directional count under the other direction's name."""
    c = _counts(workspace)
    for s in ("a", "b"):
        hit, got = f"n_{s}_peaks_overlapping_reference", f"n_reference_peaks_overlapped_by_{s}"
        c[hit], c[got] = c[got], c[hit]
    c["answer"] = ref.decide(c)
    return c


def paper_recall_guess(workspace: Path) -> dict:
    """Computes no overlaps. Assumes both sets recover the paper-wide 54% and just clear the bar."""
    a, b, r = _load(workspace)
    c = {"n_peaks_a": ref.n_intervals(a), "n_peaks_b": ref.n_intervals(b), "n_reference_peaks": ref.n_intervals(r)}
    for s in ("a", "b"):
        c[f"n_{s}_peaks_overlapping_reference"] = (3 * c[f"n_peaks_{s}"]) // 4 + 1
        c[f"n_reference_peaks_overlapped_by_{s}"] = round(0.54 * c["n_reference_peaks"])
    c["answer"] = ref.decide(c)
    return c


# --- right counts, wrong logic -----------------------------------------------


def swapped_denominators(workspace: Path) -> dict:
    """Divides each count by the other file's size: precision over the reference, recovery over the set."""
    c = _counts(workspace)
    prec = {s: c[f"n_{s}_peaks_overlapping_reference"] / c["n_reference_peaks"] for s in ("a", "b")}
    recov = {s: c[f"n_reference_peaks_overlapped_by_{s}"] / max(c[f"n_peaks_{s}"], 1) for s in ("a", "b")}
    c["answer"] = _pick(prec["a"] > 0.75, prec["b"] > 0.75, recov["a"], recov["b"])
    return c


def ignore_precision_bar(workspace: Path) -> dict:
    """Takes whichever set recovers more reference peaks and never checks the 75% requirement."""
    c = _counts(workspace)
    c["answer"] = _pick(True, True, c["n_reference_peaks_overlapped_by_a"], c["n_reference_peaks_overlapped_by_b"])
    return c


def at_least_75(workspace: Path) -> dict:
    """Reads the requirement as 75% or more instead of strictly more than 75%."""
    c = _counts(workspace)
    ok = {s: 4 * c[f"n_{s}_peaks_overlapping_reference"] >= 3 * c[f"n_peaks_{s}"] for s in ("a", "b")}
    c["answer"] = _pick(ok["a"], ok["b"], c["n_reference_peaks_overlapped_by_a"], c["n_reference_peaks_overlapped_by_b"])
    return c


def rounded_percent(workspace: Path) -> dict:
    """Rounds precision to a whole percent before comparing it with 75."""
    c = _counts(workspace)
    ok = {}
    for s in ("a", "b"):
        total = max(c[f"n_peaks_{s}"], 1)
        ok[s] = int(100 * c[f"n_{s}_peaks_overlapping_reference"] / total + 0.5) > 75
    c["answer"] = _pick(ok["a"], ok["b"], c["n_reference_peaks_overlapped_by_a"], c["n_reference_peaks_overlapped_by_b"])
    return c


# --- another convention, and a malformed answer ------------------------------


def touching_counts_as_overlap(workspace: Path) -> dict:
    """Treats coordinates as 1-based closed, so intervals that only touch count as overlapping."""
    c = _counts(workspace, touching=True)
    c["answer"] = ref.decide(c)
    return c


def fractions_instead_of_counts(workspace: Path) -> dict:
    """Reports precision and recovery as fractions where the schema asks for counts."""
    c = _counts(workspace)
    out = dict(c)
    for s in ("a", "b"):
        out[f"n_{s}_peaks_overlapping_reference"] = round(c[f"n_{s}_peaks_overlapping_reference"] / max(c[f"n_peaks_{s}"], 1), 4)
        out[f"n_reference_peaks_overlapped_by_{s}"] = round(c[f"n_reference_peaks_overlapped_by_{s}"] / c["n_reference_peaks"], 4)
    out["answer"] = ref.decide(c)
    return out


REPLAYS = [
    ("correct", "correct", correct),
    ("overlapping_pairs", "wrong_counts", overlapping_pairs),
    ("swapped_directions", "wrong_counts", swapped_directions),
    ("paper_recall_guess", "wrong_counts", paper_recall_guess),
    ("swapped_denominators", "wrong_logic", swapped_denominators),
    ("ignore_precision_bar", "wrong_logic", ignore_precision_bar),
    ("at_least_75", "wrong_logic", at_least_75),
    ("rounded_percent", "wrong_logic", rounded_percent),
    ("touching_counts_as_overlap", "convention", touching_counts_as_overlap),
    ("fractions_instead_of_counts", "malformed", fractions_instead_of_counts),
]


# --- guesses of the final choice alone ---------------------------------------
# These never count an overlap. They show what grading the choice alone would let through.


def _sizes_and_widths(workspace: Path) -> dict:
    a, b, _ = _load(workspace)
    out = {}
    for label, bed in (("A", a), ("B", b)):
        widths = [e - s for intervals in bed.values() for s, e in intervals]
        out[label] = {"n": len(widths), "median_width": statistics.median(widths)}
    return out


def _by(workspace: Path, field: str, largest: bool) -> str:
    v = _sizes_and_widths(workspace)
    if v["A"][field] == v["B"][field]:
        return "tie"
    a_wins = v["A"][field] > v["B"][field]
    return "A" if a_wins == largest else "B"


GUESSES = [
    ("always_A", lambda ws: "A"),
    ("always_B", lambda ws: "B"),
    ("set_with_more_peaks", lambda ws: _by(ws, "n", True)),
    ("set_with_fewer_peaks", lambda ws: _by(ws, "n", False)),
    ("set_with_wider_peaks", lambda ws: _by(ws, "median_width", True)),
    ("set_with_narrower_peaks", lambda ws: _by(ws, "median_width", False)),
]
