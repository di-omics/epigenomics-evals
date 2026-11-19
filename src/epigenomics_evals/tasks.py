"""Reference solutions and scripted controls for every task.

Each task function returns a list of Control objects:

    reference  the analysis the task is built around. Must pass.
    variant    a defensible alternative (another convention, another filter).
               Must pass too: the tolerance is set so that it does.
    trap       a plausible wrong analysis named in the task notes. Must fail.

That is the acceptance rule for a task in this repository. A grader window is
only as good as the answers it has been shown to separate, so every window
here is tested against both sides before the task ships.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import data
from .intervals import (
    Intervals,
    clean_peaks,
    count_overlaps,
    drop_blacklisted,
    intersect_bp,
    keep_standard,
    merged_bp,
    overlaps_any,
    subtract,
)


@dataclass
class Control:
    name: str
    kind: str  # reference | variant | trap
    description: str
    answer: dict
    task_id: str = ""
    expect_pass: bool = field(init=False)

    def __post_init__(self) -> None:
        self.expect_pass = self.kind != "trap"


def _frac(mask: np.ndarray) -> float:
    return float(np.mean(mask)) if mask.size else 0.0


def _col(iv: Intervals, bed_column: int) -> np.ndarray:
    """A numeric column by its 1-based BED column number."""
    return iv.extra[bed_column - 1].to_numpy(dtype=float)


# --- CT01: ENCODE capture in both directions ---------------------------------


def ct01() -> list[Control]:
    bl = data.load("blacklist")
    enc = data.load("encode_h3k27ac")
    raw = data.load("ac_macs2")

    def answer(ct: Intervals, ref: Intervals = enc, closed: bool = False) -> dict:
        return {
            "n_library_peaks": len(ct),
            "frac_library_in_encode": _frac(overlaps_any(ct, ref, closed)),
            "frac_encode_captured": _frac(overlaps_any(ref, ct, closed)),
        }

    clean = clean_peaks(raw, bl)
    ref = answer(clean)
    swapped = dict(ref)
    swapped["frac_library_in_encode"], swapped["frac_encode_captured"] = (
        ref["frac_encode_captured"],
        ref["frac_library_in_encode"],
    )
    shared_bp = intersect_bp(clean, enc)
    return [
        Control("reference", "reference",
                "Drop peaks in blacklisted regions and on non-standard chromosomes, then count "
                "peak overlaps in each direction.", ref),
        Control("blacklist_only", "variant",
                "Drop blacklisted peaks but keep unplaced and mitochondrial contigs.",
                answer(drop_blacklisted(raw, bl))),
        Control("book_ended_counts", "variant",
                "Treat coordinates as 1-based closed, so touching intervals count as overlapping.",
                answer(clean, closed=True)),
        Control("encode_filtered_too", "variant",
                "Apply the same filter to the ENCODE set before computing recovery.",
                answer(clean, ref=clean_peaks(enc, bl))),
        Control("no_artifact_filter", "trap",
                "Compare the peak file as delivered, blacklisted regions included.",
                answer(raw)),
        Control("swapped_denominators", "trap",
                "Report the two fractions against the wrong denominators.", swapped),
        Control("basepair_overlap", "trap",
                "Report base-pair overlap instead of the share of peaks.",
                {"n_library_peaks": len(clean),
                 "frac_library_in_encode": shared_bp / merged_bp(clean),
                 "frac_encode_captured": shared_bp / merged_bp(enc)}),
    ]


# --- CT02: which ENCODE peaks are recovered ----------------------------------

CT02_CORRECT = 4


def ct02() -> list[Control]:
    bl = data.load("blacklist")
    enc = data.load("encode_h3k27ac")
    raw = data.load("ac_macs2")
    clean = clean_peaks(raw, bl)
    q = _col(enc, 9)

    def answer(ct: Intervals, rank_by: np.ndarray, k: int | None = None, choice: int = CT02_CORRECT) -> dict:
        captured = overlaps_any(enc, ct)
        order = np.argsort(-rank_by, kind="stable")
        k = k or len(order) // 10
        return {
            "recall_strongest_decile": _frac(captured[order[:k]]),
            "recall_weakest_decile": _frac(captured[order[-k:]]),
            "median_neglog10q_recovered": float(np.median(q[captured])),
            "median_neglog10q_missed": float(np.median(q[~captured])),
            "interpretation": choice,
        }

    ref = answer(clean, q)
    overall = _frac(overlaps_any(enc, clean))
    # Trap: stratify the CUT&Tag peaks by their own q-value instead.
    ct_q = _col(clean, 9)
    in_enc = overlaps_any(clean, enc)
    order = np.argsort(-ct_q, kind="stable")
    k = len(order) // 10
    wrong_axis = {
        "recall_strongest_decile": _frac(in_enc[order[:k]]),
        "recall_weakest_decile": _frac(in_enc[order[-k:]]),
        "median_neglog10q_recovered": float(np.median(ct_q[in_enc])),
        "median_neglog10q_missed": float(np.median(ct_q[~in_enc])),
        "interpretation": CT02_CORRECT,
    }
    absent = dict(ref)
    absent["interpretation"] = 1
    return [
        Control("reference", "reference",
                "Rank ENCODE peaks by ENCODE q-value and compute recovery in the top and bottom tenth.",
                ref),
        Control("rank_by_signal_value", "variant",
                "Rank ENCODE peaks by signalValue (column 7) instead of q-value.",
                answer(clean, _col(enc, 7))),
        Control("rank_by_p_value", "variant",
                "Rank ENCODE peaks by p-value (column 8).", answer(clean, _col(enc, 8))),
        Control("unfiltered_peaks", "variant",
                "Use the CUT&Tag peaks as delivered.", answer(raw, q)),
        Control("decile_rounded_up", "variant",
                "Take ceil(n/10) peaks per decile.", answer(clean, q, k=-(-len(q) // 10))),
        Control("uniform_recovery_assumed", "trap",
                "Report the overall recovery for both strata, with no stratification.",
                {"recall_strongest_decile": overall, "recall_weakest_decile": overall,
                 "median_neglog10q_recovered": float(np.median(q)),
                 "median_neglog10q_missed": float(np.median(q)), "interpretation": CT02_CORRECT}),
        Control("stratified_wrong_peak_set", "trap",
                "Stratify the CUT&Tag peaks by their own q-value and report their overlap with ENCODE.",
                wrong_axis),
        Control("absence_read_as_absent", "trap",
                "Correct numbers, but concludes the mark is absent where no peak was called.", absent),
    ]


# --- CT03: duplicates --------------------------------------------------------

CT03_CORRECT = 6


def ct03() -> list[Control]:
    bl = data.load("blacklist")
    enc = data.load("encode_h3k27ac")
    names = ["ac_macs2", "ac_macs2_withdup", "ac_seacr", "ac_seacr_withdup"]
    keys = ["macs2_dedup", "macs2_withdup", "seacr_dedup", "seacr_withdup"]
    dup_rate = data.picard_duplication_rate()

    def answer(prep, closed: bool = False, choice: int = CT03_CORRECT, metric: str = "precision") -> dict:
        sets = {k: prep(data.load(n)) for k, n in zip(keys, names)}
        out = {"duplication_rate": dup_rate}
        for k, pk in sets.items():
            if metric == "precision":
                out[f"frac_in_encode_{k}"] = _frac(overlaps_any(pk, enc, closed))
            else:
                out[f"frac_in_encode_{k}"] = _frac(overlaps_any(enc, pk, closed))
        dup_only = subtract(sets["macs2_withdup"], sets["macs2_dedup"])
        out["n_macs2_duplicate_only_peaks"] = len(dup_only)
        out["frac_duplicate_only_in_encode"] = _frac(overlaps_any(dup_only, enc, closed))
        out["recommendation"] = choice
        return out

    def clean(iv: Intervals) -> Intervals:
        return clean_peaks(iv, bl)

    ref = answer(clean)
    # Trap: duplicate-only defined by exact coordinates rather than by overlap.
    wd, dd = clean(data.load("ac_macs2_withdup")), clean(data.load("ac_macs2"))
    dedup_keys = set(zip(dd.chrom.tolist(), dd.start.tolist(), dd.end.tolist()))
    exact = np.array([(c, s, e) not in dedup_keys
                      for c, s, e in zip(wd.chrom.tolist(), wd.start.tolist(), wd.end.tolist())])
    exact_match = dict(ref)
    exact_match["n_macs2_duplicate_only_peaks"] = int(exact.sum())
    exact_match["frac_duplicate_only_in_encode"] = _frac(overlaps_any(wd.take(exact), enc))
    keep = dict(ref)
    keep["recommendation"] = 1
    stringency = dict(ref)
    stringency["recommendation"] = 3
    return [
        Control("reference", "reference",
                "Filter each peak set, compute the share of peaks in ENCODE, isolate MACS2 peaks "
                "that exist only when duplicates are kept.", ref),
        Control("unfiltered_peaks", "variant",
                "Use the four peak files as delivered.", answer(lambda iv: iv)),
        Control("blacklist_only", "variant",
                "Drop blacklisted peaks only.", answer(lambda iv: drop_blacklisted(iv, bl))),
        Control("book_ended_counts", "variant",
                "1-based closed coordinates.", answer(clean, closed=True)),
        Control("keep_duplicates_by_default", "trap",
                "Correct numbers, but follows the usual CUT&Tag advice to keep duplicates.", keep),
        Control("raise_stringency_instead", "trap",
                "Correct numbers, but proposes filtering the extra peaks by significance.", stringency),
        Control("judged_on_recovery", "trap",
                "Reports the share of ENCODE peaks recovered by each set instead of the share of "
                "each set in ENCODE, and keeps duplicates because recovery rises.",
                answer(clean, choice=1, metric="recall")),
        Control("duplicate_only_by_exact_match", "trap",
                "Counts a peak as duplicate-only unless its exact coordinates exist in the "
                "duplicate-removed set.", exact_match),
    ]


# --- CT04: antibody panel, each library against its own mark -----------------


def ct04() -> list[Control]:
    bl = data.load("blacklist")
    refs = {"H3K27ac": data.load("encode_h3k27ac"), "H3K27me3": data.load("encode_h3k27me3")}
    sheet = data.panel_sheet()
    pooled = Intervals(
        np.r_[refs["H3K27ac"].chrom, refs["H3K27me3"].chrom],
        np.r_[refs["H3K27ac"].start, refs["H3K27me3"].start],
        np.r_[refs["H3K27ac"].end, refs["H3K27me3"].end],
    )

    def answer(prep, reference_for, closed: bool = False) -> dict:
        precision, recall = [], []
        for row in sheet.itertuples():
            pk = prep(data.load_panel_peaks(row.peak_file))
            ref = reference_for(row.target)
            precision.append(_frac(overlaps_any(pk, ref, closed)))
            recall.append(_frac(overlaps_any(ref, pk, closed)))
        precision, recall = np.array(precision), np.array(recall)
        return {
            "n_libraries_meeting_bar": int((precision >= 0.75).sum()),
            "lowest_library_precision": float(precision.min()),
            "mean_encode_recovery": float(recall.mean()),
        }

    def clean(iv: Intervals) -> Intervals:
        return clean_peaks(iv, bl)

    def own(target: str) -> Intervals:
        return refs[target]

    return [
        Control("reference", "reference",
                "Score every library against the ENCODE set for the mark its antibody targets.",
                answer(clean, own)),
        Control("unfiltered_peaks", "variant", "Use the peak files as delivered.",
                answer(lambda iv: iv, own)),
        Control("blacklist_only", "variant", "Drop blacklisted peaks only.",
                answer(lambda iv: drop_blacklisted(iv, bl), own)),
        Control("book_ended_counts", "variant", "1-based closed coordinates.",
                answer(clean, own, closed=True)),
        Control("one_reference_for_all", "trap",
                "Score all 20 libraries against the H3K27ac set, the mark most of the panel targets.",
                answer(clean, lambda t: refs["H3K27ac"])),
        Control("pooled_reference", "trap",
                "Score all libraries against both ENCODE sets pooled into one reference.",
                answer(clean, lambda t: pooled)),
    ]


# --- CT05: promoters, per peak or per promoter -------------------------------

CT05_CORRECT = 3


def ct05() -> list[Control]:
    bl = data.load("blacklist")
    raw = {"macs2": data.load("ac_macs2"), "seacr": data.load("ac_seacr")}

    def answer(prep, tss=None, half_width: int = 1500, inclusive: bool = True,
               unit: str = "tss", metric: str = "promoters", choice: int = CT05_CORRECT) -> dict:
        tss = data.refseq_tss() if tss is None else tss
        win = data.promoter_windows(tss, half_width, inclusive)
        out = {}
        for k, iv in raw.items():
            pk = prep(iv)
            if metric == "peaks":
                out[f"frac_promoters_with_peak_{k}"] = _frac(overlaps_any(pk, win))
            elif unit == "gene":
                hit = overlaps_any(win, pk)
                genes = tss["gene"].to_numpy(dtype=object)
                out[f"frac_promoters_with_peak_{k}"] = len(set(genes[hit])) / len(set(genes))
            else:
                out[f"frac_promoters_with_peak_{k}"] = _frac(overlaps_any(win, pk))
            out[f"median_peak_width_{k}"] = float(np.median(pk.width))
        out["conclusion"] = choice
        return out

    def clean(iv: Intervals) -> Intervals:
        return clean_peaks(iv, bl)

    tss_no_strand = data.refseq_tss().drop_duplicates(["chrom", "pos"]).reset_index(drop=True)
    return [
        Control("reference", "reference",
                "Build one window per unique RefSeq transcript start and ask what share of "
                "windows carries a peak from each caller.", answer(clean)),
        Control("gene_level", "variant",
                "Count a gene as recovered if any of its transcript starts carries a peak.",
                answer(clean, unit="gene")),
        Control("window_3000bp_half_open", "variant",
                "Windows of exactly 3,000 bp.", answer(clean, inclusive=False)),
        Control("tss_unique_by_position", "variant",
                "Collapse start sites that share a position on opposite strands.",
                answer(clean, tss=tss_no_strand)),
        Control("unfiltered_peaks", "variant", "Use the peak files as delivered.",
                answer(lambda iv: iv)),
        Control("share_of_peaks_at_promoters", "trap",
                "Report the share of each caller's peaks that touch a promoter and conclude the "
                "caller with the higher share recovers more promoters.",
                answer(clean, metric="peaks", choice=1)),
        Control("strand_ignored", "trap",
                "Take the GTF start coordinate as the start site for every transcript, minus "
                "strand included.", answer(clean, tss=data.refseq_tss(strand_aware=False))),
    ]


# --- CT06: accessibility or the mark, with a second library as the control ---

CT06_CORRECT = 5


def ct06() -> list[Control]:
    enc = data.load("encode_h3k27ac")
    me3_domains = data.load("encode_h3k27me3")
    atac = data.load("atac")
    frags = {"h3k27ac": data.load("ac_fragments"), "h3k27me3": data.load("me3_fragments")}

    def answer(prep=lambda f: f, closed: bool = False, split_domains: bool = True,
               partition: bool = True, choice: int = CT06_CORRECT) -> dict:
        enc_in_atac = overlaps_any(enc, atac, closed)
        marked_accessible = enc.take(enc_in_atac)
        marked_closed = enc.take(~enc_in_atac)
        unmarked = ~overlaps_any(atac, enc, closed)
        if split_domains:
            unmarked &= ~overlaps_any(atac, me3_domains, closed)
        accessible_unmarked = atac.take(unmarked) if partition else atac
        if not partition:
            marked_accessible = atac
        ac, me = prep(frags["h3k27ac"]), prep(frags["h3k27me3"])
        me_long = me.take(me.width >= 100)
        return {
            "h3k27ac_lib_marked_accessible": _frac(overlaps_any(ac, marked_accessible, closed)),
            "h3k27ac_lib_marked_not_accessible": _frac(overlaps_any(ac, marked_closed, closed)),
            "h3k27ac_lib_accessible_unmarked": _frac(overlaps_any(ac, accessible_unmarked, closed)),
            "h3k27me3_lib_marked_accessible": _frac(overlaps_any(me, marked_accessible, closed)),
            "h3k27me3_lib_accessible_unmarked": _frac(overlaps_any(me, accessible_unmarked, closed)),
            "h3k27me3_lib_any_atac": _frac(overlaps_any(me, atac, closed)),
            "h3k27me3_lib_any_atac_min100bp": _frac(overlaps_any(me_long, atac, closed)),
            "interpretation": choice,
        }

    ref = answer()
    just_atac = dict(ref)
    just_atac["interpretation"] = 2
    inflated = dict(ref)
    inflated["interpretation"] = 3
    return [
        Control("reference", "reference",
                "Split accessible sites by which ENCODE mark they carry, and read the H3K27me3 "
                "library at sites without H3K27me3 as the antibody-independent background.", ref),
        Control("standard_chromosomes_only", "variant",
                "Count fragments on chr1-22, X and Y only.", answer(prep=keep_standard)),
        Control("book_ended_counts", "variant", "1-based closed coordinates.", answer(closed=True)),
        Control("h3k27me3_domains_ignored", "trap",
                "Define unmarked accessible sites against H3K27ac only, so accessible sites inside "
                "H3K27me3 domains count as background for the H3K27me3 library.",
                answer(split_domains=False, choice=3)),
        Control("accessible_sites_not_partitioned", "trap",
                "Use all ATAC peaks for every accessible class and conclude the library behaves "
                "like ATAC-seq.", answer(partition=False, choice=2)),
        Control("called_atac_seq", "trap",
                "Correct numbers, but concludes the profile is explained by accessibility.",
                just_atac),
        Control("background_overestimated", "trap",
                "Correct numbers, but takes the H3K27me3 library's share in ATAC peaks as the background rate.",
                inflated),
    ]


TASKS = {
    "ct01_encode_capture": ct01,
    "ct02_strong_vs_weak_recovery": ct02,
    "ct03_duplicates_decision": ct03,
    "ct04_antibody_panel_reference": ct04,
    "ct05_promoter_denominator": ct05,
    "ct06_open_chromatin_control": ct06,
}


def controls(task_id: str) -> list[Control]:
    out = TASKS[task_id]()
    for c in out:
        c.task_id = task_id
    return out
