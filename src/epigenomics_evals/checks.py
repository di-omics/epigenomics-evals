"""Evidence behind the answer keys, and anchors back to the paper.

key_checks(task_id)  sentences with computed numbers that back each clause of
                     a multiple-choice key. They go into the task notes.
paper_anchors()      quantities the paper states, recomputed here from the
                     public files, so a reader can see the reference pipeline
                     lands where the paper did.
"""
from __future__ import annotations

import re

import numpy as np

from . import data
from .intervals import (
    STANDARD_CHROMS,
    clean_peaks,
    count_overlaps,
    merged_bp,
    overlaps_any,
    read_bed,
)


def genome_bp() -> int:
    total = 0
    for line in data.path("chrom_sizes").read_text().splitlines():
        chrom, size = line.split("\t")[:2]
        if chrom in STANDARD_CHROMS:
            total += int(size)
    return total


def _f(mask) -> float:
    return float(np.mean(mask))


def ct02() -> list[str]:
    bl, enc = data.load("blacklist"), data.load("encode_h3k27ac")
    ct = clean_peaks(data.load("ac_macs2"), bl)
    hit = overlaps_any(enc, ct)
    q = enc.extra[8].to_numpy(dtype=float)
    order = np.argsort(-q, kind="stable")
    k = len(q) // 10
    top, bottom = _f(hit[order[:k]]), _f(hit[order[-k:]])
    return [
        f"Overall, {_f(hit):.1%} of ENCODE peaks are recovered; {top:.1%} in the strongest tenth "
        f"and {bottom:.1%} in the weakest tenth, a {top / bottom:.0f}-fold difference.",
        f"So for a peak in the weakest tenth, no CUT&Tag peak is the outcome in about "
        f"{1 - bottom:.0%} of cases where ENCODE does report the mark.",
    ]


def ct03() -> list[str]:
    bl, enc = data.load("blacklist"), data.load("encode_h3k27ac")
    dd = clean_peaks(data.load("ac_macs2"), bl)
    wd = clean_peaks(data.load("ac_macs2_withdup"), bl)
    sd = clean_peaks(data.load("ac_seacr"), bl)
    sw = clean_peaks(data.load("ac_seacr_withdup"), bl)
    q = wd.extra[8].to_numpy(dtype=float)
    hit = overlaps_any(wd, enc)
    dup_only = ~overlaps_any(wd, dd)
    order = np.argsort(-q, kind="stable")
    running = np.cumsum(hit[order]) / np.arange(1, len(q) + 1)
    target = _f(overlaps_any(dd, enc))
    k = int(np.flatnonzero(running >= target).max()) + 1
    top = wd.take(np.isin(np.arange(len(wd)), order[:k]))
    at100 = q > 100
    return [
        f"MACS2 calls {len(dd):,} peaks with duplicates removed and {len(wd):,} with duplicates "
        f"kept ({len(wd) / len(dd):.1f}-fold). SEACR calls {len(sd):,} and {len(sw):,}.",
        f"ENCODE peaks recovered: MACS2 {_f(overlaps_any(enc, dd)):.1%} without duplicates and "
        f"{_f(overlaps_any(enc, wd)):.1%} with them. That gain is real, and it is why judging on "
        f"recovery alone points the wrong way.",
        f"The duplicate-only peaks are not weak calls: their median -log10 q in the duplicates-kept "
        f"file is {np.median(q[dup_only]):.0f}. At -log10 q above 100 the duplicates-kept set "
        f"still holds {int(at100.sum()):,} peaks, {int((dup_only & at100).sum()):,} of them "
        f"duplicate-only, with {_f(hit[at100]):.1%} in ENCODE.",
        f"To match the {target:.1%} of the duplicates-removed set, the duplicates-kept set has to "
        f"be cut to its top {k:,} peaks (-log10 q of {q[order[k - 1]]:.0f} or more), which recover "
        f"{_f(overlaps_any(enc, top)):.1%} of ENCODE against {_f(overlaps_any(enc, dd)):.1%} for "
        f"the duplicates-removed set. A stricter cutoff does not rescue it.",
    ]


def ct05() -> list[str]:
    bl = data.load("blacklist")
    macs2 = clean_peaks(data.load("ac_macs2"), bl)
    seacr = clean_peaks(data.load("ac_seacr"), bl)
    win = data.promoter_windows(data.refseq_tss())
    g = genome_bp()
    out = [f"{len(win):,} promoter windows cover {merged_bp(win) / 1e6:.1f} Mb, "
           f"{merged_bp(win) / g:.1%} of the genome."]
    for name, pk in (("MACS2", macs2), ("SEACR", seacr)):
        out.append(
            f"{name}: {len(pk):,} peaks, median width {np.median(pk.width):.0f} bp; "
            f"{_f(overlaps_any(pk, win)):.1%} of peaks touch a promoter; "
            f"{int(overlaps_any(win, pk).sum()):,} promoters ({_f(overlaps_any(win, pk)):.1%}) carry a peak."
        )
    return out


def ct06() -> list[str]:
    enc, me3, atac = data.load("encode_h3k27ac"), data.load("encode_h3k27me3"), data.load("atac")
    ac, me = data.load("ac_fragments"), data.load("me3_fragments")
    g = genome_bp()
    in_atac = overlaps_any(enc, atac)
    marked_accessible, marked_closed = enc.take(in_atac), enc.take(~in_atac)
    no_ac = ~overlaps_any(atac, enc)
    unmarked = atac.take(no_ac & ~overlaps_any(atac, me3))

    def share(frags, region) -> float:
        return _f(overlaps_any(frags, region))

    def fold(frags, region) -> float:
        return share(frags, region) / (merged_bp(region) / g)

    me_in_atac = overlaps_any(me, atac)
    me_short, me_long = me.take(me.width < 100), me.take(me.width >= 100)
    return [
        f"{_f(in_atac):.1%} of ENCODE H3K27ac peaks overlap an ATAC peak.",
        f"Marked and accessible sites: H3K27ac library {share(ac, marked_accessible):.1%} of "
        f"fragments ({fold(ac, marked_accessible):.1f}-fold over a uniform spread), H3K27me3 "
        f"library {share(me, marked_accessible):.1%} ({fold(me, marked_accessible):.1f}-fold). "
        f"Ratio {share(ac, marked_accessible) / share(me, marked_accessible):.1f}.",
        f"Marked, not accessible: H3K27ac library {fold(ac, marked_closed):.1f}-fold. The signal "
        f"is there without an ATAC peak.",
        f"Accessible, unmarked: H3K27ac library {fold(ac, unmarked):.1f}-fold, H3K27me3 library "
        f"{fold(me, unmarked):.1f}-fold.",
        f"Of the H3K27me3-library fragments in ATAC peaks, "
        f"{_f(overlaps_any(me.take(me_in_atac), me3)):.1%} also overlap an ENCODE H3K27me3 peak.",
        f"At accessible, unmarked sites the H3K27me3 library is {fold(me_short, unmarked):.1f}-fold "
        f"for fragments under 100 bp and {fold(me_long, unmarked):.1f}-fold for longer ones.",
    ]


KEY_CHECKS = {
    "ct02_strong_vs_weak_recovery": ct02,
    "ct03_duplicates_decision": ct03,
    "ct05_promoter_denominator": ct05,
    "ct06_open_chromatin_control": ct06,
}


def key_checks(task_id: str) -> list[str]:
    fn = KEY_CHECKS.get(task_id)
    return fn() if fn else []


# --- figures quoted in the task prose, computed here so none is typed by hand ---


def facts(task_id: str) -> dict:
    """Formatted figures for the intent and trap sentences of one task."""
    bl = data.load("blacklist")
    enc = data.load("encode_h3k27ac")
    if task_id == "ct01_encode_capture":
        raw = data.load("ac_macs2")
        kept = clean_peaks(raw, bl)
        removed = raw.take(overlaps_any(raw, bl) | ~np.isin(raw.chrom, STANDARD_CHROMS))
        return {"one_in": f"{len(raw) / len(removed):.0f}",
                "removed_supported": f"{_f(overlaps_any(removed, enc)):.0%}",
                "points": f"{100 * (_f(overlaps_any(kept, enc)) - _f(overlaps_any(raw, enc))):.1f}"}
    if task_id == "ct02_strong_vs_weak_recovery":
        hit = overlaps_any(enc, clean_peaks(data.load("ac_macs2"), bl))
        order = np.argsort(-enc.extra[8].to_numpy(dtype=float), kind="stable")
        k = len(order) // 10
        return {"fold": f"{_f(hit[order[:k]]) / _f(hit[order[-k:]]):.0f}"}
    if task_id == "ct03_duplicates_decision":
        dd = clean_peaks(data.load("ac_macs2"), bl)
        wd = clean_peaks(data.load("ac_macs2_withdup"), bl)
        added = ~overlaps_any(wd, dd)
        return {"dup_pct": f"{data.picard_duplication_rate():.1%}",
                "fold": f"{len(wd) / len(dd):.1f}",
                "added_supported": f"{_f(overlaps_any(wd.take(added), enc)):.1%}",
                "median_q": f"{np.median(wd.extra[8].to_numpy(dtype=float)[added]):.0f}"}
    if task_id == "ct04_antibody_panel_reference":
        sheet = data.panel_sheet()
        me3_ref = data.load("encode_h3k27me3")
        me3 = sheet[sheet["target"] == "H3K27me3"]
        wrong, own = [], []
        for row in me3.itertuples():
            pk = clean_peaks(data.load_panel_peaks(row.peak_file), bl)
            wrong.append(_f(overlaps_any(pk, enc)))
            own.append(_f(overlaps_any(pk, me3_ref)))
        return {"n_libraries": str(len(sheet)), "n_me3": str(len(me3)),
                "wrong_max": f"{max(wrong):.1%}", "me3_clearing": str(sum(v >= 0.75 for v in own))}
    if task_id == "ct05_promoter_denominator":
        macs2 = clean_peaks(data.load("ac_macs2"), bl)
        seacr = clean_peaks(data.load("ac_seacr"), bl)
        win = data.promoter_windows(data.refseq_tss())
        per_peak = _f(overlaps_any(seacr, win)) - _f(overlaps_any(macs2, win))
        per_promoter = abs(_f(overlaps_any(win, macs2)) - _f(overlaps_any(win, seacr)))
        return {"points": f"{100 * per_peak:.1f}", "pct_as_many": f"{len(seacr) / len(macs2):.0%}",
                "width_ratio": f"{np.median(seacr.width) / np.median(macs2.width):.1f}",
                "promoter_gap": f"{100 * per_promoter:.1f}"}
    if task_id == "ct06_open_chromatin_control":
        me3_ref, atac = data.load("encode_h3k27me3"), data.load("atac")
        me = data.load("me3_fragments")
        in_atac = overlaps_any(me, atac)
        marked_accessible = enc.take(overlaps_any(enc, atac))
        return {"any_pct": f"{_f(in_atac):.1%}",
                "in_domains": f"{_f(overlaps_any(me.take(in_atac), me3_ref)):.1%}",
                "control_marked": f"{_f(overlaps_any(me, marked_accessible)):.1%}"}
    return {}


def ct04() -> list[str]:
    g = genome_bp()
    ac, me3 = data.load("encode_h3k27ac"), data.load("encode_h3k27me3")
    return [
        f"The H3K27ac reference covers {merged_bp(ac) / 1e6:.0f} Mb ({merged_bp(ac) / g:.1%} of the genome) and the "
        f"H3K27me3 reference {merged_bp(me3) / 1e6:.0f} Mb ({merged_bp(me3) / g:.1%}). A share of peaks inside a "
        f"reference that broad is far easier to reach, so the 75% bar is not a like-for-like comparison across the "
        f"2 marks.",
        "The task therefore asks only which reference is the right one for each library. It does not rank "
        "antibodies for one mark against antibodies for the other, and mean_encode_recovery, which averages over "
        "libraries scored against 2 different references, is there to expose a wrong reference, not to measure quality.",
    ]


KEY_CHECKS["ct04_antibody_panel_reference"] = ct04


def option_claims() -> list[dict]:
    """Quantities stated in fixed answer-option text, each checked against the data."""
    bl = data.load("blacklist")
    enc, atac = data.load("encode_h3k27ac"), data.load("atac")
    macs2 = clean_peaks(data.load("ac_macs2"), bl)
    seacr = clean_peaks(data.load("ac_seacr"), bl)
    win = data.promoter_windows(data.refseq_tss())
    gap = abs(_f(overlaps_any(win, macs2)) - _f(overlaps_any(win, seacr)))
    width_ratio = float(np.median(seacr.width) / np.median(macs2.width))
    hit = overlaps_any(enc, macs2)
    order = np.argsort(-enc.extra[8].to_numpy(dtype=float), kind="stable")
    weakest = _f(hit[order[-(len(order) // 10):]])
    ac, me = data.load("ac_fragments"), data.load("me3_fragments")
    in_atac = overlaps_any(enc, atac)
    marked = enc.take(in_atac)
    unmarked = atac.take(~overlaps_any(atac, enc) & ~overlaps_any(atac, data.load("encode_h3k27me3")))
    ratio_control = _f(overlaps_any(ac, marked)) / _f(overlaps_any(me, marked))
    ratio_unmarked = _f(overlaps_any(ac, marked)) / _f(overlaps_any(ac, unmarked))
    return [
        {"task": "ct02_strong_vs_weak_recovery", "claim": "only a small minority of ENCODE peaks of that strength",
         "measured": f"{weakest:.1%} of the weakest tenth recovered", "holds": weakest < 0.25},
        {"task": "ct05_promoter_denominator", "claim": "a similar share of promoters (within 5 percentage points)",
         "measured": f"{100 * gap:.1f} points apart", "holds": gap <= 0.05},
        {"task": "ct05_promoter_denominator", "claim": "several-fold wider peaks",
         "measured": f"median width ratio {width_ratio:.1f}", "holds": width_ratio >= 3},
        {"task": "ct06_open_chromatin_control", "claim": "about twenty times the H3K27me3 library's share",
         "measured": f"ratio {ratio_control:.1f}", "holds": 15 <= ratio_control <= 25},
        {"task": "ct06_open_chromatin_control", "claim": "several times its own share at unmarked accessible sites",
         "measured": f"ratio {ratio_unmarked:.1f}", "holds": ratio_unmarked >= 3},
    ]


def paper_anchors() -> list[dict]:
    """Statements in the paper, recomputed from its public files."""
    bl, enc = data.load("blacklist"), data.load("encode_h3k27ac")
    peaks_dir = data.data_root() / "geo" / "peaks"
    recall, precision = {}, {}
    for f in sorted(peaks_dir.glob("*.bed.gz")):
        m = re.match(r"GSM\d+_(ab177|ab4729|diag)-(\d)_rmDup_(?:individual_seacr|q1e-5_macs2)", f.name)
        if not m:
            continue
        caller = "SEACR" if "seacr" in f.name else "MACS2"
        pk = clean_peaks(read_bed(f), bl)
        key = (f"{m.group(1)}-{m.group(2)}", caller)
        recall[key] = _f(overlaps_any(enc, pk))
        precision[key] = _f(overlaps_any(pk, enc))
    values = np.array(list(recall.values()))
    ab177 = [precision[k] for k in precision if k[0].startswith("ab177") and k[1] == "MACS2"]
    macs2 = clean_peaks(data.load("ac_macs2"), bl)
    seacr = clean_peaks(data.load("ac_seacr"), bl)
    per_seacr = count_overlaps(seacr, macs2)
    me = data.load("me3_fragments")
    atac = data.load("atac")
    return [
        {"statement": "Average ENCODE recall of 54% across CUT&Tag experiments (sd 4.99)",
         "where": "Abstract; Results, Fig. 5b",
         "recomputed": f"{values.mean():.1%} (sd {values.std(ddof=1) * 100:.2f}) over the 6 H3K27ac "
                       f"benchmarking libraries x 2 peak callers in GEO",
         "value": float(values.mean())},
        {"statement": "Antibody ab177178 did not reach 75% precision with MACS2",
         "where": "Results, Fig. 4a",
         "recomputed": "MACS2 precision " + " and ".join(f"{v:.1%}" for v in ab177)
                       + " for the 2 ab177178 libraries",
         "value": float(max(ab177))},
        {"statement": "Nearly 70% of ENCODE H3K27ac peaks overlapped with ATAC-seq peaks",
         "where": "Results, Fig. 3c",
         "recomputed": f"{_f(overlaps_any(enc, atac)):.1%}",
         "value": _f(overlaps_any(enc, atac))},
        {"statement": "An average of 1.35 to 1.68 MACS2 peaks overlapping a SEACR peak per sample",
         "where": "Results, peak calling",
         "recomputed": f"{per_seacr[per_seacr > 0].mean():.2f} for the library used in the tasks",
         "value": float(per_seacr[per_seacr > 0].mean())},
        {"statement": "Duplication rates between 55.49% and 98.45% across samples",
         "where": "Results; Supplementary Table 2",
         "recomputed": f"{data.picard_duplication_rate():.2%} for the library used in the tasks",
         "value": data.picard_duplication_rate()},
        {"statement": "About 8% of H3K27me3 CUT&Tag reads fell into ATAC peaks (8.29%), 5.88% after "
                      "removing fragments under 100 bp",
         "where": "Results, Fig. 3c",
         "recomputed": f"{_f(overlaps_any(me, atac)):.2%} and "
                       f"{_f(overlaps_any(me.take(me.width >= 100), atac)):.2%} on the first-round "
                       f"H3K27me3 library in the code repository (the paper's figures are for 2 "
                       f"later libraries)",
         "value": _f(overlaps_any(me, atac))},
    ]
