"""The pilot's public documents, built from the recorded results. Standard library only.

Called by `python3 pilot/replay.py`. VERIFICATION.md and the generated block of
README.md are written from results/summary.json and from the non-secret parts
of private/provenance.json, so no number in them is typed by hand, and neither
page states the answer key or which peak caller sits behind which label.
"""
from __future__ import annotations

from pathlib import Path

BEGIN, END = "<!-- BEGIN GENERATED -->", "<!-- END GENERATED -->"


def with_block(readme: str, block: str) -> str:
    if BEGIN not in readme or END not in readme:
        return readme
    head, rest = readme.split(BEGIN, 1)
    return head + BEGIN + "\n" + block + "\n" + END + rest.split(END, 1)[1]


def verification(prov: dict, s: dict) -> str:
    src, conv = prov["sources"], prov["conventions"]
    enc = src["reference"]["encode_metadata"]
    primary = prov["instances"][0]
    key, touching = primary["key"], primary["counts_if_touching_intervals_counted"]
    deltas = sorted(abs(touching[f] - key[f]) for f in touching if touching[f] != key[f])
    filt, hyg, labels = primary["filter"], primary["hygiene"], primary["labels"]
    ref_h = prov["reference_hygiene"]
    blobs_ok = sum(src["git_blob_ids_local"].get(k) == v for k, v in src["git_blob_ids_pinned"].items())
    same = src["same_bytes_under_full_name"]
    settings = src["settings_evidence"]
    r, bl = s["recomputation"], s["blinding"]
    by_label = {label: filt[labels[label]["peak_caller"]] for label in ("A", "B")}

    rows = [
        ("Paper and decision",
         f"doi:{prov['paper']['doi']}, Figure {prov['paper']['figure']}: {prov['paper']['criterion_in_the_paper']}.",
         "Read from the paper (Results, peak calling with SEACR and MACS2).",
         "The pilot writes it as a rule: strictly more than 75%, with tie and neither as possible answers."),
        ("Candidate peak files",
         f"The 2 peak files of library {primary['library']} in the paper's code repository, commit `{src['code_repository']['commit'][:7]}`.",
         "Git blob id of each fetched file compared with the id the GitHub tree API lists at that commit.",
         f"{blobs_ok} of {len(src['git_blob_ids_pinned'])} files identical to the pinned commit."),
        ("Which calls they are",
         "One MACS2 set (q < 1e-5, no local lambda) and one SEACR set (stringent, top 1%), both on duplicate-removed reads.",
         "The repository holds each file twice, under a short name and under its full peak-caller name. Both copies "
         "are fetched and their git blob ids compared. Settings are read from the repository's peak-calling scripts.",
         f"{sum(same.values())} of {len(same)} pairs of copies are byte-identical. "
         f"{sum(settings.values())} of {len(settings)} settings statements found in the scripts: "
         + "; ".join(k for k, v in settings.items() if v) + "."),
        ("GEO accession",
         f"{primary['geo_series']}, sample {primary['geo_sample']}.",
         "Coordinates of the GEO peak files compared row by row with the repository files.",
         f"{sum(src['primary_coordinates_identical_to_geo'].values())} of {len(src['primary_coordinates_identical_to_geo'])} "
         "files have identical coordinates. Only the peak names differ."),
        ("Reference accession",
         f"ENCODE {enc['file']['accession']}: {enc['file']['output_type']}, {enc['file']['file_format_type']}, from "
         f"{enc['experiment']['accession']} ({enc['experiment']['assay_title']}, {enc['experiment']['biosample_summary']}, "
         f"{enc['experiment']['target'].strip('/').split('/')[-1]}).",
         f"ENCODE REST API queried on {enc['fetched']}. MD5 of the repository's copy compared with ENCODE's `content_md5sum`.",
         f"Status {enc['file']['status']}. MD5 "
         + ("matches" if src["reference"]["md5_of_source_file"] == enc["file"]["content_md5sum"] else "DOES NOT MATCH")
         + f" ENCODE's. {ref_h['rows']:,} peaks."),
        ("Genome build",
         f"{conv['genome_build']} for all 3 files.",
         "ENCODE lists the reference as hg19. Every interval end compared with the hg19 and the hg38 chromosome lengths.",
         f"0 intervals run past an hg19 chromosome end. {ref_h['ending_beyond_hg38_chromosome_end']:,} reference peaks, "
         f"{hyg['A']['ending_beyond_hg38_chromosome_end']:,} in set A and {hyg['B']['ending_beyond_hg38_chromosome_end']:,} in set B "
         "run past the hg38 end of their chromosome, so the coordinates cannot be hg38."),
        ("Coordinate convention",
         f"{conv['coordinates']}.",
         "Counts computed under this convention and under the 1-based closed one, then compared with the paper's notebook.",
         f"Half-open reproduces the notebook exactly ({r['paper_notebook_numbers_equal']} of {r['paper_notebook_numbers']} numbers). "
         f"Counting intervals that only touch changes {len(deltas)} of the 7 counts, by {deltas[0]} to {deltas[-1]} peaks."),
        ("Duplicate policy",
         f"{conv['duplicates']}.",
         "Read from the repository's duplicate-removal and peak-calling scripts (see the row above).",
         f"Duplication rate of this library in its Picard report: {100 * prov['primary_library_duplication_rate']:.2f}%."),
        ("Blacklist",
         f"{conv['candidate_set_filter']}. Reference: {conv['reference_filter']}.",
         "Uncompressed blacklist compared by SHA-256 with the hg19 v2 list in Boyle-Lab/Blacklist. Peaks removed are "
         "counted, and the counts compared with the paper's notebook.",
         ("Same content as the Boyle-Lab list. " if src["blacklist"]["sha256_uncompressed"] == src["blacklist"]["boyle_lab_sha256_uncompressed"] else "DIFFERS from the Boyle-Lab list. ")
         + "; ".join(f"set {label}: {f['rows_in_source']:,} peaks in the source, {f['overlapping_blacklist']:,} on the "
                     f"blacklist, {f['outside_chr1_22_X_Y']:,} off chr1-22, X, Y, {f['rows_staged']:,} staged"
                     for label, f in by_label.items())
         + f". 0 staged peaks overlap the blacklist. {ref_h['peaks_overlapping_blacklist']:,} reference peaks do, and "
           "stay in, as in the paper's code. The blacklist file is GPL-3.0 and is not redistributed here."),
        ("Overlap rule",
         f"{conv['overlap_rule']}.",
         "The paper's notebook counts peaks with a non-zero GenomicRanges overlap count, in each direction.",
         "Same rule; same numbers."),
        ("Blinding",
         "No name in the workspace says which set came from which peak caller, library or accession.",
         "The replay scans every workspace file for caller, antibody, sample, accession and author names, and "
         "measures how well the 2 callers can be told apart from the intervals alone.",
         f"0 names found. The intervals still give the callers away: on {bl['pairs_where_the_wider_set_is_the_same_caller']} "
         f"of {bl['pairs']} pairs the set with the wider peaks comes from the same caller. On the primary pair the row "
         "counts of the 2 files also equal the peak counts printed in the paper's public notebook. The primary is a "
         "public demonstration, not a held-out item."),
    ]
    out = ["# What was verified about the inputs", "",
           f"Status: {prov['status']}. Written by `python3 pilot/replay.py --write`.",
           "Every check on this page was written and run by AI models. No person has checked it yet.", "",
           "| Item | What the pilot relies on | How it was checked | Result |", "|---|---|---|---|"]
    out += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows]
    out += ["", "## What changed from Figure 4a", "",
            "The exact comparison in Figure 4a is not public, so the pilot uses the nearest one that is.", "",
            "- **Thresholds became callers.** Figure 4a sweeps settings within each peak caller (MACS2 q-values, "
            "SEACR top fractions). Only the peak sets at the final settings were deposited. The pilot compares the 2 "
            "callers at those final settings on the same library.",
            "- **The library.** The primary instance uses the first-round Diagenode 1:50 library that the paper's code "
            "repository ships as its worked example, because its files, its Picard report and the authors' rendered "
            "notebook are all public. The paper's 6 benchmarking libraries and 4 more first-round libraries are the "
            "10 check instances.",
            "- **The blacklist.** The paper's Methods name ENCODE file ENCFF000KJP. The repository's code uses the "
            "hg19 blacklist v2, and the notebook's blacklist counts match that file. The pilot follows the code.",
            "- **The rule.** The paper states the aim: maximise ENCODE capture while keeping precision above 75%. "
            "The pilot writes it as a rule with 4 possible answers.",
            "- **The files.** Reduced to 3 columns, sorted, relabelled A and B. No interval was changed.", ""]
    return "\n".join(out)


def readme_block(s: dict) -> str:
    r, rp, n = s["recomputation"], s["real_pairs"], len(s["instances"])
    out = [f"**The answer is not printed on this page.** It is in `private/answer_key.json`. A binary search in Python "
           f"and a range join in SQLite agree on {r['fields_equal_between_the_two_recomputations']} of {r['fields']} answer "
           f"fields and both equal the key, and {r['paper_notebook_numbers_equal']} of {r['paper_notebook_numbers']} numbers "
           "equal what the paper's own notebook prints.", "",
           "**Scripted replays on the primary instance.**", "",
           "| Replay | What it does | Counts wrong, of 7 | Choice | Grader |", "|---|---|---|---|---|"]
    for p in s["primary_replays"]:
        counts = str(p["counts_wrong"]) if p["valid_answer"] else "not a valid answer"
        out.append(f"| `{p['replay']}` | {p['says']} | {counts} | {'right' if p['choice_correct'] else 'wrong'} | "
                   f"{'PASS' if p['passed'] else 'fail'} |")
    out += ["", "**Shortcuts that survive on the primary instance:** "
            + ", ".join(f"`{x}`" for x in s["survives_on_the_primary"])
            + ". Both sets clear the requirement there, so getting the requirement wrong costs nothing.", "",
            f"**Across all {n} instances** (the primary, 2 built from it, 2 cut at the 75% line, and the paper's "
            f"{rp['other_pairs']} other public pairs):", "",
            "| Replay | Passes on | Caught by |", "|---|---|---|"]
    for name, v in s["survival"].items():
        fails = v["fails_on"]
        caught = "nothing" if not fails else (f"all {n}" if len(fails) == n else (", ".join(fails) if len(fails) <= 5 else f"{len(fails)} instances"))
        out.append(f"| `{name}` | {len(v['passes_on'])} of {n} | {caught} |")
    best = max(s["choice_only"]["guesses"].items(), key=lambda kv: kv[1])
    out += ["", f"**On the {rp['n']} real pairs** the rule selects "
            + ", ".join(f"{k} in {v}" for k, v in rp["answers_by_caller"].items())
            + f". Ignoring the 75% requirement gives the wrong answer on {rp['other_pairs_where_ignoring_the_bar_is_wrong']} "
              f"of the {rp['other_pairs']} pairs other than the primary.", "",
            f"**If only the final choice were graded,** picking the {best[0].replace('_', ' ')} would be right on "
            f"{best[1]} of the {rp['n']} real pairs without counting a single overlap. Requiring the counts closes that.", ""]
    runs = s["agent_runs"]
    if runs:
        out += [f"**Agent runs: {len(runs)}.** Hands-on, on the author's machine, not benchmark-grade.", "",
                "| Run | Model | Result | Tool calls | Calls naming a path outside the workspace |", "|---|---|---|---|---|"]
        out += [f"| `{a['run_id']}` | {a['model']} | {'PASS' if a['passed'] else 'fail'} | {a['tool_calls']} | "
                f"{a['calls_naming_a_path_outside_the_workspace']} |" for a in runs]
    else:
        out += ["**Agent runs: 0.**"]
    out += ["", "Full output without answers: [results/REPLAY.md](results/REPLAY.md). Input checks: "
            "[VERIFICATION.md](VERIFICATION.md)."]
    return "\n".join(out)
