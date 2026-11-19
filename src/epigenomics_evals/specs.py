"""Task specifications: prompt, answer fields, tolerances, notes.

The numbers in a published eval.json are never typed by hand. build.py
computes the reference answer from the data, rounds it, and writes it next to
the tolerances declared here.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .tasks import CT02_CORRECT, CT03_CORRECT, CT05_CORRECT, CT06_CORRECT

PAPER = (
    "Abbasova L, Urbanaviciute P, Hu D, et al. CUT&Tag recovers up to half of ENCODE "
    "ChIP-seq histone acetylation peaks. Nat Commun 16, 2993 (2025). "
    "doi:10.1038/s41467-025-58137-2"
)
TAXONOMY_VERSION = "v1.1.0"

AC_MACS2 = "Diagenode_C15410196_1_50_MACS2.bed"
AC_MACS2_DUP = "Diagenode_C15410196_1_50_MACS2_withDup.bed"
AC_SEACR = "Diagenode_C15410196_1_50_SEACR.bed"
AC_SEACR_DUP = "Diagenode_C15410196_1_50_SEACR_withDup.bed"
AC_PICARD = "Diagenode_C15410196_1_50_picard.dupMark.txt"
AC_FRAGS = "Diagenode_C15410196_1_50_2M_rmDup_bowtie2.fragments.bed"
ME3_FRAGS = "CST9733_1_100_H3K27me3_2M_rmDup_bowtie2.fragments.bed"
ENC_AC = "H3K27ac_ENCFF044JNJ.bed.narrowPeak"
ENC_ME3 = "H3K27me3_ENCFF000BXB.bed.broadPeak"
ATAC = "K562_ATACseq.mRp.clN_peaks.narrowPeak"
GTF = "hg19.ncbiRefSeq.gtf.gz"
BLACKLIST = "hg19-blacklist.v2.bed.gz"
CHROM_SIZES = "hg19.chrom.sizes"

OVERLAP_RULE = "Two intervals overlap if they share at least 1 bp."
COMMON_FILES = (
    f"The hg19 blacklist (`{BLACKLIST}`) and chromosome sizes (`{CHROM_SIZES}`) are also in "
    "the workspace."
)


@dataclass
class Field:
    name: str
    kind: str  # int | float | choice
    tolerance: dict | None = None  # {"value": x} or {"lower": a, "upper": b}; None for choice
    digits: int = 4
    meaning: str = ""


@dataclass
class Spec:
    id: str
    title: str
    paper_anchor: str
    category: str
    failure_modes: list[str]
    data: list[str]  # logical names from data.FILES, or "panel"
    prompt: str
    fields: list[Field]
    choice_question: str = ""
    choices: list[str] = field(default_factory=list)
    correct_choice: int | None = None
    intent: str = ""
    trap: str = ""
    tolerance_rationale: str = ""
    not_said: str = ""
    provenance: str = ""
    key_review: str = "Numeric key, computed from the data. No decision field."


def _answer_block(fields: list[Field]) -> str:
    kinds = {"int": "<integer>", "float": "<float>", "choice": "<integer>"}
    body = ", ".join(f'"{f.name}": {kinds[f.kind]}' for f in fields)
    return "Return EXACTLY:\n\n<EVAL_ANSWER>\n{" + body + "}\n</EVAL_ANSWER>"


def full_prompt(spec: Spec) -> str:
    parts = [spec.prompt.strip(), _answer_block(spec.fields)]
    defs = [f"`{f.name}`: {f.meaning}" for f in spec.fields if f.meaning]
    if defs:
        parts.append("\n".join(defs))
    if spec.choices:
        choice_field = next(f.name for f in spec.fields if f.kind == "choice")
        lines = [f"`{choice_field}`: {spec.choice_question}"]
        lines += [f"{i}. {text}" for i, text in enumerate(spec.choices, 1)]
        parts.append("\n".join(lines))
    return "\n\n".join(parts)


FRACTION = "a fraction between 0 and 1"

CT01 = Spec(
    id="ct01_encode_capture",
    title="How much of the ENCODE reference does one library recover",
    paper_anchor="Fig. 5b (precision and recall of ENCODE capture); Methods, peak filtering",
    category="qc",
    failure_modes=[
        "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        "metric_statistic_mismatch.wrong_metric_summary_rank_or_association",
    ],
    data=["ac_macs2", "encode_h3k27ac", "blacklist", "chrom_sizes"],
    prompt=f"""You are given peak calls from an H3K27ac CUT&Tag library made in human K562 cells (`{AC_MACS2}`: MACS2 narrow peaks at q < 1e-5, PCR duplicates removed, hg19) and the ENCODE H3K27ac ChIP-seq replicated peak set for K562 (`{ENC_AC}`, hg19). {COMMON_FILES}

The lab wants to know how far this CUT&Tag library can stand in for the ENCODE ChIP-seq reference before it is used to annotate regulatory elements. Compare the two peak sets. Report how many CUT&Tag peaks you carried into the comparison, what fraction of those peaks overlap an ENCODE peak, and what fraction of ENCODE peaks are overlapped by a CUT&Tag peak. {OVERLAP_RULE}""",
    fields=[
        Field("n_library_peaks", "int", {"value": 300}, 0,
              "number of CUT&Tag peaks in your comparison"),
        Field("frac_library_in_encode", "float", {"value": 0.02}, 4,
              f"share of those CUT&Tag peaks that overlap an ENCODE peak, {FRACTION}"),
        Field("frac_encode_captured", "float", {"value": 0.01}, 4,
              f"share of ENCODE peaks overlapped by a CUT&Tag peak, {FRACTION}"),
    ],
    intent="Test whether the agent removes known artifact regions before comparing two peak sets, "
           "and keeps the two directions of the comparison apart.",
    trap="The peak file is compared as delivered. About 1 in {one_in} of this library's peaks sits in "
         "a blacklisted region or on an unplaced contig, and only {removed_supported} of those have an "
         "ENCODE partner, so the share of CUT&Tag peaks in ENCODE comes out {points} points too low "
         "and the peak count too high.",
    tolerance_rationale="The windows accept blacklist-only filtering, 1-based closed coordinates, "
                        "and filtering the ENCODE set as well. They do not accept the unfiltered "
                        "file, swapped denominators, or base-pair overlap.",
    not_said="The prompt lists the blacklist among the files but does not say to use it, and does "
             "not use the words precision or recall.",
    provenance="Same filter as the paper (blacklist v2, standard chromosomes) and the same "
               "overlap rule as its GenomicRanges code.",
)

CT02 = Spec(
    id="ct02_strong_vs_weak_recovery",
    title="Which ENCODE peaks are recovered, and what a missing peak means",
    paper_anchor="Fig. 5c (ENCODE q-values of captured and missed peaks)",
    category="qc",
    failure_modes=[
        "marker_reference_label_scope.dropout_panel_limit_or_absence_overclaim",
        "biological_interpretation.prior_driven_or_overclaimed_biological_interpretation",
    ],
    data=["ac_macs2", "encode_h3k27ac", "blacklist", "chrom_sizes"],
    prompt=f"""You are given peak calls from an H3K27ac CUT&Tag library made in human K562 cells (`{AC_MACS2}`: MACS2 narrow peaks at q < 1e-5, PCR duplicates removed, hg19) and the ENCODE H3K27ac ChIP-seq replicated peak set for K562 (`{ENC_AC}`, hg19). {COMMON_FILES}

A collaborator wants to use this CUT&Tag peak set to decide which K562 regulatory elements carry H3K27ac, including some weak enhancers. Before they do, establish whether the library recovers the ENCODE reference evenly or whether recovery depends on how strong the ENCODE peak is.

Rank the ENCODE peaks by ENCODE's own q-value (column 9 of the narrowPeak file, -log10 q). Report the fraction of ENCODE peaks overlapped by a CUT&Tag peak among the strongest 10% and among the weakest 10% of the ENCODE set, and the median -log10 q of the ENCODE peaks that are recovered and of those that are missed. {OVERLAP_RULE} Then answer the collaborator's question below.""",
    fields=[
        Field("recall_strongest_decile", "float", {"value": 0.03}, 4,
              f"share of the strongest 10% of ENCODE peaks overlapped by a CUT&Tag peak, {FRACTION}"),
        Field("recall_weakest_decile", "float", {"value": 0.04}, 4,
              f"the same for the weakest 10% of ENCODE peaks, {FRACTION}"),
        Field("median_neglog10q_recovered", "float", {"value": 5.0}, 2,
              "median ENCODE -log10 q of the ENCODE peaks overlapped by a CUT&Tag peak"),
        Field("median_neglog10q_missed", "float", {"value": 1.5}, 2,
              "median ENCODE -log10 q of the ENCODE peaks with no CUT&Tag peak"),
        Field("interpretation", "choice"),
    ],
    choice_question="A weak ENCODE H3K27ac peak (in the weakest tenth of the ENCODE set) at an "
                    "enhancer of interest has no CUT&Tag peak in this library. Which conclusion is "
                    "best supported? Answer with the option number.",
    choices=[
        "The enhancer does not carry H3K27ac in K562. A site that this library misses is most likely a ChIP-seq false positive, so the ENCODE peak should be dropped from the reference before the library is used.",
        "CUT&Tag reports H3K27ac at promoters but not at enhancers, so the result is uninformative for this site and for every other enhancer in the reference, whatever the strength of the ENCODE peak.",
        "The library recovers most ENCODE peaks at every strength, so a missing peak is a reliable sign that the site lacks the mark, and the collaborator can treat absent peaks as negative calls.",
        "The missing peak says little either way. The library recovers only a small minority of ENCODE peaks of that strength, so a miss is expected whether or not the mark is there, and it needs another check.",
        "The library has failed quality control. Its recovery of the ENCODE reference is too low for any of its peaks to be used, strong or weak, and the experiment should be repeated before any annotation is attempted.",
        "Recovery is about the same across the ENCODE set, so a missing weak peak is as informative as a missing strong peak and counts against the mark at this enhancer just as strongly.",
        "A weak ENCODE peak is expected to be recovered by a library of this quality, so its absence points to a real difference between this K562 stock and the one that ENCODE profiled.",
        "Weak ENCODE peaks are recovered more often than strong ones in this library, so the absence of a peak at this enhancer is strong evidence against the mark being present there.",
    ],
    correct_choice=CT02_CORRECT,
    key_review="Decision key checked against the data by script. Not yet reviewed by a person.",
    intent="Test whether the agent treats recovery as a function of signal strength and refuses "
           "to read a missing peak as a missing mark.",
    trap="One overall recovery figure is carried into the interpretation, and the absence of a "
         "peak is read as absence of H3K27ac. Recovery here is {fold} times higher in the "
         "strongest tenth of the reference than in the weakest tenth.",
    tolerance_rationale="The windows accept ranking by signal value or p-value instead of q-value, "
                        "unfiltered CUT&Tag peaks, and either rounding of the decile size. They do "
                        "not accept an unstratified answer or stratifying the CUT&Tag peaks instead "
                        "of the ENCODE peaks.",
    not_said="The prompt asks for the numbers but does not say what they imply for a missing peak.",
    provenance="The paper reports the same contrast across 8 samples (captured peaks are the more "
               "significant ENCODE peaks). The decile summary is this benchmark's own readout.",
)

CT03 = Spec(
    id="ct03_duplicates_decision",
    title="Keep or remove duplicates in a heavily duplicated library",
    paper_anchor="Fig. 5b, Supplementary Fig. 3e and 5a (duplicates, by peak caller)",
    category="peak_calling",
    failure_modes=[
        "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        "qc_control_artifact.context_inappropriate_qc_threshold",
        "metric_statistic_mismatch.wrong_metric_summary_rank_or_association",
    ],
    data=["ac_macs2", "ac_macs2_withdup", "ac_seacr", "ac_seacr_withdup", "ac_picard",
          "encode_h3k27ac", "blacklist", "chrom_sizes"],
    prompt=f"""An H3K27ac CUT&Tag library from human K562 cells was sequenced deeply. Its Picard MarkDuplicates report is in the workspace (`{AC_PICARD}`). Two analysts disagree about whether duplicate fragments should be kept before peak calling, so peaks were called four ways (hg19):

- `{AC_MACS2}`: MACS2 narrow peaks at q < 1e-5, duplicates removed
- `{AC_MACS2_DUP}`: MACS2, same settings, duplicates kept
- `{AC_SEACR}`: SEACR stringent, top 1%, duplicates removed
- `{AC_SEACR_DUP}`: SEACR, same settings, duplicates kept

The ENCODE H3K27ac ChIP-seq replicated peak set for K562 (`{ENC_AC}`) is provided as the reference. {COMMON_FILES}

Work out what keeping duplicates does to each caller's peak set and recommend how duplicates should be handled for this library. Report the duplication rate from the Picard report; for each of the four peak sets, the fraction of its peaks that overlap an ENCODE peak; the number of MACS2 peaks that appear only when duplicates are kept (peaks in the duplicates-kept set that overlap no peak in the duplicates-removed set); and the fraction of those duplicate-only peaks that overlap an ENCODE peak. {OVERLAP_RULE}""",
    fields=[
        Field("duplication_rate", "float", {"value": 0.01}, 4,
              f"fraction of read pairs marked as duplicates, {FRACTION}"),
        Field("frac_in_encode_macs2_dedup", "float", {"value": 0.08}, 4,
              f"share of MACS2 duplicates-removed peaks that overlap an ENCODE peak, {FRACTION}"),
        Field("frac_in_encode_macs2_withdup", "float", {"value": 0.05}, 4,
              "the same for MACS2 with duplicates kept"),
        Field("frac_in_encode_seacr_dedup", "float", {"value": 0.08}, 4,
              "the same for SEACR with duplicates removed"),
        Field("frac_in_encode_seacr_withdup", "float", {"value": 0.08}, 4,
              "the same for SEACR with duplicates kept"),
        Field("n_macs2_duplicate_only_peaks", "int", {"value": 8000}, 0,
              "MACS2 duplicates-kept peaks that overlap no MACS2 duplicates-removed peak"),
        Field("frac_duplicate_only_in_encode", "float", {"value": 0.03}, 4,
              f"share of those duplicate-only peaks that overlap an ENCODE peak, {FRACTION}"),
        Field("recommendation", "choice"),
    ],
    choice_question="Which handling of duplicates is best supported for this library? Answer with "
                    "the option number.",
    choices=[
        "Keep duplicates for both callers. In CUT&Tag, fragments with identical ends mostly come from Tn5 cutting the same site in different cells, and the gain in ENCODE peaks recovered when duplicates are kept confirms they are real signal.",
        "Keep duplicates for MACS2 and remove them for SEACR, because SEACR loses ENCODE support when duplicates are kept and MACS2 does not.",
        "Keep duplicates and raise the MACS2 significance threshold instead. The peaks that appear only with duplicates are low-confidence calls, so a stricter cutoff removes them and keeps the extra sensitivity.",
        "It makes no practical difference for either caller at this duplication rate, so the choice can follow whatever the rest of the pipeline expects.",
        "Remove duplicates for both callers, because both lose most of their ENCODE-supported peaks when duplicates are kept.",
        "Remove duplicates before calling peaks with MACS2: keeping them multiplies the peak count, and most of the added peaks have no ENCODE support. The SEACR peak set is nearly the same either way.",
        "Discard the library. At this duplication rate neither caller gives a usable peak set, with or without duplicates.",
        "Keep duplicates for MACS2, because a larger peak set gives more statistical power downstream and the share of peaks in ENCODE is unchanged.",
    ],
    correct_choice=CT03_CORRECT,
    key_review="Decision key checked against the data by script. Not yet reviewed by a person.",
    intent="Test whether the agent sets aside the usual CUT&Tag advice to keep duplicates when the "
           "library in front of it has a duplication rate of {dup_pct}, and whether it reads the gain "
           "in recovered ENCODE peaks next to the loss of support per peak.",
    trap="Duplicates are kept because CUT&Tag duplicates are often genuine, or because more ENCODE "
         "peaks are recovered with them. With MACS2 that multiplies the peak count {fold}-fold, and "
         "{added_supported} of the added peaks have an ENCODE partner. The added peaks are not weak "
         "calls either: their median -log10 q in the duplicates-kept file is {median_q}, so a stricter "
         "cutoff does not remove them.",
    tolerance_rationale="This task isolates the duplicates decision, so the windows are wide enough "
                        "to accept unfiltered peak files, blacklist-only filtering and 1-based "
                        "closed coordinates. They do not accept reporting ENCODE recovery in place "
                        "of the share of peaks in ENCODE, or defining duplicate-only peaks by exact "
                        "coordinates.",
    not_said="The prompt does not say that the duplication rate is high, that the callers behave "
             "differently, or which metric should drive the recommendation.",
    provenance="The paper recommends removing duplicates for MACS2 at high duplication rates and "
               "reports SEACR as robust. The four peak files are from the paper's code repository, "
               "whose MACS2 script calls the duplicates-kept set with --keep-dup all.",
)

CT04 = Spec(
    id="ct04_antibody_panel_reference",
    title="Score a 20-library antibody panel against the right reference",
    paper_anchor="Fig. 2e and Supplementary Fig. 1f (PCR cycles and DNA extraction, matched depth)",
    category="qc",
    failure_modes=[
        "subset_scope_compartment.wrong_comparator_contrast_or_baseline",
        "assay_intervention_feature_linkage.barcode_feature_observation_or_metadata_linkage_error",
    ],
    data=["panel", "encode_h3k27ac", "encode_h3k27me3", "blacklist", "chrom_sizes"],
    prompt=f"""A lab optimised its CUT&Tag protocol in human K562 cells by testing five antibodies under four library preparation conditions (11 or 13 PCR cycles, SDS-based or column-based DNA extraction): 20 libraries in all. Every library was downsampled to a common read depth before peak calling, and peaks were called with SEACR (hg19). `samplesheet.csv` describes each library and names its peak file. Two ENCODE K562 ChIP-seq peak sets are provided as references: `{ENC_AC}` and `{ENC_ME3}`. {COMMON_FILES}

The lab's acceptance bar for a library is that at least 75% of its peaks overlap an ENCODE peak. Evaluate the panel against ENCODE. Report how many of the 20 libraries meet the bar, the lowest share of peaks overlapping ENCODE seen in any library, and the mean share of ENCODE peaks recovered across the 20 libraries. {OVERLAP_RULE}""",
    fields=[
        Field("n_libraries_meeting_bar", "int", {"lower": 2, "upper": 1}, 0,
              "number of libraries with at least 75% of peaks overlapping an ENCODE peak"),
        Field("lowest_library_precision", "float", {"lower": 0.10, "upper": 0.05}, 4,
              f"the smallest share of peaks overlapping ENCODE among the 20 libraries, {FRACTION}"),
        Field("mean_encode_recovery", "float", {"value": 0.015}, 4,
              f"share of ENCODE peaks overlapped by a library's peaks, averaged over the 20 libraries, {FRACTION}"),
    ],
    intent="Test whether the agent reads the sample sheet and scores each library against the "
           "ENCODE set for the mark its antibody targets. {n_me3} of the {n_libraries} libraries were "
           "made with an H3K27me3 antibody.",
    trap="All {n_libraries} peak files are looped over one reference. Against the H3K27ac set, the "
         "{n_me3} H3K27me3 libraries score {wrong_max} or less and look like failed experiments. "
         "Against the reference for their own mark, {me3_clearing} of the {n_me3} clear the bar.",
    tolerance_rationale="The windows accept unfiltered peak files (which move one library across "
                        "the bar), blacklist-only filtering and 1-based closed coordinates. They do "
                        "not accept one reference for all libraries, or both references pooled.",
    not_said="The prompt never says that the panel contains two marks. The only place that says so "
             "is the target column of the sample sheet.",
    provenance="Peak files are the 20 depth-matched SEACR sets in the paper's GEO deposit "
               "(GSE286492). The 75% bar is the precision threshold the paper used when choosing "
               "peak-calling settings. The sample sheet was written for this benchmark from the "
               "GEO sample names and the antibody targets given in the paper.",
)

CT05 = Spec(
    id="ct05_promoter_denominator",
    title="Does the wider-peak caller recover more promoters",
    paper_anchor="Fig. 6a and Supplementary Fig. 6c-e (promoter overlap by peak caller)",
    category="annotation",
    failure_modes=[
        "metric_statistic_mismatch.wrong_metric_summary_rank_or_association",
        "biological_interpretation.sign_direction_temporal_or_trajectory_orientation_error",
    ],
    data=["ac_macs2", "ac_seacr", "refseq_gtf", "blacklist", "chrom_sizes"],
    prompt=f"""Peaks were called on one H3K27ac CUT&Tag library from human K562 cells with two peak callers (hg19, PCR duplicates removed): `{AC_MACS2}` (MACS2 narrow peaks at q < 1e-5) and `{AC_SEACR}` (SEACR stringent, top 1%). A RefSeq annotation (`{GTF}`) is provided. {COMMON_FILES}

A colleague has read that SEACR peaks are more promoter-focused than MACS2 peaks and wants to switch callers in order to recover more promoters. Check whether this library supports that. Define a promoter as the window from 1,500 bp upstream to 1,500 bp downstream of a transcription start site, and use every unique transcript start site in the RefSeq annotation on chr1-22, X and Y. For each caller, report the fraction of promoters that are overlapped by at least one peak, and the median peak width in bp. {OVERLAP_RULE}""",
    fields=[
        Field("frac_promoters_with_peak_macs2", "float", {"value": 0.025}, 4,
              f"share of promoters overlapped by at least one MACS2 peak, {FRACTION}"),
        Field("median_peak_width_macs2", "float", {"value": 40}, 0, "median MACS2 peak width, bp"),
        Field("frac_promoters_with_peak_seacr", "float", {"value": 0.025}, 4,
              f"share of promoters overlapped by at least one SEACR peak, {FRACTION}"),
        Field("median_peak_width_seacr", "float", {"value": 150}, 0, "median SEACR peak width, bp"),
        Field("conclusion", "choice"),
    ],
    choice_question="Which statement is best supported? Answer with the option number.",
    choices=[
        "SEACR recovers clearly more promoters than MACS2: the difference is more than 10 percentage points of all promoters, in line with its more promoter-focused peaks, so switching callers would recover more promoters.",
        "MACS2 recovers clearly more promoters than SEACR: the difference is more than 10 percentage points of all promoters, because it calls nearly twice as many peaks and each peak marks a different promoter.",
        "The two callers recover a similar share of promoters (within 5 percentage points). SEACR calls fewer, several-fold wider peaks, so more of its peaks touch a promoter without more promoters being recovered.",
        "The two callers recover a similar share of promoters (within 5 percentage points) and call peaks of similar width, so the choice of caller changes neither promoter-level nor peak-level summaries of this library.",
        "Neither caller recovers promoters at a rate above what the same number of randomly placed peaks of the same widths would give, so neither peak set carries promoter information in this library.",
        "SEACR recovers fewer than half as many promoters as MACS2, a difference of more than 10 percentage points, because it calls about half as many peaks and each peak can cover only one promoter.",
    ],
    correct_choice=CT05_CORRECT,
    key_review="Decision key checked against the data by script. Not yet reviewed by a person.",
    intent="Test whether the agent answers a question about promoters with promoters in the "
           "denominator, and whether it reads strand correctly when placing start sites.",
    trap="The share of each caller's peaks that touch a promoter is reported as promoter recovery. "
         "SEACR wins that comparison by {points} points, because it calls {pct_as_many} as many peaks "
         "as MACS2 and its peaks are {width_ratio} times wider. Counted per promoter, the two callers "
         "are {promoter_gap} points apart.",
    tolerance_rationale="The windows accept gene-level counting, 3,000 bp half-open windows, start "
                        "sites collapsed across strands, and unfiltered peak files. They do not "
                        "accept the per-peak share, or start sites placed at the GTF start "
                        "coordinate regardless of strand.",
    not_said="The prompt repeats the colleague's claim without saying what is wrong with it.",
    provenance="The paper reports a higher share of SEACR peaks at promoters with comparable total "
               "promoter capture between callers. It defined promoters with a UCSC knownGene "
               "database; this task uses the RefSeq GTF shipped in the paper's code repository.",
)

CT06 = Spec(
    id="ct06_open_chromatin_control",
    title="Accessibility or the mark: using a second library as the control",
    paper_anchor="Fig. 3c and the short-fragment analysis in Results (open chromatin bias)",
    category="qc",
    failure_modes=[
        "metric_statistic_mismatch.null_background_or_control_calibration_error",
        "qc_control_artifact.technical_artifact_control_or_detection_bias_as_biology",
        "biological_interpretation.prior_driven_or_overclaimed_biological_interpretation",
    ],
    data=["ac_fragments", "me3_fragments", "atac", "encode_h3k27ac", "encode_h3k27me3",
          "blacklist", "chrom_sizes"],
    prompt=f"""CUT&Tag tethers the Tn5 transposase to an antibody, and Tn5 on its own prefers accessible chromatin. A reviewer asks how much of an H3K27ac CUT&Tag profile from human K562 cells could reflect accessibility rather than H3K27ac.

You have 2 million deduplicated fragments from the H3K27ac library (`{AC_FRAGS}`) and 2 million from an H3K27me3 CUT&Tag library made in the same cell line with the same protocol (`{ME3_FRAGS}`); K562 ATAC-seq peaks (`{ATAC}`); and ENCODE K562 ChIP-seq peak sets for H3K27ac (`{ENC_AC}`) and H3K27me3 (`{ENC_ME3}`). All files are hg19. {COMMON_FILES}

Use the ATAC peaks as the set of accessible sites, and define three classes of region:

- marked and accessible: ENCODE H3K27ac peaks that overlap an ATAC peak
- marked, not accessible: ENCODE H3K27ac peaks that overlap no ATAC peak
- accessible, unmarked: ATAC peaks that overlap neither an ENCODE H3K27ac peak nor an ENCODE H3K27me3 peak

Report the fraction of fragments from each library that overlap the classes named in the answer fields, the fraction of H3K27me3-library fragments that overlap any ATAC peak, and that same fraction after discarding fragments shorter than 100 bp. A fragment overlaps a region if they share at least 1 bp. Use every fragment in the file as the denominator. Then choose the conclusion that is best supported.""",
    fields=[
        Field("h3k27ac_lib_marked_accessible", "float", {"value": 0.02}, 4,
              "share of H3K27ac-library fragments in the marked and accessible class"),
        Field("h3k27ac_lib_marked_not_accessible", "float", {"value": 0.006}, 4,
              "share of H3K27ac-library fragments in the marked, not accessible class"),
        Field("h3k27ac_lib_accessible_unmarked", "float", {"value": 0.007}, 4,
              "share of H3K27ac-library fragments in the accessible, unmarked class"),
        Field("h3k27me3_lib_marked_accessible", "float", {"value": 0.004}, 4,
              "share of H3K27me3-library fragments in the marked and accessible class"),
        Field("h3k27me3_lib_accessible_unmarked", "float", {"value": 0.004}, 4,
              "share of H3K27me3-library fragments in the accessible, unmarked class"),
        Field("h3k27me3_lib_any_atac", "float", {"value": 0.006}, 4,
              "share of H3K27me3-library fragments overlapping any ATAC peak"),
        Field("h3k27me3_lib_any_atac_min100bp", "float", {"value": 0.006}, 4,
              "the same share among H3K27me3-library fragments of 100 bp or longer"),
        Field("interpretation", "choice"),
    ],
    choice_question="Which conclusion is best supported? Answer with the option number.",
    choices=[
        "The H3K27ac library shows no enrichment at H3K27ac sites. Its fragments are spread across the three region classes in proportion to their size, as expected for random tagmentation of the genome.",
        "The H3K27ac library should be read as an accessibility profile. About a quarter of its fragments fall in ATAC peaks, which is what untargeted tagmentation of open chromatin produces, so the mark adds little beyond accessibility in this library.",
        "About 8% of H3K27me3-library fragments fall in ATAC peaks. Taking that as the rate of untargeted tagmentation, about a third of the H3K27ac library's fragments at accessible sites should be discounted as background before any peak is interpreted.",
        "The H3K27me3 library is cross-contaminated with H3K27ac signal. That explains its fragments at accessible sites, and it makes the library unusable as a comparison for the H3K27ac profile.",
        "The signal follows the mark, not accessibility. At marked accessible sites the H3K27ac library holds about twenty times the H3K27me3 library's share of fragments, and several times its own share at unmarked accessible sites.",
        "Fragments under 100 bp carry the specific H3K27me3 signal, so discarding them removes on-target fragments, lowers the share at accessible sites for that reason, and makes the H3K27me3 library less specific.",
        "Accessibility explains the H3K27me3 library better than the H3K27ac library, because H3K27me3 domains are the most accessible part of the genome and hold most of the ATAC peaks in these cells.",
        "The ATAC peak set is not a usable reference for accessibility here, because most ENCODE H3K27ac peaks do not overlap it and most ATAC peaks fall inside H3K27me3 domains instead.",
    ],
    correct_choice=CT06_CORRECT,
    key_review="Decision key checked against the data by script. Not yet reviewed by a person. This task goes beyond the paper's own analysis, so it is first in line for the author's review.",
    intent="Test whether the agent can use a library made with a different antibody as the control "
           "for antibody-independent tagmentation, and whether it reads that control at the sites in "
           "question instead of applying one genome-wide figure.",
    trap="The H3K27me3 library's fragments in ATAC peaks ({any_pct}) are taken as the rate of "
         "untargeted tagmentation and applied to the H3K27ac library. {in_domains} of those fragments "
         "also overlap an ENCODE H3K27me3 peak, where on-target signal and accessibility cannot be "
         "separated, so that figure is an upper bound. At the accessible sites that carry H3K27ac, "
         "the H3K27me3 library holds {control_marked} of its fragments.",
    tolerance_rationale="The windows accept counting fragments on standard chromosomes only and "
                        "1-based closed coordinates. They do not accept leaving H3K27me3 domains "
                        "out of the definition of unmarked sites, or using all ATAC peaks for "
                        "every accessible class.",
    not_said="The prompt defines the region classes but does not say which comparison answers the "
             "reviewer.",
    provenance="The paper splits ENCODE H3K27ac and ATAC peaks into shared and exclusive sets and "
               "reports the short-fragment effect for H3K27me3. Splitting unmarked accessible "
               "sites by H3K27me3 domains, and reading the H3K27me3 library as a control, is this "
               "benchmark's own analysis of the paper's public files. It is not a claim made in "
               "the paper.",
)

SPECS = {s.id: s for s in (CT01, CT02, CT03, CT04, CT05, CT06)}
