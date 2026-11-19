# Pilot: which peak set would you take forward?

**Independent prototype.** One task in the style of LatchBio's EpiBench. It is
built from a published paper that the author of this repository co-wrote, as a
co-first author. It is not part of EpiBench and has not been reviewed by, or
run at, LatchBio.

Everything here was built and checked by AI models (Claude) under the author's
direction. No person has checked the code, the key or these pages yet.

The paper is Abbasova, Urbanaviciute, Hu et al., Nature Communications 2025
([doi:10.1038/s41467-025-58137-2](https://doi.org/10.1038/s41467-025-58137-2),
GEO GSE286492). Figure 4a chooses peak-calling settings with one criterion:
recover as much of the ENCODE ChIP-seq peak set as possible while keeping
precision above 75%. This pilot turns that criterion into a task.

## The task

An agent gets 2 peak sets called from the same H3K27ac CUT&Tag library,
labelled A and B, and the matched ENCODE reference. It has to say which set the
rule selects (A, B, tie or neither) and report the 7 counts behind the choice:
the size of each file, and the unique peaks that overlap in each direction.
The prompt is [workspace/task.md](workspace/task.md).

## Check it with one command

```bash
python3 pilot/replay.py
```

Standard library only: no install, no download, no model, no API key. It takes
about a minute. It verifies the staged inputs, recomputes the answer 2
independent ways, replays a correct solution and 9 scripted wrong ones through
the grader, tests the label swap and the 75% line, and compares everything
with the recorded results and documents. It ends with `PILOT REPLAY: PASS` or
`FAIL`.

To grade one answer:

```bash
python3 pilot/grader.py --answer my_answer.json --key pilot/private/answer_key.json
```

## What is where

| Path | What it is | Given to the agent |
|---|---|---|
| `workspace/task.md` | the prompt | yes |
| `workspace/manifest.json` | the 3 input files with rows, sizes and SHA-256 | yes |
| `workspace/answer_schema.json` | JSON Schema of the answer | yes |
| `workspace/inputs/` | `peaks_A.bed`, `peaks_B.bed`, `reference_peaks.bed` | yes |
| `eval.json` | the task as a public EpiBench-style entry: prompt and metadata, no key | n/a |
| `grader.py` | the deterministic grader | no |
| `replays.py` | the scripted correct and wrong analyses | no |
| `replay.py` | the one command above | no |
| `docs.py` | builds `VERIFICATION.md` and the results block of this page from the recorded results | no |
| `chromosome_lengths.json` | hg19 and hg38 chromosome lengths, for the replay's genome-build check | no |
| `results/` | the replay record with every answer taken out | no |
| `VERIFICATION.md` | what was checked about the inputs, and how | no |
| `private/` | the answer key, the labels, both solutions, the 10 check pairs, the full replay record, agent runs | no |
| `private/labels.json` | which label each peak caller got in each instance; frozen, so a rebuild cannot reshuffle them | no |
| `record_agent_run.py` | grades one agent run and files its record under `private/agent_runs/` | no |
| `build_pilot.py` | rebuilds the workspace, `eval.json`, the keys, the provenance record and the check pairs from the public data | no |

To hand the task to an agent, copy `workspace/` and nothing else.

**This repository contains the answer.** `private/` is kept out of the agent's
workspace, not out of the repository, so that the pilot can be verified from a
clone. A held-out release would leave that folder out. The pages outside it
state no key value and do not say which peak caller is A.

## What the replay recorded

<!-- BEGIN GENERATED -->
**The answer is not printed on this page.** It is in `private/answer_key.json`. A binary search in Python and a range join in SQLite agree on 8 of 8 answer fields and both equal the key, and 8 of 8 numbers equal what the paper's own notebook prints.

**Scripted replays on the primary instance.**

| Replay | What it does | Counts wrong, of 7 | Choice | Grader |
|---|---|---|---|---|
| `correct` | The reference solution: unique peaks in each direction, then the rule. | 0 | right | PASS |
| `overlapping_pairs` | Counts overlapping (peak, reference peak) pairs, as an intersect piped to a line count does. | 4 | wrong | fail |
| `swapped_directions` | Reports each directional count under the other direction's name. | 4 | wrong | fail |
| `paper_recall_guess` | Computes no overlaps. Assumes both sets recover the paper-wide 54% and just clear the bar. | 4 | wrong | fail |
| `swapped_denominators` | Divides each count by the other file's size: precision over the reference, recovery over the set. | 0 | wrong | fail |
| `ignore_precision_bar` | Takes whichever set recovers more reference peaks and never checks the 75% requirement. | 0 | right | PASS |
| `at_least_75` | Reads the requirement as 75% or more instead of strictly more than 75%. | 0 | right | PASS |
| `rounded_percent` | Rounds precision to a whole percent before comparing it with 75. | 0 | right | PASS |
| `touching_counts_as_overlap` | Treats coordinates as 1-based closed, so intervals that only touch count as overlapping. | 3 | right | fail |
| `fractions_instead_of_counts` | Reports precision and recovery as fractions where the schema asks for counts. | not a valid answer | wrong | fail |

**Shortcuts that survive on the primary instance:** `ignore_precision_bar`, `at_least_75`, `rounded_percent`. Both sets clear the requirement there, so getting the requirement wrong costs nothing.

**Across all 15 instances** (the primary, 2 built from it, 2 cut at the 75% line, and the paper's 10 other public pairs):

| Replay | Passes on | Caught by |
|---|---|---|
| `correct` | 15 of 15 | nothing |
| `overlapping_pairs` | 0 of 15 | all 15 |
| `swapped_directions` | 0 of 15 | all 15 |
| `paper_recall_guess` | 0 of 15 | all 15 |
| `swapped_denominators` | 1 of 15 | 14 instances |
| `ignore_precision_bar` | 10 of 15 | boundary_exactly_75, pair_01, pair_02, pair_04, pair_07 |
| `at_least_75` | 14 of 15 | boundary_exactly_75 |
| `rounded_percent` | 14 of 15 | boundary_one_peak_above_75 |
| `touching_counts_as_overlap` | 0 of 15 | all 15 |
| `fractions_instead_of_counts` | 0 of 15 | all 15 |

**On the 11 real pairs** the rule selects SEACR in 8, MACS2 in 2, neither in 1. Ignoring the 75% requirement gives the wrong answer on 4 of the 10 pairs other than the primary.

**If only the final choice were graded,** picking the set with wider peaks would be right on 8 of the 11 real pairs without counting a single overlap. Requiring the counts closes that.

**Agent runs: 1.** Hands-on, on the author's machine, not benchmark-grade.

| Run | Model | Result | Tool calls | Calls naming a path outside the workspace |
|---|---|---|---|---|
| `run_01` | Claude Sonnet 5.5 | PASS | 4 | 0 |

Full output without answers: [results/REPLAY.md](results/REPLAY.md). Input checks: [VERIFICATION.md](VERIFICATION.md).
<!-- END GENERATED -->

## What agreement with ENCODE establishes

It establishes concordance. A CUT&Tag peak that overlaps an ENCODE peak has
support from a second experiment with different chemistry: cross-linked,
sonicated chromatin in ChIP-seq against tagmentation of native chromatin here.
In the primary pair it is also a different antibody. The share of ENCODE
recovered says how much of that catalogue a configuration finds at this depth.
The 75% requirement is there because recovery can always be bought with more
calls: called by MACS2 with duplicates kept, the same library recovers far more
of ENCODE at a small fraction of the precision (task `ct03_duplicates_decision`
in this repository has the figures).

So the rule is a sound way to rank peak-calling configurations against one
reference.

It does not establish that ENCODE is the truth. The reference is one
thresholded peak set from 2 replicates. It misses weak sites, and some of its
own peaks lie in blacklisted regions. A CUT&Tag peak outside ENCODE has not
been shown to be false, and an ENCODE peak that was not recovered has not been
shown to be absent. The 75% line is the paper's working threshold, not a
property of the biology. And nothing here says which caller is better in
general: the rule does not select the same caller on every public pair, and on
one pair it selects neither (the split is in the results above).

These need an experiment, not a comparison of peak files:

- **Is a CUT&Tag-only peak really H3K27ac?** Measure the locus another way.
  The paper used qPCR at loci chosen from ENCODE's strongest and weakest
  peaks; a second antibody, CUT&RUN, or nucleosome spike-ins that test
  antibody specificity would each address it.
- **Is a missed ENCODE peak absent in these cells?** Deeper sequencing of the
  same library, or ChIP-seq on the same cell stock.
- **Is the signal the antibody or the transposase?** A no-antibody or IgG
  control, and a comparison with accessibility data.
- **Does the peak matter?** H3K27ac marks candidate elements. Activity needs a
  functional readout such as a reporter assay or a perturbation.

## Limits

- **The primary pair is a public demonstration, not a held-out item.** Its
  counts are printed in the paper's rendered notebook, and the row counts of
  the 2 files say which is which. The 6 example tasks in this repository use
  the same library and publish their keys. `private/` here holds the answer.
  Anyone with web access can read it.
- **The blinding hides names, not identity.** One caller's peaks are several
  times wider than the other's on every pair, so a careful agent can tell
  which caller made which set. That alone does not pass the grader, because
  the 7 counts are required.
- **The requirement does not bind on the primary pair.** Both sets clear 75%,
  so an analysis that ignores the requirement, reads it as "at least 75%", or
  rounds before comparing still passes there. Ignoring it is caught by real
  pairs. The other 2 are caught only by the 2 instances built at the 75% line:
  no real pair catches them. A single pair cannot test the rule.
- **Agent runs are few and not benchmark-grade.** The replays are scripted.
  Any agent run listed above was made on the author's own machine, where the
  answer key is on disk. The record lists the agent's tool calls so that can be
  checked, but it is not a sandbox, and one run of one model is not a result.
- **The grader checks the answer, not the method.** A solution with the right
  counts and the right choice passes whatever reasoning produced them.
- **Exact counts depend on the stated convention.** Counting intervals that
  only touch as overlapping moves the counts by a handful and fails. The
  prompt states the convention for that reason.
- **The check pairs are private, not secret.** They are kept out of the
  workspace, but they are in this repository, and anyone can rebuild them from
  GEO.
- **"Independent" means 2 code paths, not 2 people.** Both recomputations were
  written by AI models in the same working session. The agreement with the
  paper's notebook is the one check against something people produced.

## Rebuild from the public data

```bash
.venv/bin/python scripts/fetch_data.py
.venv/bin/python pilot/build_pilot.py
python3 pilot/replay.py --write
```

The inputs are derived from the paper's code repository and from ENCODE file
ENCFF044JNJ. No interval was altered: the files were filtered as in the paper,
reduced to 3 columns, sorted and relabelled. Sources, checksums and every
filter count are in `private/provenance.json`. Terms and attributions for each
source are in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
