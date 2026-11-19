# Replay results

Task `fig4a_peak_set_choice`. Status: independent prototype.
Written by `python3 pilot/replay.py --write`. The replays are scripted; the checks were written and run by AI models, and no person has checked them yet. Do not edit by hand.

No answer is printed on this page. The full record, answers included, is `private/results/REPLAY_FULL.md`.

## Verdicts

| Check | Result |
|---|---|
| inputs verified | yes |
| two recomputations equal the key on every instance | yes |
| correct solution passes everywhere | yes |
| wrong counts and malformed fail everywhere | yes |
| every logic shortcut is caught by at least one instance | yes |
| label swap flips the choice | yes |
| one peak flips the boundary choice | yes |
| recorded agent verdicts reproduce | yes |
| key equals the paper notebook | yes |

## 1. The staged inputs

| Check | Result | Detail |
|---|---|---|
| workspace holds only the prompt, the manifest, the schema and 3 inputs | ok | answer_schema.json, inputs/peaks_A.bed, inputs/peaks_B.bed, inputs/reference_peaks.bed, manifest.json, task.md |
| every input matches its SHA-256 and size in the manifest | ok | 3 files checked |
| no caller, antibody, sample, accession or paper name anywhere in the workspace | ok | 0 matches |
| peaks_A.bed: 3 columns, sorted, no duplicate rows, no zero-length interval | ok | 11,094 rows |
| peaks_A.bed: every interval fits hg19 | ok | 0 end beyond an hg19 chromosome end; 174 end beyond the hg38 end of the same chromosome, so these cannot be hg38 coordinates |
| peaks_B.bed: 3 columns, sorted, no duplicate rows, no zero-length interval | ok | 19,672 rows |
| peaks_B.bed: every interval fits hg19 | ok | 0 end beyond an hg19 chromosome end; 273 end beyond the hg38 end of the same chromosome, so these cannot be hg38 coordinates |
| reference_peaks.bed: 3 columns, sorted, no duplicate rows, no zero-length interval | ok | 51,176 rows |
| reference_peaks.bed: every interval fits hg19 | ok | 0 end beyond an hg19 chromosome end; 579 end beyond the hg38 end of the same chromosome, so these cannot be hg38 coordinates |
| peaks_A.bed: chr1-22, X, Y only | ok | 24 chromosomes |
| peaks_A.bed: no peak overlaps the hg19 blacklist v2 | ok | 0 of 11,094 peaks overlap it; counted when the pilot was built and recounted by this command whenever the blacklist file has been fetched |
| peaks_B.bed: chr1-22, X, Y only | ok | 24 chromosomes |
| peaks_B.bed: no peak overlaps the hg19 blacklist v2 | ok | 0 of 19,672 peaks overlap it; counted when the pilot was built and recounted by this command whenever the blacklist file has been fetched |

## 2. The answer, recomputed 2 ways

A binary search in Python and a range join in SQLite agree on 8 of 8 answer fields, and 8 of 8 equal the key.
8 of 8 numbers also equal what the paper's own rendered notebook prints (computed there in R with GenomicRanges). Neither recomputation reads that notebook.

## 3. Scripted replays on the primary instance

| Replay | Kind | What it does | Counts wrong, of 7 | Choice | Grader |
|---|---|---|---|---|---|
| correct | correct | The reference solution: unique peaks in each direction, then the rule. | 0 | right | PASS |
| overlapping_pairs | wrong_counts | Counts overlapping (peak, reference peak) pairs, as an intersect piped to a line count does. | 4 | wrong | fail |
| swapped_directions | wrong_counts | Reports each directional count under the other direction's name. | 4 | wrong | fail |
| paper_recall_guess | wrong_counts | Computes no overlaps. Assumes both sets recover the paper-wide 54% and just clear the bar. | 4 | wrong | fail |
| swapped_denominators | wrong_logic | Divides each count by the other file's size: precision over the reference, recovery over the set. | 0 | wrong | fail |
| ignore_precision_bar | wrong_logic | Takes whichever set recovers more reference peaks and never checks the 75% requirement. | 0 | right | PASS |
| at_least_75 | wrong_logic | Reads the requirement as 75% or more instead of strictly more than 75%. | 0 | right | PASS |
| rounded_percent | wrong_logic | Rounds precision to a whole percent before comparing it with 75. | 0 | right | PASS |
| touching_counts_as_overlap | convention | Treats coordinates as 1-based closed, so intervals that only touch count as overlapping. | 3 | right | fail |
| fractions_instead_of_counts | malformed | Reports precision and recovery as fractions where the schema asks for counts. | not a valid answer | wrong | fail |

Shortcuts that survive on the primary instance: `ignore_precision_bar`, `at_least_75`, `rounded_percent`.

## 4. Label swap, the 75% line, and a tie

- **Label swap.** With the 2 files exchanged the correct choice changes, and the primary's correct answer fails against the swapped key.
- **Exactly 75%.** Built from a real pair whose higher-recovery set starts at 74.49% precision, with 189 of its peaks removed so that it sits exactly on 75%. Precision 75.0000%: the set is still excluded.
- **One peak past 75%.** The same set with 1 more unsupported peak removed. Precision 75.0027%: it qualifies, and the correct choice flips.
- **Tie.** Both labels hold the same peak set, so both qualify and recover the same reference peaks. Key is tie: yes.

## 5. Where each replay passes and fails

15 instances: the primary; 2 built from it (label swap, tie); 2 cut at the 75% line from a real pair; and the paper's 10 other public pairs.

| Replay | Kind | Passes on | Fails on |
|---|---|---|---|
| correct | correct | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 | 0: none |
| overlapping_pairs | wrong_counts | 0: none | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 |
| swapped_directions | wrong_counts | 0: none | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 |
| paper_recall_guess | wrong_counts | 0: none | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 |
| swapped_denominators | wrong_logic | 1: pair_07 | 14: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_08, pair_09, pair_10 |
| ignore_precision_bar | wrong_logic | 10: primary, label_swap, tie_identical_sets, boundary_one_peak_above_75, pair_03, pair_05, pair_06, pair_08, pair_09, pair_10 | 5: boundary_exactly_75, pair_01, pair_02, pair_04, pair_07 |
| at_least_75 | wrong_logic | 14: primary, label_swap, tie_identical_sets, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 | 1: boundary_exactly_75 |
| rounded_percent | wrong_logic | 14: primary, label_swap, tie_identical_sets, boundary_exactly_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 | 1: boundary_one_peak_above_75 |
| touching_counts_as_overlap | convention | 0: none | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 |
| fractions_instead_of_counts | malformed | 0: none | 15: primary, label_swap, tie_identical_sets, boundary_exactly_75, boundary_one_peak_above_75, pair_01, pair_02, pair_03, pair_04, pair_05, pair_06, pair_07, pair_08, pair_09, pair_10 |

## 6. The real pairs

On the 11 real pairs the rule selects A 5, B 5, neither 1 by label, and SEACR in 8, MACS2 in 2, neither in 1 by peak caller.
Ignoring the 75% requirement gives the wrong answer on 4 of the 10 pairs other than the primary.

## 7. If only the final choice were graded

| Shortcut | Right choice on | Needs to count overlaps |
|---|---|---|
| always_A | 5 of 11 | no |
| always_B | 5 of 11 | no |
| set_with_more_peaks | 4 of 11 | no |
| set_with_fewer_peaks | 6 of 11 | no |
| set_with_wider_peaks | 8 of 11 | no |
| set_with_narrower_peaks | 2 of 11 | no |
| swapped_denominators | 1 of 11 | yes |
| ignore_precision_bar | 7 of 11 | yes |
| at_least_75 | 11 of 11 | yes |
| rounded_percent | 11 of 11 | yes |

The full grader also checks the 7 counts, so none of the guesses in the first group can pass it.

## 8. How far the blinding goes

The workspace carries no names. It does not hide what the 2 sets are: on 11 of 11 pairs the set with the wider peaks comes from the same peak caller.

## 9. Agent runs

| Run | Model | Result | Tool calls | Calls naming a path outside the workspace |
|---|---|---|---|---|
| run_01 | Claude Sonnet 5.5 | PASS | 4 | 0 |

Hands-on runs on the author's machine, where the key is on disk. Not benchmark-grade.
