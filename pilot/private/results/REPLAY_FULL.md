# Replay results, full record

Held-out material: this page states the answer key, the keys of the check instances and the scripted answers.
Written by `python3 pilot/replay.py --write`. Do not edit by hand.

## The primary key, recomputed 2 ways

| Field | Binary search (Python) | Range join (SQLite) | Answer key |
|---|---|---|---|
| n_peaks_a | 11094 | 11094 | 11094 |
| n_peaks_b | 19672 | 19672 | 19672 |
| n_reference_peaks | 51176 | 51176 | 51176 |
| n_a_peaks_overlapping_reference | 10243 | 10243 | 10243 |
| n_b_peaks_overlapping_reference | 16540 | 16540 | 16540 |
| n_reference_peaks_overlapped_by_a | 19190 | 19190 | 19190 |
| n_reference_peaks_overlapped_by_b | 17365 | 17365 | 17365 |
| answer | A | A | A |

Set A: 10,243 of 11,094 (92.3292%) of its peaks overlap the reference; it recovers 19,190 of 51,176 (37.4980%) reference peaks.
Set B: 16,540 of 19,672 (84.0789%) of its peaks overlap the reference; it recovers 17,365 of 51,176 (33.9319%) reference peaks.
Private labels: {'A': 'SEACR', 'B': 'MACS2'}.

| Peak caller | Quantity | Paper's notebook | Recomputed here |
|---|---|---|---|
| MACS2 | peaks after filter | 19,672 | 19,672 |
| MACS2 | peaks overlapping encode | 16,540 | 16,540 |
| MACS2 | encode peaks recovered | 17,365 | 17,365 |
| MACS2 | blacklisted peaks in source | 1,548 | 1,548 |
| SEACR | peaks after filter | 11,094 | 11,094 |
| SEACR | peaks overlapping encode | 10,243 | 10,243 |
| SEACR | encode peaks recovered | 19,190 | 19,190 |
| SEACR | blacklisted peaks in source | 789 | 789 |

## Scripted answers on the primary instance

| Replay | Its choice | Grader | Wrong fields |
|---|---|---|---|
| correct | A | PASS | none |
| overlapping_pairs | B | fail | n_a_peaks_overlapping_reference=19,273; n_b_peaks_overlapping_reference=20,102; n_reference_peaks_overlapped_by_a=19,273; n_reference_peaks_overlapped_by_b=20,102; the choice |
| swapped_directions | B | fail | n_a_peaks_overlapping_reference=19,190; n_b_peaks_overlapping_reference=17,365; n_reference_peaks_overlapped_by_a=10,243; n_reference_peaks_overlapped_by_b=16,540; the choice |
| paper_recall_guess | tie | fail | n_a_peaks_overlapping_reference=8,321; n_b_peaks_overlapping_reference=14,755; n_reference_peaks_overlapped_by_a=27,635; n_reference_peaks_overlapped_by_b=27,635; the choice |
| swapped_denominators | neither | fail | the choice |
| ignore_precision_bar | A | PASS | none |
| at_least_75 | A | PASS | none |
| rounded_percent | A | PASS | none |
| touching_counts_as_overlap | A | fail | n_b_peaks_overlapping_reference=16,541; n_reference_peaks_overlapped_by_a=19,194; n_reference_peaks_overlapped_by_b=17,372 |
| fractions_instead_of_counts | A | fail | not a valid answer: n_a_peaks_overlapping_reference: must be a whole number |

## Every instance

| Instance | Precision A | Precision B | Recovered by A | Recovered by B | Key | Ignoring the bar picks |
|---|---|---|---|---|---|---|
| primary | 92.3292% | 84.0789% | 19,190 | 17,365 | A | A |
| label_swap | 84.0789% | 92.3292% | 17,365 | 19,190 | B | B |
| tie_identical_sets | 92.3292% | 92.3292% | 19,190 | 19,190 | tie | tie |
| boundary_exactly_75 | 75.0000% | 95.0441% | 28,899 | 24,351 | B | A (wrong) |
| boundary_one_peak_above_75 | 75.0027% | 95.0441% | 28,899 | 24,351 | A | A |
| pair_01 | 84.2766% | 63.0692% | 29,059 | 29,660 | A | B (wrong) |
| pair_02 | 64.5766% | 84.6174% | 29,732 | 29,663 | B | A (wrong) |
| pair_03 | 80.1019% | 86.6362% | 20,017 | 25,208 | B | B |
| pair_04 | 74.4916% | 95.0441% | 28,899 | 24,351 | B | A (wrong) |
| pair_05 | 84.9072% | 95.8463% | 25,669 | 25,315 | A | A |
| pair_06 | 94.4679% | 77.8458% | 28,436 | 30,479 | B | B |
| pair_07 | 69.9112% | 61.7674% | 6,310 | 12,848 | neither | B (wrong) |
| pair_08 | 89.7211% | 78.7506% | 23,473 | 22,517 | A | A |
| pair_09 | 77.1726% | 63.1669% | 2,356 | 1,135 | A | A |
| pair_10 | 88.9385% | 92.7530% | 5,924 | 9,280 | B | B |
