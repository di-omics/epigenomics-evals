# Third-party notices

The code in this repository is MIT licensed (see LICENSE). The data it uses
comes from the sources below. Each keeps its own terms. Only the files named
under "Redistributed here" are stored in this repository; everything else is
downloaded by `scripts/fetch_data.py`.

## The paper

Abbasova L, Urbanaviciute P, Hu D, Ismail JN, Schilder BM, Nott A, Skene NG,
Marzi SJ. CUT&Tag recovers up to half of ENCODE ChIP-seq histone acetylation
peaks. Nature Communications 16, 2993 (2025). doi:10.1038/s41467-025-58137-2

## Marzi-lab/CUTnTag-benchmarking-analysis (the paper's code repository)

https://github.com/Marzi-lab/CUTnTag-benchmarking-analysis , commit
52270f5bdf80fe868b74cc07917840929160eab9. MIT licence.

Redistributed here: `pilot/workspace/inputs/peaks_A.bed` and
`pilot/workspace/inputs/peaks_B.bed`. They are the 2 peak files of one library
from that repository, filtered as in the paper, reduced to 3 columns, sorted
and relabelled. No interval was changed.

The licence notice of that repository:

```text
MIT License

Copyright (c) 2022 leylaabbasova

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The repository also holds copies of files that come from other sources (the
ENCODE peak sets, ATAC-seq peaks, the hg19 blacklist, a RefSeq annotation,
chromosome sizes). Those are covered by the entries below, not by this notice.

## GEO GSE286492 (the paper's sequencing deposit)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286492

Redistributed here: `pilot/private/siblings/*/peaks_A.bed.gz` and
`peaks_B.bed.gz`. They are the MACS2 and SEACR peak files of 10 further
libraries from that deposit, treated the same way as the 2 files above. Cite
the paper when using them.

## ENCODE

Redistributed here: `pilot/workspace/inputs/reference_peaks.bed`, the 3
coordinate columns of ENCODE file ENCFF044JNJ (H3K27ac ChIP-seq replicated
peaks, K562, hg19, experiment ENCSR000AKP), sorted.

Fetched, not redistributed: ENCFF000BXB (H3K27me3 ChIP-seq peaks, K562,
experiment ENCSR000EWB), and K562 ATAC-seq peaks that the paper's authors
processed from ENCODE libraries ENCLB918NXF and ENCLB758GEG.

ENCODE data are released for use without restriction. ENCODE asks users to
cite the consortium and the accessions used:
https://www.encodeproject.org/help/citing-encode/

The ENCODE Project Consortium. An integrated encyclopedia of DNA elements in
the human genome. Nature 489, 57-74 (2012).

## ENCODE blacklist, hg19 v2 (Boyle Lab)

https://github.com/Boyle-Lab/Blacklist , GPL-3.0. Fetched and used to filter
the candidate peak sets. Not redistributed here.

Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: Identification of
Problematic Regions of the Genome. Scientific Reports 9, 9354 (2019).

## UCSC Genome Browser

Chromosome sizes for hg19 and hg38, and the hg19 NCBI RefSeq annotation
(`hg19.ncbiRefSeq.gtf.gz`). Fetched. `pilot/chromosome_lengths.json` stores the
lengths of chr1-22, X and Y for both assemblies, taken from the 2 chromosome
size files. The annotation is not redistributed.

## latch-eval-tools (LatchBio)

Installed from PyPI as an optional dependency and called to cross-check the
graders. None of its code is copied into this repository. EpiBench, scBench
and their public example tasks belong to LatchBio; this repository reuses the
JSON keys and grader type names of the public examples and nothing else.
