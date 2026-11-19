"""Where the input files live and how each task reads them."""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .intervals import STANDARD_CHROMS, Intervals, read_bed

REPO = Path(__file__).resolve().parents[2]
DEMO = "marzi_demo/demo"

# Logical name -> path under the data root. Tasks refer to logical names only.
FILES = {
    "blacklist": f"{DEMO}/resources/hg19/hg19-blacklist.v2.bed.gz",
    "chrom_sizes": f"{DEMO}/resources/hg19/hg19.chrom.sizes",
    "refseq_gtf": f"{DEMO}/resources/hg19/hg19.ncbiRefSeq.gtf.gz",
    "encode_h3k27ac": f"{DEMO}/resources/ENCODE/peaks/H3K27ac_ENCFF044JNJ.bed.narrowPeak",
    "encode_h3k27me3": f"{DEMO}/resources/ENCODE/peaks/H3K27me3_ENCFF000BXB.bed.broadPeak",
    "atac": f"{DEMO}/resources/ATAC/peaks/K562_ATACseq.mRp.clN_peaks.narrowPeak",
    "ac_macs2": f"{DEMO}/peakCalling/all_peaks/Diagenode_C15410196_1_50_MACS2.bed",
    "ac_macs2_withdup": f"{DEMO}/peakCalling/all_peaks/Diagenode_C15410196_1_50_MACS2_withDup.bed",
    "ac_seacr": f"{DEMO}/peakCalling/all_peaks/Diagenode_C15410196_1_50_SEACR.bed",
    "ac_seacr_withdup": f"{DEMO}/peakCalling/all_peaks/Diagenode_C15410196_1_50_SEACR_withDup.bed",
    "ac_picard": f"{DEMO}/removeDuplicate/picard_summary/Diagenode_C15410196_1_50_picard.dupMark.txt",
    "ac_fragments": f"{DEMO}/alignment/bed/Diagenode_C15410196_1_50_2M_rmDup_bowtie2.fragments.bed",
    "me3_fragments": f"{DEMO}/alignment/bed/CST9733_1_100_H3K27me3_2M_rmDup_bowtie2.fragments.bed",
}


def data_root() -> Path:
    return Path(os.environ.get("EPIGENOMICS_EVALS_DATA", REPO / "data"))


def path(name: str) -> Path:
    return data_root() / FILES[name]


def available() -> bool:
    return all((data_root() / rel).exists() for rel in FILES.values())


@lru_cache(maxsize=None)
def load(name: str) -> Intervals:
    return read_bed(path(name))


def picard_duplication_rate(name: str = "ac_picard") -> float:
    """PERCENT_DUPLICATION from a Picard MarkDuplicates metrics file."""
    lines = [ln for ln in path(name).read_text().splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    return float(values[header.index("PERCENT_DUPLICATION")])


@lru_cache(maxsize=None)
def refseq_tss(strand_aware: bool = True) -> pd.DataFrame:
    """Unique transcript start sites on chr1-22, X, Y from the RefSeq GTF.

    Columns: chrom, pos (0-based position of the TSS base), strand, gene.
    strand_aware=False reproduces a common bug: taking the GTF start
    coordinate as the TSS for every transcript, minus strand included.
    """
    gtf = pd.read_csv(
        path("refseq_gtf"), sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 6, 8]
    )
    tx = gtf[gtf[2] == "transcript"]
    tx = tx[np.isin(tx[0], STANDARD_CHROMS)]
    gene = tx[8].str.extract(r'gene_name "([^"]+)"')[0]
    gene = gene.fillna(tx[8].str.extract(r'gene_id "([^"]+)"')[0])
    if strand_aware:
        pos = np.where(tx[6] == "+", tx[3] - 1, tx[4] - 1)
    else:
        pos = (tx[3] - 1).to_numpy()
    out = pd.DataFrame(
        {"chrom": tx[0].to_numpy(), "pos": pos, "strand": tx[6].to_numpy(), "gene": gene.to_numpy()}
    )
    return out.drop_duplicates(["chrom", "pos", "strand"]).reset_index(drop=True)


def promoter_windows(tss: pd.DataFrame, half_width: int = 1500, inclusive: bool = True) -> Intervals:
    pos = tss["pos"].to_numpy(dtype=np.int64)
    start = np.maximum(pos - half_width, 0)
    end = pos + half_width + (1 if inclusive else 0)
    return Intervals(tss["chrom"].to_numpy(dtype=object), start, end)


# --- GEO panel used by the antibody panel task -------------------------------

PANEL_SHEET = REPO / "evals" / "ct04_antibody_panel_reference" / "samplesheet.csv"


def panel_sheet() -> pd.DataFrame:
    return pd.read_csv(PANEL_SHEET)


@lru_cache(maxsize=None)
def load_panel_peaks(filename: str) -> Intervals:
    return read_bed(data_root() / "geo" / "peaks" / filename)


def panel_available() -> bool:
    if not PANEL_SHEET.exists():
        return False
    return all((data_root() / "geo" / "peaks" / f).exists() for f in panel_sheet()["peak_file"])
