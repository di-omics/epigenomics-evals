#!/usr/bin/env python3
"""Fetch the public input data for epigenomics-evals.

Two sources, both from the paper's own data availability statement:

1. The paper's code repository (MIT), pinned to one commit:
   https://github.com/Marzi-lab/CUTnTag-benchmarking-analysis
2. The paper's GEO deposit: GSE286492 (peak files for all 38 libraries).

Nothing is executed from the downloaded files. Every file is verified against
the SHA-256 recorded in scripts/data.sha256.

Usage:
    python scripts/fetch_data.py            # fetch what is missing, verify
    python scripts/fetch_data.py --verify   # verify only
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
COMMIT = "52270f5bdf80fe868b74cc07917840929160eab9"
RAW = "https://raw.githubusercontent.com/Marzi-lab/CUTnTag-benchmarking-analysis/" + COMMIT + "/"
GEO_TAR = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286492/suppl/GSE286492_RAW.tar"
# Small files used only to verify the inputs of the pilot (blacklist identity, genome build).
VERIFICATION = {
    "verification/boyle_blacklist/hg19-blacklist.v2.bed.gz":
        "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
    "verification/ucsc/hg38.chrom.sizes":
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, attempts: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "epigenomics-evals-fetch"})
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(dest)
            return
        except OSError:
            if attempt == attempts:
                raise
            time.sleep(5 * attempt)


def repo_files() -> list[str]:
    listing = ROOT / "scripts" / "marzi_demo.files"
    return [ln.strip() for ln in listing.read_text().splitlines() if ln.strip()]


def fetch_all() -> None:
    jobs = []
    for rel in repo_files():
        dest = DATA / "marzi_demo" / rel
        if not dest.exists():
            jobs.append((RAW + rel, dest))
    geo = DATA / "geo" / "GSE286492_RAW.tar"
    if not geo.exists():
        jobs.append((GEO_TAR, geo))
    for rel, url in VERIFICATION.items():
        if not (DATA / rel).exists():
            jobs.append((url, DATA / rel))
    if jobs:
        print(f"fetching {len(jobs)} files")
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(download, url, dest) for url, dest in jobs]
        for future in futures:
            future.result()  # re-raise any download error
    peaks = DATA / "geo" / "peaks"
    if not peaks.exists():
        peaks.mkdir(parents=True)
        with tarfile.open(geo) as tar:
            for member in tar.getmembers():
                name = Path(member.name).name
                if not member.isfile() or not name.endswith(".bed.gz"):
                    continue
                src = tar.extractfile(member)
                if src is not None:
                    (peaks / name).write_bytes(src.read())


def verify() -> int:
    sums = ROOT / "scripts" / "data.sha256"
    files = sorted(p for p in DATA.rglob("*") if p.is_file() and p.name != ".DS_Store")
    current = {str(p.relative_to(DATA)): sha256(p) for p in files}
    if not sums.exists():
        sums.write_text("".join(f"{v}  {k}\n" for k, v in sorted(current.items())))
        print(f"recorded {len(current)} checksums in {sums}")
        return 0
    recorded = {}
    for line in sums.read_text().splitlines():
        digest, name = line.split("  ", 1)
        recorded[name] = digest
    bad = [k for k, v in recorded.items() if current.get(k) != v]
    for k in bad:
        print("MISMATCH or missing:", k)
    print(f"{len(recorded) - len(bad)} of {len(recorded)} files match")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    if not args.verify:
        fetch_all()
    return verify()


if __name__ == "__main__":
    sys.exit(main())
