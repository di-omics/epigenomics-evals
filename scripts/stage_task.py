#!/usr/bin/env python3
"""Stage one task as an agent would receive it: input files and the prompt.

    python scripts/stage_task.py ct03_duplicates_decision /tmp/ct03

Writes the task's input files and task.md into the destination folder. The
answer key and the notes stay behind.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def locate(name: str, task_dir: Path) -> Path:
    for candidate in (task_dir / name, DATA / "geo" / "peaks" / name):
        if candidate.exists():
            return candidate
    hits = list((DATA / "marzi_demo").rglob(name))
    if len(hits) != 1:
        raise FileNotFoundError(f"{name}: found {len(hits)} candidates; run scripts/fetch_data.py")
    return hits[0]


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    task_id, dest = sys.argv[1], Path(sys.argv[2])
    task_dir = ROOT / "evals" / task_id
    spec = json.loads((task_dir / "eval.json").read_text())
    dest.mkdir(parents=True, exist_ok=True)
    for entry in spec["data_files"]:
        shutil.copy2(locate(entry["name"], task_dir), dest / entry["name"])
    (dest / "task.md").write_text(spec["task"] + "\n")
    print(f"staged {len(spec['data_files'])} files and task.md in {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
