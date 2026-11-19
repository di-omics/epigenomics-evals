#!/usr/bin/env python3
"""Grade one agent run on the pilot and file a record of it. Standard library only.

    python3 pilot/record_agent_run.py --run-id run_01 --model "..." --harness "..." \
        --workspace /path/to/the/copy/the/agent/worked/in --final-message reply.txt \
        [--transcript session.jsonl] [--started ...] [--finished ...] [--note "..."]

It reads the agent's final reply, extracts the answer, grades it with grader.py,
checks that the staged inputs were not modified, keeps the scripts the agent
left in work/, and, when a session transcript is given, lists every tool call
and flags any path outside the workspace. The raw transcript is never copied:
only tool names, commands and file paths are kept.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

PILOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PILOT))
import grader  # noqa: E402

RUNS = PILOT / "private" / "agent_runs"
KEEP_SUFFIXES = {".py", ".sh", ".sql", ".r", ".R", ".awk", ".txt", ".json", ".md", ".tsv", ".csv"}
MAX_KEEP_BYTES = 200_000


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tool_calls(transcript: Path) -> list[dict]:
    """Every tool call in a session transcript: the tool and its command or path."""
    calls = []

    def walk(node) -> None:
        if isinstance(node, dict):
            if node.get("type") == "tool_use" and "name" in node:
                given = node.get("input") or {}
                detail = given.get("command") or given.get("file_path") or given.get("path") \
                    or given.get("pattern") or given.get("url") or given.get("query") or ""
                calls.append({"tool": node["name"], "detail": str(detail)})
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    for line in transcript.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            walk(json.loads(line))
        except json.JSONDecodeError:
            continue
    return calls


def outside_paths(calls: list[dict], workspace: Path) -> list[dict]:
    """Tool calls that name an absolute path outside the workspace."""
    roots = {str(workspace), str(workspace.resolve())}
    flagged = []
    for call in calls:
        for path in re.findall(r"(?<![\w.])(/(?:Users|private|tmp|var|Volumes|opt|usr|etc|Library|Applications)[^\s'\";|&)>]*)", call["detail"]):
            if any(path.startswith(root) for root in roots):
                continue
            if path.startswith(("/usr/bin", "/usr/local/bin", "/opt/homebrew", "/tmp/uv", "/private/tmp/uv")) \
                    or "/.local/bin/" in path:
                continue  # running an installed program is not reading user files
            flagged.append({"tool": call["tool"], "path": path, "detail": call["detail"][:300]})
    return flagged


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--harness", required=True)
    ap.add_argument("--workspace", required=True, type=Path)
    ap.add_argument("--final-message", required=True, type=Path)
    ap.add_argument("--transcript", type=Path)
    ap.add_argument("--started", default="")
    ap.add_argument("--finished", default="")
    ap.add_argument("--note", action="append", default=[])
    args = ap.parse_args()

    out = RUNS / args.run_id
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    reply = args.final_message.read_text()
    key = json.loads((PILOT / "private" / "answer_key.json").read_text())["key"]
    try:
        answer = grader.extract_answer(reply)
        verdict = grader.grade(answer, key)
    except (ValueError, json.JSONDecodeError) as exc:
        answer, verdict = None, {"passed": False, "schema_errors": [f"no readable answer: {exc}"],
                                 "fields": {}, "wrong_fields": [], "choice_correct": False,
                                 "choice_follows_from_counts": None}

    manifest = json.loads((PILOT / "workspace" / "manifest.json").read_text())
    untouched = all((args.workspace / f["path"]).exists() and sha256(args.workspace / f["path"]) == f["sha256"]
                    for f in manifest["files"])
    prompt_untouched = (args.workspace / "task.md").read_text() == (PILOT / "workspace" / "task.md").read_text()

    kept = []
    work = args.workspace / "work"
    if work.is_dir():
        for path in sorted(work.rglob("*")):
            rel = path.relative_to(work)
            if (not path.is_file() or path.suffix not in KEEP_SUFFIXES or path.stat().st_size > MAX_KEEP_BYTES
                    or any(part.startswith(".") or part in ("__pycache__", "venv", "site-packages") for part in rel.parts)):
                continue
            (out / "work" / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, out / "work" / rel)
            kept.append(str(rel))
    created = sorted(str(p.relative_to(args.workspace)) for p in args.workspace.rglob("*")
                     if p.is_file() and not any(part.startswith(".") or part in ("__pycache__", "site-packages")
                                                for part in p.relative_to(args.workspace).parts)
                     and str(p.relative_to(args.workspace)) not in {f["path"] for f in manifest["files"]}
                     | {"task.md", "manifest.json", "answer_schema.json"})

    record = {
        "run_id": args.run_id, "task_id": manifest["task_id"], "instance": "primary",
        "model": args.model, "harness": args.harness, "started": args.started, "finished": args.finished,
        "answer": answer, "verdict": verdict,
        "inputs_unmodified_after_the_run": untouched, "prompt_unmodified_after_the_run": prompt_untouched,
        "files_the_agent_created": created, "scripts_kept": kept, "notes": args.note,
    }
    if args.transcript and args.transcript.exists():
        calls = tool_calls(args.transcript)
        flagged = outside_paths(calls, args.workspace)
        by_tool: dict[str, int] = {}
        for call in calls:
            by_tool[call["tool"]] = by_tool.get(call["tool"], 0) + 1
        record["tool_calls"] = {"total": len(calls), "by_tool": by_tool,
                                "naming_a_path_outside_the_workspace": flagged}
        def tidy(text: str) -> str:
            """Keep the record portable: the agent's folder is written as <workspace>."""
            for root in sorted({str(args.workspace), str(args.workspace.resolve())}, key=len, reverse=True):
                text = text.replace(root, "<workspace>")
            return text[:2000]

        (out / "tool_calls.json").write_text(json.dumps(
            [{"tool": c["tool"], "detail": tidy(c["detail"])} for c in calls], indent=1) + "\n")
    shown = reply
    for root in sorted({str(args.workspace), str(args.workspace.resolve())}, key=len, reverse=True):
        shown = shown.replace(root, "<workspace>")
    record["notes"].append("The agent's reply is kept verbatim, except that the absolute path of its folder is written as <workspace>.")
    (out / "final_message.md").write_text(shown if shown.endswith("\n") else shown + "\n")
    (out / "run.json").write_text(json.dumps(record, indent=1) + "\n")

    print(f"{args.run_id}: {'PASS' if verdict['passed'] else 'fail'}; wrong fields: {verdict['wrong_fields'] or 'none'}; "
          f"choice correct: {verdict['choice_correct']}; schema errors: {verdict['schema_errors'] or 'none'}")
    print(f"inputs unmodified: {untouched}; scripts kept: {len(kept)}; files created: {len(created)}")
    if "tool_calls" in record:
        print(f"tool calls: {record['tool_calls']['total']} {record['tool_calls']['by_tool']}; "
              f"naming a path outside the workspace: {len(record['tool_calls']['naming_a_path_outside_the_workspace'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
