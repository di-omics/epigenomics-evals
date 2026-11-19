#!/usr/bin/env python3
"""Check the whole pilot with one command and no model run. Standard library only.

    python3 pilot/replay.py            run every check and compare with the recorded results
    python3 pilot/replay.py --write    run every check and record the results

What it does, in order:
  1. verifies the staged inputs against the manifest and checks them for leaks
  2. recomputes the answer 2 independent ways and compares both with the key
  3. replays one correct solution and the scripted wrong analyses through the grader
  4. tests the label swap, the 75% line and a tie on instances built for that
  5. repeats the replays on the paper's 10 other public pairs
  6. reports every shortcut that survives, and where

Two records are written. private/results/ holds the full record, answers included.
results/ holds a summary with no answer in it, which is what the public pages are built from.

Exit code 0 means every check held. The checks were written and run by AI models.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sqlite3
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

PILOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PILOT))

import grader  # noqa: E402
import replays  # noqa: E402

WORKSPACE = PILOT / "workspace"
PRIVATE = PILOT / "private"
RESULTS = PILOT / "results"
FULL = PRIVATE / "results"
LOCAL_BLACKLIST = PILOT.parent / "data" / "marzi_demo" / "demo" / "resources" / "hg19" / "hg19-blacklist.v2.bed.gz"
STANDARD = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ORDER = {c: i for i, c in enumerate(STANDARD)}
COUNT_FIELDS = grader.COUNT_FIELDS
LABELS = re.compile(
    r"MACS|SEACR|Diagenode|C15410196|Abcam|ab177|ab4729|Active.?Motif|GSM\d|GSE\d|ENCFF|ENCSR|rmDup|10\.1038|Marzi|Abbasova",
    re.I,
)
EXPECTED_WORKSPACE = {"task.md", "manifest.json", "answer_schema.json", "inputs/peaks_A.bed",
                      "inputs/peaks_B.bed", "inputs/reference_peaks.bed"}

Row = tuple


# --- small helpers -----------------------------------------------------------


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[Row]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [(p[0], int(p[1]), int(p[2])) for p in (ln.rstrip("\n").split("\t") for ln in handle if ln.strip())]


def write_rows(rows: list[Row], path: Path) -> None:
    path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))


def sort_rows(rows: list[Row]) -> list[Row]:
    return sorted(rows, key=lambda r: (ORDER[r[0]], r[1], r[2]))


def stage(dest: Path, a: list[Row], b: list[Row], ref: list[Row]) -> Path:
    """Write a workspace that holds the 3 input files and nothing else."""
    (dest / "inputs").mkdir(parents=True)
    write_rows(a, dest / "inputs" / "peaks_A.bed")
    write_rows(b, dest / "inputs" / "peaks_B.bed")
    write_rows(ref, dest / "inputs" / "reference_peaks.bed")
    return dest


def run_script(script: Path, workspace: Path) -> dict:
    """Run a solution in its own process. It is given the workspace path and nothing else."""
    done = subprocess.run([sys.executable, "-I", str(script), str(workspace)],
                          capture_output=True, text=True, check=True)
    return json.loads(done.stdout)


def sql_overlap_flags(query: list[Row], subject: list[Row]) -> list[bool]:
    """Which query intervals share at least 1 bp with a subject interval. SQL, not the reference code."""
    con = sqlite3.connect(":memory:")
    for name, rows in (("q", query), ("r", subject)):
        con.execute(f"CREATE TABLE {name} (id INTEGER PRIMARY KEY, chrom TEXT, s INTEGER, e INTEGER)")
        con.executemany(f"INSERT INTO {name} VALUES (?, ?, ?, ?)", [(i, c, s, e) for i, (c, s, e) in enumerate(rows)])
    con.execute("CREATE INDEX r_pos ON r (chrom, s)")
    widest = max((e - s for _, s, e in subject), default=0)
    hits = {row[0] for row in con.execute(
        "SELECT q.id FROM q WHERE EXISTS (SELECT 1 FROM r WHERE r.chrom = q.chrom "
        "AND r.s < q.e AND r.e > q.s AND r.s > q.s - ?)", (widest,))}
    con.close()
    return [i in hits for i in range(len(query))]


def swap_key(key: dict) -> dict:
    out = dict(key)
    for field in COUNT_FIELDS:
        twin = field.replace("_a_", "_x_").replace("_b_", "_a_").replace("_x_", "_b_")
        if twin.endswith("_a"):
            twin = twin[:-2] + "_b"
        elif twin.endswith("_b"):
            twin = twin[:-2] + "_a"
        out[twin] = key[field]
    out["answer"] = {"A": "B", "B": "A"}.get(key["answer"], key["answer"])
    return out


def ratio(hit: int, total: int) -> str:
    return f"{hit:,} of {total:,} ({100 * hit / total:.4f}%)" if total else "0 of 0"


# --- 1. inputs ---------------------------------------------------------------


def check_inputs(provenance: dict) -> list[dict]:
    checks = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    manifest = json.loads((WORKSPACE / "manifest.json").read_text())
    present = {str(p.relative_to(WORKSPACE)) for p in WORKSPACE.rglob("*") if p.is_file()}
    add("workspace holds only the prompt, the manifest, the schema and 3 inputs",
        present == EXPECTED_WORKSPACE, ", ".join(sorted(present)))
    bad = [f["path"] for f in manifest["files"]
           if sha256(WORKSPACE / f["path"]) != f["sha256"] or (WORKSPACE / f["path"]).stat().st_size != f["bytes"]]
    add("every input matches its SHA-256 and size in the manifest", not bad,
        f"{len(manifest['files'])} files checked" + (f"; mismatch: {bad}" if bad else ""))

    leaks = []
    for rel in sorted(present):
        for i, line in enumerate((WORKSPACE / rel).read_text().splitlines(), 1):
            if LABELS.search(line):
                leaks.append(f"{rel}:{i}")
    add("no caller, antibody, sample, accession or paper name anywhere in the workspace", not leaks,
        f"{len(leaks)} matches" + (f" at {leaks[:5]}" if leaks else ""))

    lengths = json.loads((PILOT / "chromosome_lengths.json").read_text())["lengths"]
    recorded = provenance["instances"][0]["hygiene"]
    rows = {}
    for name in ("peaks_A.bed", "peaks_B.bed", "reference_peaks.bed"):
        path = WORKSPACE / "inputs" / name
        columns = {ln.count("\t") + 1 for ln in path.read_text().splitlines()}
        r = read_rows(path)
        rows[name] = r
        add(f"{name}: 3 columns, sorted, no duplicate rows, no zero-length interval",
            columns == {3} and r == sort_rows(r) and len(set(r)) == len(r) and all(s < e for _, s, e in r),
            f"{len(r):,} rows")
        beyond19 = sum(1 for c, _, e in r if e > lengths["hg19"][c])
        beyond38 = sum(1 for c, _, e in r if e > lengths["hg38"][c])
        add(f"{name}: every interval fits hg19", beyond19 == 0,
            f"{beyond19} end beyond an hg19 chromosome end; {beyond38:,} end beyond the hg38 end of the same "
            f"chromosome, so these cannot be hg38 coordinates")
    for label, name in (("A", "peaks_A.bed"), ("B", "peaks_B.bed")):
        chroms = {c for c, _, _ in rows[name]}
        add(f"{name}: chr1-22, X, Y only", chroms <= set(STANDARD), f"{len(chroms)} chromosomes")
        at_build = recorded[label]["overlapping_blacklist_after_filter"]
        if LOCAL_BLACKLIST.exists():  # recount when the blacklist has been fetched; it is not redistributed here
            now = sum(sql_overlap_flags(rows[name], read_rows(LOCAL_BLACKLIST)))
            assert now == at_build, f"{name}: {now} peaks overlap the blacklist, {at_build} recorded at build"
        add(f"{name}: no peak overlaps the hg19 blacklist v2", at_build == 0,
            f"{at_build} of {len(rows[name]):,} peaks overlap it; counted when the pilot was built and recounted "
            "by this command whenever the blacklist file has been fetched")
    return checks


# --- 2 and 3. recompute and replay -------------------------------------------


def recompute(workspace: Path) -> dict:
    first = run_script(PRIVATE / "reference_solution.py", workspace)
    second = run_script(PRIVATE / "recompute_sqlite.py", workspace)
    return {"reference_solution": first, "sqlite": second["answer"],
            "overlapping_pairs": second["overlapping_pairs"], "agree": first == second["answer"]}


def run_replays(workspace: Path, key: dict) -> list[dict]:
    out = []
    for name, kind, fn in replays.REPLAYS:
        answer = fn(workspace)
        verdict = grader.grade(answer, key)
        out.append({"replay": name, "kind": kind, "says": (fn.__doc__ or "").strip().split("\n")[0],
                    "answer": answer, "passed": verdict["passed"], "schema_errors": verdict["schema_errors"],
                    "wrong_fields": verdict["wrong_fields"], "choice_correct": verdict["choice_correct"],
                    "choice_follows_from_counts": verdict["choice_follows_from_counts"]})
    return out


def describe(key: dict) -> dict:
    return {
        "precision_a": ratio(key["n_a_peaks_overlapping_reference"], key["n_peaks_a"]),
        "precision_b": ratio(key["n_b_peaks_overlapping_reference"], key["n_peaks_b"]),
        "recovered_by_a": ratio(key["n_reference_peaks_overlapped_by_a"], key["n_reference_peaks"]),
        "recovered_by_b": ratio(key["n_reference_peaks_overlapped_by_b"], key["n_reference_peaks"]),
        "a_meets_bar": grader.qualifies(key["n_a_peaks_overlapping_reference"], key["n_peaks_a"]),
        "b_meets_bar": grader.qualifies(key["n_b_peaks_overlapping_reference"], key["n_peaks_b"]),
    }


# --- 4. instances built to test the label swap, the boundary and a tie -------


def build_fixtures(a: list[Row], b: list[Row], key: dict) -> dict:
    """Instances derived from the primary inputs. Their keys follow from how they are built."""
    winner = key["answer"]
    assert winner in ("A", "B"), "these fixtures need a primary instance with a winner"
    loser = "B" if winner == "A" else "A"
    sets = {"A": a, "B": b}
    w, lo = winner.lower(), loser.lower()

    fixtures = {"label_swap": {
        "what": "the primary instance with the 2 labels exchanged",
        "a": b, "b": a, "key": swap_key(key)}}

    tie_key = dict(key)
    for field in COUNT_FIELDS:
        if f"_{lo}_" in field:
            tie_key[field] = key[field.replace(f"_{lo}_", f"_{w}_")]
        elif field.endswith(f"_{lo}"):
            tie_key[field] = key[field[:-1] + w]
    tie_key["answer"] = "tie"
    fixtures["tie_identical_sets"] = {
        "what": "both labels hold the same peak set, so both qualify and recover the same reference peaks",
        "a": sets[winner], "b": sets[winner], "key": tie_key}
    return fixtures


def build_boundary(real: list[tuple], ref: list[Row]) -> dict:
    """2 instances on either side of the 75% line, cut from the real pair that sits closest below it.

    Take the real pair in which one set recovers more reference peaks but misses
    the requirement. Remove peaks that have no reference partner from that set
    until exactly 75% of what is left overlaps the reference: still excluded.
    Remove one more: now it qualifies, and because it recovers more, the answer flips.
    """
    best = None
    for name, _, _, ia, ib, ikey in real:
        if ikey["answer"] not in ("A", "B"):
            continue
        other = "B" if ikey["answer"] == "A" else "A"
        o = other.lower()
        more = ikey[f"n_reference_peaks_overlapped_by_{o}"] > ikey[f"n_reference_peaks_overlapped_by_{ikey['answer'].lower()}"]
        if more:
            precision = ikey[f"n_{o}_peaks_overlapping_reference"] / ikey[f"n_peaks_{o}"]
            if best is None or precision > best[0]:
                best = (precision, name, other, ia, ib, ikey)
    assert best is not None, "no real pair has a set that recovers more but misses the requirement"
    start, name, x_label, ia, ib, ikey = best
    y_label = "B" if x_label == "A" else "A"
    x, y = x_label.lower(), y_label.lower()
    sets = {"A": ia, "B": ib}
    flags = sql_overlap_flags(sets[x_label], ref)
    hits = [r for r, f in zip(sets[x_label], flags) if f]
    misses = [r for r, f in zip(sets[x_label], flags) if not f]
    k = len(hits) - len(hits) % 3
    d = k // 3
    assert 0 < d <= len(misses) and 4 * k == 3 * (k + d)
    fixtures = {}
    for fixture, n_miss in (("boundary_exactly_75", d), ("boundary_one_peak_above_75", d - 1)):
        trimmed = sort_rows(hits[:k] + misses[:n_miss])
        recovered = sum(sql_overlap_flags(ref, trimmed))
        assert recovered > ikey[f"n_reference_peaks_overlapped_by_{y}"], "the trimmed set must still recover more"
        fkey = dict(ikey)
        fkey[f"n_peaks_{x}"] = len(trimmed)
        fkey[f"n_{x}_peaks_overlapping_reference"] = k
        fkey[f"n_reference_peaks_overlapped_by_{x}"] = recovered
        fkey["answer"] = y_label if n_miss == d else x_label
        staged = {x_label: trimmed, y_label: sets[y_label]}
        fixtures[fixture] = {
            "what": (f"a real pair whose higher-recovery set starts at {100 * start:.2f}% precision, with "
                     f"{len(sets[x_label]) - len(trimmed):,} of its peaks removed so that it sits "
                     + ("exactly on 75%" if n_miss == d else "one peak past 75%")),
            "built_from": name, "trimmed_label": x_label,
            "a": staged["A"], "b": staged["B"], "key": fkey}
    return fixtures


# --- the whole run -----------------------------------------------------------


def run() -> dict:
    provenance = json.loads((PRIVATE / "provenance.json").read_text())
    key = json.loads((PRIVATE / "answer_key.json").read_text())["key"]
    a = read_rows(WORKSPACE / "inputs" / "peaks_A.bed")
    b = read_rows(WORKSPACE / "inputs" / "peaks_B.bed")
    ref = read_rows(WORKSPACE / "inputs" / "reference_peaks.bed")
    results: dict = {"task_id": provenance["task_id"], "status": provenance["status"]}
    results["inputs"] = check_inputs(provenance)

    labels = {r["instance"]: r for r in provenance["instances"]}
    real = [("primary", "real", "the pilot instance an agent receives", a, b, key)]
    for folder in sorted((PRIVATE / "siblings").iterdir()):
        sibling_key = json.loads((folder / "key.json").read_text())["key"]
        real.append((folder.name, "real", "another public pair from the same paper, blinded the same way",
                     read_rows(folder / "peaks_A.bed.gz"), read_rows(folder / "peaks_B.bed.gz"), sibling_key))
    built = dict(build_fixtures(a, b, key))
    built.update(build_boundary(real, ref))
    instances = [real[0]]  # (name, kind, what, a rows, b rows, key)
    instances += [(name, "built", fx["what"], fx["a"], fx["b"], fx["key"]) for name, fx in built.items()]
    instances += real[1:]

    per_instance = []
    with tempfile.TemporaryDirectory() as tmp:
        for name, kind, what, ia, ib, ikey in instances:
            workspace = stage(Path(tmp) / name, ia, ib, ref)
            twice = recompute(workspace)
            runs = run_replays(workspace, ikey)
            guesses = {g: fn(workspace) for g, fn in replays.GUESSES}
            per_instance.append({
                "instance": name, "kind": kind, "what": what, "key": ikey, "derived": describe(ikey),
                "recomputation": twice,
                "both_recomputations_equal_the_key": twice["agree"] and twice["reference_solution"] == ikey,
                "replays": runs,
                "choice_only_guesses": {g: {"guess": v, "correct": v == ikey["answer"]} for g, v in guesses.items()},
            })
    results["instances"] = per_instance
    by_name = {p["instance"]: p for p in per_instance}

    primary = by_name["primary"]
    results["paper_notebook_check"] = provenance.get("paper_notebook_check")
    results["private_labels_for_the_primary"] = {l: v["peak_caller"] for l, v in labels["primary"]["labels"].items()}

    # Label swap: the twin's key is the primary key mirrored, and the primary's correct answer must fail on it.
    correct_primary = next(r for r in primary["replays"] if r["replay"] == "correct")["answer"]
    results["label_swap"] = {
        "primary_choice": key["answer"], "twin_choice": by_name["label_swap"]["key"]["answer"],
        "twin_solved_correctly": by_name["label_swap"]["both_recomputations_equal_the_key"],
        "primary_answer_graded_against_twin_key": grader.grade(correct_primary, by_name["label_swap"]["key"])["passed"],
    }
    exact, above = by_name["boundary_exactly_75"], by_name["boundary_one_peak_above_75"]
    trimmed = built["boundary_exactly_75"]["trimmed_label"].lower()
    results["boundary"] = {
        "built_from": built["boundary_exactly_75"]["built_from"],
        "exactly_75": {"key": exact["key"], "derived": exact["derived"]},
        "one_peak_above": {"key": above["key"], "derived": above["derived"]},
        "precision_exact_percent": f"{100 * exact['key'][f'n_{trimmed}_peaks_overlapping_reference'] / exact['key'][f'n_peaks_{trimmed}']:.4f}",
        "precision_above_percent": f"{100 * above['key'][f'n_{trimmed}_peaks_overlapping_reference'] / above['key'][f'n_peaks_{trimmed}']:.4f}",
        "the_two_instances_differ_by_rows": abs(
            (exact["key"]["n_peaks_a"] + exact["key"]["n_peaks_b"]) - (above["key"]["n_peaks_a"] + above["key"]["n_peaks_b"])),
        "choice_flips": exact["key"]["answer"] != above["key"]["answer"],
    }

    # Where each replay passes and fails.
    survival = {}
    for name, kind, _ in replays.REPLAYS:
        outcome = {p["instance"]: next(r for r in p["replays"] if r["replay"] == name)["passed"] for p in per_instance}
        survival[name] = {"kind": kind,
                          "passes_on": [i for i, ok in outcome.items() if ok],
                          "fails_on": [i for i, ok in outcome.items() if not ok]}
    results["survival"] = survival

    real_results = [p for p in per_instance if p["kind"] == "real"]
    results["choice_only"] = {
        "note": "what would pass if only the final choice were graded, on the real pairs",
        "answers": {v: sum(p["key"]["answer"] == v for p in real_results) for v in ("A", "B", "tie", "neither")},
        "guesses": {g: sum(p["choice_only_guesses"][g]["correct"] for p in real_results) for g, _ in replays.GUESSES},
        "wrong_logic": {name: sum(next(r for r in p["replays"] if r["replay"] == name)["choice_correct"] for p in real_results)
                        for name, kind, _ in replays.REPLAYS if kind == "wrong_logic"},
        "n_real_pairs": len(real_results),
    }

    # How well the blinding holds: can the peak caller be told from the intervals alone?
    wider_caller = []
    for name, _, _, ia, ib, _ in real:
        med = {"A": statistics.median(e - s for _, s, e in ia), "B": statistics.median(e - s for _, s, e in ib)}
        wider = "A" if med["A"] > med["B"] else "B"
        wider_caller.append(labels[name]["labels"][wider]["peak_caller"])
    by_caller: dict[str, int] = {}
    for rec in provenance["instances"]:
        ans = rec["key"]["answer"]
        who = rec["labels"][ans]["peak_caller"] if ans in ("A", "B") else ans
        by_caller[who] = by_caller.get(who, 0) + 1
    results["blinding"] = {
        "pairs": len(real),
        "pairs_where_the_wider_set_is_the_same_caller": max(wider_caller.count(c) for c in set(wider_caller)),
        "rule_selects_by_caller": dict(sorted(by_caller.items(), key=lambda kv: (-kv[1], kv[0]))),
    }

    # Agent runs on record, regraded here.
    runs = []
    for path in sorted((PRIVATE / "agent_runs").glob("*/run.json")):
        run_record = json.loads(path.read_text())
        regraded = grader.grade(run_record["answer"], key)["passed"] if run_record["answer"] is not None else False
        calls = run_record.get("tool_calls", {})
        runs.append({"run_id": run_record["run_id"], "model": run_record["model"], "passed": regraded,
                     "recorded_verdict_reproduces": regraded == run_record["verdict"]["passed"],
                     "tool_calls": calls.get("total"),
                     "calls_naming_a_path_outside_the_workspace": len(calls.get("naming_a_path_outside_the_workspace", []))})
    results["agent_runs"] = runs

    # The checks that decide the exit code.
    n_instances = len(per_instance)
    verdicts = {
        "inputs_verified": all(c["ok"] for c in results["inputs"]),
        "two_recomputations_equal_the_key_on_every_instance": all(p["both_recomputations_equal_the_key"] for p in per_instance),
        "correct_solution_passes_everywhere": len(survival["correct"]["passes_on"]) == n_instances,
        "wrong_counts_and_malformed_fail_everywhere": all(
            not s["passes_on"] for s in survival.values() if s["kind"] in ("wrong_counts", "malformed")),
        "every_logic_shortcut_is_caught_by_at_least_one_instance": all(
            s["fails_on"] for s in survival.values() if s["kind"] in ("wrong_logic", "convention")),
        "label_swap_flips_the_choice": results["label_swap"]["primary_choice"] != results["label_swap"]["twin_choice"]
        and results["label_swap"]["twin_solved_correctly"]
        and not results["label_swap"]["primary_answer_graded_against_twin_key"],
        "one_peak_flips_the_boundary_choice": results["boundary"]["choice_flips"]
        and results["boundary"]["the_two_instances_differ_by_rows"] == 1,
        "recorded_agent_verdicts_reproduce": all(r["recorded_verdict_reproduces"] for r in runs),
    }
    nb = results["paper_notebook_check"]
    if nb:
        verdicts["key_equals_the_paper_notebook"] = all(
            v["paper_notebook"] == v["recomputed"] for caller in nb.values() for v in caller.values() if isinstance(v, dict))
    results["verdicts"] = verdicts
    results["survives_on_the_primary"] = [n for n, s in survival.items()
                                          if n != "correct" and "primary" in s["passes_on"]]
    return results


# --- the public summary: everything above with the answers taken out ----------


def summarise(results: dict) -> dict:
    """A record that can be published: verdicts and counts of agreement, no key value and no label mapping."""
    by_name = {p["instance"]: p for p in results["instances"]}
    primary = by_name["primary"]
    real = [p for p in results["instances"] if p["kind"] == "real"]
    nb = results.get("paper_notebook_check") or {}
    nb_values = [v for caller in nb.values() for v in caller.values() if isinstance(v, dict)]
    twice = primary["recomputation"]
    fields = list(COUNT_FIELDS) + ["answer"]

    def naive_wrong(p: dict) -> bool:
        return next(r for r in p["replays"] if r["replay"] == "ignore_precision_bar")["answer"]["answer"] != p["key"]["answer"]

    return {
        "task_id": results["task_id"],
        "status": results["status"],
        "verdicts": results["verdicts"],
        "inputs": results["inputs"],
        "recomputation": {
            "fields": len(fields),
            "fields_equal_between_the_two_recomputations": sum(twice["reference_solution"][f] == twice["sqlite"][f] for f in fields),
            "fields_equal_to_the_key": sum(twice["reference_solution"][f] == primary["key"][f] for f in fields),
            "paper_notebook_numbers": len(nb_values),
            "paper_notebook_numbers_equal": sum(v["paper_notebook"] == v["recomputed"] for v in nb_values),
        },
        "primary_replays": [
            {"replay": r["replay"], "kind": r["kind"], "says": r["says"], "passed": r["passed"],
             "counts_wrong": len(r["wrong_fields"]), "choice_correct": r["choice_correct"],
             "valid_answer": not r["schema_errors"]}
            for r in primary["replays"]],
        "survives_on_the_primary": results["survives_on_the_primary"],
        "label_swap": {"choice_changes": results["label_swap"]["primary_choice"] != results["label_swap"]["twin_choice"],
                       "primary_answer_fails_on_the_twin": not results["label_swap"]["primary_answer_graded_against_twin_key"]},
        "boundary": {"what": by_name["boundary_exactly_75"]["what"],
                     "precision_exact_percent": results["boundary"]["precision_exact_percent"],
                     "precision_above_percent": results["boundary"]["precision_above_percent"],
                     "differ_by_rows": results["boundary"]["the_two_instances_differ_by_rows"],
                     "choice_flips": results["boundary"]["choice_flips"]},
        "tie": {"what": by_name["tie_identical_sets"]["what"], "key_is_tie": by_name["tie_identical_sets"]["key"]["answer"] == "tie"},
        "instances": [{"instance": p["instance"], "kind": p["kind"], "what": p["what"]} for p in results["instances"]],
        "survival": results["survival"],
        "real_pairs": {
            "n": len(real),
            "answers_by_label": {k: v for k, v in results["choice_only"]["answers"].items() if v},
            "answers_by_caller": results["blinding"]["rule_selects_by_caller"],
            "other_pairs": len(real) - 1,
            "other_pairs_where_ignoring_the_bar_is_wrong": sum(naive_wrong(p) for p in real if p["instance"] != "primary"),
        },
        "choice_only": {"guesses": results["choice_only"]["guesses"], "wrong_logic": results["choice_only"]["wrong_logic"]},
        "blinding": {"pairs": results["blinding"]["pairs"],
                     "pairs_where_the_wider_set_is_the_same_caller": results["blinding"]["pairs_where_the_wider_set_is_the_same_caller"]},
        "agent_runs": results["agent_runs"],
    }


def public_report(s: dict) -> str:
    out = ["# Replay results", "",
           f"Task `{s['task_id']}`. Status: {s['status']}.",
           "Written by `python3 pilot/replay.py --write`. The replays are scripted; the checks were written and run "
           "by AI models, and no person has checked them yet. Do not edit by hand.", "",
           "No answer is printed on this page. The full record, answers included, is `private/results/REPLAY_FULL.md`.", "",
           "## Verdicts", "", "| Check | Result |", "|---|---|"]
    out += [f"| {name.replace('_', ' ')} | {'yes' if ok else 'NO'} |" for name, ok in s["verdicts"].items()]
    out += ["", "## 1. The staged inputs", "", "| Check | Result | Detail |", "|---|---|---|"]
    out += [f"| {c['check']} | {'ok' if c['ok'] else 'FAILED'} | {c['detail']} |" for c in s["inputs"]]
    r = s["recomputation"]
    out += ["", "## 2. The answer, recomputed 2 ways", "",
            f"A binary search in Python and a range join in SQLite agree on {r['fields_equal_between_the_two_recomputations']} "
            f"of {r['fields']} answer fields, and {r['fields_equal_to_the_key']} of {r['fields']} equal the key.",
            f"{r['paper_notebook_numbers_equal']} of {r['paper_notebook_numbers']} numbers also equal what the paper's own "
            "rendered notebook prints (computed there in R with GenomicRanges). Neither recomputation reads that notebook."]
    out += ["", "## 3. Scripted replays on the primary instance", "",
            "| Replay | Kind | What it does | Counts wrong, of 7 | Choice | Grader |", "|---|---|---|---|---|---|"]
    for p in s["primary_replays"]:
        counts = str(p["counts_wrong"]) if p["valid_answer"] else "not a valid answer"
        out.append(f"| {p['replay']} | {p['kind']} | {p['says']} | {counts} | "
                   f"{'right' if p['choice_correct'] else 'wrong'} | {'PASS' if p['passed'] else 'fail'} |")
    surv = s["survives_on_the_primary"]
    out += ["", ("Shortcuts that survive on the primary instance: " + ", ".join(f"`{x}`" for x in surv) + "."
                 if surv else "No shortcut survives on the primary instance.")]
    b = s["boundary"]
    out += ["", "## 4. Label swap, the 75% line, and a tie", "",
            f"- **Label swap.** With the 2 files exchanged the correct choice "
            f"{'changes' if s['label_swap']['choice_changes'] else 'DOES NOT change'}, and the primary's correct answer "
            f"{'fails' if s['label_swap']['primary_answer_fails_on_the_twin'] else 'PASSES'} against the swapped key.",
            f"- **Exactly 75%.** Built from {b['what']}. Precision {b['precision_exact_percent']}%: the set is still excluded.",
            f"- **One peak past 75%.** The same set with {b['differ_by_rows']} more unsupported peak removed. Precision "
            f"{b['precision_above_percent']}%: it qualifies, and the correct choice "
            f"{'flips' if b['choice_flips'] else 'DOES NOT flip'}.",
            f"- **Tie.** {s['tie']['what'][0].upper() + s['tie']['what'][1:]}. Key is tie: {'yes' if s['tie']['key_is_tie'] else 'NO'}."]
    n = len(s["instances"])
    out += ["", "## 5. Where each replay passes and fails", "",
            f"{n} instances: the primary; 2 built from it (label swap, tie); 2 cut at the 75% line from a real pair; "
            f"and the paper's {s['real_pairs']['other_pairs']} other public pairs.", "",
            "| Replay | Kind | Passes on | Fails on |", "|---|---|---|---|"]
    for name, v in s["survival"].items():
        out.append(f"| {name} | {v['kind']} | {len(v['passes_on'])}: {', '.join(v['passes_on']) or 'none'} | "
                   f"{len(v['fails_on'])}: {', '.join(v['fails_on']) or 'none'} |")
    rp = s["real_pairs"]
    out += ["", "## 6. The real pairs", "",
            f"On the {rp['n']} real pairs the rule selects " + ", ".join(f"{k} {v}" for k, v in rp["answers_by_label"].items())
            + " by label, and " + ", ".join(f"{k} in {v}" for k, v in rp["answers_by_caller"].items()) + " by peak caller.",
            f"Ignoring the 75% requirement gives the wrong answer on {rp['other_pairs_where_ignoring_the_bar_is_wrong']} of the "
            f"{rp['other_pairs']} pairs other than the primary."]
    out += ["", "## 7. If only the final choice were graded", "",
            "| Shortcut | Right choice on | Needs to count overlaps |", "|---|---|---|"]
    out += [f"| {g} | {c} of {rp['n']} | no |" for g, c in s["choice_only"]["guesses"].items()]
    out += [f"| {g} | {c} of {rp['n']} | yes |" for g, c in s["choice_only"]["wrong_logic"].items()]
    out += ["", "The full grader also checks the 7 counts, so none of the guesses in the first group can pass it."]
    bl = s["blinding"]
    out += ["", "## 8. How far the blinding goes", "",
            f"The workspace carries no names. It does not hide what the 2 sets are: on "
            f"{bl['pairs_where_the_wider_set_is_the_same_caller']} of {bl['pairs']} pairs the set with the wider peaks comes "
            "from the same peak caller."]
    out += ["", "## 9. Agent runs", ""]
    if s["agent_runs"]:
        out += ["| Run | Model | Result | Tool calls | Calls naming a path outside the workspace |", "|---|---|---|---|---|"]
        out += [f"| {a['run_id']} | {a['model']} | {'PASS' if a['passed'] else 'fail'} | {a['tool_calls']} | "
                f"{a['calls_naming_a_path_outside_the_workspace']} |" for a in s["agent_runs"]]
        out += ["", "Hands-on runs on the author's machine, where the key is on disk. Not benchmark-grade."]
    else:
        out += ["None recorded."]
    return "\n".join(out) + "\n"


# --- the full record, kept with the held-out material -------------------------


def full_report(results: dict) -> str:
    by_name = {p["instance"]: p for p in results["instances"]}
    primary = by_name["primary"]
    key = primary["key"]
    twice = primary["recomputation"]
    out = ["# Replay results, full record", "",
           "Held-out material: this page states the answer key, the keys of the check instances and the scripted answers.",
           "Written by `python3 pilot/replay.py --write`. Do not edit by hand.", "",
           "## The primary key, recomputed 2 ways", "",
           "| Field | Binary search (Python) | Range join (SQLite) | Answer key |", "|---|---|---|---|"]
    for field in list(COUNT_FIELDS) + ["answer"]:
        out.append(f"| {field} | {twice['reference_solution'][field]} | {twice['sqlite'][field]} | {key[field]} |")
    d = primary["derived"]
    out += ["", f"Set A: {d['precision_a']} of its peaks overlap the reference; it recovers {d['recovered_by_a']} reference peaks.",
            f"Set B: {d['precision_b']} of its peaks overlap the reference; it recovers {d['recovered_by_b']} reference peaks.",
            f"Private labels: {results['private_labels_for_the_primary']}."]
    nb = results.get("paper_notebook_check")
    if nb:
        out += ["", "| Peak caller | Quantity | Paper's notebook | Recomputed here |", "|---|---|---|---|"]
        for caller, values in nb.items():
            for name, v in values.items():
                if isinstance(v, dict):
                    out.append(f"| {caller} | {name.replace('_', ' ')} | {v['paper_notebook']:,} | {v['recomputed']:,} |")
    out += ["", "## Scripted answers on the primary instance", "",
            "| Replay | Its choice | Grader | Wrong fields |", "|---|---|---|---|"]
    for r in primary["replays"]:
        if r["schema_errors"]:
            wrong = f"not a valid answer: {r['schema_errors'][0]}"
        else:
            wrong = "; ".join([f"{f}={r['answer'][f]:,}" for f in r["wrong_fields"]]
                              + ([] if r["choice_correct"] else ["the choice"])) or "none"
        out.append(f"| {r['replay']} | {r['answer'].get('answer')} | {'PASS' if r['passed'] else 'fail'} | {wrong} |")
    out += ["", "## Every instance", "",
            "| Instance | Precision A | Precision B | Recovered by A | Recovered by B | Key | Ignoring the bar picks |",
            "|---|---|---|---|---|---|---|"]
    for p in results["instances"]:
        k = p["key"]
        naive = next(r for r in p["replays"] if r["replay"] == "ignore_precision_bar")["answer"]["answer"]
        out.append(f"| {p['instance']} | {100 * k['n_a_peaks_overlapping_reference'] / k['n_peaks_a']:.4f}% | "
                   f"{100 * k['n_b_peaks_overlapping_reference'] / k['n_peaks_b']:.4f}% | "
                   f"{k['n_reference_peaks_overlapped_by_a']:,} | {k['n_reference_peaks_overlapped_by_b']:,} | "
                   f"{k['answer']} | {naive}{'' if naive == k['answer'] else ' (wrong)'} |")
    return "\n".join(out) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Replay the pilot without a model run.")
    parser.add_argument("--write", action="store_true", help="record the results instead of comparing with them")
    args = parser.parse_args()
    results = run()
    summary = summarise(results)

    import docs
    provenance = json.loads((PRIVATE / "provenance.json").read_text())
    readme = (PILOT / "README.md").read_text()
    wanted = {
        FULL / "replay.json": json.dumps(results, indent=1) + "\n",
        FULL / "REPLAY_FULL.md": full_report(results),
        RESULTS / "summary.json": json.dumps(summary, indent=1) + "\n",
        RESULTS / "REPLAY.md": public_report(summary),
        PILOT / "VERIFICATION.md": docs.verification(provenance, summary),
        PILOT / "README.md": docs.with_block(readme, docs.readme_block(summary)),
    }
    if args.write:
        for path, text in wanted.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    stale = [str(path.relative_to(PILOT)) for path, text in wanted.items() if not path.exists() or path.read_text() != text]

    for name, ok in results["verdicts"].items():
        print(f"{'ok  ' if ok else 'FAIL'}  {name.replace('_', ' ')}")
    survivors = results["survives_on_the_primary"]
    print("survives on the primary instance:", ", ".join(survivors) if survivors else "nothing")
    print("matches the recorded results and documents:", "yes" if not stale else "NO: " + ", ".join(stale))
    good = all(results["verdicts"].values()) and not stale
    print("PILOT REPLAY:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
