"""The Figure 4a pilot: replay, grader edges, independent recomputation, leak checks.

These tests need nothing but the repository: the pilot ships its inputs.
"""
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PILOT = REPO / "pilot"
WORKSPACE = PILOT / "workspace"
PRIVATE = PILOT / "private"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


grader = load("pilot_grader", PILOT / "grader.py")
KEY = json.loads((PRIVATE / "answer_key.json").read_text())["key"]
RESULTS = json.loads((PRIVATE / "results" / "replay.json").read_text())
SUMMARY = json.loads((PILOT / "results" / "summary.json").read_text())
SECRET = {KEY[f] for f in grader.COUNT_FIELDS if "overlap" in f}
LABELS = re.compile(
    r"MACS|SEACR|Diagenode|C15410196|Abcam|ab177|ab4729|Active.?Motif|GSM\d|GSE\d|ENCFF|ENCSR|rmDup|10\.1038|Marzi|Abbasova",
    re.I,
)


def public_text_files() -> list[Path]:
    """Every text file of the repository outside pilot/private (tracked files when this is a git checkout)."""
    listed = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                            cwd=REPO, capture_output=True, text=True)
    names = [n for n in listed.stdout.split("\n") if n] if listed.returncode == 0 and listed.stdout else [
        str(p.relative_to(REPO)) for p in REPO.rglob("*")
        if p.is_file() and not {".venv", "data", ".git", "__pycache__", ".pytest_cache"} & set(p.relative_to(REPO).parts)]
    out = []
    for name in names:
        path = REPO / name
        if name.startswith("pilot/private/") or path.suffix in (".gz", ".bed") or not path.is_file():
            continue
        try:
            path.read_text()
        except UnicodeDecodeError:
            continue
        out.append(path)
    return out


def test_one_command_replay_passes_and_reproduces_the_recorded_results():
    done = subprocess.run([sys.executable, "-I", str(PILOT / "replay.py")], capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    assert "PILOT REPLAY: PASS" in done.stdout
    assert "matches the recorded results and documents: yes" in done.stdout


def test_workspace_holds_only_what_an_agent_should_see():
    files = {str(p.relative_to(WORKSPACE)) for p in WORKSPACE.rglob("*") if p.is_file()}
    assert files == {"task.md", "manifest.json", "answer_schema.json", "inputs/peaks_A.bed",
                     "inputs/peaks_B.bed", "inputs/reference_peaks.bed"}
    for rel in files:
        assert not LABELS.search((WORKSPACE / rel).read_text()), f"label leak in {rel}"
    for rel in ("task.md", "manifest.json", "answer_schema.json"):
        text = (WORKSPACE / rel).read_text()
        assert not any(str(value) in text for value in SECRET), f"key value in {rel}"


def test_no_page_outside_the_private_folder_states_the_key_or_the_label_mapping():
    counts = [re.compile(r"(?<![\d,.])" + re.escape(form) + r"(?![\d,])")
              for value in SECRET for form in (str(value), f"{value:,}")]
    mapping = re.compile(r"(set [AB]\b|peaks_[AB]\b).{0,80}(SEACR|MACS2)|(SEACR|MACS2).{0,80}(set [AB]\b|peaks_[AB]\b)", re.I)
    files = public_text_files()
    assert len(files) > 40
    for path in files:
        in_pilot = path.is_relative_to(PILOT)
        for number, line in enumerate(path.read_text().splitlines(), 1):
            assert not any(c.search(line) for c in counts), f"key count at {path.relative_to(REPO)}:{number}"
            assert not (in_pilot and mapping.search(line)), f"label mapping at {path.relative_to(REPO)}:{number}"
    published = json.dumps(SUMMARY)
    assert not any(str(value) in published for value in SECRET)
    assert "private_labels" not in published


def test_public_eval_entry_carries_no_key():
    entry = json.loads((PILOT / "eval.json").read_text())
    assert set(entry) == {"id", "canary", "task", "data_node", "data_availability", "metadata"}
    assert entry["task"] == (WORKSPACE / "task.md").read_text()
    assert entry["metadata"]["status"] == "independent prototype"
    assert "<EVAL_ANSWER>" in entry["task"] and not LABELS.search(entry["task"])


def test_schema_file_and_grader_ask_for_the_same_fields():
    schema = json.loads((WORKSPACE / "answer_schema.json").read_text())
    assert schema["required"] == list(grader.ALL_FIELDS)
    assert set(schema["properties"]) == set(grader.ALL_FIELDS)
    assert schema["additionalProperties"] is False
    assert schema["properties"]["answer"]["enum"] == ["A", "B", "tie", "neither"]


@pytest.mark.parametrize("overlapping,total,expected", [
    (3, 4, False), (75, 100, False), (76, 100, True), (301, 400, True), (0, 0, False), (0, 10, False),
    (10, 10, True), (3 * 10**12, 4 * 10**12, False), (3 * 10**12 + 1, 4 * 10**12, True),
])
def test_the_75_percent_requirement_is_strict_and_exact(overlapping, total, expected):
    assert grader.qualifies(overlapping, total) is expected


def counts(n_a, hit_a, got_a, n_b, hit_b, got_b, n_ref=1000):
    return {"n_peaks_a": n_a, "n_peaks_b": n_b, "n_reference_peaks": n_ref,
            "n_a_peaks_overlapping_reference": hit_a, "n_b_peaks_overlapping_reference": hit_b,
            "n_reference_peaks_overlapped_by_a": got_a, "n_reference_peaks_overlapped_by_b": got_b}


@pytest.mark.parametrize("c,expected", [
    (counts(100, 90, 500, 100, 80, 400), "A"),       # both qualify, A recovers more
    (counts(100, 90, 400, 100, 80, 500), "B"),       # both qualify, B recovers more
    (counts(100, 90, 400, 100, 80, 400), "tie"),     # both qualify, same recovery
    (counts(100, 75, 900, 100, 80, 400), "B"),       # A recovers more but sits exactly on 75%
    (counts(100, 76, 900, 100, 80, 400), "A"),       # one peak more and A qualifies
    (counts(100, 70, 900, 100, 60, 400), "neither"),
    (counts(100, 70, 900, 100, 80, 1), "B"),         # only B qualifies
])
def test_the_rule(c, expected):
    assert grader.implied_choice(c) == expected


def test_the_key_passes_and_any_single_change_fails():
    assert grader.grade(dict(KEY), KEY)["passed"]
    for field in grader.COUNT_FIELDS:
        for delta in (-1, 1):
            wrong = dict(KEY)
            wrong[field] += delta
            verdict = grader.grade(wrong, KEY)
            assert not verdict["passed"] and verdict["wrong_fields"] == [field]
    for choice in ("A", "B", "tie", "neither"):
        if choice != KEY["answer"]:
            assert not grader.grade({**KEY, "answer": choice}, KEY)["passed"]


def test_grader_accepts_harmless_formatting_and_rejects_malformed_answers():
    assert grader.grade({**KEY, "answer": KEY["answer"].lower() + " "}, KEY)["passed"]
    assert grader.grade({k: (float(v) if isinstance(v, int) else v) for k, v in KEY.items()}, KEY)["passed"]
    text = "Here is my answer.\n<EVAL_ANSWER>\n" + json.dumps(KEY) + "\n</EVAL_ANSWER>\nDone."
    assert grader.grade(grader.extract_answer(text), KEY)["passed"]
    for bad in (
        {k: v for k, v in KEY.items() if k != "n_reference_peaks"},        # a field is missing
        {**KEY, "precision_a": 0.92},                                       # an extra field
        {**KEY, "n_peaks_a": str(KEY["n_peaks_a"])},                        # a count as text
        {**KEY, "n_peaks_a": True},
        {**KEY, "n_peaks_a": KEY["n_peaks_a"] + 0.5},
        {**KEY, "n_peaks_a": -1},
        {**KEY, "answer": "C"},
    ):
        verdict = grader.grade(bad, KEY)
        assert not verdict["passed"] and verdict["schema_errors"]


def test_counts_follow_from_the_definition_itself():
    """A third computation that uses no trick: compare every pair of intervals on each chromosome."""
    np = pytest.importorskip("numpy")

    def read(name):
        by_chrom = {}
        for line in (WORKSPACE / "inputs" / name).read_text().splitlines():
            chrom, start, end = line.split("\t")
            by_chrom.setdefault(chrom, []).append((int(start), int(end)))
        return {c: np.array(v) for c, v in by_chrom.items()}

    ref = read("reference_peaks.bed")
    for label in ("a", "b"):
        peaks = read(f"peaks_{label.upper()}.bed")
        hit = got = 0
        for chrom, p in peaks.items():
            r = ref[chrom]
            overlap = (p[:, None, 0] < r[None, :, 1]) & (p[:, None, 1] > r[None, :, 0])
            hit += int(overlap.any(axis=1).sum())
            got += int(overlap.any(axis=0).sum())
        assert hit == KEY[f"n_{label}_peaks_overlapping_reference"]
        assert got == KEY[f"n_reference_peaks_overlapped_by_{label}"]


def test_latch_grader_agrees_on_every_recorded_replay():
    graders = pytest.importorskip("latch_eval_tools.graders")
    all_of = graders.get_grader("all_of")
    compared = 0
    for instance in RESULTS["instances"]:
        key = instance["key"]
        config = {"children": [
            {"type": "numeric_tolerance",
             "config": {"ground_truth": {f: key[f] for f in grader.COUNT_FIELDS},
                        "tolerances": {f: {"type": "absolute", "value": 0} for f in grader.COUNT_FIELDS}}},
            {"type": "multiple_choice", "config": {"correct_answer": key["answer"]}},
        ]}
        for replay in instance["replays"]:
            theirs = all_of.evaluate(replay["answer"], config).passed
            assert theirs == replay["passed"], (instance["instance"], replay["replay"])
            compared += 1
    assert compared == len(RESULTS["instances"]) * 10
    shipped = json.loads((PRIVATE / "grader_config.latch.json").read_text())
    assert all_of.evaluate(dict(KEY), shipped["config"]).passed


def test_public_documents_are_the_ones_the_results_generate():
    docs = load("pilot_docs", PILOT / "docs.py")
    prov = json.loads((PRIVATE / "provenance.json").read_text())
    readme = (PILOT / "README.md").read_text()
    assert (PILOT / "VERIFICATION.md").read_text() == docs.verification(prov, SUMMARY)
    assert docs.readme_block(SUMMARY) in readme
    assert "independent prototype" in readme.lower()
    assert "AI models" in readme and "No person has checked" in readme


def test_recorded_verdicts_all_hold():
    assert all(RESULTS["verdicts"].values()) and RESULTS["verdicts"] == SUMMARY["verdicts"]
    by_name = {p["instance"]: p for p in RESULTS["instances"]}
    assert len(by_name) == 15
    assert all(p["both_recomputations_equal_the_key"] for p in by_name.values())
    answers = {p["key"]["answer"] for p in by_name.values()}
    assert answers == {"A", "B", "tie", "neither"}, "every possible answer occurs in the suite"


def test_recorded_agent_runs_regrade_to_the_recorded_verdict():
    for path in sorted((PRIVATE / "agent_runs").glob("*/run.json")):
        run = json.loads(path.read_text())
        assert run["inputs_unmodified_after_the_run"] and run["prompt_unmodified_after_the_run"]
        if run["answer"] is not None:
            assert grader.grade(run["answer"], KEY)["passed"] == run["verdict"]["passed"]
        text = (path.parent / "final_message.md").read_text()
        assert not LABELS.search(text), "the agent's reply names a paper label it was never given"


def test_sources_are_credited_and_the_gpl_blacklist_is_not_shipped():
    notices = (REPO / "THIRD_PARTY_NOTICES.md").read_text()
    for needle in ("Copyright (c) 2022 leylaabbasova", "ENCFF044JNJ", "GSE286492", "Amemiya", "GPL-3.0"):
        assert needle in notices
    shipped = [p for p in PILOT.rglob("*blacklist*")]
    assert not shipped, "the hg19 blacklist is GPL-3.0 and must not be redistributed here"
