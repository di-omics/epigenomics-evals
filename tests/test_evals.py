"""The published task files: format, internal consistency, house rules.

These tests read evals/ and results/ only, so they run without the data.
"""
import json
import re
from pathlib import Path

import pytest

from epigenomics_evals import grading
from epigenomics_evals.specs import SPECS

REPO = Path(__file__).resolve().parents[1]
EVALS = sorted((REPO / "evals").glob("*/eval.json"))
LATCH_GRADER_TYPES = {"numeric_tolerance", "all_of"}


def load(path):
    return json.loads(path.read_text())


def test_manifest_lists_every_task():
    manifest = load(REPO / "evals" / "manifest.json")
    assert [e["id"] for e in manifest["evals"]] == list(SPECS)
    assert len(EVALS) == len(SPECS) == 6


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_eval_has_latch_keys_and_a_checkable_answer(path):
    ev = load(path)
    for key in ("id", "canary", "task", "data_node", "grader", "metadata", "notes"):
        assert key in ev
    assert ev["id"] == path.parent.name
    assert "canary GUID" in ev["canary"]
    assert "<EVAL_ANSWER>" in ev["task"] and "</EVAL_ANSWER>" in ev["task"]
    assert ev["metadata"]["eval_type"] == "scientific"
    assert ev["metadata"]["failure_modes"]
    assert ev["grader"]["type"] in LATCH_GRADER_TYPES
    block = re.search(r"<EVAL_ANSWER>\s*(\{.*?\})\s*</EVAL_ANSWER>", ev["task"], re.S).group(1)
    promised = set(re.findall(r'"([a-z0-9_]+)":', block))
    graded = set(grading.windows(ev["grader"]))
    assert promised == graded, "every field the prompt asks for is graded, and no other"


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_windows_are_well_formed_and_contain_the_key(path):
    ev = load(path)
    children = ev["grader"]["config"]["children"] if ev["grader"]["type"] == "all_of" else [ev["grader"]]
    truth = {k: v for child in children for k, v in child["config"]["ground_truth"].items()}
    for name, (lo, hi) in grading.windows(ev["grader"]).items():
        assert lo <= truth[name] <= hi
    assert grading.own_check(truth, ev["grader"])["passed"]


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_choice_key_points_at_a_listed_option(path):
    ev = load(path)
    spec = SPECS[ev["id"]]
    if not spec.choices:
        return
    options = re.findall(r"^(\d+)\. ", ev["task"], re.M)
    assert [int(o) for o in options] == list(range(1, len(spec.choices) + 1))
    assert 1 <= spec.correct_choice <= len(spec.choices)


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_notes_name_the_trap_and_the_measurements(path):
    notes = load(path)["notes"]
    for section in ("1) Scientific intent", "2) Primary trap", "3) Accepted analyses",
                    "4) Scripted wrong analyses", "5) Tolerance", "6) Data",
                    "7) What the prompt does not say", "8) Source"):
        assert section in notes
    assert "outside by" in notes or "the key is" in notes


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_loads_in_latch_eval_tools(path):
    types = pytest.importorskip("latch_eval_tools.types")
    ev = types.Eval(**load(path))
    assert ev.id == path.parent.name and ev.grader is not None


def test_results_table_matches_the_rule():
    results = load(REPO / "results" / "controls.json")
    controls = results["controls"]
    assert len(controls) >= 40
    for c in controls:
        assert c["own"]["passed"] == c["expect_pass"], (c["task"], c["name"])
        if c["latch"] is not None:
            assert c["latch"]["passed"] == c["expect_pass"], (c["task"], c["name"])
    for task_id in SPECS:
        kinds = {c["kind"] for c in controls if c["task"] == task_id}
        assert kinds == {"reference", "variant", "trap"}, "each task has all three kinds"


def test_readme_numbers_come_from_the_results():
    results = load(REPO / "results" / "controls.json")["controls"]
    readme = (REPO / "README.md").read_text()
    n = len(results)
    n_pass = sum(c["expect_pass"] for c in results)
    assert f"{n} scripted analyses ({n_pass} that should pass, {n - n_pass} that should fail)" in readme
    assert f"{n} of {n} behave as expected" in readme


def test_tracked_text_is_ascii_without_em_dashes():
    paths = [p for pattern in ("*.md", "*.toml", "docs/*.md", "docs/*.svg", ".github/workflows/*", "results/*", "evals/**/*",
                               "src/**/*.py", "scripts/*", "tests/*.py", "pilot/*.py", "pilot/*.md",
                               "pilot/*.json", "pilot/workspace/*.md", "pilot/workspace/*.json",
                               "pilot/private/*.py", "pilot/private/*.json", "pilot/private/*.md", "pilot/private/results/*",
                               "pilot/results/*")
             for p in REPO.glob(pattern) if p.is_file()]
    assert paths
    for p in paths:
        data = p.read_bytes()
        assert all(b < 128 for b in data), f"non-ASCII byte in {p}"


@pytest.mark.parametrize("path", EVALS, ids=lambda p: p.parent.name)
def test_keyed_option_does_not_stand_out_and_status_is_stated(path):
    ev = load(path)
    spec = SPECS[ev["id"]]
    assert ev["metadata"]["kit"] == "hts"
    assert ev["metadata"]["key_review"] == spec.key_review
    assert "10) Review status" in ev["notes"] and "{" not in ev["notes"], "no unfilled figure in the notes"
    if spec.choices:
        lengths = [len(c) for c in spec.choices]
        assert min(lengths) < lengths[spec.correct_choice - 1] < max(lengths)
        assert "Not yet reviewed by a person" in ev["metadata"]["key_review"]


def test_readme_says_who_checked_and_what_is_pending():
    readme = (REPO / "README.md").read_text()
    assert "No person has" in readme and "AI models" in readme
    assert "Not affiliated with LatchBio" in readme
    for task_id, spec in SPECS.items():
        if spec.choices:
            assert f"`{task_id}`" in readme, "tasks whose decision key awaits review are named"
    assert (REPO / "THIRD_PARTY_NOTICES.md").exists()
