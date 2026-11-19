"""Recompute every control from the data and hold the tasks to the rule.

Skipped when the input data has not been fetched (scripts/fetch_data.py).
"""
import json
from pathlib import Path

import pytest

from epigenomics_evals import data, grading, tasks
from epigenomics_evals.build import grade_controls
from epigenomics_evals.specs import SPECS

REPO = Path(__file__).resolve().parents[1]
needs_data = pytest.mark.skipif(
    not (data.available() and data.panel_available()), reason="input data not fetched"
)


@pytest.fixture(scope="module")
def graded():
    return {task_id: grade_controls(spec) for task_id, spec in SPECS.items()}


@needs_data
@pytest.mark.parametrize("task_id", list(SPECS))
def test_published_key_is_the_recomputed_key(task_id, graded):
    grader, _ = graded[task_id]
    published = json.loads((REPO / "evals" / task_id / "eval.json").read_text())["grader"]
    assert grader == published


@needs_data
@pytest.mark.parametrize("task_id", list(SPECS))
def test_reference_and_variants_pass_and_traps_fail(task_id, graded):
    _, rows = graded[task_id]
    assert sum(r["kind"] == "reference" for r in rows) == 1
    assert any(r["kind"] == "variant" for r in rows) and any(r["kind"] == "trap" for r in rows)
    for r in rows:
        assert r["own"]["passed"] == r["expect_pass"], r["name"]


@needs_data
@pytest.mark.parametrize("task_id", list(SPECS))
def test_latch_grader_agrees_with_own_grader(task_id, graded):
    pytest.importorskip("latch_eval_tools.graders")
    _, rows = graded[task_id]
    for r in rows:
        assert r["latch"]["passed"] == r["own"]["passed"], r["name"]


@needs_data
@pytest.mark.parametrize("task_id", list(SPECS))
def test_no_accepted_analysis_sits_on_a_window_edge(task_id, graded):
    """A variant that passes by a hair makes the verdict a matter of rounding."""
    grader, rows = graded[task_id]
    for r in rows:
        if not r["expect_pass"]:
            continue
        for name, (lo, hi) in grading.windows(grader).items():
            if hi == lo:
                continue
            value = r["answer"][name]
            margin = min(value - lo, hi - value) / (hi - lo)
            assert margin >= 0.03, f"{r['name']}.{name} is within 3% of the window edge"


@needs_data
def test_every_trap_is_described_in_the_task_notes():
    for task_id in SPECS:
        notes = json.loads((REPO / "evals" / task_id / "eval.json").read_text())["notes"]
        for control in tasks.controls(task_id):
            assert control.name in notes, (task_id, control.name)


@needs_data
def test_every_quantity_stated_in_an_answer_option_holds_in_the_data():
    from epigenomics_evals import checks

    claims = checks.option_claims()
    assert claims and all(c["holds"] for c in claims), [c for c in claims if not c["holds"]]
    for task_id in SPECS:
        assert all(isinstance(v, str) and v for v in checks.facts(task_id).values())
