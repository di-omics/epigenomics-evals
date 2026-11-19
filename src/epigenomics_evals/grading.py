"""Grader configs in LatchBio's public format, and two ways to apply them.

own_check     a few lines of arithmetic written here, with no dependencies.
latch_check   the published latch-eval-tools package, when it is installed.

Every scripted control is graded both ways and the two verdicts must agree.
"""
from __future__ import annotations

from .specs import Field, Spec


def _round(value: float, f: Field):
    if f.kind in ("int", "choice") or f.digits == 0:
        return int(round(value))
    return round(float(value), f.digits)


def _tolerance(f: Field) -> dict:
    tol = {"type": "absolute"}
    tol.update(f.tolerance or {"value": 0})
    return tol


def grader_config(spec: Spec, reference: dict) -> dict:
    """Build the eval.json grader block from a spec and its reference answer."""
    numbers = [f for f in spec.fields if f.kind != "choice"]
    number_child = {
        "type": "numeric_tolerance",
        "config": {
            "ground_truth": {f.name: _round(reference[f.name], f) for f in numbers},
            "tolerances": {f.name: _tolerance(f) for f in numbers},
        },
    }
    choices = [f for f in spec.fields if f.kind == "choice"]
    if not choices:
        return number_child
    choice_child = {
        "type": "numeric_tolerance",
        "config": {
            "ground_truth": {f.name: spec.correct_choice for f in choices},
            "tolerances": {f.name: {"type": "absolute", "value": 0} for f in choices},
        },
    }
    return {"type": "all_of", "config": {"children": [number_child, choice_child]}}


def windows(grader: dict) -> dict[str, tuple[float, float]]:
    """Accepted closed interval per answer field."""
    children = grader["config"]["children"] if grader["type"] == "all_of" else [grader]
    out = {}
    for child in children:
        cfg = child["config"]
        for name, truth in cfg["ground_truth"].items():
            tol = cfg["tolerances"][name]
            value = tol.get("value", 0)
            out[name] = (truth - tol.get("lower", value), truth + tol.get("upper", value))
    return out


def own_check(answer: dict, grader: dict) -> dict:
    """Independent verdict: every field present, numeric and inside its window."""
    detail = {}
    for name, (lo, hi) in windows(grader).items():
        value = answer.get(name)
        ok = isinstance(value, (int, float)) and not isinstance(value, bool) and lo <= value <= hi
        gap = 0.0
        if isinstance(value, (int, float)) and not ok:
            gap = lo - value if value < lo else value - hi
        detail[name] = {"value": value, "window": [lo, hi], "ok": ok, "gap": gap}
    return {"passed": all(d["ok"] for d in detail.values()), "fields": detail}


def latch_available() -> bool:
    try:
        import latch_eval_tools.graders  # noqa: F401
    except Exception:
        return False
    return True


def latch_version() -> str | None:
    try:
        from importlib.metadata import version

        return version("latch-eval-tools")
    except Exception:
        return None


def latch_check(answer: dict, grader: dict) -> dict | None:
    """Verdict from LatchBio's own grader. None when the package is absent."""
    if not latch_available():
        return None
    from latch_eval_tools.graders import get_grader

    result = get_grader(grader["type"]).evaluate(answer, grader["config"])
    return {"passed": bool(result.passed), "score": float(result.score)}
