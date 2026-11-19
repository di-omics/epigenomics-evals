#!/usr/bin/env python3
"""Deterministic grader for the peak-set choice task. Standard library only.

    python3 grader.py --answer answer.json --key private/answer_key.json

The answer may be a JSON file or any text that contains the JSON between
<EVAL_ANSWER> and </EVAL_ANSWER>. No model is involved and nothing is random:
the same answer and key always give the same verdict.

An answer passes when it is well formed, all 7 counts equal the key exactly,
and the final choice equals the key. The grader also reports whether the
choice follows from the submitted counts under the task's rule. That does not
change the verdict; it tells a wrong count apart from a wrong decision.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

COUNT_FIELDS = (
    "n_peaks_a",
    "n_peaks_b",
    "n_reference_peaks",
    "n_a_peaks_overlapping_reference",
    "n_b_peaks_overlapping_reference",
    "n_reference_peaks_overlapped_by_a",
    "n_reference_peaks_overlapped_by_b",
)
CHOICES = {"a": "A", "b": "B", "tie": "tie", "neither": "neither"}
ALL_FIELDS = COUNT_FIELDS + ("answer",)


def extract_answer(text: str) -> dict:
    """The answer object from raw agent output, or from a plain JSON file."""
    blocks = re.findall(r"<EVAL_ANSWER>(.*?)</EVAL_ANSWER>", text, flags=re.S)
    payload = blocks[-1] if blocks else text
    value = json.loads(payload)
    if not isinstance(value, dict):
        raise ValueError("answer must be a JSON object")
    return value


def normalise(answer: dict) -> tuple[dict, list[str]]:
    """Return (clean answer, schema errors). Whole-number floats count as integers."""
    errors: list[str] = []
    clean: dict = {}
    for name in sorted(set(answer) - set(ALL_FIELDS)):
        errors.append(f"unexpected field: {name}")
    for name in COUNT_FIELDS:
        if name not in answer:
            errors.append(f"missing field: {name}")
            continue
        value = answer[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            errors.append(f"{name}: must be an integer")
        elif isinstance(value, float) and value != int(value):
            errors.append(f"{name}: must be a whole number")
        elif value < 0:
            errors.append(f"{name}: must not be negative")
        else:
            clean[name] = int(value)
    if "answer" not in answer:
        errors.append("missing field: answer")
    elif not isinstance(answer["answer"], str) or answer["answer"].strip().lower() not in CHOICES:
        errors.append('answer: must be one of "A", "B", "tie", "neither"')
    else:
        clean["answer"] = CHOICES[answer["answer"].strip().lower()]
    return clean, errors


def qualifies(overlapping: int, total: int) -> bool:
    """Strictly more than 75% of the set's peaks overlap a reference peak."""
    return 4 * overlapping > 3 * total


def implied_choice(counts: dict) -> str:
    """The choice the task's rule gives for a set of counts."""
    ok_a = qualifies(counts["n_a_peaks_overlapping_reference"], counts["n_peaks_a"])
    ok_b = qualifies(counts["n_b_peaks_overlapping_reference"], counts["n_peaks_b"])
    if not ok_a and not ok_b:
        return "neither"
    if ok_a != ok_b:
        return "A" if ok_a else "B"
    got_a = counts["n_reference_peaks_overlapped_by_a"]
    got_b = counts["n_reference_peaks_overlapped_by_b"]
    if got_a == got_b:
        return "tie"
    return "A" if got_a > got_b else "B"


def grade(answer: dict, key: dict) -> dict:
    clean, errors = normalise(answer)
    if errors:
        return {"passed": False, "schema_errors": errors, "fields": {}, "wrong_fields": [],
                "choice_correct": False, "choice_follows_from_counts": None}
    fields = {
        name: {"submitted": clean[name], "expected": key[name], "ok": clean[name] == key[name]}
        for name in COUNT_FIELDS
    }
    wrong = [name for name in COUNT_FIELDS if not fields[name]["ok"]]
    choice_correct = clean["answer"] == key["answer"]
    return {
        "passed": not wrong and choice_correct,
        "schema_errors": [],
        "fields": fields,
        "wrong_fields": wrong,
        "choice_correct": choice_correct,
        "choice_follows_from_counts": implied_choice(clean) == clean["answer"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--answer", required=True, type=Path)
    parser.add_argument("--key", required=True, type=Path)
    args = parser.parse_args()
    key = json.loads(args.key.read_text())["key"]
    try:
        answer = extract_answer(args.answer.read_text())
    except (ValueError, json.JSONDecodeError) as exc:
        print(json.dumps({"passed": False, "schema_errors": [f"could not read an answer: {exc}"]}, indent=1))
        return 1
    verdict = grade(answer, key)
    print(json.dumps(verdict, indent=1))
    return 0 if verdict["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
