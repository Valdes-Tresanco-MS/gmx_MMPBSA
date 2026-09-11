#!/usr/bin/env python3
"""Compare isolated current and 1.6.5 combination-matrix runs."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

try:
    from .common import compare_csv
except ImportError:
    from common import compare_csv


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _first_diagnostic(log_path: Path) -> str:
    if not log_path.is_file():
        return "baseline log is missing"
    lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        if "[ERROR" in line or "error:" in line.lower() or "Error:" in line:
            return line.strip()
    return "baseline calculation did not complete successfully"


def compare_runs(current: dict, baseline: dict, *, atol: float, rtol: float) -> dict:
    current_by_id = {record["id"]: record for record in current["cases"]}
    baseline_by_id = {record["id"]: record for record in baseline["cases"]}
    comparisons = []
    for row_id in sorted(set(current_by_id) | set(baseline_by_id)):
        left = current_by_id.get(row_id)
        right = baseline_by_id.get(row_id)
        record = {"id": row_id, "status": "SKIP"}
        if left is None or right is None:
            record["status"] = "FAIL"
            record["reason"] = "row missing from one run"
            comparisons.append(record)
            continue
        if left.get("status") != "PASS" or right.get("status") != "PASS":
            if left.get("status") == "PASS" and right.get("status") != "PASS":
                record["status"] = "BASELINE-UNSUPPORTED"
                record["reason"] = _first_diagnostic(Path(right["log"]))
                record["baseline_status"] = right.get("status")
            else:
                record["status"] = "FAIL"
                record["reason"] = (
                    f"run status: {left.get('status')} vs {right.get('status')}"
                )
            comparisons.append(record)
            continue

        differences = []
        for output in left.get("expected_outputs", []):
            if not output.endswith(".csv"):
                continue
            errors = compare_csv(
                Path(left["workdir"]) / output,
                Path(right["workdir"]) / output,
                atol=atol,
                rtol=rtol,
            )
            differences.extend(f"{output}: {error}" for error in errors)
        record["status"] = "PASS" if not differences else "FAIL"
        record["differences"] = differences[:50]
        comparisons.append(record)
    return {
        "current_run": current,
        "baseline_run": baseline,
        "atol": atol,
        "rtol": rtol,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "comparisons": comparisons,
        "summary": {
            status: sum(item["status"] == status for item in comparisons)
            for status in sorted({item["status"] for item in comparisons})
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current_run", type=Path)
    parser.add_argument("baseline_run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atol", type=float, default=1e-3)
    parser.add_argument("--rtol", type=float, default=1e-6)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = compare_runs(
        _load(args.current_run),
        _load(args.baseline_run),
        atol=args.atol,
        rtol=args.rtol,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Comparison: {args.output}")
    print(json.dumps(report["summary"], sort_keys=True))
    return 0 if all(status in {"PASS", "BASELINE-UNSUPPORTED"} for status in report["summary"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
