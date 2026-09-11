#!/usr/bin/env python3
"""Print a calculation-combination plan without running any calculations.

This is intentionally a planning tool.  It resolves matrix rows to the
manifest command and records input overrides, but it never copies examples,
edits inputs, launches gmx_MMPBSA, or creates result directories.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

try:
    from .common import MANIFEST_PATH, REPO_ROOT, load_cases
except ImportError:
    from common import MANIFEST_PATH, REPO_ROOT, load_cases


MATRIX_PATH = Path(__file__).with_name("calculation_combination_matrix.json")


def load_matrix(matrix_path: Path = MATRIX_PATH) -> dict[str, Any]:
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    if matrix.get("version") != 1:
        raise ValueError(f"Unsupported matrix version: {matrix.get('version')!r}")
    cases = matrix.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("The combination matrix must define a non-empty cases list")
    ids = [entry.get("id") for entry in cases]
    if any(not isinstance(case_id, str) or not case_id for case_id in ids):
        raise ValueError("Every matrix case must have a non-empty string id")
    if len(ids) != len(set(ids)):
        raise ValueError("Matrix case IDs must be unique")
    return matrix


def resolve_rows(matrix: dict[str, Any], selectors: list[str]) -> list[dict[str, Any]]:
    rows = {entry["id"]: entry for entry in matrix["cases"]}
    tiers = matrix.get("tiers", {})
    selected_ids: list[str] = []
    for selector in selectors or ["core"]:
        names = tiers.get(selector, [selector])
        if isinstance(names, str):
            names = [names]
        for case_id in names:
            if case_id not in rows:
                raise ValueError(f"Unknown combination-matrix case or tier: {case_id!r}")
            if case_id not in selected_ids:
                selected_ids.append(case_id)
    return [rows[case_id] for case_id in selected_ids]


def build_plan(matrix: dict[str, Any], selectors: list[str]) -> list[dict[str, Any]]:
    _, manifest_cases = load_cases(MANIFEST_PATH)
    plan = []
    for row in resolve_rows(matrix, selectors):
        base_case_id = int(row["base_case"])
        try:
            base_case = manifest_cases[base_case_id]
        except KeyError as exc:
            raise ValueError(f"Matrix row {row['id']!r} references unknown manifest case {base_case_id}") from exc
        plan.append({
            "id": row["id"],
            "description": row["description"],
            "base_case": base_case_id,
            "base_name": base_case.name,
            "workdir": str(REPO_ROOT / "examples" / base_case.workdir),
            "input": base_case.input,
            "command_args": list(base_case.command_args),
            "overrides": row.get("overrides", {}),
            "legacy_candidate": bool(row.get("legacy_candidate", False)),
            "slow": bool(row.get("slow", base_case.slow)),
            "requires": list(row.get("requires", [])),
            "expected_log": list(row.get("expected_log", [])),
        })
    return plan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", help="Matrix tiers or case IDs (default: core)")
    parser.add_argument("--matrix", type=Path, default=MATRIX_PATH,
                        help="Combination-matrix metadata file")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    matrix = load_matrix(args.matrix.expanduser().resolve())
    plan = build_plan(matrix, args.selectors)
    if args.json:
        print(json.dumps({"matrix": matrix["name"], "plan": plan}, indent=2, sort_keys=True))
    else:
        print(f"Matrix: {matrix['name']}")
        print("Planning only: no calculations, input edits, or result directories will be created.")
        for row in plan:
            flags = []
            if row["slow"]:
                flags.append("slow")
            if row["legacy_candidate"]:
                flags.append("legacy-candidate")
            suffix = f" [{', '.join(flags)}]" if flags else ""
            print(f"{row['id']}: case {row['base_case']} ({row['base_name']}){suffix}")
            print(f"  {row['description']}")
            print(f"  overrides={json.dumps(row['overrides'], sort_keys=True)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
