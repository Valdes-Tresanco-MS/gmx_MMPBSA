#!/usr/bin/env python3
"""Exercise expected CLI failures and diagnostic-bundle policy."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

try:
    from .common import DEFAULT_EXAMPLES, DEFAULT_RESULTS_ROOT, copy_examples, entrypoint_command, environment_for, load_cases, timestamped_results_root
    from .run_combination_matrix import apply_overrides
except ImportError:
    from common import DEFAULT_EXAMPLES, DEFAULT_RESULTS_ROOT, copy_examples, entrypoint_command, environment_for, load_cases, timestamped_results_root
    from run_combination_matrix import apply_overrides


def _run_case(root: Path, examples: Path, case, python: Path, source: Path,
              environment: dict[str, str], no_bundle: bool) -> dict:
    copied = copy_examples(examples, root)
    apply_overrides(copied, case.workdir, case.input, {"gb": {"igb": 999}})
    workdir = copied / case.workdir
    command = entrypoint_command(case, python, source)
    if no_bundle:
        command.append("--no-error-bundle")
    completed = subprocess.run(
        command, cwd=workdir, env=environment, capture_output=True, text=True, check=False
    )
    bundles = sorted(path.name for path in workdir.glob("gmx_MMPBSA_error_bundle_*.zip"))
    return {
        "no_error_bundle": no_bundle,
        "return_code": completed.returncode,
        "stderr": completed.stderr[-2000:],
        "bundles": bundles,
        "status": "PASS" if completed.returncode != 0 and bool(bundles) == (not no_bundle) else "FAIL",
    }


def validate(python: Path, source: Path, examples: Path, results_root: Path | None) -> dict:
    _, cases = load_cases()
    root = timestamped_results_root(results_root, "negative")
    environment = environment_for(python, source)
    records = []
    for no_bundle in (False, True):
        scenario_root = root / ("no-bundle" if no_bundle else "bundle")
        scenario_root.mkdir()
        records.append(_run_case(scenario_root, examples, cases[3], python, source, environment, no_bundle))
    return {"root": str(root), "cases": records,
            "status": "PASS" if all(item["status"] == "PASS" for item in records) else "FAIL"}


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", required=True,
                        help="Required acknowledgement that failing CLI cases will be launched")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.python.expanduser().resolve(), args.source, args.examples, args.results_root)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
