#!/usr/bin/env python3
"""Run independent calculation rows concurrently in isolated copies."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from .common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        _fingerprint,
        copy_examples,
        entrypoint_command,
        environment_for,
        launched_command,
        load_cases,
        metadata_exit_code,
        timestamped_results_root,
    )
    from .run_combination_matrix import apply_overrides
    from .plan_combination_matrix import build_plan, load_matrix
except ImportError:
    from common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        _fingerprint,
        copy_examples,
        entrypoint_command,
        environment_for,
        launched_command,
        load_cases,
        metadata_exit_code,
        timestamped_results_root,
    )
    from run_combination_matrix import apply_overrides
    from plan_combination_matrix import build_plan, load_matrix


def _run_row(row, result_root, source, examples, python, env, ranks, timeout):
    _, manifest_cases = load_cases()
    case = manifest_cases[row["base_case"]]
    row_root = result_root / "cases" / row["id"]
    row_root.mkdir(parents=True, exist_ok=False)
    copied = copy_examples(examples, row_root)
    apply_overrides(copied, case.workdir, case.input, row["overrides"])
    workdir = copied / case.workdir
    command = launched_command(entrypoint_command(case, python, source), ranks, env)
    log_path = row_root / "cases.log"
    before = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
    record = {
        "id": row["id"],
        "description": row["description"],
        "workdir": str(workdir),
        "command": command,
        "log": str(log_path),
        "expected_outputs": list(case.expected_outputs),
        "before": before,
    }
    started = time.monotonic()
    try:
        with log_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                command, cwd=workdir, env=env, stdout=handle, stderr=subprocess.STDOUT,
                timeout=timeout, check=False,
            )
        record["return_code"] = completed.returncode
    except subprocess.TimeoutExpired:
        record["return_code"] = None
        record["status"] = "TIMEOUT"
    except Exception as exc:
        record["return_code"] = None
        record["status"] = "LAUNCH-ERROR"
        record["error"] = repr(exc)
    record["duration_seconds"] = round(time.monotonic() - started, 3)
    after = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
    record["after"] = after
    missing = [name for name, fingerprint in after.items() if not fingerprint["exists"]]
    if "status" not in record:
        record["status"] = "PASS" if record["return_code"] == 0 and not missing else "FAIL"
    return record


def run_validation(matrix, selectors, source, examples, python, results_root, label,
                   ranks, timeout, workers):
    _, cases = load_cases()
    plan = build_plan(matrix, selectors)
    if len(plan) < 2:
        raise ValueError("Concurrency validation requires at least two matrix rows")
    result_root = timestamped_results_root(results_root, label)
    env = environment_for(python, source)
    metadata = {
        "selectors": selectors,
        "workers": workers,
        "ranks": ranks,
        "source": str(source.expanduser().resolve()),
        "examples": str(examples.expanduser().resolve()),
        "python": str(python.expanduser().resolve()),
        "cases": [],
    }
    with ThreadPoolExecutor(max_workers=min(workers, len(plan))) as pool:
        futures = {
            pool.submit(_run_row, row, result_root, source, examples, python, env, ranks, timeout): row
            for row in plan
        }
        for future in as_completed(futures):
            record = future.result()
            metadata["cases"].append(record)
            print(f"[{record['id']}] {record['status']} ({record['duration_seconds']:.1f}s)")
    metadata["cases"].sort(key=lambda record: record["id"])
    metadata["summary"] = {
        status: sum(record["status"] == status for record in metadata["cases"])
        for status in sorted({record["status"] for record in metadata["cases"]})
    }
    (result_root / "run.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Results: {result_root}")
    return metadata


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", default=["GB-01", "GB-02"])
    parser.add_argument("--matrix", type=Path, default=Path(__file__).with_name("calculation_combination_matrix.json"))
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--label", default="concurrency-validation")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--execute", action="store_true", required=True,
                        help="Required acknowledgement that calculations will be launched")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    metadata = run_validation(
        load_matrix(args.matrix.expanduser().resolve()), args.selectors,
        args.source, args.examples, args.python, args.results_root, args.label,
        args.ranks, args.timeout, args.workers,
    )
    return metadata_exit_code(metadata)


if __name__ == "__main__":
    raise SystemExit(main())
