#!/usr/bin/env python3
"""Execute the calculation-combination matrix in isolated external copies."""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from .common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        _fingerprint,
        _set_namelist_value,
        apply_legacy_settings,
        copy_examples,
        entrypoint_command,
        environment_for,
        launched_command,
        load_cases,
        metadata_exit_code,
        timestamped_results_root,
    )
    from .plan_combination_matrix import build_plan, load_matrix
except ImportError:
    from common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        _fingerprint,
        _set_namelist_value,
        apply_legacy_settings,
        copy_examples,
        entrypoint_command,
        environment_for,
        launched_command,
        load_cases,
        metadata_exit_code,
        timestamped_results_root,
    )
    from plan_combination_matrix import build_plan, load_matrix


def _format_input_value(value: Any) -> str:
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, list):
        return ",".join(_format_input_value(item) for item in value)
    return str(value)


def apply_overrides(examples_dir: Path, workdir: str, input_name: str,
                    overrides: dict[str, dict[str, Any]]) -> None:
    input_path = examples_dir / workdir / input_name
    text = input_path.read_text(encoding="utf-8")
    for section, values in overrides.items():
        for key, value in values.items():
            text = _set_namelist_value(text, section, key, _format_input_value(value))
    input_path.write_text(text, encoding="utf-8")


def _expected_log_tokens(log_path: Path, tokens: list[str], *extra_paths: Path) -> list[str]:
    if not tokens:
        return []
    paths = (log_path,) + extra_paths
    text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in paths
        if path.is_file()
    ).lower()
    return [token for token in tokens if token.lower() not in text]


def run_combination_matrix(
    *,
    matrix: dict[str, Any],
    selectors: list[str],
    source: Path,
    examples: Path,
    python: Path,
    results_root: Path | None,
    label: str,
    ranks: int,
    legacy_settings: bool,
    timeout: int | None,
) -> dict[str, Any]:
    _, manifest_cases = load_cases()
    plan = build_plan(matrix, selectors)
    result_root = timestamped_results_root(results_root, label)
    env = environment_for(python, source)
    metadata: dict[str, Any] = {
        "matrix": matrix["name"],
        "selectors": selectors or ["core"],
        "source": str(source.expanduser().resolve()),
        "examples": str(examples.expanduser().resolve()),
        "python": str(python.expanduser().resolve()),
        "ranks": ranks,
        "legacy_settings": legacy_settings,
        "timeout": timeout,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [],
    }

    for row in plan:
        row_root = result_root / "cases" / row["id"]
        row_root.mkdir(parents=True, exist_ok=False)
        copied_examples = copy_examples(examples, row_root)
        case = manifest_cases[row["base_case"]]
        if legacy_settings:
            apply_legacy_settings(copied_examples, case)
        apply_overrides(copied_examples, case.workdir, case.input, row["overrides"])
        workdir = copied_examples / case.workdir
        command = launched_command(
            entrypoint_command(case, python, source), ranks, env
        )
        log_file = row_root / "cases.log"
        before = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
        record: dict[str, Any] = {
            "id": row["id"],
            "base_case": row["base_case"],
            "description": row["description"],
            "workdir": str(workdir),
            "input": str(workdir / case.input),
            "overrides": row["overrides"],
            "legacy_candidate": row["legacy_candidate"],
            "slow": row["slow"],
            "requires": row["requires"],
            "expected_log": row["expected_log"],
            "command": command,
            "log": str(log_file),
            "expected_outputs": list(case.expected_outputs),
            "before": before,
        }
        print(f"[{row['id']}] {row['description']}")
        print("     " + " ".join(command))
        started = time.monotonic()
        try:
            with log_file.open("w", encoding="utf-8") as log_handle:
                import subprocess
                completed = subprocess.run(
                    command,
                    cwd=workdir,
                    env=env,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
            record["return_code"] = completed.returncode
        except subprocess.TimeoutExpired:
            record["return_code"] = None
            record["status"] = "TIMEOUT"
        except Exception as exc:  # capture launch failures per row
            record["return_code"] = None
            record["status"] = "LAUNCH-ERROR"
            record["error"] = repr(exc)

        record["duration_seconds"] = round(time.monotonic() - started, 3)
        after = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
        record["after"] = after
        missing = [name for name, fingerprint in after.items() if not fingerprint["exists"]]
        unchanged = [
            name for name in case.expected_outputs
            if before[name].get("exists") and before[name] == after[name]
        ]
        record["missing_outputs"] = missing
        record["unchanged_outputs"] = unchanged
        missing_log_tokens = _expected_log_tokens(
            log_file,
            row["expected_log"],
            workdir / "gmx_MMPBSA.log",
            workdir / "_GMXMMPBSA_info",
        )
        record["missing_expected_log"] = missing_log_tokens
        if "status" not in record:
            if record["return_code"] != 0:
                record["status"] = "FAIL"
            elif missing:
                record["status"] = "FAIL-MISSING-OUTPUT"
            elif unchanged:
                record["status"] = "FAIL-STALE-OUTPUT"
            elif missing_log_tokens:
                record["status"] = "FAIL-MISSING-EXPECTED-LOG"
            else:
                record["status"] = "PASS"
        metadata["cases"].append(record)
        print(f"     {record['status']} ({record.get('duration_seconds', 0):.1f}s)")

    metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["summary"] = {
        status: sum(record.get("status") == status for record in metadata["cases"])
        for status in sorted({record.get("status") for record in metadata["cases"]})
    }
    metadata_path = result_root / "run.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Results: {result_root}")
    return metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", help="Matrix tiers or case IDs (default: core)")
    parser.add_argument("--matrix", type=Path,
                        default=Path(__file__).with_name("calculation_combination_matrix.json"))
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--label", default="combination-matrix")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--legacy-settings", action="store_true")
    parser.add_argument("--execute", action="store_true",
                        help="Required acknowledgement that calculations will be launched")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.execute:
        raise SystemExit("Refusing to launch: add --execute to run the matrix")
    matrix = load_matrix(args.matrix.expanduser().resolve())
    metadata = run_combination_matrix(
        matrix=matrix,
        selectors=args.selectors,
        source=args.source,
        examples=args.examples,
        python=args.python,
        results_root=args.results_root,
        label=args.label,
        ranks=args.ranks,
        legacy_settings=args.legacy_settings,
        timeout=args.timeout,
    )
    return metadata_exit_code(metadata)


if __name__ == "__main__":
    raise SystemExit(main())
