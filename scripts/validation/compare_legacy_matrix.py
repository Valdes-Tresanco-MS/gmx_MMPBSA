#!/usr/bin/env python3
"""Run and numerically compare legacy-compatible calculations.

The current implementation and the 1.6.5 source are run independently against
copied examples with explicit legacy defaults. CSV values are compared with
configurable tolerances; calculation failures are never hidden as parity.
Documented comparison exceptions remain non-zero unless explicitly allowed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

try:
    from .common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        REPO_ROOT,
        csv_comparison_for_cases,
        load_cases,
        metadata_exit_code,
        resolve_selectors,
        run_matrix,
    )
except ImportError:
    from common import (
        DEFAULT_EXAMPLES,
        DEFAULT_RESULTS_ROOT,
        REPO_ROOT,
        csv_comparison_for_cases,
        load_cases,
        metadata_exit_code,
        resolve_selectors,
        run_matrix,
    )


DEFAULT_LEGACY_SELECTORS = ["3", "4", "5", "7", "12", "13", "14", "15", "16", "19", "20", "21", "22", "23", "24"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", default=DEFAULT_LEGACY_SELECTORS,
                        help="Legacy-compatible manifest IDs (default: conservative compatibility set)")
    parser.add_argument("--baseline-python", type=Path, required=True,
                        help="Python interpreter in the isolated official 1.6.5 environment")
    parser.add_argument("--current-python", type=Path, required=True,
                        help="Python interpreter in the current development environment")
    parser.add_argument("--current-source", type=Path, default=REPO_ROOT)
    parser.add_argument("--baseline-source", type=Path, default=None,
                        help="Existing 1.6.5 source checkout; otherwise archive --baseline-ref")
    parser.add_argument("--baseline-ref", default="1.6.5",
                        help="Git ref to archive outside the checkout when --baseline-source is omitted")
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--atol", type=float, default=1e-3,
                        help="Absolute CSV comparison tolerance in kcal/mol")
    parser.add_argument("--rtol", type=float, default=1e-6,
                        help="Relative CSV comparison tolerance")
    parser.add_argument("--allow-known-exceptions", action="store_true",
                        help="Allow documented baseline gaps and expected output changes")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def archive_ref(ref: str, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=False)
    command = ["git", "archive", ref]
    archive = subprocess.run(command, cwd=REPO_ROOT, stdout=subprocess.PIPE, check=True)
    import tarfile
    import io
    with tarfile.open(fileobj=io.BytesIO(archive.stdout), mode="r:") as handle:
        destination = destination.resolve()
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if destination != target and destination not in target.parents:
                raise ValueError(f"Unsafe path in git archive: {member.name!r}")
            try:
                handle.extract(member, destination, filter="data")
            except TypeError:  # Python 3.11 has no extraction filter argument.
                handle.extract(member, destination)
    return destination


def main() -> int:
    args = build_parser().parse_args()
    raw, cases = load_cases()
    selected_ids = resolve_selectors(raw, args.selectors)
    root = args.results_root.expanduser().resolve() / "legacy-comparison"
    root.mkdir(parents=True, exist_ok=True)
    temporary_baseline = None
    baseline_source = args.baseline_source
    try:
        if baseline_source is None:
            temporary_baseline = Path(tempfile.mkdtemp(prefix="gmxMMPBSA-1.6.5-source-", dir=root))
            baseline_source = archive_ref(args.baseline_ref, temporary_baseline / "source")
        current = run_matrix(
            label="current-legacy-settings",
            source_root=args.current_source,
            examples=args.examples,
            python=args.current_python,
            cases=cases,
            selected_ids=selected_ids,
            results_root=root,
            ranks=args.ranks,
            legacy_settings=True,
            legacy_cli_compat=True,
            timeout=args.timeout,
            dry_run=args.dry_run,
        )
        baseline = run_matrix(
            label="1.6.5-legacy-settings",
            source_root=baseline_source,
            examples=args.examples,
            python=args.baseline_python,
            cases=cases,
            selected_ids=selected_ids,
            results_root=root,
            ranks=args.ranks,
            legacy_settings=True,
            legacy_cli_compat=True,
            timeout=args.timeout,
            dry_run=args.dry_run,
        )
        comparison = csv_comparison_for_cases(current, baseline, atol=args.atol, rtol=args.rtol)
        comparison.update({
            "selected_ids": selected_ids,
            "atol": args.atol,
            "rtol": args.rtol,
            "current_run": current,
            "baseline_run": baseline,
        })
        current_root = Path(current["cases"][0]["log"]).parents[1]
        comparison_path = current_root / "comparison.json"
        comparison_path.write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Comparison: {comparison_path}")
        statuses = [record["status"] for record in comparison["comparisons"]]
        accepted_statuses = {"PASS"}
        if args.allow_known_exceptions:
            accepted_statuses.update({"BASELINE-UNSUPPORTED", "EXPECTED-DIFFERENCE"})
        baseline_status_by_id = {
            record["id"]: record.get("status") for record in comparison["comparisons"]
        }
        baseline_ok = metadata_exit_code(baseline) == 0
        if args.allow_known_exceptions and not baseline_ok:
            baseline_ok = all(
                record.get("status") == "PASS"
                or baseline_status_by_id.get(record["id"]) == "BASELINE-UNSUPPORTED"
                for record in baseline["cases"]
            )
        return 0 if (
            metadata_exit_code(current) == 0
            and baseline_ok
            and all(status in accepted_statuses for status in statuses)
        ) else 1
    finally:
        if temporary_baseline is not None:
            shutil.rmtree(temporary_baseline, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
