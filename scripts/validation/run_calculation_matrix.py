#!/usr/bin/env python3
"""Run selected bundled calculation examples in an external result directory."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    from .common import DEFAULT_EXAMPLES, REPO_ROOT, load_cases, metadata_exit_code, resolve_selectors, run_matrix
except ImportError:
    from common import DEFAULT_EXAMPLES, REPO_ROOT, load_cases, metadata_exit_code, resolve_selectors, run_matrix


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selectors", nargs="*", default=["101"],
                        help="Manifest IDs, suite IDs, or aliases (default: 101/all)")
    parser.add_argument("--source", type=Path, default=REPO_ROOT,
                        help="Source checkout whose GMXMMPBSA code is executed")
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES,
                        help="Canonical examples directory to copy")
    parser.add_argument("--python", type=Path, required=True,
                        help="Explicit Python interpreter for the selected environment")
    parser.add_argument("--results-root", type=Path, default=None,
                        help="External parent directory for timestamped results")
    parser.add_argument("--label", default="current-calculations",
                        help="Result-directory label")
    parser.add_argument("--ranks", type=int, default=1,
                        help="MPI ranks; 1 runs directly, values >1 use mpirun")
    parser.add_argument("--timeout", type=int, default=None,
                        help="Per-example timeout in seconds")
    parser.add_argument("--legacy-settings", action="store_true",
                        help="Patch copied inputs to explicit 1.6.5 defaults")
    parser.add_argument("--dry-run", action="store_true",
                        help="Copy inputs and print commands without executing calculations")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    raw, cases = load_cases()
    try:
        selected_ids = resolve_selectors(raw, args.selectors)
        metadata = run_matrix(
            label=args.label,
            source_root=args.source,
            examples=args.examples,
            python=args.python,
            cases=cases,
            selected_ids=selected_ids,
            results_root=args.results_root,
            ranks=args.ranks,
            legacy_settings=args.legacy_settings,
            timeout=args.timeout,
            dry_run=args.dry_run,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        parser.error(str(exc))
    return metadata_exit_code(metadata)


if __name__ == "__main__":
    raise SystemExit(main())
