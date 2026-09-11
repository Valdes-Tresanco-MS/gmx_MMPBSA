#!/usr/bin/env python3
"""Exercise cleanup levels in a persistent external validation directory."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

try:
    from .common import DEFAULT_RESULTS_ROOT, timestamped_results_root
except ImportError:
    from common import DEFAULT_RESULTS_ROOT, timestamped_results_root

from GMXMMPBSA.utils import remove


MINIMAL_FILES = {
    "_GMXMMPBSA_info": "metadata\n",
    "_GMXMMPBSA_pb.mdin": "temporary mdin\n",
    "COM.prmtop": "temporary topology\n",
    "COMPACT_MMXSA_RESULTS.mmxsa": "compact result\n",
    "GMXMMPBSA_membrane_parameters.csv": "diagnostic\n",
    "GMXMMPBSA_membrane_parameters.png": b"diagnostic",
}
FULL_FILES = {
    **MINIMAL_FILES,
    "FINAL_RESULTS_MMPBSA.dat": "summary\n",
    "FINAL_RESULTS_MMPBSA.csv": "summary\n",
    "FINAL_DECOMP_MMPBSA.dat": "decomposition\n",
    "FINAL_DECOMP_MMPBSA.csv": "decomposition\n",
}


def _write_files(root: Path, files: dict[str, str | bytes]) -> None:
    for name, contents in files.items():
        path = root / name
        if isinstance(contents, bytes):
            path.write_bytes(contents)
        else:
            path.write_text(contents, encoding="utf-8")


def _names(root: Path) -> set[str]:
    return {path.name for path in root.iterdir()}


def validate(results_root: Path | None = None) -> dict:
    root = timestamped_results_root(results_root, "cleanup")
    minimal = root / "minimal"
    full = root / "full"
    minimal.mkdir()
    full.mkdir()
    _write_files(minimal, MINIMAL_FILES)
    _write_files(full, FULL_FILES)

    old_cwd = Path.cwd()
    try:
        os.chdir(minimal)
        remove(0)
        minimal_after = _names(minimal)
        os.chdir(full)
        remove(-1)
        full_after = _names(full)
    finally:
        os.chdir(old_cwd)

    minimal_expected = {
        "_GMXMMPBSA_info",
        "COMPACT_MMXSA_RESULTS.mmxsa",
        "GMXMMPBSA_membrane_parameters.csv",
        "GMXMMPBSA_membrane_parameters.png",
    }
    full_expected = {
        "GMXMMPBSA_membrane_parameters.csv",
        "GMXMMPBSA_membrane_parameters.png",
    }
    record = {
        "root": str(root),
        "minimal": {"remaining": sorted(minimal_after), "expected": sorted(minimal_expected)},
        "full": {"remaining": sorted(full_after), "expected": sorted(full_expected)},
    }
    record["status"] = "PASS" if minimal_after == minimal_expected and full_after == full_expected else "FAIL"
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.results_root)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
