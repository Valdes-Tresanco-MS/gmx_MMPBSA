#!/usr/bin/env python3
"""Validate disk-backed analyzer success or its explicit dependency guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from GMXMMPBSA import API


def validate(result: Path) -> dict:
    record = {"result": str(result.expanduser().resolve()), "mode": None, "error": None}
    try:
        api = API.load(result)
        data = api.get_ana_data(
            performance_options={"energy_memory": False, "decomp_memory": False},
            verbose=False,
        )
        record["mode"] = "disk-backed"
        record["sections"] = sorted(data)
        record["status"] = "PASS"
    except ImportError as exc:
        message = str(exc)
        record["mode"] = "dependency-guard"
        record["error"] = message
        record["status"] = "PASS" if "pyarrow" in message and "fastparquet" in message else "FAIL"
    except Exception as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
        record["status"] = "FAIL"
    return record


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path, help="External _GMXMMPBSA_info or compact result")
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.result)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
