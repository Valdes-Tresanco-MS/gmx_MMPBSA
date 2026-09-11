#!/usr/bin/env python3
"""Compare paired gmx_MMPBSA API results from two external result trees.

The comparison is deliberately independent of CSV formatting.  It loads each
paired ``_GMXMMPBSA_info`` or ``COMPACT_MMXSA_RESULTS.mmxsa`` file through the
current public API and compares numeric cells in raw data, summaries,
correlations, binding data, and decomposition data.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    from GMXMMPBSA import API
except ModuleNotFoundError:  # pragma: no cover - supports direct script use
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from GMXMMPBSA import API


ENTROPY_MODELS = {"nmode", "qh", "ie", "c2"}
DEFAULT_POLICY = Path(__file__).with_name("api_comparison_policy.json")


def _discover(root: Path, filename: str) -> dict[Path, Path]:
    root = root.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"API result directory does not exist: {root}")
    return {path.relative_to(root): path for path in root.rglob(filename)}


def _policy_for(relative: Path, policy: dict[str, Any]) -> dict[str, Any] | None:
    text = relative.as_posix()
    matches = [
        (prefix, value)
        for prefix, value in policy.get("cases", {}).items()
        if text == prefix or text.startswith(prefix + "/")
    ]
    return max(matches, key=lambda item: len(item[0]))[1] if matches else None


def _key(value: Any) -> str:
    return repr(value)


def _numeric(value: Any) -> bool:
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)


def _flatten(value: Any, prefix: tuple[str, ...] = ()) -> dict[tuple[str, ...], float]:
    """Flatten numeric Series/DataFrame leaves with their labels included."""

    result: dict[tuple[str, ...], float] = {}
    if isinstance(value, dict):
        for name, child in value.items():
            result.update(_flatten(child, prefix + (_key(name),)))
    elif isinstance(value, pd.Series):
        for index, item in value.items():
            if _numeric(item):
                result[prefix + ("index", _key(index))] = float(item)
    elif isinstance(value, pd.DataFrame):
        for index, row in value.iterrows():
            for column, item in row.items():
                if _numeric(item):
                    result[prefix + ("cell", _key(index), _key(column))] = float(item)
    return result


def _api_payload(api: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    data = api.data.get("normal", {})
    energy_models = tuple(name for name in data if name not in ENTROPY_MODELS)
    entropy_models = tuple(name for name in data if name in ENTROPY_MODELS)
    energy = api.get_energy(verbose=False)
    entropy = api.get_entropy(verbose=False) if entropy_models else {"data": {}, "summary": {}, "correlation": {}}
    decomp_present = bool(api.data.get("decomp_normal") or api.data.get("decomp_mutant"))
    decomp = api.get_decomp_energy(verbose=False) if decomp_present else {"data": {}}
    payload = {
        "energy.data": energy["data"],
        "energy.summary": energy["summary"],
        "energy.correlation": energy["correlation"],
        "entropy.data": entropy["data"],
        "entropy.summary": entropy["summary"],
        "decomposition.data": decomp["data"],
    }
    if energy_models and entropy_models:
        binding = api.get_binding(energy["summary"], entropy["summary"], verbose=False)
        payload.update({
            "binding.data": binding["data"],
            "binding.correlation": binding["correlation"],
        })
    metadata = {
        "frames": int(api.get_info()["numframes"]),
        "nmode_frames": int(api.get_info()["numframes_nmode"]),
        "stability": bool(api.get_files().stability),
    }
    return payload, metadata


def compare_pair(current_path: Path, baseline_path: Path, atol: float, rtol: float) -> dict[str, Any]:
    result: dict[str, Any] = {
        "current": str(current_path),
        "baseline": str(baseline_path),
        "differences": [],
        "errors": [],
    }
    try:
        current_payload, current_meta = _api_payload(API.load(current_path))
        baseline_payload, baseline_meta = _api_payload(API.load(baseline_path))
        result["metadata"] = {"current": current_meta, "baseline": baseline_meta}
        if current_meta != baseline_meta:
            result["differences"].append({
                "field": "metadata",
                "current": current_meta,
                "baseline": baseline_meta,
            })
        current_flat = _flatten(current_payload)
        baseline_flat = _flatten(baseline_payload)
        for key in sorted(set(current_flat) | set(baseline_flat)):
            if key not in current_flat or key not in baseline_flat:
                result["differences"].append({
                    "field": ".".join(key),
                    "kind": "missing-field",
                    "current": current_flat.get(key),
                    "baseline": baseline_flat.get(key),
                })
                continue
            left, right = current_flat[key], baseline_flat[key]
            if math.isnan(left) and math.isnan(right):
                continue
            if not np.isclose(left, right, atol=atol, rtol=rtol, equal_nan=True):
                result["differences"].append({
                    "field": ".".join(key),
                    "kind": "numeric",
                    "current": left,
                    "baseline": right,
                    "absolute_difference": abs(left - right),
                })
        result["numeric_fields"] = len(set(current_flat) | set(baseline_flat))
    except Exception as exc:  # retain all pair outcomes in the report
        result["errors"].append({"type": type(exc).__name__, "message": str(exc)})
    result["status"] = "PASS" if not result["errors"] and not result["differences"] else "FAIL"
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atol", type=float, default=1e-3)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--formats", nargs="+", choices=("info", "compact"), default=("info", "compact"))
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    policy = json.loads(args.policy.expanduser().resolve().read_text(encoding="utf-8"))
    pairs: list[tuple[Path, Path, Path]] = []
    for fmt in args.formats:
        filename = "_GMXMMPBSA_info" if fmt == "info" else "COMPACT_MMXSA_RESULTS.mmxsa"
        current = _discover(args.current_root, filename)
        baseline = _discover(args.baseline_root, filename)
        for relative in sorted(current.keys() & baseline.keys()):
            pairs.append((relative, current[relative], baseline[relative]))
    if not pairs:
        raise SystemExit("No paired API result files found")

    comparisons = []
    for relative, current, baseline in pairs:
        print(f"[{relative}]")
        comparison = compare_pair(current, baseline, args.atol, args.rtol)
        comparison["relative"] = str(relative)
        expected = _policy_for(relative, policy)
        if expected and comparison["status"] == "FAIL" and expected.get("outcome") == "expected-difference":
            comparison["status"] = "EXPECTED-DIFFERENCE"
            comparison["expected_difference"] = expected["reason"]
        comparisons.append(comparison)
        print(f"     {comparison['status']} ({len(comparison['differences'])} differences)")
    summary = {
        status: sum(item["status"] == status for item in comparisons)
        for status in sorted({item["status"] for item in comparisons})
    }
    report = {
        "current_root": str(args.current_root.expanduser().resolve()),
        "baseline_root": str(args.baseline_root.expanduser().resolve()),
        "atol": args.atol,
        "rtol": args.rtol,
        "summary": summary,
        "comparisons": comparisons,
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Summary: {summary}")
    print(f"Report: {output}")
    return 0 if summary.get("FAIL", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
