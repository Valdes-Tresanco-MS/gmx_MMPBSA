#!/usr/bin/env python3
"""Validate the public API against external calculation and compact outputs.

The inputs are files or directories outside the checkout.  Directories are
searched recursively for ``_GMXMMPBSA_info`` and
``COMPACT_MMXSA_RESULTS.mmxsa``.  The runner exercises loading, public data
accessors, binding/decomposition routing, and in-memory analyzer preparation;
it does not launch a calculation.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

try:
    from GMXMMPBSA import API
except ModuleNotFoundError:  # pragma: no cover - supports direct script use
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from GMXMMPBSA import API


ENTROPY_MODELS = {"nmode", "qh", "ie", "c2"}
SUMMARY_LABELS = ("Average", "SD", "SEM", "Block SD", "Block SEM")


def discover(paths: Iterable[Path]) -> list[Path]:
    """Return unique API result files from the supplied paths."""

    found: set[Path] = set()
    for supplied in paths:
        path = supplied.expanduser().resolve()
        if path.is_file():
            found.add(path)
            continue
        if not path.is_dir():
            raise FileNotFoundError(f"API input does not exist: {path}")
        found.update(path.rglob("_GMXMMPBSA_info"))
        found.update(path.rglob("COMPACT_MMXSA_RESULTS.mmxsa"))
    return sorted(found)


def _walk_summary(value: Any):
    if isinstance(value, dict):
        for child in value.values():
            yield from _walk_summary(child)
    elif isinstance(value, (pd.Series, pd.DataFrame)):
        yield value


def _summary_statistics(summary: Any) -> dict[str, int]:
    leaves = list(_walk_summary(summary))
    missing = 0
    for leaf in leaves:
        labels = set(leaf.index) if isinstance(leaf, (pd.Series, pd.DataFrame)) else set()
        if not set(SUMMARY_LABELS).issubset(labels):
            missing += 1
    return {"leaves": len(leaves), "missing_statistics": missing}


def _data_models(api: Any, entropy: bool = False) -> tuple[str, ...]:
    models = tuple(api.data.get("normal", {}).keys())
    return tuple(model for model in models if (model in ENTROPY_MODELS) == entropy)


def validate_result(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": str(path),
        "format": "compact" if path.suffix == ".mmxsa" else "info",
        "checks": {},
        "errors": [],
    }
    try:
        api = API.load(path)
        info = api.get_info()
        inputs = api.get_input()
        files = api.get_files()
        record["checks"]["load"] = True
        record["frames"] = int(info["numframes"])
        record["stability"] = bool(files.stability)
        record["input_sections"] = sorted(inputs)
        if record["frames"] < 1:
            raise ValueError("loaded result reports no energy frames")
        for accessor, value in (("get_info", info), ("get_input", inputs), ("get_files", files)):
            if value is None:
                raise ValueError(f"{accessor} returned None")
            record["checks"][accessor] = True

        energy_models = _data_models(api)
        energy = api.get_energy(verbose=False)
        record["energy"] = _summary_statistics(energy["summary"])
        record["checks"]["get_energy"] = True
        if energy_models and record["energy"]["missing_statistics"]:
            raise ValueError("energy summary is missing one or more statistics")

        entropy_models = _data_models(api, entropy=True)
        entropy = api.get_entropy(verbose=False) if entropy_models else {"summary": {}}
        record["entropy"] = _summary_statistics(entropy["summary"])
        record["checks"]["get_entropy"] = True

        decomp_present = bool(api.data.get("decomp_normal") or api.data.get("decomp_mutant"))
        decomp = api.get_decomp_energy(verbose=False) if decomp_present else {"map": {}, "data": {}}
        record["decomposition"] = {"present": decomp_present, "sections": sorted(decomp["data"])}
        record["checks"]["get_decomp_energy"] = True

        if energy_models and entropy_models:
            binding = api.get_binding(energy["summary"], entropy["summary"], verbose=False)
            record["binding"] = {"sections": sorted(binding["data"])}
        else:
            record["binding"] = {"sections": []}
        record["checks"]["get_binding"] = True

        energy_options = {"model": energy_models} if energy_models else None
        entropy_options = {"model": entropy_models} if entropy_models else None
        decomp_options = {"res_threshold": 0} if decomp_present else None
        analyzer = api.get_ana_data(
            energy_options=energy_options,
            entropy_options=entropy_options,
            decomp_options=decomp_options,
            performance_options={"energy_memory": True, "decomp_memory": True},
            verbose=False,
        )
        record["analyzer"] = {
            section: len(value.get("keys", {}))
            for section, value in analyzer.items()
            if isinstance(value, dict) and "keys" in value
        }
        record["checks"]["get_ana_data_inmemory"] = True

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            alias = API.load_gmxmmpbsa_info(path)
        if not isinstance(alias, API.MMPBSA_API):
            raise TypeError("legacy loader did not return MMPBSA_API")
        record["checks"]["legacy_loader"] = True
        record["legacy_loader_warnings"] = sorted({type(item.message).__name__ for item in caught})
    except Exception as exc:  # keep the full suite running and report the case
        record["errors"].append({"type": type(exc).__name__, "message": str(exc)})

    record["status"] = "PASS" if not record["errors"] else "FAIL"
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="API result files or external result directories")
    parser.add_argument("--output", type=Path, required=True, help="External JSON report path")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    inputs = discover(args.paths)
    if not inputs:
        raise SystemExit("No API result files found")
    records = []
    for path in inputs:
        print(f"[{path.name}] {path}")
        record = validate_result(path)
        records.append(record)
        print(f"     {record['status']}")

    summary = {
        status: sum(record["status"] == status for record in records)
        for status in sorted({record["status"] for record in records})
    }
    report = {"inputs": [str(path) for path in inputs], "summary": summary, "results": records}
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Summary: {summary}")
    print(f"Report: {output}")
    return 0 if summary.get("FAIL", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
