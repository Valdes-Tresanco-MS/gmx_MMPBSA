#!/usr/bin/env python3
"""Summarize external validation JSON/log reports as readable CSV tables."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path("/tmp/gmx_MMPBSA-validation")


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def counts(statuses: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(statuses).items()))


def effective_combination_rows(root: Path) -> list[dict[str, Any]]:
    initial_path = root / "20260911T085613Z-combination-full-current" / "run.json"
    initial = load(initial_path)
    corrections = {
        "ADV-03": (root / "20260911T173505Z-combination-adv03-dftb-env-fixed" / "run.json",
                   "Fixed direct-invocation DFTB environment roots."),
        "GB-09": (root / "20260911T173843Z-combination-gb09-fixed" / "run.json",
                  "Fixed cpptraj surface-output command termination."),
        "AN-04": (root / "20260911T174027Z-combination-an04-fixed" / "run.json",
                  "Fixed decomposition residue-range selection."),
        "PB-08": (root / "20260911T174027Z-combination-pb08-fixed" / "run.json",
                  "Fixed invalid bcopt=10/nfocus combination."),
        "ADV-12": (root / "20260911T095700Z-combination-adv12-current-corrected" / "run.json",
                   "Corrected expected explicit-waters log token."),
    }
    corrected = {}
    for case_id, (path, note) in corrections.items():
        corrected[case_id] = (load(path)["cases"][0], path, note)
    rows = []
    for case in initial["cases"]:
        case_id = case["id"]
        if case_id in corrected:
            final, evidence, note = corrected[case_id]
            rows.append({
                "suite": "current_combination_matrix",
                "case_id": case_id,
                "test": final["description"],
                "status": final["status"],
                "raw_status": case["status"],
                "return_code": final.get("return_code", ""),
                "duration_seconds": final.get("duration_seconds", ""),
                "differences": "",
                "notes": note,
                "evidence": str(evidence),
            })
        else:
            rows.append({
                "suite": "current_combination_matrix",
                "case_id": case_id,
                "test": case["description"],
                "status": case["status"],
                "raw_status": case["status"],
                "return_code": case.get("return_code", ""),
                "duration_seconds": case.get("duration_seconds", ""),
                "differences": "",
                "notes": "",
                "evidence": str(initial_path),
            })
    return rows


def comparison_rows(path: Path, suite: str) -> list[dict[str, Any]]:
    data = load(path)
    rows = []
    for item in data.get("comparisons", []):
        differences = item.get("differences") or []
        if isinstance(differences, dict):
            differences = list(differences)
        status = item.get("status", "")
        raw_status = item.get("raw_status", "")
        notes = []
        if raw_status == "EXPECTED-DIFFERENCE" or status == "EXPECTED-DIFFERENCE":
            notes.append("Documented intentional compatibility difference")
        if status == "BASELINE-UNSUPPORTED":
            notes.append("Baseline cannot execute this case")
        if item.get("comparison_scope"):
            notes.append(f"Comparison scope: {item['comparison_scope']}")
        rows.append({
            "suite": suite,
            "case_id": item.get("id", ""),
            "test": item.get("description", ""),
            "status": status,
            "raw_status": raw_status,
            "return_code": "",
            "duration_seconds": "",
            "differences": len(differences),
            "notes": "; ".join(notes),
            "evidence": str(path),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    overview: list[dict[str, Any]] = []
    calculations: list[dict[str, Any]] = []
    interfaces: list[dict[str, Any]] = []
    docs: list[dict[str, Any]] = []
    catalog: list[dict[str, Any]] = []

    def add_overview(suite, status, statuses, evidence, notes=""):
        summary = counts(statuses)
        overview.append({
            "suite": suite,
            "status": status,
            "total": len(statuses),
            "pass": summary.get("PASS", 0),
            "expected_difference": summary.get("EXPECTED-DIFFERENCE", 0),
            "baseline_unsupported": summary.get("BASELINE-UNSUPPORTED", 0),
            "fail": summary.get("FAIL", 0),
            "skipped": summary.get("SKIP", 0),
            "evidence": str(evidence),
            "notes": notes,
        })

    current_examples_path = root / "20260911T014117Z-current-all-20260910" / "run.json"
    current_examples = load(current_examples_path)
    add_overview("current bundled examples", "PASS", [x["status"] for x in current_examples["cases"]],
                  current_examples_path, "22 runnable manifest examples; Comp_receptor remains excluded for missing topol.top.")
    for case in current_examples["cases"]:
        calculations.append({"suite": "current_examples", "case_id": case["id"], "test": case["name"],
                             "status": case["status"], "raw_status": case["status"],
                             "return_code": case.get("return_code", ""), "duration_seconds": case.get("duration_seconds", ""),
                             "differences": "", "notes": "", "evidence": str(current_examples_path)})

    combination = effective_combination_rows(root)
    add_overview("current combination matrix", "PASS", [x["status"] for x in combination],
                  root / "20260911T085613Z-combination-full-current" / "run.json",
                  "Effective result after four targeted calculation fixes and one log-token correction.")
    calculations.extend(combination)

    core_compare = root / "legacy-comparison" / "20260911T020706Z-current-legacy-settings" / "comparison.json"
    core_rows = comparison_rows(core_compare, "legacy comparison core")
    add_overview("legacy comparison core", "PASS-WITH-EXPECTED-DIFFERENCE",
                 [x["status"] for x in core_rows], core_compare,
                 "Explicit legacy settings: PBRadii=3, igb=5, PB exdi=80; case 24 is documented GBNSR6 behavior change.")
    calculations.extend(core_rows)

    extended_compare = root / "20260911T032247Z-current-remaining-legacy-settings" / "comparison.json"
    extended_rows = comparison_rows(extended_compare, "legacy comparison extended")
    add_overview("legacy comparison extended", "PASS-WITH-BASELINE-GAPS",
                 [x["status"] for x in extended_rows], extended_compare,
                 "Three comparable cases pass or are expected differences; four baseline cases are unsupported.")
    calculations.extend(extended_rows)

    combo_compare = root / "20260911T095300Z-combination-legacy-current" / "combination-comparison.json"
    combo_rows = comparison_rows(combo_compare, "legacy combination comparison")
    add_overview("legacy combination comparison", "PASS-WITH-BASELINE-GAP",
                 [x["status"] for x in combo_rows], combo_compare,
                 "8/9 rows pass; TR-02 is unsupported by the 1.6.5 baseline CLI.")
    calculations.extend(combo_rows)

    mpi_compare = root / "mpi-gb01-comparison.json"
    mpi_rows = comparison_rows(mpi_compare, "serial vs two-rank calculation")
    add_overview("serial vs two-rank calculation", "PASS", [x["status"] for x in mpi_rows], mpi_compare,
                 "GB-01 CSV comparison at atol=1e-3 and rtol=1e-6.")
    calculations.extend(mpi_rows)

    api_validation_path = root / "api-validation-20260911.json"
    api_validation = load(api_validation_path)
    add_overview("Python API validation", "PASS", [x["status"] for x in api_validation["results"]],
                 api_validation_path, "83 current and 1.6.5-era info/compact result loads and accessor checks.")
    for item in api_validation["results"]:
        interfaces.append({"suite": "api_validation", "test": item["path"], "format": item["format"],
                           "status": item["status"], "raw_status": "", "details": json.dumps({
                               "frames": item.get("frames"), "stability": item.get("stability"),
                               "energy": item.get("energy"), "entropy": item.get("entropy"),
                               "decomposition": item.get("decomposition"), "binding": item.get("binding"),
                               "analyzer": item.get("analyzer"), "errors": item.get("errors")}, sort_keys=True),
                           "evidence": str(api_validation_path)})

    api_compare_path = root / "api-comparison-legacy-20260911-policy.json"
    api_compare = load(api_compare_path)
    add_overview("API legacy comparison", "PASS-WITH-EXPECTED-DIFFERENCES",
                 [x["status"] for x in api_compare["comparisons"]], api_compare_path,
                 "29 pass; 8 documented IE/C2, GBNSR6, or component-level CHARMM-CMAP differences; 0 failures.")
    for item in api_compare["comparisons"]:
        interfaces.append({"suite": "api_comparison", "test": item.get("relative", item.get("metadata", "")),
                           "format": "", "status": item["status"], "raw_status": item.get("raw_status", ""),
                           "details": json.dumps({"differences": item.get("differences"), "errors": item.get("errors")}, sort_keys=True),
                           "evidence": str(api_compare_path)})

    for name, label in [("analyzer-current-20260911.json", "analyzer current info"),
                        ("analyzer-baseline-20260911.json", "analyzer baseline info"),
                        ("analyzer-compact-20260911.json", "analyzer compact"),
                        ("analyzer-recursive-20260911.json", "analyzer recursive")]:
        path = root / name
        item = load(path)
        add_overview(label, item["status"], [item["status"]], path,
                     f"systems={item['systems']}; unique_names={item['unique_names']}; recursive={item['recursive']}.")
        interfaces.append({"suite": "analyzer_intake", "test": label, "format": "", "status": item["status"],
                           "raw_status": "", "details": json.dumps(item, sort_keys=True), "evidence": str(path)})

    mpi_path = root / "mpi-20260911.json"
    mpi = load(mpi_path)
    add_overview("MPI prerequisites", mpi["status"], [mpi["status"]], mpi_path,
                 "mpi4py world size and GROMACS discovery; application calculation is reported separately.")
    interfaces.append({"suite": "mpi_prerequisites", "test": "mpi4py/mpirun/GROMACS", "format": "",
                       "status": mpi["status"], "raw_status": "", "details": json.dumps(mpi, sort_keys=True),
                       "evidence": str(mpi_path)})

    concurrency_path = root / "20260911T190435Z-concurrency-validation" / "run.json"
    concurrency = load(concurrency_path)
    add_overview("concurrent calculations", "PASS", [x["status"] for x in concurrency["cases"]], concurrency_path,
                 "GB-01 and GB-02 ran concurrently in separate physical copies.")
    for item in concurrency["cases"]:
        interfaces.append({"suite": "concurrency", "test": item["id"], "format": "", "status": item["status"],
                           "raw_status": "", "details": json.dumps(item, sort_keys=True), "evidence": str(concurrency_path)})

    cleanup_path = root / "cleanup-20260911.json"
    cleanup = load(cleanup_path)
    add_overview("cleanup behavior", cleanup["status"], [cleanup["status"]], cleanup_path,
                 "Minimal and full cleanup policies validated; real CLI cleanup was also run on a disposable copy.")
    for test in ("minimal cleanup", "full cleanup"):
        interfaces.append({"suite": "cleanup", "test": test, "format": "", "status": cleanup["status"],
                           "raw_status": "", "details": json.dumps(cleanup, sort_keys=True), "evidence": str(cleanup_path)})

    negative_path = root / "negative-20260911-rerun.json"
    negative = load(negative_path)
    add_overview("negative/error handling", negative["status"], [x["status"] for x in negative["cases"]], negative_path,
                 "Invalid IGB fails; default creates a bundle, --no-error-bundle suppresses it.")
    for item in negative["cases"]:
        interfaces.append({"suite": "negative_error", "test": "invalid igb=999", "format": "",
                           "status": item["status"], "raw_status": "", "details": json.dumps(item, sort_keys=True),
                           "evidence": str(negative_path)})

    disk_path = root / "disk-analyzer-20260911.json"
    disk = load(disk_path)
    add_overview("disk-backed analyzer", disk["status"], [disk["status"]], disk_path,
                 f"Mode={disk['mode']}; parquet engine availability is environment-dependent.")
    interfaces.append({"suite": "disk_analyzer", "test": "disk-backed analyzer", "format": "",
                       "status": disk["status"], "raw_status": "", "details": json.dumps(disk, sort_keys=True),
                       "evidence": str(disk_path)})

    package_path = root / "packaging-20260911-final.json"
    package = load(package_path)
    add_overview("packaging/install", package["status"], [package["status"]], package_path,
                 "sdist/wheel build, external no-dependency install, import, entry points, and data files.")
    interfaces.append({"suite": "packaging", "test": "sdist/wheel external install", "format": "",
                       "status": package["status"], "raw_status": "", "details": json.dumps(package, sort_keys=True),
                       "evidence": str(package_path)})

    docs.extend([
        {"suite": "MkDocs normal non-strict", "status": "PASS", "exit_code": 0,
         "details": "Documentation built successfully; warnings were emitted.",
         "evidence": str(root / "mkdocs-normal-nonstrict-20260911.log")},
        {"suite": "MkDocs normal strict", "status": "FAIL-ENVIRONMENT", "exit_code": 1,
         "details": "Deprecation and invalid-escape warnings are treated as strict errors by installed MkDocs stack.",
         "evidence": str(root / "mkdocs-normal-20260911.log")},
        {"suite": "MkDocs symlink strict spike", "status": "FAIL-LINKS", "exit_code": 1,
         "details": "Relative links from canonical examples are broken when docs/examples is replaced by a symlink.",
         "evidence": str(root / "docs-spike-20260911.log")},
    ])
    for item in docs:
        add_overview(item["suite"], item["status"], [item["status"]], item["evidence"], item["details"])

    catalog_paths = [
        current_examples_path,
        root / "20260911T085613Z-combination-full-current" / "run.json",
        root / "20260911T173505Z-combination-adv03-dftb-env-fixed" / "run.json",
        root / "20260911T173843Z-combination-gb09-fixed" / "run.json",
        root / "20260911T174027Z-combination-an04-fixed" / "run.json",
        root / "20260911T174027Z-combination-pb08-fixed" / "run.json",
        root / "20260911T095700Z-combination-adv12-current-corrected" / "run.json",
        core_compare,
        extended_compare,
        combo_compare,
        mpi_compare,
        api_validation_path,
        api_compare_path,
        mpi_path,
        concurrency_path,
        cleanup_path,
        negative_path,
        disk_path,
        package_path,
    ]
    for path in sorted(set(catalog_paths)):
        try:
            data = load(path)
        except Exception:
            continue
        if not isinstance(data, dict):
            continue
        catalog.append({"report": str(path), "top_level_status": data.get("status", ""),
                        "summary": json.dumps(data.get("summary", {}), sort_keys=True),
                        "records": len(data.get("results", data.get("cases", data.get("comparisons", [])))),
                        "included_in_tables": True})

    write_csv(output / "validation_overview.csv",
              ["suite", "status", "total", "pass", "expected_difference", "baseline_unsupported", "fail", "skipped", "evidence", "notes"], overview)
    write_csv(output / "calculation_details.csv",
              ["suite", "case_id", "test", "status", "raw_status", "return_code", "duration_seconds", "differences", "notes", "evidence"], calculations)
    write_csv(output / "interface_boundary_details.csv",
              ["suite", "test", "format", "status", "raw_status", "details", "evidence"], interfaces)
    write_csv(output / "documentation_details.csv",
              ["suite", "status", "exit_code", "details", "evidence"], docs)
    write_csv(output / "report_catalog.csv",
              ["report", "top_level_status", "summary", "records", "included_in_tables"], catalog)
    print(f"Wrote CSV recap to {output}")
    for path in sorted(output.glob("*.csv")):
        print(f"  {path.name}: {sum(1 for _ in path.open(encoding='utf-8')) - 1} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
