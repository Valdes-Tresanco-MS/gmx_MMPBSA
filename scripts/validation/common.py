"""Shared machinery for isolated calculation validation.

The runners in this package deliberately execute copied examples under an
external results directory.  They do not modify the checkout or reuse output
files from a previous run.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO_ROOT / "GMXMMPBSA" / "data" / "gmx_MMPBSA_test_manifest.json"
COMPARISON_POLICY_PATH = Path(__file__).with_name("legacy_comparison_policy.json")
DEFAULT_EXAMPLES = REPO_ROOT / "examples"
DEFAULT_RESULTS_ROOT = Path(tempfile.gettempdir()) / "gmx_MMPBSA-validation"


@dataclass(frozen=True)
class Case:
    id: int
    name: str
    path: str
    workdir: str
    input: str
    executable: str
    command_args: tuple[str, ...]
    expected_outputs: tuple[str, ...]
    slow: bool
    requires: tuple[str, ...]


def load_cases(manifest_path: Path = MANIFEST_PATH) -> tuple[dict[str, Any], dict[int, Case]]:
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = {
        int(case_id): Case(
            id=int(case_id),
            name=entry["name"],
            path=entry["path"],
            workdir=entry["workdir"],
            input=entry["input"],
            executable=entry["executable"],
            command_args=tuple(entry["command_args"]),
            expected_outputs=tuple(entry["expected_outputs"]),
            slow=bool(entry.get("slow", False)),
            requires=tuple(entry.get("requires", ())),
        )
        for case_id, entry in raw["tests"].items()
    }
    return raw, cases


def load_comparison_policies(policy_path: Path = COMPARISON_POLICY_PATH) -> dict[int, dict[str, Any]]:
    """Load explicit current-vs-1.6.5 comparison exceptions."""

    raw = json.loads(policy_path.read_text(encoding="utf-8"))
    return {int(case_id): dict(policy) for case_id, policy in raw.get("cases", {}).items()}


def resolve_selectors(raw: dict[str, Any], selectors: Iterable[str]) -> list[int]:
    """Resolve manifest IDs, suite IDs, aliases, and the all-suite alias."""

    tests = {int(key) for key in raw["tests"]}
    suites = raw["suites"]
    suite_by_id = {str(value["id"]): name for name, value in suites.items()}
    aliases = raw.get("aliases", {})
    resolved: list[int] = []
    seen: set[int] = set()

    def expand(token: str) -> list[int]:
        token = str(token).strip()
        if token in suite_by_id:
            return list(suites[suite_by_id[token]]["tests"])
        if token == "101":
            return list(suites["all"]["tests"])
        if token in aliases:
            target = aliases[token]
            if isinstance(target, int):
                return [target]
            if target in suites:
                return list(suites[target]["tests"])
            raise ValueError(f"Unknown manifest alias target: {target!r}")
        if token.isdigit() and int(token) in tests:
            return [int(token)]
        raise ValueError(f"Invalid calculation selector: {token!r}")

    for selector in selectors:
        for case_id in expand(selector):
            if case_id not in tests:
                raise ValueError(f"Manifest selector resolved to unknown test {case_id}")
            if case_id not in seen:
                seen.add(case_id)
                resolved.append(case_id)
    if not resolved:
        raise ValueError("At least one calculation selector is required")
    return resolved


def timestamped_results_root(root: Path | None, label: str) -> Path:
    base = (root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    result = base / f"{stamp}-{label}"
    result.mkdir(parents=True, exist_ok=False)
    return result


def environment_for(python: Path, source_root: Path) -> dict[str, str]:
    """Build an environment tied to the requested Python environment/source."""

    python = python.expanduser().resolve()
    source_root = source_root.expanduser().resolve()
    env = os.environ.copy()
    env_bin = python.parent
    env["PATH"] = os.pathsep.join((str(env_bin), env.get("PATH", "")))
    # Direct invocation of an environment's interpreter does not run conda's
    # activation hooks. AmberTools therefore needs its data roots made
    # explicit, especially for QM/MM DFTB Slater-Koster files. Preserve an
    # intentionally configured custom DFTB_PREFIX.
    env_root = python.parent.parent
    env["AMBERHOME"] = env.get("AMBERHOME") or str(env_root)
    env["DFTB_PREFIX"] = env.get("DFTB_PREFIX") or str(env_root / "dat" / "slko")
    existing_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = os.pathsep.join(
        value for value in (str(source_root), existing_pythonpath) if value
    )
    return env


def entrypoint_command(
    case: Case,
    python: Path,
    source_root: Path,
    *,
    legacy_cli_compat: bool = False,
) -> list[str]:
    if case.executable == "gmx_MMPBSA":
        function = "gmxmmpbsa"
    elif case.executable == "amber_MMPBSA":
        function = "gmxmmpbsa_amber"
    else:
        raise ValueError(f"Unsupported manifest executable: {case.executable}")

    bootstrap = f"from GMXMMPBSA.app import {function}; {function}()"
    args = list(case.command_args)
    if legacy_cli_compat:
        for option in ("-rg", "-lg"):
            try:
                index = args.index(option)
            except ValueError:
                continue
            if index + 1 < len(args) and args[index + 1] == "System":
                args[index + 1] = "0"
    command = [str(python), "-c", bootstrap, *args]
    if "-nogui" not in command and "--nogui" not in command:
        command.append("-nogui")
    return command


def launched_command(command: list[str], ranks: int, env: dict[str, str]) -> list[str]:
    if ranks < 1:
        raise ValueError("MPI ranks must be at least 1")
    if ranks == 1:
        return command
    mpirun = shutil.which("mpirun", path=env.get("PATH"))
    if not mpirun:
        raise FileNotFoundError("mpirun is required when --ranks is greater than 1")
    return [mpirun, "-np", str(ranks), *command]


def copy_examples(examples: Path, destination: Path) -> Path:
    examples = examples.expanduser().resolve()
    if not examples.is_dir():
        raise FileNotFoundError(f"Examples directory does not exist: {examples}")
    copied = destination / "examples"
    shutil.copytree(examples, copied, symlinks=True, ignore=_generated_example_files)
    return copied


def _generated_example_files(directory: str, names: list[str]) -> set[str]:
    """Exclude ignored calculation products that may exist in a dirty example tree."""

    ignored = {
        ".gmx_mmpbsa_temp",
        "FINAL_RESULTS_MMPBSA.dat",
        "FINAL_RESULTS_MMPBSA.csv",
        "FINAL_DECOMP_MMPBSA.dat",
        "FINAL_DECOMP_MMPBSA.csv",
        "gmx_MMPBSA.log",
        "COM.prmtop",
        "REC.prmtop",
        "LIG.prmtop",
        "COM_traj_0.xtc",
    }
    return {
        name for name in names
        if name in ignored
        or name.endswith(".log")
        or (name.startswith("_GMXMMPBSA_") and name != "_GMXMMPBSA_COM_FIXED.pdb")
    }


def _fingerprint(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.is_file():
        return {"exists": False}
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    return {"exists": True, "size": stat.st_size, "sha256": digest.hexdigest()}


def apply_legacy_settings(examples_dir: Path, case: Case) -> None:
    """Make a copied input explicit with the 1.6.5 scientific defaults.

    This is only for compatibility runs.  It intentionally changes the copied
    input, never the canonical example, and is recorded in the run metadata.
    """

    input_path = examples_dir / case.workdir / case.input
    text = input_path.read_text(encoding="utf-8")
    text = _set_namelist_value(text, "general", "PBRadii", "3")
    if _has_namelist(text, "gb"):
        text = _set_namelist_value(text, "gb", "igb", "5")
    if _has_namelist(text, "pb"):
        text = _set_namelist_value(text, "pb", "exdi", "80")
    input_path.write_text(text, encoding="utf-8")


def _section_bounds(text: str, section: str) -> tuple[int, int] | None:
    match = __import__("re").search(rf"(?im)^\s*&{section}\b", text)
    if not match:
        return None
    end = __import__("re").search(r"(?m)^\s*/\s*(?:#.*)?$", text[match.end():])
    if not end:
        raise ValueError(f"Unterminated &{section} namelist")
    return match.end(), match.end() + end.start()


def _has_namelist(text: str, section: str) -> bool:
    return _section_bounds(text, section) is not None


def _set_namelist_value(text: str, section: str, key: str, value: str) -> str:
    import re

    bounds = _section_bounds(text, section)
    if bounds is None:
        return text
    start, end = bounds
    body = text[start:end]
    pattern = re.compile(rf"(\b{re.escape(key)}\s*=\s*)[^,\n/]+", re.IGNORECASE)
    if pattern.search(body):
        body = pattern.sub(rf"\g<1>{value}", body)
    else:
        body = body.rstrip() + f"\n{key}={value},\n"
    return text[:start] + body + text[end:]


def _read_csv_rows(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.reader(handle))


def _select_csv_section(rows: list[list[str]], section: str) -> list[list[str]]:
    """Return one named result section, including its header row."""

    try:
        section_start = next(
            index for index, row in enumerate(rows)
            if row and row[0].strip() == section
        )
    except StopIteration:
        return []

    data_start = section_start + 1
    while data_start < len(rows) and not any(value.strip() for value in rows[data_start]):
        data_start += 1
    data_end = data_start
    while data_end < len(rows) and any(value.strip() for value in rows[data_end]):
        data_end += 1
    return rows[data_start:data_end]


def compare_csv(
    left: Path,
    right: Path,
    *,
    atol: float,
    rtol: float,
    section: str | None = None,
) -> list[str]:
    errors: list[str] = []
    left_rows = _read_csv_rows(left)
    right_rows = _read_csv_rows(right)
    if section is not None:
        left_rows = _select_csv_section(left_rows, section)
        right_rows = _select_csv_section(right_rows, section)
        if not left_rows or not right_rows:
            return [f"section {section!r} is missing from one or both CSV files"]
    if len(left_rows) != len(right_rows):
        errors.append(f"row count differs: {len(left_rows)} != {len(right_rows)}")
        return errors
    for row_number, (left_row, right_row) in enumerate(zip(left_rows, right_rows), start=1):
        if len(left_row) != len(right_row):
            errors.append(f"row {row_number} column count differs")
            continue
        for column, (left_value, right_value) in enumerate(zip(left_row, right_row), start=1):
            try:
                left_number = float(left_value)
                right_number = float(right_value)
            except ValueError:
                if left_value.strip() != right_value.strip():
                    errors.append(f"row {row_number}, column {column}: {left_value!r} != {right_value!r}")
            else:
                if not math.isclose(left_number, right_number, rel_tol=rtol, abs_tol=atol):
                    errors.append(
                        f"row {row_number}, column {column}: {left_number} != {right_number}"
                    )
    return errors


def run_matrix(
    *,
    label: str,
    source_root: Path,
    examples: Path,
    python: Path,
    cases: dict[int, Case],
    selected_ids: list[int],
    results_root: Path | None,
    ranks: int,
    legacy_settings: bool = False,
    legacy_cli_compat: bool = False,
    timeout: int | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    result_root = timestamped_results_root(results_root, label)
    copied_examples = copy_examples(examples, result_root)
    env = environment_for(python, source_root)
    metadata: dict[str, Any] = {
        "label": label,
        "source_root": str(source_root.expanduser().resolve()),
        "examples": str(examples.expanduser().resolve()),
        "python": str(python.expanduser().resolve()),
        "ranks": ranks,
        "legacy_settings": legacy_settings,
        "legacy_cli_compat": legacy_cli_compat,
        "dry_run": dry_run,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [],
    }

    for case_id in selected_ids:
        case = cases[case_id]
        workdir = copied_examples / case.workdir
        if legacy_settings:
            apply_legacy_settings(copied_examples, case)
        command = entrypoint_command(
            case, python, source_root, legacy_cli_compat=legacy_cli_compat
        )
        command = launched_command(command, ranks, env)
        log_file = result_root / "cases" / f"{case_id:02d}.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        before = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
        case_record: dict[str, Any] = {
            "id": case.id,
            "name": case.name,
            "workdir": case.workdir,
            "command": command,
            "log": str(log_file),
            "expected_outputs": list(case.expected_outputs),
            "before": before,
        }
        print(f"[{case.id:02d}] {case.name}")
        print("     " + " ".join(command))
        if dry_run:
            case_record["status"] = "DRY-RUN"
            metadata["cases"].append(case_record)
            continue

        started = time.monotonic()
        try:
            with log_file.open("w", encoding="utf-8") as log_handle:
                completed = subprocess.run(
                    command,
                    cwd=workdir,
                    env=env,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
            case_record["return_code"] = completed.returncode
        except subprocess.TimeoutExpired:
            case_record["return_code"] = None
            case_record["status"] = "TIMEOUT"
        except Exception as exc:  # capture launch failures per case
            case_record["return_code"] = None
            case_record["status"] = "LAUNCH-ERROR"
            case_record["error"] = repr(exc)

        case_record["duration_seconds"] = round(time.monotonic() - started, 3)
        after = {name: _fingerprint(workdir / name) for name in case.expected_outputs}
        case_record["after"] = after
        missing = [name for name, fingerprint in after.items() if not fingerprint["exists"]]
        unchanged = [
            name for name in case.expected_outputs
            if before[name].get("exists") and before[name] == after[name]
        ]
        case_record["missing_outputs"] = missing
        case_record["unchanged_outputs"] = unchanged
        if "status" not in case_record:
            if case_record["return_code"] != 0:
                case_record["status"] = "FAIL"
            elif missing:
                case_record["status"] = "FAIL-MISSING-OUTPUT"
            elif unchanged:
                case_record["status"] = "FAIL-STALE-OUTPUT"
            else:
                case_record["status"] = "PASS"
        metadata["cases"].append(case_record)
        print(f"     {case_record['status']} ({case_record.get('duration_seconds', 0):.1f}s)")

    metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["summary"] = {
        status: sum(record.get("status") == status for record in metadata["cases"])
        for status in sorted({record.get("status") for record in metadata["cases"]})
    }
    metadata_path = result_root / "run.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Results: {result_root}")
    return metadata


def metadata_exit_code(metadata: dict[str, Any]) -> int:
    return 0 if all(record.get("status") in {"PASS", "DRY-RUN"} for record in metadata["cases"]) else 1


def csv_comparison_for_cases(
    left_metadata: dict[str, Any],
    right_metadata: dict[str, Any],
    *,
    atol: float,
    rtol: float,
    policies: dict[int, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    policies = load_comparison_policies() if policies is None else policies
    left_by_id = {record["id"]: record for record in left_metadata["cases"]}
    right_by_id = {record["id"]: record for record in right_metadata["cases"]}
    comparisons = []
    for case_id in sorted(set(left_by_id) | set(right_by_id)):
        left = left_by_id.get(case_id)
        right = right_by_id.get(case_id)
        policy = policies.get(case_id, {})
        record: dict[str, Any] = {"id": case_id, "status": "SKIP"}
        if not left or not right:
            record["reason"] = "case missing from one run"
        elif left.get("status") != "PASS" or right.get("status") != "PASS":
            baseline_policy = policy.get("baseline", {})
            baseline_log = Path(right["log"])
            patterns = baseline_policy.get("log_patterns", [])
            baseline_text = baseline_log.read_text(encoding="utf-8", errors="replace") if baseline_log.is_file() else ""
            if (
                left.get("status") == "PASS"
                and right.get("status") != "PASS"
                and baseline_policy.get("status") == "unsupported"
                and patterns
                and any(pattern in baseline_text for pattern in patterns)
            ):
                record["status"] = "BASELINE-UNSUPPORTED"
                record["reason"] = baseline_policy.get("reason", "baseline does not support this case")
                record["baseline_status"] = right.get("status")
            else:
                record["status"] = "FAIL"
                record["reason"] = f"run status: {left.get('status')} vs {right.get('status')}"
        else:
            left_root = Path(left_metadata["cases"][0]["log"]).parents[1]
            right_root = Path(right_metadata["cases"][0]["log"]).parents[1]
            raw_errors = []
            for output in left["expected_outputs"]:
                left_path = left_root / "examples" / left["workdir"] / output
                right_path = right_root / "examples" / right["workdir"] / output
                if output.endswith(".csv"):
                    raw_errors.extend(
                        f"{output}: {error}"
                        for error in compare_csv(left_path, right_path, atol=atol, rtol=rtol)
                    )

            scope = policy.get("csv_scope")
            scoped_errors = raw_errors
            if scope == "delta":
                scoped_errors = []
                for output in left["expected_outputs"]:
                    left_path = left_root / "examples" / left["workdir"] / output
                    right_path = right_root / "examples" / right["workdir"] / output
                    if output.endswith(".csv"):
                        scoped_errors.extend(
                            f"{output}: {error}"
                            for error in compare_csv(
                                left_path,
                                right_path,
                                atol=atol,
                                rtol=rtol,
                                section="Delta Energy Terms",
                            )
                        )
                record["comparison_scope"] = "Delta Energy Terms"
                record["raw_status"] = "PASS" if not raw_errors else "EXPECTED-DIFFERENCE"
                record["raw_differences"] = raw_errors[:50]

            if policy.get("outcome") == "expected-difference" and raw_errors:
                record["status"] = "EXPECTED-DIFFERENCE"
                record["reason"] = policy.get("reason", "documented post-1.6.5 behavior change")
                record["differences"] = raw_errors[:50]
            elif scoped_errors:
                record["status"] = "FAIL"
                record["differences"] = scoped_errors[:50]
            else:
                record["status"] = "PASS"
                record["differences"] = []
            if policy.get("reason") and "reason" not in record:
                record["reason"] = policy["reason"]
        comparisons.append(record)
    return {"comparisons": comparisons}
