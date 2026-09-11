#!/usr/bin/env python3
"""Check MPI prerequisites without confusing them with an app MPI gate.

The probe launches the selected interpreter with ``mpirun`` and verifies that
the requested world size is visible to ``mpi4py``.  It also records the
GROMACS build's MPI mode.  A successful probe is runtime evidence for the
launcher and Python MPI layer only; an application calculation should be run
through ``run_combination_matrix.py --ranks N`` separately.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

try:
    from .common import REPO_ROOT, environment_for
except ImportError:
    from common import REPO_ROOT, environment_for


PROBE = (
    "from mpi4py import MPI; "
    "print(f'rank={MPI.COMM_WORLD.Get_rank()} size={MPI.COMM_WORLD.Get_size()}', flush=True)"
)


def validate(python: Path, source: Path, ranks: int) -> dict:
    if ranks < 2:
        raise ValueError("MPI validation requires at least two ranks")
    python = python.expanduser().resolve()
    env = environment_for(python, source)
    record = {
        "python": str(python),
        "source": str(source.expanduser().resolve()),
        "requested_ranks": ranks,
        "checks": {},
        "errors": [],
        "application_calculation": "NOT_RUN",
    }

    mpirun = shutil.which("mpirun", path=env.get("PATH"))
    record["mpirun"] = mpirun
    if not mpirun:
        record["errors"].append("mpirun was not found on the selected environment PATH")
    else:
        command = [mpirun, "-np", str(ranks), str(python), "-c", PROBE]
        record["probe_command"] = command
        completed = subprocess.run(
            command, env=env, cwd=source, capture_output=True, text=True, check=False
        )
        record["probe_return_code"] = completed.returncode
        record["probe_stdout"] = completed.stdout.strip()
        record["probe_stderr"] = completed.stderr.strip()
        observed = sorted(
            {int(size) for size in re.findall(r"\bsize=(\d+)\b", completed.stdout)}
        )
        record["observed_world_sizes"] = observed
        record["checks"]["mpi4py_world"] = completed.returncode == 0 and observed == [ranks]
        if not record["checks"]["mpi4py_world"]:
            record["errors"].append(
                f"mpi4py did not report the requested world size {ranks}: {observed}"
            )

    gmx = shutil.which("gmx", path=env.get("PATH"))
    record["gmx"] = gmx
    if not gmx:
        record["errors"].append("gmx was not found on the selected environment PATH")
    else:
        version = subprocess.run(
            [gmx, "--version"], env=env, cwd=source, capture_output=True, text=True, check=False
        )
        record["gmx_return_code"] = version.returncode
        version_text = version.stdout + version.stderr
        record["gmx_version_lines"] = [
            line.strip() for line in version_text.splitlines()
            if line.strip().startswith(("GROMACS version:", "MPI library:", "MPI version:"))
        ]
        record["checks"]["gromacs"] = version.returncode == 0
        if version.returncode:
            record["errors"].append("gmx --version failed")

    record["status"] = "PASS" if not record["errors"] else "FAIL"
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=REPO_ROOT)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.python, args.source, args.ranks)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
