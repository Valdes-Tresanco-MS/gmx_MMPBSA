#!/usr/bin/env python3
"""Build and inspect wheel/sdist artifacts from an external checkout copy."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def validate(source: Path, python: Path, results_root: Path) -> dict:
    root = results_root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    build_root = Path(tempfile.mkdtemp(prefix="gmx_MMPBSA-package-", dir=root))
    checkout = build_root / "checkout"
    target = build_root / "installed"
    shutil.copytree(source.expanduser().resolve(), checkout, symlinks=True,
                    ignore=shutil.ignore_patterns("build", "*.egg-info"))
    # Keep the Git metadata in this disposable copy so Versioneer resolves the
    # same valid PEP 440 version as the working tree without writing artifacts
    # into the user's checkout.
    build = subprocess.run(
        [str(python), "setup.py", "sdist", "bdist_wheel"],
        cwd=checkout, capture_output=True, text=True, check=False,
    )
    dist = checkout / "dist"
    artifacts = sorted(path.name for path in dist.glob("*") if path.is_file())
    wheel = next((path for path in dist.glob("*.whl")), None)
    install_return_code = None
    import_return_code = None
    import_output = ""
    if wheel:
        install = subprocess.run(
            [str(python), "-m", "pip", "install", "--no-deps", "--target", str(target), str(wheel)],
            cwd=checkout, capture_output=True, text=True, check=False,
        )
        install_return_code = install.returncode
        if install.returncode == 0:
            probe = subprocess.run(
                [str(python), "-I", "-c",
                 "import sys; sys.path.insert(0, %r); import GMXMMPBSA, importlib.metadata as m; "
                 "print(m.version('gmx_MMPBSA')); print(GMXMMPBSA.__file__)" % str(target)],
                cwd=build_root, capture_output=True, text=True, check=False,
            )
            import_return_code = probe.returncode
            import_output = (probe.stdout + probe.stderr).strip()
    required_paths = [
        target / "GMXMMPBSA" / "data" / "gmx_MMPBSA_test_manifest.json",
        target / "GMXMMPBSA" / "data" / "xvv_files" / "tip3p.xvv",
    ]
    record = {
        "root": str(build_root),
        "artifacts": artifacts,
        "build_return_code": build.returncode,
        "build_output_tail": (build.stdout + build.stderr)[-2000:],
        "wheel": str(wheel) if wheel else None,
        "install_return_code": install_return_code,
        "import_return_code": import_return_code,
        "import_output": import_output,
        "required_paths": {str(path.relative_to(target)): path.is_file() for path in required_paths},
    }
    record["status"] = "PASS" if (
        build.returncode == 0 and bool(wheel) and install_return_code == 0
        and import_return_code == 0 and all(record["required_paths"].values())
    ) else "FAIL"
    return record


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, default=Path("/tmp/gmx_MMPBSA-validation"))
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.source, args.python.expanduser().resolve(), args.results_root)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
