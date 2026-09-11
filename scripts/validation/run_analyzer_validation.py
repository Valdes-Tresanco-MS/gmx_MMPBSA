#!/usr/bin/env python3
"""Exercise non-interactive analyzer discovery and Qt initialization.

This runner uses the real analyzer file-discovery code and Qt initialization
dialog in offscreen mode. It does not enter the GUI event loop, open a window,
or launch calculations.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from PyQt6.QtWidgets import QApplication
except ImportError:  # pragma: no cover - exercised in PyQt5 environments
    from PyQt5.QtWidgets import QApplication

from GMXMMPBSA.analyzer.dialogs import InitDialog
from GMXMMPBSA.analyzer.utils import get_files
from GMXMMPBSA.commandlineparser import anaparser


def discover(paths: list[Path], recursive: bool) -> list[Path]:
    arguments = ["-f", *(str(path.expanduser().resolve()) for path in paths)]
    if recursive:
        arguments.append("--recursive")
    parser_args = anaparser.parse_args(arguments)
    return get_files(parser_args)


def validate(paths: list[Path], recursive: bool) -> dict:
    files = discover(paths, recursive)
    qt_app = QApplication.instance() or QApplication([])
    dialog = InitDialog(None)
    try:
        dialog.get_files_info(files)
        names = [dialog.f_item.child(index).text(2) for index in range(dialog.f_item.childCount())]
        result = {
            "requested": [str(path.expanduser().resolve()) for path in paths],
            "recursive": recursive,
            "files": [str(path) for path in files],
            "systems": len(names),
            "names": names,
            "unique_names": len(names) == len(set(names)),
        }
    finally:
        dialog.close()
        qt_app.processEvents()
    result["status"] = "PASS" if result["systems"] and result["unique_names"] else "FAIL"
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate(args.paths, args.recursive)
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Report: {output}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
