#!/usr/bin/env python3
"""Extract selected values from a real gmx_MMPBSA result with the Python API."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from GMXMMPBSA import API  # noqa: E402


RESULT_FILE = Path(__file__).with_name("COMPACT_MMXSA_RESULTS.mmxsa")
OUTPUT_CSV = Path(__file__).with_name("gb_delta_total.csv")
GENERATED_PDB = Path(__file__).with_name("_GMXMMPBSA_COM_FIXED.pdb")


def main() -> None:
    api = API.load(RESULT_FILE)

    info = api.get_info()
    inputs = api.get_input()
    energy = api.get_energy(model=("gb",), mol=("delta",), term=("TOTAL",))

    summary = energy["summary"]["normal"]["gb"]["delta"]["TOTAL"]
    gb_delta = energy["data"]["normal"]["gb"]["delta"]
    total = gb_delta["TOTAL"]

    total.to_frame("GB delta TOTAL").to_csv(OUTPUT_CSV)

    print("Loaded result:", RESULT_FILE.relative_to(ROOT))
    print("Frames:", info["numframes"])
    print("Temperature:", inputs["general"]["temperature"])
    print("GB delta TOTAL average:", f"{summary['Average']:.6f}")
    print("GB delta TOTAL SD:", f"{summary['SD']:.6f}")
    print("GB delta TOTAL SEM:", f"{summary['SEM']:.6f}")
    print("Wrote:", OUTPUT_CSV.relative_to(ROOT))

    if GENERATED_PDB.exists():
        GENERATED_PDB.unlink()


if __name__ == "__main__":
    main()
