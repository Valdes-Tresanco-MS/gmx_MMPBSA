---
template: main.html
title: Python API example
---

# Python API example

This example shows how to use the Python API on a real `gmx_MMPBSA` result. It
loads the compact output included in the API example folder and extracts the GB
binding `TOTAL` energy for the `delta` system.

## Input result

The script reads:

```text
examples/API/COMPACT_MMXSA_RESULTS.mmxsa
```

This compact `.mmxsa` file was generated from the protein-ligand single
trajectory example. It contains the parsed result payload without requiring the
calculation to be rerun.

## Script

The runnable script is available at `examples/API/extract_api_data.py`:

```python
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
```

## Run it

From the repository root, run:

```bash
python examples/API/extract_api_data.py
```

In the supported project environment, this prints:

```text
Loaded result: examples/API/COMPACT_MMXSA_RESULTS.mmxsa
Frames: 10
Temperature: 298.15
GB delta TOTAL average: -15.017773
GB delta TOTAL SD: 1.686023
GB delta TOTAL SEM: 0.533167
Wrote: examples/API/gb_delta_total.csv
```

## Extracted CSV

The run writes `examples/API/gb_delta_total.csv`:

```csv
Frames,GB delta TOTAL
1,-17.07779711999989
2,-14.646319519999896
3,-15.423147679999886
4,-17.63201023999989
5,-15.716009839999954
6,-11.549836800000074
7,-13.420182960000057
8,-14.013321119999999
9,-15.962333759999742
10,-14.736771999999775
Average,-15.017773103999914
SD,1.6860229037427301
SEM,0.5331672563037857
```
