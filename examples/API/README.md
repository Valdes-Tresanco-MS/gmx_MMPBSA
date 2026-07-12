# Python API extraction example

This example includes a compact result file generated from the protein-ligand
single trajectory example:

```text
examples/API/COMPACT_MMXSA_RESULTS.mmxsa
```

Run it from the repository root in an environment where `gmx_MMPBSA` and its
Python dependencies are available:

```bash
python examples/API/extract_api_data.py
```

The script loads the result with `GMXMMPBSA.API.load()`, extracts metadata and
the GB delta `TOTAL` energy, then writes:

```text
examples/API/gb_delta_total.csv
```

Expected terminal output:

```text
Loaded result: examples/API/COMPACT_MMXSA_RESULTS.mmxsa
Frames: 10
Temperature: 298.15
GB delta TOTAL average: -15.017773
GB delta TOTAL SD: 1.686023
GB delta TOTAL SEM: 0.533167
Wrote: examples/API/gb_delta_total.csv
```
