---
template: main.html
title: amber_MMPBSA
---

# amber_MMPBSA

`amber_MMPBSA` is an independent command-line module included with `gmx_MMPBSA` for running end-state free energy
calculations directly from native AMBER input files.

Unlike `gmx_MMPBSA`, which prepares and processes GROMACS files before building the AMBER-compatible calculation
workflow, `amber_MMPBSA` starts from AMBER topology, coordinate, trajectory, and mask inputs. This makes it useful
when the system has already been prepared in the AMBER ecosystem and no GROMACS files are needed.

Although `amber_MMPBSA` is an independent module, the calculation and analysis features available in `gmx_MMPBSA`
also work with `amber_MMPBSA` workflows. This includes the same input file options for supported calculation types,
output generation, decomposition analysis, entropy calculations, CSV export, result rewriting, and post-processing
with `gmx_MMPBSA_ana`.

## Input files

`amber_MMPBSA` uses native AMBER files:

| Input | Option | Description |
|:------|:------:|:------------|
| Complex topology | `-cp` | AMBER topology file for the complex |
| Complex structure | `-cs` | AMBER coordinate file for the complex. Supported formats include `*.pdb`, `*.inpcrd`, and `*.rst7` |
| Complex trajectory | `-ct` | AMBER trajectory file. Supported format: `*.mdcrd` |
| Complex masks | `-cm` | Receptor and ligand masks from the complex |
| Receptor topology | `-rp` | AMBER topology file for the receptor |
| Ligand topology | `-lp` | AMBER topology file for the ligand |

## Example

A basic single-trajectory command is:

``` bash
amber_MMPBSA -O -i mmpbsa.in \
  -cp ras-raf_complex.prmtop \
  -cs ras-raf_complex.inpcrd \
  -ct prod_complex.mdcrd \
  -cm ":1-166" ":167-242" \
  -o FINAL_RESULTS_MMPBSA.dat \
  -eo FINAL_RESULTS_MMPBSA.csv
```

See the [AMBER input files example](examples/AMBER/README.md) for a complete runnable example.
