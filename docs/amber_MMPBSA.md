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

`amber_MMPBSA` reuses the same calculation engine, input namelists (where supported), output writers, and
`gmx_MMPBSA_ana` post-processing as `gmx_MMPBSA`. It is **not** a full feature-parity replacement: several GROMACS-only
workflows are unsupported. See [Supported and unsupported features](#supported-and-unsupported-features) below.

## Input files

`amber_MMPBSA` uses native AMBER files:

| Input | Option | Required | Description |
|:------|:------:|:--------:|:------------|
| Complex topology | `-cp` | yes | AMBER topology file for the complex |
| Complex trajectory | `-ct` | yes | Trajectory file readable by cpptraj with the supplied AMBER topology. Supported formats include `*.mdcrd`, `*.nc`, `*.crd`, `*.rst7`, `*.inpcrd`, `*.xtc`, `*.trr`, `*.pdb`, `*.gro`, and `*.dcd` |
| Complex masks | `-cm` | yes | Receptor and ligand masks from the complex (Amber residue-number masks, including non-contiguous ranges) |
| Receptor topology | `-rp` | no | AMBER topology file for the receptor (otherwise built from the complex) |
| Ligand topology | `-lp` | no | AMBER topology file for the ligand (otherwise built from the complex) |

For explicit-water calculations, `-cp` must be the original solvated AMBER topology and `-ct` must contain the
matching solvated trajectory. The receptor and ligand masks in `-cm` must select solute residues only; solvent and ions
are removed during setup.

### Native AMBER GB radii

Native AMBER topologies already contain the per-atom `RADII` and `SCREEN` values, together with the `RADIUS_SET`
metadata written during `tleap` preparation. `amber_MMPBSA` preserves these values for the generated working
topologies; the input-file `PBRadii` setting does not replace them. Choose the radius set when building the AMBER
topology, for example:

```text
set default PBradii mbondi3
saveamberparm complex complex.prmtop
```

For GB calculations, the conventional pairings are `igb=1`/`mbondi`, `igb=2` or `5`/`mbondi2`, `igb=7`/`bondi`, and
`igb=8`/`mbondi3`. If a recognized `RADIUS_SET` differs from that pairing, `amber_MMPBSA` prints a warning but still
uses the radii in the topology. This is intentional: unusual combinations may be valid user choices, so they are
not changed or rejected automatically. If the warning is not intentional, rebuild the topology with the desired
`PBradii` in `tleap`.

## Supported and unsupported features

**Supported (single-trajectory focus)**

- ST GB / PB / GBNSR6 / RISM / decomposition / IE / C2 / alanine scanning (with the same input namelists as
  `gmx_MMPBSA`, where applicable)
- ST GB / GBNSR6 / PB / RISM / NMODE / QM/MMGBSA with `explicit_waters > 0`; selected waters are assigned to the receptor and selected with the same
  `explicit_waters_mask` and `cpptraj closest` workflow used by `gmx_MMPBSA`
- Optional separate receptor/ligand topologies (`-rp` / `-lp`)
- Non-contiguous Amber residue masks
- Radii/`SCREEN` values preserved from the input prmtop (input `PBRadii` does not rebuild normal topologies); an advisory warning is emitted when the topology radius set does not match the conventional choice for the selected `igb`
- Result rewriting and analysis with `gmx_MMPBSA_ana`

**Not supported / different behavior**

- Explicit receptor waters are restricted to ST GB/GBNSR6/PB/RISM/NMODE/QM/MMGBSA and require `SOLVATED_TRAJECTORY=1`
- Ligand multiple-trajectory (`-lt`) — rejected
- Receptor multiple-trajectory (`-rt`) — incomplete; prefer ST
- Solvent or ions are not retained in the final working topologies; the original solvated `-cp` is accepted only for
  explicit-water preprocessing
- QM/MM + explicit waters — one-frame native-AMBER smoke-tested with PM6-DH+; the existing 1–4 EEL consistency warning
  remains and requires scientific review before production interpretation

## Example

A basic single-trajectory command is:

``` bash
amber_MMPBSA -O -i mmpbsa.in \
  -cp ras-raf_complex.prmtop \
  -ct prod_complex.mdcrd \
  -cm ":1-166" ":167-242" \
  -o FINAL_RESULTS_MMPBSA.dat \
  -eo FINAL_RESULTS_MMPBSA.csv
```

No complex structure flag is required. Setup extracts the reference structure from frame 1 of `-ct`; use `-cr` when
explicit chain and residue mapping should follow a separate reference PDB.

AMBER masks can also select non-contiguous residue ranges. For example, the receptor can be residues 1-120 and
181-260 while the ligand is residues 121-180:

``` bash
amber_MMPBSA -O -i mmpbsa.in \
  -cp complex.prmtop \
  -ct prod.mdcrd \
  -cm ":1-120,181-260" ":121-180" \
  -o FINAL_RESULTS_MMPBSA.dat \
  -eo FINAL_RESULTS_MMPBSA.csv
```

See the [AMBER input files example](examples/AMBER/README.md) for a complete runnable example.

An explicit-water calculation uses the solvated topology directly:

``` bash
amber_MMPBSA -O -i mmpbsa_explicit.in \
  -cp solvated.prmtop \
  -ct production.mdcrd \
  -cm ":1-166" ":167-242" \
  -o FINAL_RESULTS_MMPBSA.dat \
  -eo FINAL_RESULTS_MMPBSA.csv
```

The corresponding input contains, for example:

``` text
&general
  explicit_waters=10,
  explicit_waters_mask="within 4",
  solvated_trajectory=1,
/
&gb
  igb=2,
/
```

The interface reference mask is static, while `cpptraj closest` selects the nearest water molecules for each frame.
The selected waters are kept in `COM.prmtop` and `REC.prmtop`; `LIG.prmtop` remains dry. OPC/TIP4P-style extra-point
waters are rejected by default and can only be used with the existing explicit-point stripping approximation.
For AMBER systems with nonstandard water residue names, provide their comma-separated names through
`explicit_waters_group`.

## Progress display

`amber_MMPBSA` and `gmx_MMPBSA` accept
`--progress-style {auto,rich,classic,plain,none}`. The default `auto` mode uses an adaptive Rich progress display in
an interactive terminal and falls back to the classic bar when MPI or output forwarding hides the terminal. Use
`rich` to force the richer renderer through MPI, `plain` for milestone log messages, or `none` to disable progress.
Rich and classic displays also record clean 10% checkpoints in `gmx_MMPBSA.log`, including frame count, processing
rate, elapsed time, ETA, and MPI ranks, so cluster jobs can be followed with `tail -f gmx_MMPBSA.log`.
