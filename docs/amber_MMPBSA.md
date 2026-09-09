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
| Receptor trajectory | `-rt` | no | Unbound receptor trajectory for multiple-trajectory calculations; requires `-rp` and `-rm` |
| Receptor mask | `-rm` | no | Residue mask in the unbound receptor topology when `-rt` is used |
| Ligand topology | `-lp` | no | AMBER topology file for the ligand (otherwise built from the complex) |
| Ligand trajectory | `-lt` | no | Unbound ligand trajectory for multiple-trajectory calculations; requires `-lp` and `-lm` |
| Ligand mask | `-lm` | no | Residue mask in the unbound ligand topology when `-lt` is used |

For explicit-water calculations, `-cp` must be the original solvated AMBER topology and `-ct` must contain the
matching solvated trajectory. The receptor and ligand masks in `-cm` must select solute residues only; solvent and ions
are removed during setup.

### Native AMBER GB radii

Native AMBER topologies already contain the per-atom `RADII` and `SCREEN` values, together with the `RADIUS_SET`
metadata written during `tleap` preparation. `amber_MMPBSA` preserves these values for the generated working
topologies; the input-file `PBRadii` setting is validated but does not replace them in the normal native-AMBER
workflow. Choose the radius set when building the AMBER topology, for example:

```text
set default PBradii mbondi3
saveamberparm complex complex.prmtop
```

For GB calculations, the conventional pairings are `igb=1`/`mbondi`, `igb=2` or `5`/`mbondi2`, `igb=7`/`bondi`, and
`igb=8`/`mbondi3`. If a recognized `RADIUS_SET` differs from that pairing, `amber_MMPBSA` prints a warning but still
uses the radii in the topology. This is intentional: unusual combinations may be valid user choices, so they are
not changed or rejected automatically. If the warning is not intentional, rebuild the topology with the desired
`PBradii` in `tleap`. During alanine-scanning mutant construction, the inherited topology radius name is reused; the
input `PBRadii` value is only a fallback if the source `RADIUS_SET` cannot be identified.

Every calculation writes `GMXMMPBSA_radii.json`. It records requested versus effective radius sets separately for the
complex, receptor, ligand, and any mutant topologies, the `RADIUS_SET` label, ParmEd version, assignment route,
source force-field family, atom representation, active GB/PB/GBNSR6 models, and SHA-256 checksums of the final
`RADII` and `SCREEN` arrays. The record explicitly states that native topology arrays take precedence over the input
`PBRadii`; generated GROMACS topologies record when `PBRadii` was applied through ParmEd. Set `radii_audit=1` to
additionally write one per-atom CSV file per topology.

The provenance warnings are advisory. CHARMM with AMBER `mbondi*` radii for GB is a cross-parameterization protocol;
`charmm_radii` is suggested for CHARMM PB but is never selected automatically. OPLS with AMBER GB radii is labeled
empirically unvalidated, and GROMOS or united-atom inputs receive a strong experimental-support warning. GBNSR6 is
reported independently rather than being judged by the pairwise-GB `igb`/radius compatibility map. No warning silently
changes a force field, radius set, or topology.

## Supported and unsupported features

**Supported**

- ST and MT GB / PB / GBNSR6 / RISM / NMODE / alanine scanning (with the same input namelists as `gmx_MMPBSA`,
  where applicable)
- IE/C2 can be requested with MT, but this is experimental: the estimators use frame-indexed ΔGGAS values from
  independently sampled bound and unbound trajectories and are not validated as independent-ensemble estimators;
  prefer ST or NMODE for production entropy interpretation
- ST GB / GBNSR6 / PB / RISM / NMODE / QM/MMGBSA with `explicit_waters > 0`; selected waters are assigned to the receptor and selected with the same
  `explicit_waters_mask` and `cpptraj closest` workflow used by `gmx_MMPBSA`
- Optional separate receptor/ligand topologies (`-rp` / `-lp`)
- Non-contiguous Amber residue masks
- Radii/`SCREEN` values preserved from the input prmtop (input `PBRadii` does not rebuild normal topologies); an advisory warning is emitted when the topology radius set does not match the conventional choice for the selected `igb`
- Result rewriting and analysis with `gmx_MMPBSA_ana`

**Not supported / different behavior**

- Explicit receptor waters are restricted to ST GB/GBNSR6/PB/RISM/NMODE/QM/MMGBSA and require `SOLVATED_TRAJECTORY=1`
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

For a multiple-trajectory calculation, provide unbound topologies, masks, and trajectories with `-rp/-rm/-rt` and
`-lp/-lm/-lt`:

``` bash
amber_MMPBSA -O -i mmpbsa.in \
  -cp complex.prmtop -ct complex.mdcrd \
  -cm ":1-166" ":167-242" \
  -rp receptor.prmtop -rm ":1-166" -rt receptor.mdcrd \
  -lp ligand.prmtop -lm ":1-76" -lt ligand.mdcrd
```

The three selected trajectories must have the same number of frames after `startframe`, `endframe`, and `interval`
are applied. If more than one file is supplied to any trajectory option, files are concatenated in command-line order
and pooled as one trajectory; they are not treated as independent replicas.

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

## Diagnostic error bundles

Failed calculations create a diagnostic zip bundle by default. Use `--no-error-bundle` to disable that behavior. This
option does not suppress logging or change the exit status. See [Logging and progress](logging.md#warnings-and-errors)
for the bundle contents and data-sharing considerations.

### Automatic CSV filenames

Per-frame CSV output is generated automatically. The default is the summary filename with its suffix
replaced by `.csv`. If that would name the summary itself, `.frames.csv` is used instead: `-o results.csv`
produces the text summary `results.csv` and the energy vectors `results.frames.csv`. The same rule applies
to decomposition output (`-do`/`-deo`). Explicit `-eo` and `-deo` values are preserved. Active output paths
must refer to distinct files; collisions are rejected before opening the output files.
