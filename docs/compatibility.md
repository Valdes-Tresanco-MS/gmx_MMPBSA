---
template: main.html
title: Compatibility and upgrades
---

# Compatibility and upgrades

This page preserves the migration guidance for older gmx_MMPBSA releases. For the 1.7.0 release,
follow the migration section below and then consult the [changelog](changelog.md) and current
[installation requirements](installation.md#requirements).

## Upgrade the installed package

In an activated conda environment, upgrade with:

```bash
python -m pip install --upgrade gmx_MMPBSA
```

For an installation tied to a compiled AmberTools environment, use its interpreter instead:

```bash
amber.python -m pip install --upgrade gmx_MMPBSA
```

Confirm which installation is active before migrating calculations or results:

```bash
gmx_MMPBSA -v
python -m pip show gmx_MMPBSA
```

Keep a copy of the original input, result, and `_GMXMMPBSA_info` files. Major releases can change accepted input
variables, generated intermediates, and analyzer compatibility; preserving the original environment is the safest way
to reproduce an older calculation.

## Migrating from 1.6.5 to 1.7.0

!!! note "Release migration"
    This section describes the 1.7.0 release workflow. Keep the 1.6.5 environment available until the new workflow
    has been validated for your system.

The changes below can affect installation, accepted inputs, output files, uncertainty estimates, or the interpretation of
results. Do not replace an existing 1.6.5 working directory in place.

### 1. Preserve the 1.6.5 workflow first

Before testing the 1.7.0 workflow, preserve:

- the complete 1.6.5 input file and command line;
- the original topology, coordinate, index, and trajectory files;
- the text result files, per-frame CSV files, `_GMXMMPBSA_info`, and any fixed/reference PDB files; and
- the 1.6.5 environment used to create the result.

Use the old environment to read or rewrite an archived result when reproducibility matters. A result rewritten or
recalculated with 1.7.0 is not numerically equivalent merely because it uses the same input
files. For a migration check, copy the inputs to a new directory, run a short calculation there, and compare the exit
status, result schema, warnings, and provenance files before starting a production calculation.

### 2. Create a separate environment and verify dependencies

The 1.7.0 Python `>=3.11,<3.13` requirement is part of the supported environment. Release validation used AmberTools `>=24.8,<27` and GROMACS
`>=2022,<2027` as the tested and recommended conda ranges; these ranges are not hard runtime version checks for the
external programs. At runtime, AmberTools executables are resolved from the active environment/PATH, while GROMACS
executables are resolved from `PATH` or the `gmx_path` input variable. Other installed versions can work when their
executables and the selected input/topology workflow are compatible. The Python package also requires NumPy `>=1.26.4,<2`, pandas `>=2.2,<3`, Matplotlib `>=3.8,<4`,
SciPy `>=1.14.1,<2`, Seaborn `>=0.13,<0.14`, `mpi4py >=4.0.1,<5`, ParmEd `>=4.2.2,<5`, and Rich `>=13,<16`.
Use the matching [conda environment file](env.yml) or the [installation instructions](installation.md#requirements),
then run `python -m pip check` before testing a calculation.

Keep MPI and GROMACS from the same environment visible to the selected Python interpreter. In particular, verify
`python -c "import mpi4py; print(mpi4py.__version__)"`, `mpirun --version`, and `gmx --version` from the environment
that will run the calculation.

### 3. Review workflow and input changes

#### Implicit-solvent defaults

The 1.7.0 defaults for omitted implicit-solvent settings are not numerically equivalent to the 1.6.5 defaults:

| Setting | 1.6.5 default | 1.7.0 default | Effect |
| --- | ---: | ---: | --- |
| `&gb igb` | `5` | `8` | Uses the GB-Neck2 model by default |
| `&general PBRadii` | `3` | `4` | Uses `mbondi3`, paired with the default `igb=8` |
| `&pb exdi` | `80.0` | `78.5` | Changes the PB external dielectric constant |

Inputs that explicitly set these variables retain their selected values. To reproduce the implicit-solvent defaults of
an older 1.6.5 calculation, preserve the original input or specify the corresponding values explicitly. Do not compare
1.6.5 and 1.7.0 results until the GB model, topology radii, PB dielectric, topology route, and trajectory frames have
been confirmed to match.

#### Required GROMACS topology (`-cp`)

GROMACS calculations now **require** a complex topology (`-cp`). The legacy path that rebuilt Amber topologies with
tleap from extracted PDBs (optionally with `-lm` mol2) is removed. Always pass the GROMACS `*.top` (and referenced
`*.itp` files) from the MD setup. For unbound MT trajectories, also supply `-rp` / `-lp` with matching tops.
Small-molecule ligands must already be included in that topology tree.

#### Native AMBER workflows

The 1.7.0 release includes the separate `amber_MMPBSA` entry point for native AMBER topology workflows. It is not a
full feature-parity replacement for `gmx_MMPBSA`; use the [native AMBER guide](amber_MMPBSA.md) as the support matrix.
The native interface uses `-cp`, `-ct`, and `-cm`; it does not use the GROMACS-only `-cs` option. Separate receptor
and ligand topology/trajectory options are available where the guide marks the multiple-trajectory path as supported,
but explicit-water, QM/MM, entropy, and ligand-MT restrictions still apply.

Native AMBER and GROMACS-derived calculations also differ in radius handling. Native AMBER working topologies preserve
the input `RADII`, `SCREEN`, and `RADIUS_SET` values. GROMACS conversion applies the selected `PBRadii` through ParmEd.
Do not assume that setting the same `PBRadii` value makes those two workflows equivalent.

#### Composite alanine/glycine mutations

`mutant_res` can now select one or several residues for one composite mutation. All selected residues must belong to the
same component and use the same `mutant` target (`ALA` or `GLY`). The result is one combined mutant effect; it is not an
independent per-residue scan. Composite mutations require `cas_intdiel=0`, with the dielectric chosen explicitly in the
GB/PB input. Use `-cr` when chain IDs, residue numbering, or insertion codes are significant.

#### Quasi-harmonic entropy

`qh_entropy=1` is rejected for new calculations in 1.7.0. Set `qh_entropy=0` and use a supported entropy
method for new work. Historical QH result files remain readable by the analyzer during this compatibility window only;
keep the old environment if an archived result cannot be read reliably. QH support is scheduled for removal after 1.7.0.

#### Interaction Entropy and C2

The Interaction Entropy estimator now makes the full selected ensemble the primary IE estimate. `ie_segment` retains
the final percentage of the cumulative curve as a tail-convergence diagnostic; it no longer replaces the full-ensemble
estimate. The estimator also uses a shared running ensemble mean and a numerically stable log-sum-exp evaluation.

IE and C2 now report deterministic non-overlapping block diagnostics. Main summaries and API tables can contain
`Average`, `SD`, `SEM`, `Block SD`, and `Block SEM`; `Block SEM` is the preferred uncertainty for correlated
trajectories when enough blocks are available. The legacy frame-based SD/SEM values remain for compatibility, but they
should not be relabeled as block estimates. Short trajectories provide weak block evidence, so report the block size
and number of blocks when interpreting the uncertainty.

#### Deprecated `sander_apbs`

New calculations reject `sander_apbs=1`. Use the built-in PBSA solver. The legacy value remains recognized only so
archived results can still be read; do not carry `sander_apbs=1` forward into a 1.7.0 input file.

#### Expected numeric differences vs 1.6.5

Even when the GB model, `PBRadii`, PB dielectric, topology route, and trajectory frames are matched, some quantities
are intentionally not numerically identical to 1.6.5. Treat the following as accepted known differences, not as
regression failures:

- **IE and C2.** The corrected full-ensemble IE estimator, shared running mean, stable log-sum-exp evaluation, and
  deterministic block diagnostics change primary IE/C2 values and uncertainty columns relative to 1.6.5. Do not expect
  bit-for-bit or kcal/mol parity with historical IE or C2 summaries.
- **GBNSR6.** The 1.7.0 parser and frame-ID/term merge correct post-1.6.5 frame and term handling. GBNSR6 totals,
  per-frame series, and related API fields can differ from 1.6.5 even with matched inputs.
- **CHARMM CMAP components.** GROMACS topology conversion omits component-level CHARMM CMAP terms that cancel in the
  single-trajectory binding delta (`C − R − L`). Complex/receptor/ligand component totals and API component fields may
  therefore differ from 1.6.5 while the STP binding Δ remains the compatibility quantity. Multiple-trajectory workflows
  do not get that cancellation, so MT CHARMM CMAP cases can differ in both components and Δ.

Parity checks that pin legacy `igb` / `PBRadii` / `exdi` values still allow these documented exceptions. Compare 1.6.5
and 1.7.0 results only after confirming the model settings above and after accounting for this list.

### 4. Check output paths and result consumers

Per-frame CSV output is generated automatically from the summary output name. For example:

```text
-o results.dat       -> results.csv
-o results.csv       -> results.frames.csv
-do decomposition.dat -> decomposition.csv
```

Explicit `-eo` and `-deo` paths take precedence. Output-path collisions are rejected before files are opened, so update
scripts that assumed the summary and per-frame data could share a filename. Decomposition vector CSV files contain
per-frame values; summary tables additionally expose the block statistics described above.

Every calculation also writes `GMXMMPBSA_radii.json`. Use it to record the requested versus effective radius set,
assignment route, force-field context, model context, and checksums of the final `RADII` and `SCREEN` arrays. For an
additional per-atom audit, set `radii_audit=1`. Preserve this file with the result when comparing 1.6.5 and 1.7.0
calculations.

Failed calculations create a diagnostic zip bundle by default. It can contain logs, input/setup files, generated
intermediates, and up to five trajectory frames. Review it before sharing because it may contain scientific inputs.
Use `--no-error-bundle` when coordinates or other inputs must not be copied into an archive; this changes bundle
creation only, not logging or the calculation exit status.

### 5. Migration checklist

1. Keep the original 1.6.5 environment, inputs, results, and `_GMXMMPBSA_info` files untouched.
2. Build the 1.7.0 environment from the bounded requirements and run `python -m pip check`.
3. Run a short copy of the calculation and confirm the input parser, topology route, MPI/GROMACS discovery, and exit
   status.
4. Inspect the output schema, automatic CSV names, `GMXMMPBSA_radii.json`, and `gmx_MMPBSA.log`.
5. If entropy is enabled, record whether the result is full-ensemble IE, a tail diagnostic, C2, or historical QH, and
   report block size and block count with block uncertainties.
6. Compare the 1.7.0 result with the preserved 1.6.5 result only after confirming that the model, topology route,
   frames, and uncertainty convention are the same, and after accounting for the
   [expected numeric differences](#expected-numeric-differences-vs-165) above.

## Upgrading from 1.4.3 to 1.5.x

Version 1.5.0 substantially changed the calculation, result-processing, and analyzer layers. Treat it as a breaking
migration rather than reusing a 1.4.x working directory in place.

### Input variables

The 1.5 series introduced or exposed additional controls, including:

- `c2_entropy` in [`&general`](input_file.md#general-namelist-variables)
- `extdiel` in [`&gb`](input_file.md#gb-namelist-variables)
- PB controls such as `smoothopt`, `iprob`, `arcres`, `mprob`, `npbopt`, `accept`, `nbuffer`, `fscale`, `npbgrid`,
  `scalec`, `nsnba`, `decompopt`, `use_rmin`, `sprob`, `vprob`, `rhow_effect`, `use_sav`, `maxsph`, and `npbverb`
- RISM controls such as `noasympcorr`, `ljTolerance`, `asympKSpaceTolerance`, the `tree*` settings, `mdiis_del`,
  `mdiis_nvec`, `mdiis_restart`, `maxstep`, `npropagate`, and `entropicDecomp`

The behavior or interpretation of `PBRadii`, `interaction_entropy`, `assign_chainID`, `solvated_trajectory`, `verbose`,
and `temperature` also changed. Review their current definitions before reusing an old input file.

The old topology-preparation settings `protein_forcefield`, `ligand_forcefield`, `forcefields`, and `use_sander` are
not part of the current topology-based workflow. Bonded, nonbonded, charge, ligand, and ion parameters must already
be present in the supplied topology files.

### Calculations and analyzer

The 1.5 series added C2 entropy and nonlinear-PB workflows and exposed more PB and RISM controls. The analyzer was
also reworked and did not accept every result produced by earlier versions. Retain the older installation when an
archived result cannot be read reliably; rerunning or rewriting a result is not scientifically equivalent unless the
same inputs, models, and trajectory frames are preserved.

## Upgrading from 1.3.x to 1.4.x

The 1.4 series renamed several input variables:

| 1.3.x setting | 1.4.x replacement | Introduced |
|---|---|---|
| `entropy = 1` | `qh_entropy = 1` | 1.4.2 |
| `entropy = 2` | `interaction_entropy = 1` | 1.4.2 |
| `entropy_seg` | `ie_segment` | 1.4.2 |
| `protein_forcefield` and `ligand_forcefield` | `forcefields` | 1.4.1 |
| `entropy_temp` | `temperature` | 1.4.1 |

The `sys_name` and `exp_ki` variables were added in 1.4.0, `print_res` changed, and `complex_fixed` became an internal
info-file value. The legacy entropy and force-field variables were deprecated during 1.4.x and removed in 1.5.0.

The analyzer command changed from `-p` to the more flexible `-f` input option. See the
[`gmx_MMPBSA_ana` command-line reference](gmx_MMPBSA_ana_command-line.md).

### Working with archived results

For an existing calculation, first try loading the untouched result with the matching historical gmx_MMPBSA version.
If a copied `_GMXMMPBSA_info` file must be adapted for a 1.4.x analyzer, back it up before adding values such as:

```python
INPUT['temperature'] = 298.15
INPUT['exp_ki'] = 0.0
INPUT['sys_name'] = 'Protein-Ligand'
```

`exp_ki` is needed only for correlation analysis. Do not renumber residues or change topology/trajectory identity merely
to make an old result load. If the archived info file refers to `FILES.complex_fixed`, preserve the corresponding fixed
PDB with its chain IDs and residue numbering.
