---
template: main.html
title: Compatibility and upgrades
---

# Compatibility and upgrades

This page preserves the migration guidance for older gmx_MMPBSA releases. For changes after the 1.5 transition,
consult the [changelog](changelog.md) and the current [installation requirements](installation.md#requirements).

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

The old `protein_forcefield`, `ligand_forcefield`, and `use_sander` variables were removed. Use the consolidated
`forcefields` setting where force-field selection is required.

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
