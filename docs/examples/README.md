---
template: main.html
title: Examples
---

These examples cover most of the calculations and analyses available in `gmx_MMPBSA`. Although each example focuses on
a specific case, gmx_MMPBSA can process systems containing several components (_e.g._, metalloprotein-ligand or
protein-DNA-ligand complexes). A single run can also combine multiple calculation types, such as GB with alanine
scanning and per-residue decomposition, or PB with interaction entropy and pairwise decomposition.

## Maintaining example documentation

Example README files under `examples/` are the **canonical** copies used by
`gmx_MMPBSA_test` and GitHub browsing. The MkDocs site reads published copies
under `docs/examples/`.

When you edit an example README:

1. Change the file under `examples/` only.
2. Run from the repository root:

   ```bash
   python scripts/sync_example_docs.py
   ```

3. Commit both the `examples/` change and the synced `docs/examples/` copy.

CI runs `python scripts/sync_example_docs.py --check` and fails if the docs
copies are stale.

The examples available through `gmx_MMPBSA_test` are defined in
`GMXMMPBSA/data/gmx_MMPBSA_test_manifest.json`. When adding, removing, or
renaming a testable example, update the manifest and run:

```bash
python scripts/validate_gmx_MMPBSA_test_docs.py
```

## Jupyter notebooks

Two Jupyter notebooks are available for interactive use:

* [Google Colab notebook](https://colab.research.google.com/github/Valdes-Tresanco-MS/gmx_MMPBSA/blob/colab-notebook/notebooks/gmx_MMPBSA_Colab.ipynb):
  installs a conda-based CPU environment, runs bundled examples, supports uploaded user files, and displays results
  through the Python API.
* [Local notebook](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/blob/colab-notebook/notebooks/gmx_MMPBSA_Local.ipynb):
  runs or loads local results, extracts data with the Python API, and plots energy terms with seaborn.

## Systems

The following examples represent systems that gmx_MMPBSA can process and analyze. The list is not exhaustive;
gmx_MMPBSA can also process other systems with compatible input structures.

* [Protein-protein](Protein_protein/README.md)[^1][^2][^3]
* [Protein-ligand](Protein_ligand/ST/README.md)[^1][^2]
* [Protein-DNA](Protein_DNA/README.md)[^1][^2][^3]
* [Protein-glycan](Protein_glycan/README.md)[^1][^2][^3]
* [MMPBSA with membrane proteins](Protein_membrane/README.md)[^1][^2]
* [Metalloprotein-ligand](Metalloprotein_ligand/README.md)[^1][^2]
* [Multicomponent system (Comp_receptor)](Comp_receptor/README.md)[^1][^2][^3]
<!--
* COVID-19 related proteins
    * [Info](COVID-19_related_proteins/README.md)
    * [Main protease](COVID-19_related_proteins/Main_protease_7l5d/README.md)
    * [Papain-like protease](COVID-19_related_proteins/Papain-like_protease_7koj/README.md)
    * [S1-ACE2 complex](COVID-19_related_proteins/S1-ACE2_complex_7dmu/README.md)
    * [S1 RBD with antibody](COVID-19_related_proteins/S1_RBD_with_antibody_6zlr/README.md)
-->

## Analyses

This section covers the analyses available in gmx_MMPBSA. Although each example focuses on a specific case, one run
can combine several calculation types (_e.g._, GB with alanine scanning and per-residue decomposition, or PB with
interaction entropy and per-residue decomposition).

* [Single Trajectory Protocol](Protein_ligand/ST/README.md)[^1][^2][^3]
* [Multiple Trajectory Protocol](Protein_ligand/MT/README.md)[^1]
* Binding free energy calculations
    * [Binding free energy calculation with GB](Protein_ligand/ST/README.md)
    * [ST MM/PB(GB)SA with explicit receptor waters](Explicit_receptor_waters/README.md)[^1]
    * [Binding free energy calculation with GBNSR6](GBNSR6/README.md)[^1]
    * [Binding free energy calculation with linear PB (LPBE)](Linear_PB_solver/README.md)[^1]
    * [Binding free energy calculation with NonLinear PB (non-LPBE)](NonLinear_PB_solver/README.md)[^1]
    * [Binding free energy calculation with 3D-RISM model](3D-RISM/README.md)[^1]
* [Alanine scanning](Alanine_scanning/README.md)[^1][^2][^3]
* [Decomposition analysis](Decomposition_analysis/README.md)[^1][^2][^3]
* Entropy
    * [Interaction Entropy calculations](Entropy_calculations/Interaction_Entropy/README.md)[^1][^2][^3]
    * [NMODE Entropy calculations](Entropy_calculations/nmode/README.md)[^1]
    * [C2 Entropy calculations](Entropy_calculations/C2_Entropy/README.md)[^1]
* [Stability calculations](Stability/README.md)[^1][^2][^3]
* [QM/MMGBSA calculations](QM_MMGBSA/README.md)[^1]
* [Correlation](Correlation/README.md)
* [Python API extraction](API/README.md)
* [Local API/seaborn notebook](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/blob/colab-notebook/notebooks/gmx_MMPBSA_Local.ipynb)

## Input formats and force-field workflows

The examples below demonstrate preparation workflows for different input formats and force fields. The molecular
systems shown are representative examples, not a list of supported receptor-ligand compositions. For example, a
protein-protein tutorial under native AMBER or PSF/DCD does not mean that the workflow is restricted to
protein-protein complexes. The same preparation principles apply to other systems when the resulting topologies,
structures, trajectories, and molecular selections are compatible with `gmx_MMPBSA` or `amber_MMPBSA`.

Individual calculation models can impose narrower requirements. Review each tutorial's topology-conversion notes
and the restrictions of the selected energy, entropy, decomposition, or membrane method.

### Native AMBER inputs

Use native AMBER topologies and trajectories directly with `amber_MMPBSA`.

* [Representative protein-protein example](AMBER/README.md)[^1]

### GROMACS topologies prepared with CHARMM

These examples use supplied GROMACS topologies containing CHARMM parameters, including specialized membrane and LPH
workflows.

* [Representative protein-ligand example](Protein_ligand_CHARMMff/README.md)[^1][^2]
* [Specialized membrane protein-ligand example (CHARMM-GUI)](Protein_membrane/README.md)[^1]
* [Specialized ligand example with LPH virtual sites](Protein_ligand_LPH_atoms_CHARMMff/README.md)[^1]

### GROMACS topologies prepared with OPLS

* [Representative protein-protein example](OPLS/protein_protein/README.md)

### Converting PSF/DCD simulations

PSF and DCD files are preparation sources, not files read directly by `gmx_MMPBSA`. Convert them into a compatible
topology, structure, trajectory, index, and receptor/ligand selections before analysis.

* [Representative protein-protein conversion example](psf_dcd/protein_protein/README.md)

 [^1]: It is part of the `All` set defined with `-t 0` in `gmx_MMPBSA_test`
 [^2]: It is part of the `Minimal` set defined with `-t 1` in `gmx_MMPBSA_test`
 [^3]: It is part of the `Fast` set defined with `-t 2` in `gmx_MMPBSA_test`
