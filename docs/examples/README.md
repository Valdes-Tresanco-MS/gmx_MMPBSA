---
template: main.html
title: Examples
---

Here you can find a representation of almost all the types of calculations and analyses that you can perform 
with `gmx_MMPBSA`. Although each example focuses on specific cases, you can use `gmx_MMPBSA` on systems that combine 
a number of different components (_i.e._ metalloprotein-ligand complex, Protein-DNA-ligand, etc.). In addition, 
several types of calculations (_e.g._ GB, Alanine scanning and Per-residue decomposition; PB, Interaction Entropy, 
and Per-wise decomposition) can be also performed in the same run for a specific system.

## Systems

This is a representation of the systems that can be processed and analyzed with gmx_MMPBSA. Our program has a robust 
method to process input structures. Even if your system is not represented here, still have a go, you won't 
be disappointed! 😀

* [Protein-protein](Protein_protein/README.md)[^1][^2][^3]
* [Protein-ligand](Protein_ligand/ST/README.md)[^1][^2]  
* [Protein-DNA](Protein_DNA/README.md)[^1][^2][^3]
* [Protein-glycan](Protein_glycan/README.md)[^1][^2][^3]
* [MMPBSA with membrane proteins](Protein_membrane/README.md)[^1][^2]  
* [Metalloprotein-peptide](Metalloprotein_peptide/README.md)[^1][^2]
* [Multicomponent system (Protein-DNA-RNA-Ions-Ligand)](Comp_receptor/README.md)[^1][^2][^3]
* COVID-19 related proteins
    * [Info](COVID-19_related_proteins/README.md)
    * [Main protease](COVID-19_related_proteins/Main_protease_7l5d/README.md)
    * [Papain-like protease](COVID-19_related_proteins/Papain-like_protease_7koj/README.md)
    * [S1-ACE2 complex](COVID-19_related_proteins/S1-ACE2_complex_7dmu/README.md)
    * [S1 RBD with antibody](COVID-19_related_proteins/S1_RBD_with_antibody_6zlr/README.md)

## CHARMMff support

This section focuses more on how to work with systems prepared with CHARMM force fields. We only show few examples 
for better clarity.

* [Protein-Ligand](Protein_ligand_CHARMMff/README.md)[^1][^2]
* [Protein-Ligand complex embedded in membrane](Protein_membrane_CHARMMff/README.md)[^1]
* [Protein-Ligand with LPH atoms](Protein_ligand_LPH_atoms_CHARMMff/README.md)

## Analysis

This section focuses on the analysis that can be performed with gmx_MMPBSA. Although each example focuses on specific 
cases, you can use gmx_MMPBSA to perform several types of calculations (_e.g._ GB, Alanine scanning and Per-residue 
decomposition; PB, Interaction Entropy, and Per-residue decomposition) in the same run for a specific system.

* [Single Trajectory Protocol](Protein_ligand/ST/README.md)[^1][^2][^3]
* [Multiple Trajectory Protocol](Protein_ligand/MT/README.md)[^1]
* Binding free energy calculations
    * [Binding free energy calculation with GB](Protein_ligand/ST/README.md)
    * [Binding free energy calculation with linear PB (LPBE)](Linear_PB_solver/README.md)
    * [Binding free energy calculation with NonLinear PB (non-LPBE)](NonLinear_PB_solver/README.md)  
    * [Binding free energy calculation with 3D-RISM model](3D-RISM/README.md)[^1]
* [Alanine scanning](Alanine_scanning/README.md)[^1][^2][^3]
* [Decomposition analysis](Decomposition_analysis/README.md)[^1][^2][^3]
* Entropy
    * [Interaction Entropy calculations](Entropy_calculations/Interaction_Entropy/README.md)[^1][^2][^3]
    * [NMODE Entropy calculations](Entropy_calculations/nmode/README.md)[^1]
    * [C2 Entropy calculations](Entropy_calculations/C2_Entropy/README.md)
* [Stability calculations](Stability/README.md)[^1][^2][^3]
* [QM/MMGBSA calculations](QM_MMGBSA/README.md)

## Support for psf_dcd files

This section focuses on how to work with psf-dcd files. These files are used for several MD simulation 
programs such as NAMD, OpenMM or GENESIS. We plan to add more examples in the near future.

* [Protein-protein binding free energy calculations](psf_dcd/protein_protein/README.md)
* [Protein-ligand binding free energy calculations](psf_dcd/protein_ligand/README.md)

 [^1]: It is part of the `All` set defined with `-t 0` in `gmx_MMPBSA_test`
 [^2]: It is part of the `Minimal` set defined with `-t 1` in `gmx_MMPBSA_test`
 [^3]: It is part of the `Fast` set defined with `-t 2` in `gmx_MMPBSA_test`
