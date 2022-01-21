---
template: main.html
title: Examples
---

Here you can find a representation of almost all the types of calculations and analyses that you can perform with `gmx_MMPBSA`. 
Although each example focuses on specific cases, you can use `gmx_MMPBSA` on systems that combine a number of 
different components (_i.e._ metalloprotein-ligand complex, Protein-DNA-ligand, etc.). Several types of calculations 
(_e.g._ GB, Alanine scanning and Per-residue decomposition; PB, Interaction Entropy, and Per-wise decomposition) can 
be also performed in the same run for a specific system.

## Systems
* [Protein-protein binding free energy calculations](Protein_protein/README.md)[^1][^2]
* [Single Trajectory method (Protein-Ligand)](Protein_ligand/ST/README.md)[^1][^2]
* [Protein-DNA binding free energy calculations](Protein_DNA/README.md)[^1][^2]
* [Protein-glycan binding free energy calculations](Protein_glycan/README.md)[^1][^2]
* [MMPBSA with membrane proteins](Protein_membrane/README.md)[^1][^2]  
* [Metalloprotein-peptide binding free energy calculations](Metalloprotein_peptide/README.md)[^1][^2]
* [Protein-DNA-RNA-Ions-Ligand binding free energy calculations](Protein_DNA_RNA_Ion_ligand/README.md)[^1][^2]
* [Binding free energy calculations in multicomponent systems](Comp_receptor/README.md)
* COVID-19 related proteins
    * [Main protease](COVID-19_related_proteins/Main_protease_7l5d/README.md)
    * [Papain-like protease](COVID-19_related_proteins/Papain-like_protease_7koj/README.md)
    * [S1-ACE2 complex](COVID-19_related_proteins/S1-ACE2_complex_7dmu/README.md)
    * [S1 RBD with antibody](COVID-19_related_proteins/S1_RBD_with_antibody_6zlr/README.md)

## CHARMMff support
* [Protein-Ligand (ST)](Protein_ligand_CHARMMff/README.md)[^1][^2]
* [Protein-ligand complex embedded in membrane](Protein_membrane_CHARMMff/README.md)[^1]
* [Mycalamide A Bound to the Large Ribosomal Subunit](Ribosomal50S_Mycalamide_A/README.md)
* [Protein-ligand with LPH atoms (ST)](Protein_ligand_LPH_atoms_CHARMMff/README.md)

## Analysis
* [Alanine scanning](Alanine_scanning/README.md)[^1][^2]
* [Decomposition analysis](Decomposition_analysis/README.md)[^1][^2]
* Entropy
    * [Interaction Entropy calculations](Entropy_calculations/Interaction_Entropy/README.md)[^1][^2]
    * [nmode Entropy calculations](Entropy_calculations/nmode/README.md)[^1]
    * [C2 Entropy calculations](Entropy_calculations/C2 Entropy/README.md)
* [Stability calculations](Stability/README.md)[^1][^2]
* [MM/3D-RISM](3D-RISM/README.md)[^1]    
* [QM/MMGBSA calculations](QM_MMGBSA/README.md)
* [Multiple Trajectory method (Protein-Ligand)](Protein_ligand/MT/README.md)[^1]
* [BFE with NonLinear PB solver](NonLinear_PB_solver/README.md)



 [^1]: Can be run individually in `gmx_MMPBSA_test` and it is part of the `all` set defined with `-t` in 
 `gmx_MMPBSA_test`
 [^2]: It is part of the `minimal` (default) set defined with `-t` in `gmx_MMPBSA_test`
