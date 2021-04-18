---
template: main.html
title: Examples
---



Here you can find a representation of almost all the types of calculations that you can perform with `gmx_MMPBSA`. 
Although each example focuses on specific cases, you can use `gmx_MMPBSA` on systems that combine a number of 
characteristics (i.e. metalloprotein-ligand complex, Protein-DNA-ligand, etc.). Several types of calculations (e.g. GB/
Alanine scanning/Per-residue decomposition; PB/Interaction Entropy/Per-wise decomposition) can be also performed in 
the same run for a specific system.


## Analysis

* [MM/3D-RISM](3D-RISM/README.md)[^1]
* [Alanine scanning](Alanine_scanning/README.md)[^1][^2]
* [Decomposition analysis](Decomposition_analysis/README.md)[^1][^2]
* [Multiple Trajectory method (Protein-Ligand)](Protein_ligand/MT/README.md)[^1]
* [Stability calculations](Stability/README.md)[^1][^2]
* Entropy
    * [Interaction Entropy calculations](Entropy_calculations/Interaction_Entropy/README.md)[^1][^2]
    * [nmode Entropy calculations](Entropy_calculations/nmode/README.md)[^1]

## Systems
* [Protein-protein binding free energy calculations](Protein_protein/README.md)[^1][^2]
* [Single Trajectory method (Protein-Ligand)](Protein_ligand/ST/README.md)[^1][^2]
* [Protein-DNA binding free energy calculations](Protein_DNA/README.md)[^1][^2]
* [Metalloprotein-peptide binding free energy calculations](Metalloprotein_peptide/README.md)[^1][^2]
* [Protein-DNA-RNA-Ions-Ligand binding free energy calculations](Protein_DNA_RNA_Ion_ligand/README.md)[^1][^2]
* [Binding free energy calculations with complex receptors](Comp_receptor/README.md)
* [Protein-glycan binding free energy calculations](Protein_glycan/README.md)[^1][^2]
* [MMPBSA with membrane proteins](Protein_membrane/README.md)[^1][^2]

## CHARMMff support
* [Protein-Ligand (ST)](Protein_ligand_CHARMMff/README.md)[^1][^2]
* [Protein-ligand complex embedded in membrane](Protein_membrane_CHARMMff/README.md)[^1]



 [^1]: Can be run individually in `gmx_MMPBSA_test` and it is part of the `all` set defined with `-t` in 
 `gmx_MMPBSA_test`
 [^2]: It is part of the `minimal` (default) set defined with `-t` in `gmx_MMPBSA_test`
