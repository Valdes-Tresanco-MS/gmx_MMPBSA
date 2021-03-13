---
template: main.html
title: How gmx_MMPBSA works
---

# How gmx_MMPBSA works

gmx_MMPBSA is a new tool aid to perform end-state free energy calculations based on AMBER's MMPBSA.py with GROMACS 
files.

**But, what does that mean?**

Basically, `gmx_MMPBSA` brings in all the [MMPBSA.py][1] functionalities to GROMACS users. 

[MMPBSA.py][1] is an excellent and well-known tool to perform calculations with the PB and GB models in 
AMBER (February/2021 more than 1200 citations). On the other hand, there are tools like `g_mmpbsa` that are well 
known in the GROMACS community (February/2021 more than 1100 citations). Interestingly, [MMPBSA.py][1] is more robust 
and was published first, however, both have a similar number of citations. This is probably due to the fact that the 
GROMACS (Open source and free) community is large, while AMBER has the restriction of a paid license, or a small 
community with free academic license.

The use of [MMPBSA.py][1] for GROMACS users requires enormous effort to successfully complete the process. In that 
regard, we have decided to make our experience in this process available to the community. We have 
not limited ourselves only to the extent of utility itself, but added additional value. We have designed a tool that 
allows to perform a number of calculations in an effortless way with a graphical tool incorporated, which we believe 
has great potential for the analysis of the results obtained. We also have incorporated some functionalities that aren't
present in [MMPBSA.py][1] and made others more accessible to the beginner user in Amber.


## gmx_MMPBSA general workflow

gmx_MMPBSA can be divided into 3 parts as shown in figure 1. IN the first part, `Preparation`, the topologies 
and trajectories are generated, among other elements depending on the calculations, such as the mutants for the 
alanine or glycine scanning or the interaction residues during decomposition analysis. In the second part, `Calculation`,
the binding free energies and/or entropies are estimated using the selected models. Finally, in the last step `Analysis`,
the results can be analyzed by using the graphical user interface `gmx_MMPBSA_ana`.

![Placeholder](assets/images/workflow.svg)
**Figure 1:** gmx_MMPBSA general workflow

## Required input files

Currently, `gmx_MMPBSA` supports two families of force fields: Amber and CHARMM. Although these force fileds are quite 
similar, the generation of topologies differs. When Amber force filed was used to prepare the system(s), the topologies 
compatible with AMBER package can be generated from either GROMACS topologies or structures. Meanwhile, when CHARMM force
field is used, AMBER's topologies are generated from GROMACS topologies exclusively. (Table 1).

**Table 1:** Required input files for every force field 

| Force field |   Structure   | Index |  Trajectory   | Topology | Reference Structure | Small Molecule Mol2 |
|:-----------:|:-------------:|:-----:|:-------------:|:--------:|:-------------------:|:-------------------:|
|    AMBER    | tpr, gro, pdb |  ndx  | xtc, trr, pdb | Optional |      Optional       |   Only if not top   |
|   CHARMM    | tpr, gro, pdb |  ndx  | xtc, trr, pdb |  Always  |      Optional       |         No          |

## Topology preparation

In this section, we will go in deatil about each file and what they are used for.

**GROMACS files**

`MD Structure+mass(db) (tpr, gro, pdb):`
:   This file is used to generate the structure in pdb format of the complex with `editconf`. We recommend using the 
    *.tpr (production *.tpr) format.

`Index (ndx)` 
:   This the file that contains the index of each element contained in the *.tpr file, organized as groups. This file 
    is required for the definition of the groups corresponding to the receptor and the ligand.

`Trajectory (xtc, trr, pdb)`
:   Trajectory files.

`Topology (top)`
:   This file contains all the parameters corresponding to the force field selected during the system setup. When using 
    a GROMACS topology, `parmed` is used for converting the topologies. This method can be useful when studying a 
    complex system with many elements. However, its main strength (flexibility) is its 
    main weakness. GROMACS can parameterize your system differently than AMBER, so the conversion can result in a 
    topology with very dissimilar parameters that AMBER does not understand. If you use amber, you can consider whether 
    to use this method or not. However, when using CHARMM (any version) force field, the topology is always required. 

`Reference Structure`
:   It corresponds to a file in *.pdb format that must contain a complete structure. That is, the user must make sure 
    that this structure contains the same number of atoms and residues as the complex that he initially defined, the 
    correct residue numbering, as well as its chain ID. This structure is optional, but we recommend using it since it 
    should guarantee a smooth processing of the files in `gmx_MMPBSA`. 
    Essentially, the objective is to be able to correctly assign the mentioned parameters since internally 
    gmx_MMPBSA handles sensitive information, for example: when it extracts the receptor and the ligand from the 
    complex structure, the mutation in the alanine scan, etc.

`Small molecule parameters (mol2)`
:   This file contains a ligand parameterization with `antechamber` that is not found in the selected force field 
    (Amber family). It is only necessary to define it when studying a system that contains this type of ligand, and a
    topology has not been defined. This is used to build amber topology from structure using `tleap`.

**Here are the steps that gmx_MMPBSA follows to generate the topologies:**
{++Specific steps for the topology are highlighted like this++}

- Generates a new index and registers the groups defined by the user
- Generates the pdb of the complex, receptor, and ligand
- {++Topologies are cleaned (remove water and ions)++}
- The structures or {++parameters++} for the receptor and the ligand are generated if it is a ST approximation.
- If Alanine scanning: the mutant of the complex, and the mutant receptor or ligand are generated
- If decomposition: interaction residues are extracted
- The complex is mapped. Registers the continuity of the receptor (example: Metalloprotein structure files are 
  generally structured as follows: Receptor protein + Ligand + Receptor Ions)
- {++The PBRadii is assigned++}
- Topologies are {++converted with `parmed`++} or generated with `tleap`


**The following figure shows the process to generate AMBER topologies depending on the force field.** 

![Placeholder](assets/images/topgeneration.svg)
**Figure 2:** Topology generation workflow for Single Trajectory Approximation

[1]: https://pubs.acs.org/doi/10.1021/ct300418h