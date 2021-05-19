---
template: main.html
title: Getting started
---

[![pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![support](https://img.shields.io/badge/support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![python](https://img.shields.io/badge/python-v3.x-blue)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4569307.svg)](http://doi.org/10.5281/zenodo.4569307)

[<img src="../assets/logo.svg" height="120" width="240" align="right"/>]()

# Getting started

gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations with GROMACS 
files. **_It works with all GROMACS versions along with AmberTools20 or 21 and brings improvements in compatibility, 
versatility, analyses, and parallelization compared to existing programs (see [here](versus.md) for a detailed comparison])_**

There are many options for running `gmx_MMPBSA`. Among the types of calculations you can do are:

* Normal binding free energies
* Stability
* Alanine scanning
* Entropy corrections
* Decomposition schemes
* QM/MMGBSA

!!! note
    You can check [`gmx_MMPBSA` in a nutshell page](summary.md) for a more detailed overview of the types of calculations 
    supported in gmx_MMPBSA. Also check our [example page](examples/README.md) to see a detailed list of all the 
    examples available

In the current version, gmx_MMPBSA support a number of different systems including but not limited to:

* Protein-protein
* Protein-ligand
* Protein-DNA
* Metalloprotein-peptide
* Protein-glycan
* Membrane proteins
* Multicomponent systems (_e.g._, Protein-DNA-RNA-Ions-Ligand)

!!! note
    In the current version, gmx_MMPBSA supports Amber and CHARMM force fields. That means any system built with 
    [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) in GROMACS using 
    Amber/CHARMM force field or [CHARMM-GUI](https://www.charmm-gui.org/) is supported in gmx_MMPBSA 😀
  
## `gmx_MMPBSA` a quick overview
`gmx_MMPBSA` is a python module that contains 3 applications: 

* [gmx_MMPBSA](summary.md) is the fundamental application and carries out the calculations mentioned above
* [gmx_MMPBSA_ana](analyzer.md) provides an intuitive way to analyze the data from gmx_MMPBSA calculations and save 
  high-quality pictures
* [gmx_MMPBSA_test](command-line.md#running-gmx_mmpbsa_test) is a tool designed to test if the installation was 
  successful by running one or more available [examples](examples/README.md) in gmx_MMPBSA.
  
!!! note
    gmx_MMPBSA can run in parallel and requires just few things in order to perform any kind of calculation. That is:

    * an input parameters file (`in`, contains all the specifications regarding the type of calculation that is going to be performed)
    * a MD Structure+mass(db) file (`tpr`, `pdb`, `gro`)
    * an index file (`ndx`)
    * receptor and ligand group (group numbers in the index files)
    * a trajectory file (`xtc`, `pdb`, `gro`, `trr`)
    * In certain occasions, defining a topology file (`top`) may be required.

    Once the calculation is done, you can analyze the results in [gmx_MMPBSA_ana](analyzer.md)

    You can check [How gmx_MMPBSA works page](howworks.md) to get more details. Also check our [example page](examples/README.md)
    to see how gmx_MMPBSA works with a real example

The following video shows how a typical binding free energy calculation with GB model and Interaction entropy method 
is done in gmx_MMPBSA

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/k1aLlBhnkxo" frameborder="0" allowfullscreen></iframe>
</div>


!!! Ready
    Ready to use gmx_MMPBSA 😀? Check the [installation page](installation.md) 

!!! Citation
    Please, visit [Cite us page](cite_us.md) to see all the details. Visit [Pypi Stats](https://pypistats.org/packages/gmx-mmpbsa)
    to see how gmx_MMPBSA is doing


[<img src="../assets/images/jetbrains-variant-4.png" height="100" width="178" align="right" />][4]

## Acknowledgments
This project is possible thanks to the Open Source license of the [JetBrains][4] programs. 

  [4]: https://www.jetbrains.com/?from=gmx_MMPBSA
