---
template: main.html
title: Getting started
---

[![pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![support](https://img.shields.io/badge/support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![python](https://img.shields.io/badge/python-v3.x-blue)]()
[![Downloads](https://pepy.tech/badge/gmx-mmpbsa)](https://pepy.tech/project/gmx-mmpbsa)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.1c00645-blue)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645)

[<img src="../assets/TOC.png" height="120%" width="258" align="right"/>]()

# Getting started

gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations with GROMACS 
files. **_It works with all GROMACS versions along with AmberTools20 or 21 and brings improvements in compatibility, 
versatility, analyses, and parallelization compared to existing programs (see [here](versus.md) for a detailed comparison)_**

!!! note "Cite gmx_MMPBSA"
    <a href="https://www.scimagojr.com/journalsearch.php?q=5100155074&amp;tip=sid&amp;exact=no" title="SCImago Journal 
    &amp; Country Rank"><img border="0" align="left" width=160 src="https://www.scimagojr.com/journal_img.php?id=5100155074" 
    alt="SCImago Journal &amp; Country Rank"  /></a>
    
    `gmx_MMPBSA` official paper has been published on _Journal of Chemical Theory and Computation_ and can be accessed 
    [here](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645). If you use `gmx_MMPBSA`, please cite it as follows:
    
    Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. _gmx_MMPBSA: A New Tool to Perform 
    End-State Free Energy Calculations with GROMACS_. Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291. 
    https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645. 
    
    Download | [*.bib](gmx_MMPBSA_citation.bib) | [*.ris](gmx_MMPBSA_citation.ris)

    **Please, also consider citing MMPBSA.py's paper:**

    Bill R. Miller, T. Dwight McGee, Jason M. Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg. MMPBSA.
    py: An Efficient Program for End-State Free Energy Calculations.  Journal of Chemical Theory and Computation, 
    2012 8 (9), 3314-3321. https://pubs.acs.org/doi/10.1021/ct300418h. 

    Download | [*.bib](MMPBSA_py_citation.bib) | [*.ris](MMPBSA_py_citation.ris) | [*.xml](MMPBSA_py_citation.xml)

    Visit [Cite us page](cite_us.md#example) for more information on how to cite `gmx_MMPBSA` and the programs/methods 
    implemented in it.

Multiple calculations can be performed with `gmx_MMPBSA` such as:

* Normal binding free energies
* Alanine scanning  
* Decomposition schemes
* Entropy corrections
* Stability
* QM/MMGBSA

!!! note "There is always more..."
    You can check [`gmx_MMPBSA` in a nutshell page](summary.md) for a more detailed overview of the types of calculations 
    supported in gmx_MMPBSA. Also check our [example page](examples/README.md) to see a detailed list of all the 
    examples available

In the current version, gmx_MMPBSA supports a number of different systems including but not limited to:

* Protein-protein
* Protein-ligand
* Protein-DNA
* Metalloprotein-peptide
* Protein-glycan
* Membrane proteins
* Multicomponent systems (_e.g._, Protein-DNA-RNA-Ions-Ligand)

!!! note "Support for Amber and CHARMM force fields"
    In the current version, gmx_MMPBSA supports Amber and CHARMM force fields. That means any system built with 
    [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) in GROMACS using 
    Amber/CHARMM force field or [CHARMM-GUI](https://www.charmm-gui.org/) is supported in gmx_MMPBSA 😀

The following video shows how a typical binding free energy calculation with GB model and Interaction entropy method 
is performed in gmx_MMPBSA.

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/k1aLlBhnkxo" frameborder="0" allowfullscreen></iframe>
</div>

## `gmx_MMPBSA` a quick overview
`gmx_MMPBSA` is a python module that contains 3 applications: 

* [gmx_MMPBSA](summary.md) is the fundamental application and carries out the calculations mentioned above
* [gmx_MMPBSA_ana](analyzer.md) provides an intuitive way to analyze the data from gmx_MMPBSA calculations and save 
  high-quality pictures
* [gmx_MMPBSA_test](command-line.md#running-gmx_mmpbsa_test) is a tool designed to test if the installation was 
  successful by running one or more available [examples](examples/README.md) in gmx_MMPBSA.
  
!!! note "Easy to run"
    gmx_MMPBSA can run in parallel and requires just few things in order to perform any kind of calculation. That is:

    * an input parameters file (`*.in`, contains all the specifications regarding the type of calculation that is going to be performed)
    * a MD Structure+mass(db) file (`*.tpr`, `*.pdb`, `*.gro`)
    * an index file (`*.ndx`)
    * receptor and ligand groups (group numbers in the index file)
    * a trajectory file (`*.xtc`, `*.pdb`, `*.gro`, `*.trr`)
    * In certain occasions, defining a topology file (`*.top`) may be required.

    Once the calculation is done, you can analyze the results in [gmx_MMPBSA_ana](analyzer.md)

    You can check [How gmx_MMPBSA works page](howworks.md) to get more details. Also check our [example page](examples/README.md)
    to see how gmx_MMPBSA works with a real example

!!! Note "Ready?"
    Ready to use gmx_MMPBSA 😀? Check the [installation page](installation.md) 

!!! Note "Follow gmx_MMPBSA"
    Visit [Pypi Stats](https://pypistats.org/packages/gmx-mmpbsa) or [PePy](https://pepy.tech/project/gmx-mmpbsa)
    to see how gmx_MMPBSA is doing.


[<img src="../assets/images/jetbrains-variant-4.png" height="100" width="178" align="right" />][4]

## Acknowledgments
This project is possible thanks to the Open Source license of the [JetBrains][4] programs. 

  [4]: https://www.jetbrains.com/?from=gmx_MMPBSA
