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

gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations with GROMACS 
files. **_It works with all GROMACS versions along with AmberTools20 or 21 and brings improvements in compatibility, 
versatility, analyses, and parallelization compared to existing programs (see [here](versus.md) for a detailed comparison)_**

[comment]: <> (<div style="text-align:center">)

[comment]: <> (<a href="https://info.flagcounter.com/Zm6l"><img src="https://s11.flagcounter.)

[comment]: <> (com/map/Zm6l/size_s/txt_000000/border_CCCCCC/pageviews_1/viewers_0/flags_0/" alt="Flag Counter" border="0"></a>)

[comment]: <> (</div>)

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

## **Installation**
Ready to use gmx_MMPBSA 😀? Check the [installation page](installation.md)

## **What can be done with gmx_MMPBSA?**
Multiple calculations can be performed with `gmx_MMPBSA` such as:

* Binding free energy calculations with [PB][1], [GB][2] and/or [3D-RISM][3] models
* [Alanine scanning][4] 
* [Binding free energy decomposition][5]
* Entropy corrections ([IE][6], [C2][7], [NMODE][8])
* [Stability][9]
* [QM/MMGBSA][10]

  [1]: examples/Linear_PB_solver/README.md
  [2]: examples/Protein_ligand/ST/README.md
  [3]: examples/3D-RISM/README.md
  [4]: examples/Alanine_scanning/README.md
  [5]: examples/Decomposition_analysis/README.md
  [6]: examples/Entropy_calculations/Interaction_Entropy/README.md
  [7]: examples/Entropy_calculations/C2_Entropy/README.md
  [8]: examples/Entropy_calculations/nmode/README.md
  [9]: examples/Stability/README.md
  [10]: examples/QM_MMGBSA/README.md

!!! note "There is always more..."
    You can check [`gmx_MMPBSA` in a nutshell page](summary.md) for a more detailed overview of the types of calculations 
    supported in gmx_MMPBSA. Also check our [example page](examples/README.md) to see a detailed list of all the 
    examples available

In the current version, gmx_MMPBSA supports a number of different systems including but not limited to:

* [Protein-protein][12]
* [Protein-ligand][13]
* [Protein-DNA][14]
* [Metalloprotein-peptide][15]
* [Protein-glycan][16]
* [Membrane proteins][17]
* Multicomponent systems (_e.g._, [Protein-DNA-RNA-Ions-Ligand][18])

  [12]: examples/Protein_protein/README.md
  [13]: examples/Protein_ligand/ST/README.md
  [14]: examples/Protein_DNA/README.md
  [15]: examples/Metalloprotein_peptide/README.md
  [16]: examples/Protein_glycan/README.md
  [17]: examples/Protein_membrane/README.md
  [18]: examples/Comp_receptor/README.md

!!! note "Support for Amber and CHARMM force fields"
    In the current version, gmx_MMPBSA supports Amber and CHARMM force fields. That means any system built with 
    [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) in GROMACS using 
    Amber/CHARMM force field or [CHARMM-GUI](https://www.charmm-gui.org/) is supported in gmx_MMPBSA 😀

The following video shows how a typical binding free energy calculation with GB model and Interaction entropy method 
is performed in gmx_MMPBSA.

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/k1aLlBhnkxo" frameborder="0" allowfullscreen></iframe>
</div>

## **`gmx_MMPBSA` a quick overview**
`gmx_MMPBSA` is a python module that contains 3 applications: 

* [gmx_MMPBSA](howworks.md) is the fundamental application and carries out the calculations mentioned above
* [gmx_MMPBSA_ana](analyzer.md) provides an intuitive way to analyze the data from gmx_MMPBSA calculations and save 
  high-quality pictures
* [gmx_MMPBSA_test](examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line) is a tool designed to test if the installation 
  was 
  successful by running one or more available [examples](examples/README.md) in gmx_MMPBSA.
  
!!! note "Easy to run"
    gmx_MMPBSA can run in parallel and requires just few things in order to perform any kind of calculation. That is:

    * an input parameters file (`*.in`, contains all the specifications regarding the type of calculation that is going 
    to be performed)
    * a MD Structure+mass(db) file (`*.tpr`, `*.pdb`)
    * an index file (`*.ndx`)
    * receptor and ligand groups (group numbers in the index file)
    * a trajectory file (`*.xtc`, `*.pdb`, `*.trr`)
    * In certain occasions, defining a topology file (`*.top`) may be required.

    Once the calculation is done, you can analyze the results in [gmx_MMPBSA_ana](analyzer.md)

    You can check [How gmx_MMPBSA works page](howworks.md) to get more details. Also check our 
    [example page](examples/README.md) to see how gmx_MMPBSA works with real examples

## **Need help?**
[Help](Q&A/README.md) section contains the most frequently asked questions and errors. Also, look at our 
[Google group](https://groups.google.com/g/gmx_mmpbsa) or the [issues](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues)
section to find out about specific cases and others.

If you still have doubts or cannot solve the problem, please consider opening an 
[issue](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues) or posting in our 
[Google group](https://groups.google.com/g/gmx_mmpbsa)

## **Follow gmx_MMPBSA**
Visit [Pypi Stats](https://pypistats.org/packages/gmx-mmpbsa) or [PePy](https://pepy.tech/project/gmx-mmpbsa)
to see how gmx_MMPBSA is doing.


[<img src="../assets/images/jetbrains-variant-4.png" height="100" width="178" align="right" />][11]

## **Acknowledgments**
This project is possible thanks to the Open Source license of the [JetBrains][11] programs. 

  [11]: https://www.jetbrains.com/?from=gmx_MMPBSA
