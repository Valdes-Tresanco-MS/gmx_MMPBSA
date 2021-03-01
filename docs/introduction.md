---
template: main.html
title: Introduction
---

# Introduction

Molecular Mechanics / Poisson Boltzmann (or Generalized Born) Surface Area (MM/PB(GB)SA) is a post-processing method 
in which representative snapshots from an ensemble of conformations are used to calculate the free energy change 
between two states (typically a bound and free state of a receptor and ligand). Free energy differences are 
calculated by combining the so-called gas phase energy contributions (MM term) that are independent of the chosen 
solvent model as well as solvation free energy components (both polar and non-polar) calculated from an implicit solvent
model combination (PBSA or GBSA) for each species. Entropy contributions to the total free energy may be added as a 
further refinement.

The gas phase free energy contributions are calculated by `sander` or `mmpbsa_py_energy` within the AmberTools package 
according to the force field used in the MD simulation. The solvation free energy contributions may be further 
decomposed into a polar and non-polar contributions. The polar portion is calculated using the Poisson Boltzmann (PB) 
equation, the Generalized Born method, or the Reference Interaction Site Model (RISM). The PB equation is solved 
numerically by either the pbsa program included with AmberTools or by the Adaptive Poisson Boltzmann Solver (APBS) 
program (for more information, see http://www.poissonboltzmann.org). The non-polar contribution is approximated by 
the LCPO method implemented within sander or the molsurf method as implemented in cpptraj. The entropy calculations 
can be done in either a HCT Generalized Born solvation model or in the gas phase using a mmpbsa_py_nabnmode 
program written in the nab programming language, or via the quasi-harmonic approximation in ptraj as in the original 
[MMPBSA.py][1]. In this new module, entropy term (−TΔS) can be also estimated using the so-called [Interaction 
Entropy][2] method, which is theoretically rigorous, computationally efficient, and numerically reliable for 
calculating entropic contribution to free energy in protein–ligand binding and other interaction processes.

Usually, the Single Trajectory (ST) approximation is employed when performing MM/PB(GB)SA calculations. This 
approximation assumes that the configurational space explored by the systems are very similar between the bound and 
unbound states, so every snapshot for each species (_i.e._ complex, receptor, and ligand) is extracted from the same 
trajectory file. This approximation improves binding free energies convergence and also reduces the computing time. 
However, it only should be applied when the molecules in the unbound state present a similar behavior to that of the 
bound state. On the other hand, in the so-called Multiple Trajectory (MT) approximation, the snapshots for each one 
of the species (_i.e._ complex, receptor, and ligand) are extracted from their own trajectory file. This 
approximation, theoretically more rigorous though, leads to higher standard deviation of the binding free energies.  


## Literature
Further information can be found in [Amber manual][3]:

* [MMPBSA.py][4]
* [The Generalized Born/Surface Area Model][5]
* [PBSA][6]
* [Reference Interaction Site Model][7]
* [Generalized Born (GB) for QM/MM calculations][8]

and the foundational papers:

* [Srinivasan J. et al., 1998][9] 
* [Kollman P. A. et al., 2000][10] 
* [Gohlke H., Case D. A. 2004][11] 

as well as some reviews:

* [Genheden S., Ryde U. 2015][12] 
* [Wang et. al., 2018][13]  
* [Wang et. al., 2019][14] 

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [2]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682

  [3]: https://ambermd.org/doc12/Amber20.pdf
  [4]: https://ambermd.org/doc12/Amber20.pdf#chapter.34
  [5]: https://ambermd.org/doc12/Amber20.pdf#chapter.4
  [6]: https://ambermd.org/doc12/Amber20.pdf#chapter.6
  [7]: https://ambermd.org/doc12/Amber20.pdf#chapter.7
  [8]: https://ambermd.org/doc12/Amber20.pdf#subsection.10.1.3
  [9]: https://pubs.acs.org/doi/abs/10.1021/ja981844+
  [10]: https://pubs.acs.org/doi/abs/10.1021/ar000033j
  [11]: https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.10379
  [12]: https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1032936
  [13]: https://www.frontiersin.org/articles/10.3389/fmolb.2017.00087/full
  [14]: https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.9b00055