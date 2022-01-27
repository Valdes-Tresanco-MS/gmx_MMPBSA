---
template: main.html
title: gmx_MMPBSA in a nutshell
---

# `gmx_MMPBSA` in a nutshell
`gmx_MMPBSA` provides all the [MMPBSA.py][1] functionalities to GROMACS users. 
In addition, other functionalities have been implemented that ease a number of calculations (_e.g._ MM/PB(GB)SA 
with user-defined internal dielectric constant, interaction entropy and C2 entropy calculations). A GUI application is 
also incorporated that allows for visualizing the results and saving high-quality images.

## Types of calculations you can do
There are many options available in `gmx_MMPBSA`. These are some calculations you can perform with `gmx_MMPBSA`:

* **Normal binding free energies**, with either PB, GB or 3D-RISM solvent models. Each can be done with either
1, 2, or 3 different trajectories. The complex trajectory must always be provided. Whichever trajectories of the 
receptor and/or ligand that are NOT specified will be extracted from that of the complex. This allows a 1-, 
2-, or 3-trajectory analysis. All PB calculations and GB models are performed via the `sander` program. Calculations 
with 3D-RISM solvent model are performed with `rism3d.snglpnt` built with AmberTools.
* **Stability** calculations with any model solvent model (_i.e_ PB, GB or 3D-RISM).
* **Alanine scanning** with either PB or GB implicit solvent models. All trajectories will be mutated to match
the mutated topology files, and whichever calculations that would be carried out for the normal systems are
also carried out for the mutated systems. Note that only 1 mutation is allowed per simulation, and it must
be to an alanine or glycine. If `mutant_only` variable is not set to 1, differences resulting from the mutations are 
calculated.
* **Entropy corrections**. An entropy term can be added to the free energies calculated above using either the
quasi-harmonic, the normal mode, interaction entropy or C2 approximations. Calculations will be performed for the normal 
and mutated systems (alanine scanning) as requested. Normal mode calculations are done with the
`mmpbsa_py_nabnmode` program included with AmberTools.
* **Decomposition schemes**. The energy terms will be decomposed according to the decomposition scheme
outlined in the `idecomp` variable description. This should work with all the above, though entropy terms
cannot be decomposed.
* **QM/MMGBSA**. This is a binding free energy (or stability calculation) using the Generalized Born solvent
model allowing you to treat part of your system with a quantum mechanical Hamiltonian.
* **Support for Membrane Proteins**. Calculate the MMPBSA binding free energy for a ligand bound to a protein
that is embedded into a membrane. In this case, the membrane is implemented as a slab-like region with a uniform or 
heterogeneous dielectric constant depth profile.
  

## `gmx_MMPBSA` a technical view
`gmx_MMPBSA` is a python module that contains 3 applications: 

* [`gmx_MMPBSA`][5] is the fundamental application and carries out all the calculations mentioned above
* [`gmx_MMPBSA_ana`][6] provides an intuitive way to analyze the data from gmx_MMPBSA calculations and save 
  high-quality pictures
* [`gmx_MMPBSA_test`][7] is a tool designed to test if the installation was successful by running one or more available 
  [examples][4] in gmx_MMPBSA.


  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [2]: advanced.md#advanced-options
  [3]: #types-of-calculations-you-can-do
  [4]: examples/README.md
  [5]: gmx_MMPBSA_running.md
  [6]: gmx_MMPBSA_ana_running.md
  [7]: examples/gmx_MMPBSA_test.md#running-gmx_mmpbsa_test