---
template: main.html
title: gmx_MMPBSA in a nutshell
---

# `gmx_MMPBSA` in a nutshell
`gmx_MMPBSA` makes the capabilities of [MMPBSA.py][1] available to GROMACS users. It also supports additional
functionality, including MM/PB(GB)SA calculations with a user-defined internal dielectric constant and interaction
entropy and C2 entropy calculations. The `gmx_MMPBSA_ana` graphical application provides interactive result
visualization and can export high-quality figures.

## Types of calculations you can do
There are many options available in `gmx_MMPBSA`. These are some calculations you can perform with `gmx_MMPBSA`:

* **Standard binding free energies** with PB, GB, or 3D-RISM solvent models. These calculations can use one, two, or
three trajectories. The complex trajectory is always required. If the receptor and/or ligand trajectories are not
specified, they are extracted from the complex trajectory. PB and GB calculations are performed with `sander`, while
3D-RISM calculations are also launched through the AmberTools `sander` backend; `gmx_MMPBSA` distributes frame work
across MPI ranks and collects rank-specific RISM output.
* **Stability** calculations with any solvent model (_i.e_ PB, GB or 3D-RISM).
* **Alanine scanning** with PB or GB implicit-solvent models. The trajectories are mutated to match the mutant
  topologies, and the requested calculations are performed for both the original and mutant systems. A run can mutate
  one residue or several residues together in one composite mutant; it does not independently scan each selected
  residue. The target residues must be mutated to alanine or glycine. Unless `mutant_only` is set to `1`, the output
  also reports the differences caused by the mutation.
* **Entropy corrections**. An entropy term can be added to the free energies calculated above using the normal mode,
interaction entropy or C2 approximations. Quasi-harmonic data from historical result files can still be inspected in
this final compatibility release; all QH support will be removed afterward.
Calculations will be performed for the normal and mutated systems (alanine scanning) as requested. Normal mode calculations are done with the
`mmpbsa_py_nabnmode` program included with AmberTools.
* **Decomposition schemes**. The energy terms will be decomposed according to the decomposition scheme (per-residue or 
per-wise) outlined in the `idecomp` variable description. This should work with all the above, though entropy terms
cannot be decomposed.
* **QM/MMGBSA**. This is a binding free energy (or stability calculation) using the Generalized Born solvent
model allowing you to treat part of your system with a quantum mechanical Hamiltonian.
* **Support for Membrane Proteins**. Calculate the MMPBSA binding free energy for a ligand bound to a protein
  that is embedded into a membrane. In this case, the membrane is implemented as a slab-like region with a uniform or
  heterogeneous dielectric constant depth profile.
* **Native AMBER calculations** through [`amber_MMPBSA`][8], using AMBER topology, coordinate, trajectory, and mask
  files. Supported ST and MT methods, plus native-AMBER restrictions and experimental paths, are documented in the
  dedicated guide.
  

## A technical view of `gmx_MMPBSA`
`gmx_MMPBSA` is a Python package that contains four applications:

* [`gmx_MMPBSA`][5] is the main application and carries out the calculations described above
* [`amber_MMPBSA`][8] runs the supported native-AMBER workflows without requiring GROMACS input files
* [`gmx_MMPBSA_ana`][6] provides an intuitive way to analyze the data from gmx_MMPBSA calculations and save 
  high-quality pictures
* [`gmx_MMPBSA_test`][7] tests whether the installation was successful by running one or more available
  [examples][4] in gmx_MMPBSA.


  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [3]: #types-of-calculations-you-can-do
  [4]: examples/README.md
  [5]: howworks.md
  [6]: analyzer.md
  [7]: examples/gmx_MMPBSA_test.md#running-gmx_mmpbsa_test
  [8]: amber_MMPBSA.md
