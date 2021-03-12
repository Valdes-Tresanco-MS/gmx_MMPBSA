---
template: main.html
title: gmx_MMPBSA in a nutshell
---

# `gmx_MMPBSA` in a nutshell
`gmx_MMPBSA` brings in all the [MMPBSA.py][1] functionalities to GROMACS users. 
In addition, few other functionalities were implemented that eases a number of calculations (_e.g._ MM/PB(GB)SA 
with different internal dielectric constant, interaction entropy calculation). A GUI application is also incorporated 
that allows for visualizing the results and saving high-quality images.

## Types of calculations you can do
There are many options for running `gmx_MMPBSA`. Among the types of calculations you can do are:

* **Normal binding free energies**, with either PB or GB implicit solvent models. Each can be done with either
1, 2, or 3 different trajectories, but the complex, receptor, and ligand topology files must all be defined. The
complex trajectory must always be provided. Whichever trajectories of the receptor and/or ligand that are NOT
specified will be extracted from the complex trajectory. This allows a 1-, 2-, or 3-trajectory analysis. All PB
calculations and GB models can be performed with just AmberTools via the `mmpbsa_py_energy` program installed with 
`MMPBSA.py`.
* **Stability** calculations with any model.
* **Alanine scanning** with either PB or GB implicit solvent models. All trajectories will be mutated to match
the mutated topology files, and whichever calculations that would be carried out for the normal systems are
also carried out for the mutated systems. Note that only 1 mutation is allowed per simulation, and it must
be to an alanine. If mutant_only is not set to 1, differences resulting from the mutations are calculated. This
option is incompatible with intermediate NetCDF trajectories (see the netcdf = 1 option above). This has the
same program requirements as option 1 above.
* **Entropy corrections**. An entropy term can be added to the free energies calculated above using either the
quasi-harmonic, the normal mode or interaction entropy approximations. Calculations will be performed for the normal 
and mutated systems (alanine scanning) as requested. Normal mode calculations are done with the
mmpbsa_py_nabnmode program included with AmberTools.
* **Decomposition schemes**. The energy terms will be decomposed according to the decomposition scheme
outlined in the `idecomp` variable description. This should work with all of the above, though entropy terms
cannot be decomposed. APBS energies cannot be decomposed, either.
* **QM/MMGBSA**. This is a binding free energy (or stability calculation) using the Generalized Born solvent
model allowing you to treat part of your system with a quantum mechanical Hamiltonian. See [“Advanced
Options”][2] for tips about optimizing this option.
* **MM/3D-RISM**. This is a binding free energy (or stability calculation) using the `3D-RISM` solvation model.
This functionality is performed with `rism3d.snglpnt` built with AmberTools.
* **Membrane Protein MMPBSA**. Calculate the MMPBSA binding free energy for a ligand bound to a protein
that is embedded into a membrane.
  
  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [2]: advanced.md#advanced-options