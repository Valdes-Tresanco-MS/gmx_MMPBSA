---
template: main.html
title: Protein-ligand (Charmm)
---

!!! danger "CHARMM and MM(PB/GB)SA"
    PB model is recommended when working with CHARMMff files. Nevertheless, the combination of PB/GB models with radii 
    optimized for amber atom types (_i.e._ bondi, mbondi, mbondi2, mbondi3) and CHARMM force field hasn't been tested 
    extensively. Please, check this [thread][1] for more information and proceed with caution.

    **:material-new-box:{:.heart } in gmx_MMPBSA v1.5.x series!!!**

    In gmx_MMPBSA v1.5.0 we have added a new PB radii set named _charmm_radii_. **This radii set should be used only 
    with systems prepared with CHARMM force fields**. The atomic radii set for Poisson-Boltzmann calculations has been 
    derived from average solvent electrostatic charge distribution with explicit solvent. The accuracy has been tested 
    with free energy perturbation with explicit solvent [ref.](https://pubs.acs.org/doi/10.1021/jp970736r). Most of 
    the values were taken from a _*radii.str_ file used in PBEQ Solver 
    in [charmm-gui](https://www.charmm-gui.org/?doc=input/pbeqsolver).

    * Radii for protein atoms in 20 standard amino acids from 
    [Nina, Belogv, and Roux](https://pubs.acs.org/doi/10.1021/jp970736r)
    * Radii for nucleic acid atoms (RNA and DNA) from 
    [Banavali and Roux](https://pubs.acs.org/doi/abs/10.1021/jp025852v)
    * Halogens and other atoms from [Fortuna and Costa](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)

# Protein-ligand binding free energy calculations (Single Trajectory method) with CHARMMff files

!!! info
    This example can be found in the [docs/examples/Protein_ligand_CHARMMff][6] directory in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            |                             Required                              |          Type           | Description                                                                                                      |
|:-------------------------------|:-----------------------------------------------------------------:|:-----------------------:|:-----------------------------------------------------------------------------------------------------------------|
| Input parameters file          |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |          `in`           | input file containing all the specifications regarding the type of calculation that is going to be performed     |
| The MD Structure+mass(db) file |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |    `tpr` `pdb`    | Structure file containing the system coordinates                                                                 |
| Receptor and ligand group      |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `integers`        | Receptor and ligand group numbers in the index file                                                              |
| A trajectory file              |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     | `xtc` `pdb` `trr` | final GROMACS MD trajectory, fitted and with no pbc.                                                             |
| A topology file                |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |          `top`          | take into account that *.itp files belonging to the topology file should be also present in the folder           |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |          `pdb`          | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][2]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -cp topol.top

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -cp topol.top

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 10

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for PB calculation"
Sample input file for PB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-Lig-CHARMM",
startframe=1,
endframe=4,
# In gmx_MMPBSA v1.5.0 we have added a new PB radii set named charmm_radii. 
# This radii set should be used only with systems prepared with CHARMM force fields. 
# Uncomment the line below to use charmm_radii set
#PBRadii=7,
/
&pb
# radiopt=0 is recommended which means using radii from the prmtop file for both the PB calculation and for the NP
# calculation
istrng=0.15, fillratio=4.0, radiopt=0
/
```

!!! warning "Remember"
    `radiopt = 0` is recommended which means using radii from the `prmtop` file

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][3] as well as several [examples][4]_

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, an MD Structure+mass(db) file (`com.tpr`), an 
index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the 
index file (`1 13`) are needed. The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA 
calculation. A topology file is also needed (mandatory) in this case to generate the topology files in amber format 
with all the terms for CHARMM force field.
!!! note
    Once the calculation is done, the results can be analyzed in `gmx_MMPBSA_ana` (if `-nogui` flag was not used in the command-line). 
    Please, check the [gmx_MMPBSA_ana][5] section for more information


  [1]: http://archive.ambermd.org/201508/0382.html 
  [2]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [3]: ../../input_file.md#the-input-file
  [4]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Protein_ligand_CHARMMff
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line