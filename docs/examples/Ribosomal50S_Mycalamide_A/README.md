---
template: main.html
title: Mycalamide A Bound to the Large Ribosomal Subunit
---

!!! danger "CHARMM and MM(PB/GB)SA"
    PB model is recommended when working with CHARMMff files. Nevertheless, the combination of PB/GB models and 
    CHARMM force field hasn't been tested extensively. Please, check this [thread][1] for more information and 
    proceed with caution.

# Mycalamide A Bound to the Large Ribosomal Subunit binding free energy calculations (Single Trajectory method) with CHARMMff files

!!! info
    This example can be found in the [docs/examples/Ribosomal50S_Mycalamide_A][6] directory in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA Github repository

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`             | input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`     | Structure file containing the system coordinates|
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`          | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | final GROMACS MD trajectory, fitted and with no pbc.|
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`            | take into account that *.itp files belonging to the topology file should be also present in the folder       |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`            |  Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers       |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][2]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 41 23 -ct com_traj.xtc -cp topol.top

=== "With MPI"

        mpirun -np 6 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 41 23 -ct com_traj.xtc -cp topol.top

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for GB calculation
This input file is meant to show only that gmx_MMPBSA works. 
Although, we tried to use the input files as recommended
in the Amber manual, some parameters have been changed to 
perform more expensive calculations in a reasonable amount 
of time. Feel free to change the parameters according to 
what is better for your system.

&general
verbose=2,
/
&gb
igb=8, intdiel=10, saltcon=0.15,
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][3] as well as several [examples][4]_

## Considerations
This is an extremely complex system (PDB ID: 3I55) that contains several ions, protein chains, ribosomal RNA as 
well as a ligand (Mycalamide A) bound. As you will see, gmx_MMPBSA is able to handle successfully such a complex 
system. Of note, just a relevant part of the entire system has been considered for binding free calculations, 
since the inclusion of the rest will increase the computation time without improving the results. You can check 
the file _GMXMMPBSA_COM_FIXED.pdb during the calculation to see how the complex looks like. In this case, a single 
trajectory (ST) approximation is followed, which means the receptor and ligand structures and trajectories will 
be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), 
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`41 23`) are 
needed. The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. A topology 
file is also needed (mandatory) in this case to generate the topology files in amber format with all the terms for 
CHARMM force field.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information


  [1]: http://archive.ambermd.org/201508/0382.html 
  [2]: ../../command-line.md#gmx_mmpbsa-command-line
  [3]: ../../input_file.md#the-input-file
  [4]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Ribosomal50S_Mycalamide_A
  [7]: ../../command-line.md#gmx_mmpbsa_test-command-line