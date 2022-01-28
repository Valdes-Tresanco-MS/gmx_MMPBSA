---
template: main.html
title: Protein-ligand (MT)
---

# Protein-ligand binding free energy calculations (Multiple Trajectory method)

!!! info
    This example can be found in the [docs/examples/Protein_ligand/MT][6] directory in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA Github repository

## Requirements

!!! danger
    The ligand mol2 file must be the Antechamber output.

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | (**Complex, Receptor and Ligand**) Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | (**Complex, Receptor and Ligand**) file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | (**Complex, Receptor and Ligand**) Group numbers in the index files |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | (**Complex, Receptor and Ligand**) Final GROMACS MD trajectory, fitted and with no pbc. |
| Ligand parameters file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `mol2`         | The Antechamber output  `mol2` file of ligand parametrization|
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | (**Complex, Receptor and Ligand**) GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

!!! tip "Remember"
    When a topology file is defined, the ligand mol2 file is not needed. The ligand mol2 file only required when  
    `gmx_MMPBSA` build the amber topology from a structure  

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:
=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 16

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -rs rec.pdb -ri rec_index.ndx -rg 1 -rt rec_traj.pdb -lm ligand.mol2 -ls lig.pdb -li lig_index.ndx -lg 2 -lt lig_traj.pdb

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13\
               -ct com_traj.xtc -rs rec.pdb -ri rec_index.ndx -rg 1 -rt rec_traj.pdb \
               -lm ligand.mol2 -ls lig.pdb -li lig_index.ndx -lg 2 -lt lig_traj.pdb

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for GB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-Lig-MT",
startframe=5,
endframe=14,
verbose=2,
forcefields="oldff/leaprc.ff99SB",leaprc.gaff"
/
&gb
igb=5, saltcon=0.150,
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

## Considerations
In this case, a multiple trajectory (MT) approximation is followed, which means the receptor and ligand structures and 
trajectories are needed. For the complex, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
For the receptor, a MD Structure+mass(db) file (`rec.pdb`), an index file (`rec_index.ndx`),
a trajectory file (`rec_traj.pdb`), and the receptor group number in the rec_index file (`1`) are needed. For the ligand,
a ligand .mol2 file is needed for generating the ligand topology. Besides, a MD Structure+mass(db) file (`lig.pdb`), an 
index file (`lig_index.ndx`), a trajectory file (`lig_traj.pdb`), and ligand group number in the lig_index file (`2`) are 
needed. The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. In this case,
16 frames `(endframe-startframe)/interval = (21-5)/1 = 16` are going to be used when performing the the MM/PB(GB)SA 
calculation with the igb5 (GB-OBC2) model and a salt concentration = 0.15M.

!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information

  [1]: ../../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../input_file.md#the-input-file
  [3]: ../../../input_file.md#sample-input-files
  [4]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Protein_ligand/MT
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line