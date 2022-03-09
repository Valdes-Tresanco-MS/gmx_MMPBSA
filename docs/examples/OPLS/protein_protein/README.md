---
template: main.html
title: Protein-protein (OPLS)
---

!!! danger "OPLS and MM(PB/GB)SA"
    PB model is recommended when working with OPLSff files. Nevertheless, the combination of PB/GB models with radii 
    optimized for amber atom types (_i.e._ bondi, mbondi, mbondi2, mbondi3) and OPLS force field hasn't been tested 
    extensively. Please, proceed with caution.

# Protein-protein binding free energy calculations (Single Trajectory method) with OPLSff files

!!! info
    This example can be found in the [docs/examples/OPLS/protein_protein][6] directory in the repository folder. If you 
    didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.


## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | File containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Group numbers in the index files |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| A topology file                 | :octicons-check-circle-fill-16:{ .req .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional
 
_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ct com_traj.xtc -ci index.ndx -cg 10 11 -cp topol.top

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.pdb -ct com_traj.xtc -ci index.ndx -cg 10 11 -cp topol.top

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for PB calculation"
Sample input file for PB calculation
#This input file is meant to show only that gmx_MMPBSA works. 
Although, we tried to use the input files as recommended in the
#Amber manual, some parameters have been changed to perform 
more expensive calculations in a reasonable amount of time. 
Feel free to change the parameters #according to what is better 
for your system.

&general
sys_name="OPLS_Support",
startframe=1,
endframe=8,
/
&pb
radiopt=0, istrng=0.150,
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, an MD Structure+mass(db) file (`com.pdb`), an 
index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the 
index file (`10 11`) are needed. The `mmpbsa.in` input file will contain all the 
parameters needed for the MM/PB(GB)SA calculation. In this case, 8 frames are going to be used when performing the 
MM/PB(GB)SA calculation with the PB model (linear PB equation) and a salt concentration = 0.15M.

!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information
  
  [1]: ../../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../input_file.md#the-input-file
  [3]: ../../../input_file.md#sample-input-files
  [5]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/OPLS/protein_protein
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
