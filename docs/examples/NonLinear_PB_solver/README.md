---
template: main.html
title: Alanine scanning
---


# Binding free energy calculation with linear PB (NLPBE)

!!! info
    **:material-new-box:{:.heart } in gmx_MMPBSA v1.5.0!!!**

    This example can be found in the [docs/examples/NonLinear_PB_solver][6] directory in the repository folder. If you 
    didn't use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder 
    from gmx_MMPBSA GitHub repository.


## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 20 21 -ct com_traj.xtc

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 20 21 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for NLPB calculation"
Sample input file for NLPB
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="NonLinear_PB",
startframe=1,
endframe=10,
forcefields="leaprc.protein.ff14SB",
/
&pb
npbopt=1,
indi=1.0, istrng=0.15,                                                     
radiopt=0,                                           
eneopt=1, cutnb=8.0,
/
# check these threads 
# http://archive.ambermd.org/201203/0191.html
# http://archive.ambermd.org/201610/0114.html
# for more info on NLPB
```

!!! info "Keep in mind"
    See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]. 
    These examples are meant only to show that gmx_MMPBSA works. It is recommended to go over these variables, even 
    the ones that are not included in this input file but are available for the calculation that it's performed and
    see the values they can take (check the [input file section](../../input_file.md)). This will allow you to 
    tackle a number of potential problems or simply use fancier approximations in your calculations.

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand 
amber format topologies and trajectories will be obtained from that of the complex. To do so, an 
MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both 
the receptor and ligand group numbers in the index file (`20 21`) are needed. The `mmpbsa.in` input file will contain
all the parameters needed for the MM/PB(GB)SA calculation. In this case, 6 frames are going to be used when 
performing the MM/PB(GB)SA calculation with the Non-Linear PB solver (`npbopt=1`). The dielectric constant 
(`indi`) is set = 1. 

!!! warning
    When running a NLPB solver, `eneopt` is set = 1. That way, the total electrostatic energy and forces will be 
    computed with the particle-particle particle-mesh (P3M) procedure outlined in [Lu and Luo][8]. In doing so, 
    energy term `EPB` in the output file is set to zero, while `EEL` term includes both the reaction field 
    energy (`EPB`) and the Coulombic energy (`EEL`). The van der Waals energy is computed along with the 
    particle-particle portion of the Coulombic energy. This option requires a nonzero `cutnb` (in this 
    case, `cutnb=8.0`) and `bcopt = 5` (default option).

    It's noteworthy mentioning that `ΔGGAS` and `ΔGSOLV` as reported are no longer properly decomposed. Since 
    `EPB` and `EEL` are combined into the "gas phase" term, the gas and solvation terms can't be separated. 
    Nevertheless, the total ΔTOTAL should be perfectly fine, since everything is sum up together in the end.

!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool  
  [5]: ../../input_file.md#general-namelist-variables
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/NonLinear_PB_solver
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://aip.scitation.org/doi/10.1063/1.1622376