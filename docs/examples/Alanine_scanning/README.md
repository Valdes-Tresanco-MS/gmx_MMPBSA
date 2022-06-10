---
template: main.html
title: Alanine scanning
---


# Alanine scanning

!!! info
    This example can be found in the [docs/examples/Alanine_scanning][6] directory in the repository folder. If you 
    didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.


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

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 20 21 -ct com_traj.xtc -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 20 21 -ct com_traj.xtc -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 12

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for Alanine scanning"
Sample input file for Alanine scanning
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Alanine_Scanning",
startframe=1,
endframe=10,
forcefields="leaprc.protein.ff14SB", 
PBRadii=4,
/
&gb
igb=8, saltcon=0.150,
/
&alanine_scanning
mutant='ALA', mutant_res='A:23', cas_intdiel=1
/
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
all the parameters needed for the MM/PB(GB)SA calculation. In this case, 10 frames are going to be used when performing 
the MM/PB(GB)SA calculation with the igb8 (GB-Neck2) model and a salt concentration = 0.15M. Of note, mbondi3 
radii (`PBRadii=4`) will be used as recommended for GB-Neck2 solvation model. Also, The dielectric constant 
(`intdiel`) will be modified depending on the nature of the residue to be mutated as `cas_intdiel=1`. In this case, 
the residue `A/23` is a Tyrosine which means `intdiel = 3` will be used.

A plain text output file with all the statistics (default: `FINAL_RESULTS_MMPBSA.dat`) and a CSV-format 
output file containing all energy terms for every frame in every calculation will be saved. The file name in 
'-eo' flag will be forced to end in [.csv] (`FINAL_RESULTS_MMPBSA.csv` in this case). This file is only written when 
specified on the command-line.

!!! note
    Once the calculation is done, the results can be analyzed in `gmx_MMPBSA_ana` (if `-nogui` flag was not used in the command-line). 
    Please, check the [gmx_MMPBSA_ana][4] section for more information

### How to define properly which residue is going to be mutated?
The generated PDB files must keep the original numbering, so selection based on residue number is reliable. However, 
the chain id can vary depending on several factors. If you use the reference structure (`-cr` flag), then you don't 
have to worry about any changes. The selection will be based on this structure.

On the other hand, if this reference structure is omitted, then it will depend on:

* GROMACS version
    
    _We have noticed that in GROMACS `20xx.x` versions, `trjconv` can omit the chain IDs._

* The option assign_chainID
    
    _This option defines when chain IDs are assigned. Please see this variable in 
    [`&general` namelist variables section][5]_

!!! tip
    In any of these cases, you must verify that the selection is correct. You can see the structure of the fixed 
    Complex structure (`_GMXMMPBSA_FIXED_COM.pdb`), Receptor (`_GMXMMPBSA_REC_Fx.pdb`), and ligand 
    (`_GMXMMPBSA_LIG_Fy.pdb`) respectively. x and y represent the fragment for discontinuous molecules

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool  
  [5]: ../../input_file.md#general-namelist-variables
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Alanine_scanning
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line