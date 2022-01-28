---
template: main.html
title: Decomposition
---

# Decomposition analysis

!!! info
    This example can be found in the [docs/examples/Decomposition_analysis][6] directory in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA Github repository


## Requirements

!!! danger
    The ligand mol2 file must be the Antechamber output.

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| Ligand parameters file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `mol2`         | The Antechamber output  `mol2` file of ligand parametrization|
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

!!! tip "Remember"
    When a topology file is defined, the ligand mol2 file is not needed. The ligand mol2 file only required when  
    `gmx_MMPBSA` builds the amber topology from a structure  

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 14

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for decomposition analysis
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Decomposition",
startframe=5,
endframe=14,
/
&gb
igb=5, saltcon=0.150,
/
#make sure to include at least one residue from both the receptor
#and ligand in the print_res mask of the &decomp section.
#http://archive.ambermd.org/201308/0075.html
&decomp
idecomp=2, dec_verbose=3,
print_res="within 4"
#check _GMXMMPBSA_COM_FIXED.pdb file to select which residues are going to be printed in the output file
#print_res="40,41,44,47,78,81,82,85,88,115,118,122,215,218,219,220,232,241"
/
```
_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index 
file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index 
file (`1 13`) are needed. A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` 
input file will contain all the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames 
`(endframe-startframe)/interval = (21-5)/1 = 16` are going to be used when performing the the MM/PB(GB)SA calculation 
with the igb5 (GB-OBC2) model and a salt concentration = 0.15M.

Per-residue `decomp` with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms (`idecomp=2`) is going to be 
performed and residues within 4Å in both receptor and ligand will be printed in the output file. Please see 
`print_res` variable in [`&decomp namelist variables`][4] section 

!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information

### Alternatives for per-residue energy contribution visualization
You can also use [VMD][8] and [Chimera][9] to view the modified pdb file with the per-residue energy contributions

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../input_file.md#decomp-namelist-variables
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Decomposition_analysis
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://www.youtube.com/watch?v=PeboM8KE5SA
  [9]: https://www.youtube.com/watch?v=jKA4fuYuKps