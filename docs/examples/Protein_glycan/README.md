---
template: main.html
title: Protein-glycan
---

# Protein-glycan binding free energy calculations
## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-x-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-x-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`    | Structure file containing the system coordinates |
| An index file                  | :octicons-x-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-x-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-x-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| A topology file (not included) | :octicons-x-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-x-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `top`         | Complex reference structure file with correct assignment of chain ID and residue numbers |
              
:octicons-x-circle-fill-16:{ .req } -> Must be defined -- :octicons-x-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-x-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 26 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for GB calculation
This input file is mean to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual, 
some parameters have been changed to perform more expensive calculations.
Feel free to change the parameters according to what is better for your
system.

&general
startframe=5, endframe=20, verbose=2,
protein_forcefield="oldff/leaprc.ff99SB",
ligand_forcefield="leaprc.GLYCAM_06h-1"
/

&gb
igb=5, saltcon=0.150,
/
```


_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

  
## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and glycan structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 26`) are needed.
The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. **In this case, there
is no need to define a .mol2 for the glycan**. 12 frames are going to be used when performing the the MM/PB(GB)SA 
calculation with the igb5 (GB-OBC2) model and a salt concentration = 0.15M.

Of note, the recommended GLYCAM force fields are: * "leaprc.GLYCAM_06j-1" (Compatible with amber12SB and later), 
"leaprc.GLYCAM_06EPb" (Compatible with amber12SB and later), and "leaprc.GLYCAM_06h-1" (Compatible with amber99SB and 
earlier. It is included in `gmx_MMPBSA` package. If it is selected, it will be copied to $AMBERHOME/dat/x) Check 
[Amber manual](https://ambermd.org/doc12/Amber20.pdf#section.3.3) for more info on GLYCAM force fields.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information
  
  [1]: ../../command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana