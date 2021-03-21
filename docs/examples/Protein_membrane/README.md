---
template: main.html
title: Protein-Membrane
---

# MMPBSA with membrane proteins

!!! danger
    The ligand mol2 file must be the Antechamber output.

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| Ligand parameters file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `mol2`         | The Antechamber output  `mol2` file of ligand parametrization|
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `top`         | Complex reference structure file with correct assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

!!! tip "Remember"
    When a topology file is defined, the ligand mol2 file is not needed. The ligand mol2 file only required when  
    `gmx_MMPBSA` build the amber topology from a structure  
_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ci index.ndx -cg 1 13 -ct com_traj.pdb -lm ligand.mol2

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for MMPBSA with membrane proteins
This input file is mean to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual, 
some parameters have been changed to perform more expensive calculations.
Feel free to change the parameters according to what is better for your
system.

&general
startframe=1, endframe=4, interval=1,
debug_printlevel=2, use_sander=1,
/

&pb
radiopt=0, indi=20.0, istrng=0.150, fillratio=1.25, ipb=1, nfocus=1,
bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0, npbverb=1, solvopt=2, inp=2,
memopt=1, emem=7.0, mctrdz=-10.383, mthick=36.086, poretype=1,
maxarcdot=15000
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

  

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.pdb`), an index file (`index.ndx`),
a trajectory file (`com_traj.pdb`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` input file will contain all the
parameters needed for the MM/PB(GB)SA calculation. Of note, special parameters for MMPBSA with membrane proteins have
been included. See more [here](https://ambermd.org/doc12/Amber20.pdf#subsection.6.2.4).
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information
  
  
  [1]: ../../command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana