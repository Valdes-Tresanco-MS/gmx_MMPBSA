---
template: main.html
title: Protein-ligand (Charmm)
---

# Protein-ligand binding free energy calculations (Single Trajectory method) with CHARMMff files
## Requirements

In this case, `gmx_MMPBSA` requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb, *.gro) 
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final GROMACS MD trajectory, fitted and with no pbc.
* A topology file (*.top) -- take into account that *.itp files belonging to the topology file should be also present in the folder

!!! warning
    The topology file is always required when working with CHARMM force field

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

  [1]: ../../command-line.md#calling-gmx_mmpbsa-from-the-command-line

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -cp topol.top

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for PB calculation
This input file is mean to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations.
Feel free to change the parameters according to what is better for your
system.

&general
startframe=1, endframe=11, verbose=2,
/
&pb
# radiopt=0 is recommended which means using radii from the prmtop file
# for both the PB calculation and for the NP calculation

istrng=0.15, fillratio=4.0, radiopt=0
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][1] as well as several [examples][2]_

  [1]: ../../input_file.md#the-input-file
  [2]: ../../input_file.md#sample-input-files

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation.

PB model is recommended when working with CHARMMff files. Please, check this [thread](http://archive.ambermd.org/201508/0382.html) 
and proceed cautiously.

Once the calculation is done, the GUI app (`gmx_MMPBSA_ana`) will show up. In this app, you can visualize the results for 
the GB calculation. The results can be saved as *.csv file by clicking "File" in the upper left corner and then 
"Export GB/PB energy (csv)".
