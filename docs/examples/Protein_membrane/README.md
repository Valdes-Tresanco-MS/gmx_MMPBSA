---
template: main.html
title: Protein-Membrane
---

# MMPBSA with membrane proteins
## Requirements

In this case, `gmx_MMPBSA` requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb, *.gro)
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final GROMACS MD trajectory, fitted and with no pbc.
* A *.mol2 file of the unbound ligand used to parametrize ligand for GROMACS.


_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

  [1]: ../../command-line.md#calling-gmx_mmpbsa-from-the-command-line

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

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][1] as well as several [examples][2]_

  [1]: ../../input_file.md#the-input-file
  [2]: ../../input_file.md#sample-input-files

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.pdb`), an index file (`index.ndx`),
a trajectory file (`com_traj.pdb`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` input file will contain all the
parameters needed for the MM/PB(GB)SA calculation. Of note, special parameters for MMPBSA with membrane proteins have
been included. See more [here](https://ambermd.org/doc12/Amber20.pdf#subsection.6.2.4).

Once the calculation is done, the GUI app (`gmx_MMPBSA_ana`) will show up. In this app, you can visualize the results for 
the PB calculation. The results can be saved as *.csv file by clicking "File" in the upper left corner and then 
"Export GB/PB energy (csv)".
