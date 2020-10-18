# MMPBSA with membrane proteins
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The structure file (*.tpr, *.pdb, *.gro)
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final Gromacs MD trajectory, fitted and with no pbc.
* A *.mol2 file of the unbound ligand used to parametrize ligand for Gromacs.


_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ci index.ndx -cg 1 13 -ct com_traj.pdb -lm ligand.mol2

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for MMPBSA with membrane proteins
&general
startframe=1, endframe=4, interval=1,
debug_printlevel=2, use_sander=1,
/
&pb
radiopt=0, indi=20.0, istrng=0.150,
fillratio=1.25, ipb=1, nfocus=1,
bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0,
npbverb=1, solvopt=2, inp=2,
memopt=1, emem=7.0, mctrdz=-10.383, mthick=36.086, poretype=1,
maxarcdot=15000
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and
trajectories will be obtained from that of the complex. To do so, a MD *.tpr file (`com.pdb`), an index file (`index.ndx`),
a trajectory file (`com_traj.pdb`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` input file will contain all the
parameters needed for the MM/PB(GB)SA calculation. Of note, special parameters for MMPBSA with membrane proteins have
been included. See more [here](https://ambermd.org/doc12/Amber20.pdf#subsection.6.2.4)
