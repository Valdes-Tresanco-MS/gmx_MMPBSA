#### Protein-protein binding free energy calculations
In its simplest version, gmx_MMPBSA requires:

* The structure file (*.tpr) -- *.tpr used as input in mdrun program for MD simulation
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.trr) -- final Gromacs MD trajectory, fitted and with no pbc.
* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for GB calculation
&general
startframe=5, endframe=21, verbose=2
/
&gb
igb=2, saltcon=0.150,
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, the 
ligand is also another protein) amber format topologies will be obtained from that of the complex. To do so, a MD *.tpr 
file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand 
group numbers in the index file (`19 20`) are needed. The `mmpbsa.in` input file will contain all the parameters needed
for the MM/PB(GB)SA calculation. In this case, 16 frames `(endframe-startframe)/interval = (21-5)/1 = 16` are going to 
be used when performing the the MM/PB(GB)SA calculation with the igb5 (GB-OBC2) model and a salt concentration = 0.15M.
