# Protein-protein binding free energy calculations
In its simplest version, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb, *.gro) 
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final GROMACS MD trajectory, fitted and with no pbc.

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for GB calculation
#This input file is mean to show only that gmx_MMPBSA works. Althought, we tried to used the input files as recommended in the 
#Amber manual, some parameters have been changed to perform more expensive calculations. Feel free to change the parameters 
#according to what is better for your system.
&general
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield="oldff/leaprc.ff99SB",
/
&gb
igb=2, saltcon=0.150,
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, the 
ligand is also another protein) amber format topologies and trajectories will be obtained from that of the complex. To 
do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and 
both the receptor and ligand group numbers in the index file (`19 20`) are needed. The `mmpbsa.in` input file will 
contain all the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames `(endframe-startframe)/interval = (21-5)/1 = 16`
are going to be used when performing the the MM/PB(GB)SA calculation with the igb2 (GB-OBC1) model and a salt 
concentration = 0.15M.

Once the calculation is done, the GUI app (gmx_MMPBSA_ana) will show up. In this app, you can visualize the results for 
the GB calculation. The results can be saved as *.csv file by clicking "File" in the upper left corner and then 
"Export GB/PB energy (csv)".
