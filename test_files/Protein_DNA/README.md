# Protein-DNA binding free energy calculations
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The structure file (*.tpr, *.pdb, *.gro)
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final Gromacs MD trajectory, fitted and with no pbc.

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for GB calculation
&general
startframe=5, endframe=21, verbose=2,
protein_forcefield=3, PBRadii=4
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, the 
ligand is DNA) amber format topologies and trajectories will be obtained from that of the complex. To 
do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and
both the receptor and ligand group numbers in the index file (`1 12`) are needed. The `mmpbsa.in` input file will contain
all  the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames `(endframe-startframe)/interval = (21-5)/1 = 16`
are going to be used when performing the the MM/PB(GB)SA calculation with the igb8 (GB-Neck2) model and a salt 
concentration = 0.15M. Of note, mbondi3 radii (`PBRadii=4`) will be used as recommended for GB-Neck2 solvation model. 
Also a high dielectric constant `intdiel=10` will be used because of the high number of charged residues at the interface.