# Protein-DNA binding free energy calculations
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb) -- make sure chain labels are included
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final Gromacs MD trajectory, fitted and with no pbc.

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for Alanine scanning
&general
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield=3, PBRadii=4
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
&alanine_scanning
#make sure to change this parameter to 'ligand' is the mutation is going to be performed in the ligand
mutant='receptor'
mutant_res='B:65'
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, the 
ligand is DNA) amber format topologies and trajectories will be obtained from that of the complex. To 
do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and
both the receptor and ligand group numbers in the index file (`1 12`) are needed. The `mmpbsa.in` input file will contain
all the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames `(endframe-startframe)/interval = (21-5)/1 = 16`
are going to be used when performing the the MM/PB(GB)SA calculation with the igb8 (GB-Neck2) model and a salt 
concentration = 0.15M. Of note, mbondi3 radii (`PBRadii=4`) will be used as recommended for GB-Neck2 solvation model. 
Also a high dielectric constant `intdiel=10` will be used because of the high number of charged residues at the interface.

Alanine scanning will be performed and residue 65 (according to FIXED receptor pdb file, see note below) located in 
chain B will be mutated (of note, this is the residue number 12 in chain B in the original pdb file).

Once the calculation is done, the GUI app (gmx_MMPBSA_gui) will show up. In this app, you can visualize the results for 
the GB calculation for both the wild-type and the mutant system. The results can be saved as *.csv file by clicking 
"File" in the upper left corner and then "Export GB/PB energy (csv)".

## How to define properly which residue is going to be mutated?

* Receptor FIXED pdb files will be always renumbered starting from 1. Chain labels will be kept as they appear in the 
original pdb file. In case of error, please check `_GMXMMPBSA_REC_FIXED.pdb` file and find there the information
(_i.e._ chain and number) for the residue are you interested in.
* Ligand FIXED pdb files will be always renumbered starting from 1. Chain labels will be kept as they appear in the 
original pdb file. In case of error, please check `_GMXMMPBSA_LIG_FIXED.pdb` file and find there the information
(_i.e._ chain and number) for the residue are you interested in.
* Only one mutation is allowed on a residue different than Pro, Gly and Ala.
