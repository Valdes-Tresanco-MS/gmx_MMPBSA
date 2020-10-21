# Protein-ligand binding free energy calculations (Multiple Trajectory method)
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb, *.gro) for each one of the components (_i.e._ complex, receptor and ligand)
* An index file (*.ndx) for each one of the components (_i.e._ complex, receptor and ligand)
* Group numbers in the index files
* Trajectory files (*.xtc, *.pdb, *.gro, *.trr) -- final Gromacs MD trajectory, fitted and with no pbc -- for each one 
of the components (_i.e._ complex, receptor and ligand)
* A *.mol2 file of the unbound ligand used to parametrize ligand for Gromacs. -- Antechamber output *.mol2 is recommended

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -rs rec.pdb -ri rec_index.ndx -rg 1 -rt rec_traj.pdb -lm ligand.mol2 -ls lig.pdb -li lig_index.ndx -lg 2 -lt lig_traj.pdb

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for GB calculation
#This input file is mean to show only that gmx_MMPBSA works. Althought, we tried to used the input files as recommended in the 
#Amber manual, some parameters have been changed to perform more expensive calculations. Feel free to change the parameters 
#according to what is better for your system.
&general
startframe=5, endframe=20, verbose=2, interval=1,
protein_forcefield=3, ligand_forcefield=1
/
&gb
igb=5, saltcon=0.15,
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA#sample-input-files)_

In this case, a multiple trajectory (MT) approximation is followed, which means the receptor and ligand structures and 
trajectories are needed. For the complex, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
For the receptor, a MD Structure+mass(db) file (`rec.pdb`), an index file (`rec_index.ndx`),
a trajectory file (`rec_traj.pdb`), and the receptor group number in the rec_index file (`1`) are needed. For the ligand,
a ligand .mol2 file is needed for generating the ligand topology. Besides, a MD Structure+mass(db) file (`lig.pdb`), an 
index file (`lig_index.ndx`), a trajectory file (`lig_traj.pdb`), and ligand group number in the lig_index file (`2`) are 
needed. The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. In this case,
16 frames `(endframe-startframe)/interval = (21-5)/1 = 16` are going to be used when performing the the MM/PB(GB)SA 
calculation with the igb5 (GB-OBC2) model and a salt concentration = 0.15M.

Once the calculation is done, the GUI app (gmx_MMPBSA_gui) will show up. In this app, you can visualize the results for 
the GB calculation. The results can be saved as *.csv file by clicking "File" in the upper left corner and then 
"Export GB/PB energy (csv)".
