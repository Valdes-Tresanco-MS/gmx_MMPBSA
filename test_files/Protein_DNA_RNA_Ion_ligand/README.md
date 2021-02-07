# Protein-DNA binding free energy calculations
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb, *.gro) 
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final GROMACS MD trajectory, fitted and with no pbc.
* A *.mol2 file of the unbound ligand used to parametrize ligand for GROMACS. -- Antechamber output *.mol2 is recommended

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for GB calculation
#This input file is meant to show only that gmx_MMPBSA works. Althought, we tried to used the input files as recommended in the 
#Amber manual, some parameters have been changed to perform more expensive calculations. Feel free to change the parameters 
#according to what is better for your system.
&general
verbose=2, protein_forcefield="oldff/leaprc.ff99SBildn", ligand_forcefield="leaprc.gaff"
PBRadii=4, ions_parameters=1
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor (Protein+DNA+RNA+Ions) and 
ligand amber format topologies and trajectories will be obtained from that of the complex. To 
do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and
both the receptor and ligand group numbers in the index file (`33 14`) are needed. A ligand .mol2 file is also needed 
for generating the ligand topology.The `mmpbsa.in` input file will contain all  the parameters needed for the 
MM/PB(GB)SA calculation. In this case, 11 frames are going to be used when performing the the MM/PB(GB)SA calculation 
with the igb8 (GB-Neck2) model and a salt concentration = 0.15M. Of note, mbondi3 radii (`PBRadii=4`) will be used as 
recommended for GB-Neck2 solvation model. Also a high dielectric constant `intdiel=10` will be used because of the 
high number of charged residues at the interface.

In this case, Li/Merz ion parameters (12-6 normal usage set) for Mg ions were used. Check 
[Amber manual](https://ambermd.org/doc12/Amber20.pdf#section.3.6) for more info on ion parameters.

Once the calculation is done, the GUI app (gmx_MMPBSA_ana) will show up. In this app, you can visualize the results for 
the GB calculation. The results can be saved as *.csv file by clicking "File" in the upper left corner and then 
"Export GB/PB energy (csv)".
