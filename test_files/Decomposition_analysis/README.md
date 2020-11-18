# Decomposition analysis
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

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file with decomposition analysis
#This input file is mean to show only that gmx_MMPBSA works. Althought, we tried to used the input files as recommended in the 
#Amber manual, some parameters have been changed to perform more expensive calculations. Feel free to change the parameters 
#according to what is better for your system.
&general
startframe=5, endframe=21, interval=1,
/
&gb
igb=5, saltcon=0.150,
/
#make sure to include at least one residue from both the receptor
#and ligand in the print_res mask of the &decomp section.
#http://archive.ambermd.org/201308/0075.html
&decomp
idecomp=2, dec_verbose=3,
print_res="within 4"
#check _GMXMMPBSA_COM_FIXED.pdb file to select which residues are going to be printed in the output file
#print_res="40,41,44,47,78,81,82,85,88,115,118,122,215,218,219,220,232,241"
/
```

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files)_

In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` input file will contain all the 
parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames `(endframe-startframe)/interval = (21-5)/1 = 16`
are going to be used when performing the the MM/PB(GB)SA calculation with the igb5 (GB-OBC2) model and a 
salt concentration = 0.15M.

Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms (`idecomp=2`) is going to be 
performed and residues within 4A of the ligand will be printed in the output file. Of note, you can define the residues 
that are going to be printed in the output file in two different ways:

* print_res="within 4" -- print all the residues within 4A of the ligand
* print_res="40,41,44,47,78,81,82,85,88,115,118,122,215,218,219,220,232,241" -- print those residues in the output file

Once the calculation is done, the GUI app (gmx_MMPBSA_gui) will show up. In this app, you can visualize the results for 
the GB calculation as well as the decomposition analysis. The results can be saved as *.csv file by clicking "File" in 
the upper left corner and then "Export GB/PB energy (csv)". Besides, a *.pdb file with per-residue energies in the 
b-factor column, can be saved by clicking "File" in the upper left corner and then "Energy to bfactor". This way, 
the contribution energy of every residue with its structure can be visualized in a molecular visualization 
program.

Follow these tutorials for coloring residues in a *.pdb file according to their contributions:
* [PyMOL](https://youtu.be/ScoMb6Am6II)
* [UCSF-Chimera](https://youtu.be/jKA4fuYuKps)
* [VMD](https://youtu.be/PeboM8KE5SA)
