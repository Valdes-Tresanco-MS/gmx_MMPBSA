# Protein-DNA binding free energy calculations
In this case, gmx_MMPBSA requires:

* An input parameters file (*.in) -- input file containing all the specifications regarding the type of calculation that
is going to be performed
* The MD Structure+mass(db) file (*.tpr, *.pdb) -- make sure chain labels are included
* An index file (*.ndx) -- *.ndx file containing the receptor and ligand in separated groups
* Receptor and ligand group numbers in the index file
* A trajectory file (*.xtc, *.pdb, *.gro, *.trr) -- final GROMACS MD trajectory, fitted and with no pbc.

_See a detailed list of all the flags in gmx_MMPBSA command line [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#calling-gmx_mmpbsa-from-the-command-line)_

That being said, once you are in the folder containing all files, the command-line will be as follows:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

```
Sample input file for Alanine scanning
#This input file is mean to show only that gmx_MMPBSA works. Althought, we tried to used the input files as recommended in the 
#Amber manual, some parameters have been changed to perform more expensive calculations. Feel free to change the parameters 
#according to what is better for your system.
&general
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield="oldff/leaprc.ff99SB", PBRadii=4
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

_See a detailed list of all the options in gmx_MMPBSA input file [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#the-input-file) 
as well as several [examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files)_

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

Once the calculation is done, the GUI app (gmx_MMPBSA_ana) will show up. In this app, you can visualize the results for 
the GB calculation for both the wild-type and the mutant system. The results can be saved as *.csv file by clicking 
"File" in the upper left corner and then "Export GB/PB energy (csv)".

## How to define properly which residue is going to be mutated?
The generated PDB files must keep the original numbering, so selection based on residue number is reliable. However, the chain id can vary depending on several factors. If you use the reference structure (-cr flag), then you don't have to worry about any changes. The selection will be based on this structure.

On the other hand, if this reference structure is omitted, then it will depend on:
* The complex structure file format
    
    _The * .gro format does not contain information related to chains._

* GROMACS version
    
    _We have seen that the GROMACS 20xx.x versions, trjconv omit the chain IDs._

* Options used to generate the topology in GROMACS (referring to the -merge option)
    
    _If you use the -merge option then GROMACS will merge the chains into one._
* The option assign_chainID
    
    _This option defines when chain IDs are assigned. For the first and second option it must be assign_chainID = 1 or 2. For the 3rd it must be assign_chainID = 2. See the description here._

**Note:** In any of these cases, you must verify that the selection is correct. You can see the structure of the Complex (`_GMXMMPBSA_COM.pdb`), Receptor (`_GMXMMPBSA_REC_FIXED.pdb`), and ligand (`_GMXMMPBSA_LIG_FIXED.pdb`) respectively.
