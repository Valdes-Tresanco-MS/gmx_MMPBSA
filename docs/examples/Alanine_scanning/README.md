---
template: main.html
title: Alanine scanning
---


# Protein-DNA binding free energy calculations. Alanine scanning

!!! info
    This example can be found in the [docs/examples/Alanine_scanning][6] directory in the repository folder


## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file with correct assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t ala_scan

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for Alanine scanning
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Alanine_Scanning",
startframe=5, endframe=21, verbose=2,
forcefields="oldff/leaprc.ff99SB", PBRadii=4,
/
&gb
igb=8, saltcon=0.150, 
/
&alanine_scanning
mutant='ALA', mutant_res='B:13', cas_intdiel=1
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, 
the ligand is DNA) amber format topologies and trajectories will be obtained from that of the complex. To do so, a 
MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both 
the receptor and ligand group numbers in the index file (`1 12`) are needed. The `mmpbsa.in` input file will contain
all the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames 
`(endframe-startframe)/interval = (21-5)/1 = 16` are going to be used when performing the the MM/PB(GB)SA 
calculation with the igb8 (GB-Neck2) model and a salt concentration = 0.15M. Of note, mbondi3 radii (`PBRadii=4`) 
will be used as recommended for GB-Neck2 solvation model. Also, The dielectric constant (`intdiel`) will be modified 
depending on the nature of the residue to be mutated as `cas_intdiel=1`. In this case, the residue `B:13` is an Arginine
which means `intdiel = 5` will be used.

!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information

### How to define properly which residue is going to be mutated?
The generated PDB files must keep the original numbering, so selection based on residue number is reliable. However, 
the chain id can vary depending on several factors. If you use the reference structure (`-cr` flag), then you don't 
have to worry about any changes. The selection will be based on this structure.

On the other hand, if this reference structure is omitted, then it will depend on:

* The complex structure file format
    
    _The `*.gro` format does not contain information related to chains._

* GROMACS version
    
    _We have noticed that in GROMACS `20xx.x` versions, `trjconv` can omit the chain IDs._

* The option assign_chainID
    
    _This option defines when chain IDs are assigned. Please see this variable in 
    [`&general` namelist variables section][5]_

!!! tip
    In any of these cases, you must verify that the selection is correct. You can see the structure of the fixed 
    Complex structure (`_GMXMMPBSA_FIXED_COM.pdb`), Receptor (`_GMXMMPBSA_REC_Fx.pdb`), and ligand 
    (`_GMXMMPBSA_LIG_Fy.pdb`) respectively. x and y represent the fragment for discontinuous molecules

  [1]: ../../command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool  
  [5]: ../../input_file.md#general-namelist-variables
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Alanine_scanning
  [7]: ../../command-line.md#gmx_mmpbsa_test-command-line