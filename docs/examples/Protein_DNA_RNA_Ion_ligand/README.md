---
template: main.html
title: Protein-DNA_RNA_ION-Ligand
---

# Protein-DNA_RNA_ION-Ligand binding free energy calculations

!!! info
    This example can be found in the [docs/examples/Protein_DNA_RNA_Ion_ligand][6] directory in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.

!!! danger
    This system was also used to show the usage of `forcefields` variable in [Binding free energy calculations in multicomponent systems][8]
    example. Keep in mind that this example will be removed in version 1.5.0.

## Requirements
!!! danger
    The ligand mol2 file must be the Antechamber output.

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| Ligand parameters file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `mol2`         | The Antechamber output  `mol2` file of ligand parametrization|
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t prot_dna_rna_ions_lig


where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for GB calculation"
Sample input file for GB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-DNA-RNA-ION-Lig",
startframe=1
endframe=10
forcefields="oldff/leaprc.ff99SBildn,leaprc.gaff"
PBRadii=4, ions_parameters=1
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
```

!!! info "Keep in mind"
    See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]. 
    These examples are meant only to show that gmx_MMPBSA works. It is recommended to go over these variables, even 
    the ones that are not included in this input file but are available for the calculation that it's performed and
    see the values they can take (check the [input file section](../../input_file.md)). This will allow you to 
    tackle a number of potential problems or simply use more fancy approximations in your calculations.


## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor (Protein+DNA+RNA+Ions) and 
ligand amber format topologies and trajectories will be obtained from that of the complex. To 
do so, an MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and
both the receptor and ligand group numbers in the index file (`33 14`) are needed. A ligand .mol2 file is also needed 
for generating the ligand topology.The `mmpbsa.in` input file will contain all  the parameters needed for the 
MM/PB(GB)SA calculation. In this case, 11 frames are going to be used when performing the MM/PB(GB)SA calculation 
with the igb8 (GB-Neck2) model and a salt concentration = 0.15M. Of note, mbondi3 radii (`PBRadii=4`) will be used as 
recommended for GB-Neck2 solvation model. Also, a high dielectric constant `intdiel=10` will be used because of the 
high number of charged residues at the interface.

In this case, Li/Merz ion parameters (12-6 normal usage set) for Mg ions were used. Check 
[Amber manual](https://ambermd.org/doc12/Amber20.pdf#section.3.6) for more info on ion parameters.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][4] section for more information
  
  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Protein_DNA_RNA_Ion_ligand
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: ../Comp_receptor/README.md