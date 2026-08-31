---
template: main.html
title: Complex receptor
---

# Binding free energy calculations in multicomponent systems

!!! info
    This example can be found in the [examples/Comp_receptor][6] directory in the repository folder. If you didn't
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.

## Requirements
!!! danger
    The ligand mol2 file must be the Antechamber output.

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand groups     | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `integers` `strings` | Receptor and ligand single-token group names or zero-based group numbers in the complex index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and free of PBC artifacts. |
| Ligand parameters file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `mol2`         | The Antechamber output  `mol2` file of ligand parametrization|
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`         | GROMACS topology file; any referenced `*.itp` files must be in the same directory |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
Once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 9


where the `mmpbsa.in` input file is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for GB calculation"
Sample input file for GB calculation
This sample input is intended only to demonstrate that gmx_MMPBSA works. Although
it follows the recommendations in the Amber manual, some parameters have been adjusted
to keep computational cost reasonable. Modify them as appropriate for your system.

&general
sys_name="Complex_receptor",
forcefields="oldff/leaprc.ff99SBildn,leaprc.gaff"
PBRadii=4, ions_parameters=1
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
```

!!! info "Keep in mind"
    See all `gmx_MMPBSA` input-file options [here][2] and additional examples [here][3].
    These examples are intended only to demonstrate that gmx_MMPBSA works. Review all variables available for the
    selected calculation, including those not shown here, and confirm their accepted values in the
    [input file section](../../input_file.md). This can help you avoid problems and select suitable approximations.


## Considerations
This example uses the single-trajectory (ST) approximation, so the receptor (protein, DNA, RNA, and ions) and ligand Amber-format topologies and trajectories are generated from the complex. The `forcefields` variable
simplifies topology generation for compound receptors that contain more than one type of molecule.
To do so, an MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and
both the receptor and ligand group numbers in the index file (`33 14`) are needed. A ligand .mol2 file is also needed 
for generating the ligand topology. The `mmpbsa.in` input file contains all the parameters needed for the
MM/PB(GB)SA calculation. In this case, 11 frames are used for MM/PB(GB)SA calculation
with the igb8 (GB-Neck2) model and a salt concentration of 0.15 M. The `mbondi3` radii (`PBRadii=4`) are used, as
recommended for the GB-Neck2 solvation model. A high internal dielectric constant (`intdiel=10`) is used because of the
high number of charged residues at the interface.

In this case, Li/Merz ion parameters (12-6 normal usage set) for Mg ions were used. Check 
[Amber manual](https://ambermd.org/doc12/Amber21.pdf#section.3.6) for more info on ion parameters.

The calculation writes a plain-text statistics file (`FINAL_RESULTS_MMPBSA.dat` by default). When `-eo` is
specified, it also writes all energy terms for every frame to a CSV file. The file name supplied to `-eo` is given
the `.csv` extension (`FINAL_RESULTS_MMPBSA.csv` in this example).

!!! note
    After the calculation, the results can be analyzed with `gmx_MMPBSA_ana` unless `-nogui` was used.
    See the [gmx_MMPBSA_ana][4] section for more information
  
  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Comp_receptor
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
