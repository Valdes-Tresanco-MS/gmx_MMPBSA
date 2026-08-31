---
template: main.html
title: SARS-CoV-2 Receptor Binding Domain (RBD):CR3022 antibody complex
---

!!! danger "CHARMM and MM(PB/GB)SA"
    PB model is recommended when working with CHARMMff files. Nevertheless, the combination of PB/GB models and 
    CHARMM force field hasn't been tested extensively. See this [thread][1] for more information and
    proceed with caution.

# SARS-CoV-2 Receptor Binding Domain (RBD):CR3022 antibody complex binding free energy calculations (Single Trajectory method) with CHARMMff files

!!! info
    This example can be found in the [examples/COVID-19_related_proteins/S1_RBD_with_antibody_6zlr][6] directory
    in the repository folder. If you didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA Github repository.

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`             | input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`     | Structure file containing the system coordinates|
| Receptor and ligand groups     | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `integers` `strings` | Receptor and ligand single-token group names or zero-based group numbers in the complex index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | final GROMACS MD trajectory, fitted and free of PBC artifacts.|
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`            | GROMACS topology file; any referenced `*.itp` files must be in the same directory                           |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`            |  Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers       |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][2]_

## Command-line
Once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 21 20 -ct com_traj.xtc -cp topol.top

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 21 20 -ct com_traj.xtc -cp topol.top

where the `mmpbsa.in` input file is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for PB calculation"
Sample input file for PB calculation
This sample input is intended only to demonstrate that gmx_MMPBSA works. Although
it follows the recommendations in the Amber manual, some parameters have been adjusted
to keep the computational cost reasonable. Modify them as appropriate for your system.

&general
sys_name="SARS_CoV2_S1_Ab",
/
&pb
#radiopt=0 is recommended which means using radii from the prmtop file for both the PB calculation and for the NP calculation
istrng=0.15, fillratio=4.0, radiopt=0, inp=1,
/
```

!!! warning "Remember"
    `radiopt = 0` is recommended which means using radii from the `prmtop` file

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][3] as well as several [examples][4]_

## Considerations
This example uses the single-trajectory (ST) approximation, so the receptor and ligand structures and
trajectories are generated from the complex. To do so, an MD Structure+mass(db) file (`com.tpr`), an
index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the 
index file (`21 20`) are needed. The `mmpbsa.in` input file contains all the parameters needed for the MM/PB(GB)SA
calculation. A topology file is also needed (mandatory) in this case to generate the topology files in Amber-format
with all the terms for CHARMM force field.
!!! note
    After the calculation, the results can be analyzed with `gmx_MMPBSA_ana` unless `-nogui` was used.
    See the [gmx_MMPBSA_ana][5] section for more information


  [1]: http://archive.ambermd.org/201508/0382.html 
  [2]: ../../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [3]: ../../../input_file.md#the-input-file
  [4]: ../../../input_file.md#sample-input-files
  [5]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/COVID-19_related_proteins/S1_RBD_with_antibody_6zlr
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
