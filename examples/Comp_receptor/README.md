---
template: main.html
title: Complex receptor
---

# Binding free energy calculations in multicomponent systems

!!! info
    This example can be found in the [examples/Comp_receptor][6] directory in the repository folder. If you didn't
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.

!!! warning "GROMACS topology required"
    As of 1.7.0, `gmx_MMPBSA` **requires** a GROMACS complex topology (`-cp`). The former structure-only path that
    rebuilt Amber topologies with tleap from extracted PDBs (and `-lm` mol2) is removed. Parameters must come from
    the same topology tree used in the MD.

    This example directory does not yet ship a `topol.top` (+ referenced `*.itp` files). Supply the topology from the
    original multicomponent setup before running. Until those files are added, `gmx_MMPBSA_test -t 9` is not included
    in the default test suites.

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand groups     | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `integers` `strings` | Receptor and ligand single-token group names or zero-based group numbers in the complex index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and free of PBC artifacts. |
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`         | GROMACS topology file; any referenced `*.itp` files must be in the same directory. Must include the ligand. |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
Once you are in the folder containing all files (including `topol.top`), the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

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
    [input file section](../../docs/input_file.md). This can help you avoid problems and select suitable approximations.


## Considerations
This example uses the single-trajectory (ST) approximation, so the receptor (protein, DNA, RNA, and ions) and ligand
Amber-format topologies and trajectories are generated from the complex via **GROMACS topology conversion** (`-cp`).
An MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), receptor and
ligand group numbers in the index file (`33 14`), and the matching GROMACS topology (`topol.top`) are needed.

!!! note
    Historical versions of this example used `-lm ligand.mol2` with tleap. That workflow is deprecated: include the
    ligand in the GROMACS topology instead. The bundled `ligand.mol2` is retained only as a historical artifact.

## Final considerations
See [how gmx_MMPBSA works](../../docs/howworks.md) and the [compatibility guide](../../docs/compatibility.md) for topology
conversion limits on multicomponent systems.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Comp_receptor
