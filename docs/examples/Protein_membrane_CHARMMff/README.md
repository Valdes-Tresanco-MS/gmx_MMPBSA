---
template: main.html
title: Protein-ligand embedded in membrane (Charmm)
---

!!! danger "CHARMM and MM(PB/GB)SA"
    PB model is recommended when working with CHARMMff files. Nevertheless, the combination of PB/GB models and 
    Charmm force field hasn't been tested extensively. Please, check this [thread][1] for more information and 
    proceed with caution.

# Protein-ligand embedded in membrane binding free energy calculations (Single Trajectory method) with CHARMMff files

!!! info
    This example can be found in the [docs/examples/Protein_membrane_CHARMMff][6] directory in the repository folder

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`             | input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`     | Structure file containing the system coordinates|
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`          | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | final GROMACS MD trajectory, fitted and with no pbc.|
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`            | take into account that *.itp files belonging to the topology file should be also present in the folder       |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `top`            | take into account that *.itp files belonging to the topology file should be also present in the folder       |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][2]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t memb_charmm

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.pdb -ci index.ndx -cg 6 5 -ct md.xtc -cp topol.top

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.pdb -ci index.ndx -cg 6 5 -ct md.xtc -cp topol.top

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for MMPBSA with membrane proteins
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
interval=1,
debug_printlevel=2, use_sander=1,
/

&pb
radiopt=0, indi=20.0, istrng=0.150, fillratio=1.25, ipb=1, nfocus=1,
bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0, npbverb=1, solvopt=2, inp=2,
memopt=1, emem=7.0, mctrdz=31, mthick=31, poretype=1,
maxarcdot=15000
/
```

!!! warning "Remember"
    `radiopt = 0` is recommended which means using radii from the `prmtop` file

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][3] as well as several [examples][4]_


## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. A topology file is
also needed (mandatory) in this case to generate the topology files in amber format with all the terms for CHARMM force 
field.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information


  [1]: http://archive.ambermd.org/201508/0382.html 
  [2]: ../../command-line.md#gmx_mmpbsa-command-line
  [3]: ../../input_file.md#the-input-file
  [4]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Protein_membrane_CHARMMff
  [7]: ../../command-line.md#gmx_mmpbsa_test-command-line