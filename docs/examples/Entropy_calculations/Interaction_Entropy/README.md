---
template: main.html
title: Interaction Entropy
---

# Entropy calculations

!!! info
    This example can be found in the [docs/examples/Entropy_calculations/Interaction_Entropy][6] directory in the 
     repository folder


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
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `top`         | Complex reference structure file with correct assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t ie

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for entropy calculations (IE)
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
#
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield="oldff/leaprc.ff99SB",

#entropy variable control whether to perform a quasi-harmonic entropy (QH)
# approximation or the Interaction Entropy approximation
# (https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) 
entropy=2, entropy_seg=25, temperature=298
/

&gb
igb=2, saltcon=0.150,
/
```

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]_

  

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand (in this case, 
the ligand is also another protein) amber format topologies and trajectories will be obtained from that of the 
complex. To do so, a MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file 
(`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`19 20`) are needed. The `mmpbsa.
in` input file will contain all the parameters needed for the MM/PB(GB)SA calculation. In this case, 16 frames 
`(endframe-startframe)/interval = (21-5)/1 = 16` are going to be used when performing the the MM/PB(GB)SA 
calculation with the igb2 (GB-OBC1) model and a salt concentration = 0.15M.

[Interaction Entropy (IE)][4] will be calculated and the average for the last quartile (`entropy_seg=25`) of the 
total number of frames will be reported. Of note, two other methods (`QH` and `nmode`) can be used for estimating the 
entropic contribution, though they are way more expensive in computation as compared with IE method.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information
  
  [1]: ../../../command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../input_file.md#the-input-file
  [3]: ../../../input_file.md#sample-input-files  
  [4]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682
  [5]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Entropy_calculations/Interaction_Entropy
  [7]: ../../../command-line.md#gmx_mmpbsa_test-command-line
