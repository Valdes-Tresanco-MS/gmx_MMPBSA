---
template: main.html
title: Correlation
---


# Correlation

!!! info
    This example can be found in the [docs/examples/Correlation][6] directory in the repository folder. If you 
    didn't 
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.


## Requirements

In this case, `gmx_MMPBSA` will be run for each system and requires for each of them:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | file containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc. |
| A topology file (not included) | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }    |           `top`         | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`         | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, use `bash` to loop over the folder and run the 
calculation in each one of them:

``` bash
for i in */
>do 
>echo $i
>cd $i
>gmx_MMPBSA -O -i mmpbsa.in -cs ../com.tpr -ci ../index.ndx -cg 20 21 -ct ../com_traj.xtc -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui
>cd ..
>done
```

!!! info "Keep in mind"
    See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]. 
    These examples are meant only to show that gmx_MMPBSA works. It is recommended to go over these variables, even 
    the ones that are not included in this input file but are available for the calculation that it's performed and
    see the values they can take (check the [input file section](../../input_file.md)). This will allow you to 
    tackle a number of potential problems or simply use fancier approximations in your calculations.

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand 
amber format topologies and trajectories will be obtained from that of the complex. To do so, an 
MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`), a trajectory file (`com_traj.xtc`), and both 
the receptor and ligand group numbers in the index file (`20 21`) are needed. . As we are running a GB calculation for the wild-type and 
the mutants, we can use the same MD Structure+mass(db) (`../com.tpr`), index (`../index.ndx`), and the
trajectory file (`../com_traj.xtc`) for all systems. The `mmpbsa.in` input file will contain
all the parameters needed for the MM/PB(GB)SA calculation and will be specific for each mutant, indicating the 
mutation and the experimental Ki. In this case, 10 frames are going to be used when performing 
the MM/PB(GB)SA calculation with the igb8 (GB-Neck2) model and a salt concentration = 0.15M. Of note, mbondi3 
radii (`PBRadii=4`) will be used as recommended for GB-Neck2 solvation model.

A plain text output file with all the statistics (default: `FINAL_RESULTS_MMPBSA.dat`) and a CSV-format 
output file containing all energy terms for every frame in every calculation will be saved for each system.

## Correlation analysis with `gmx_MMPBSA_ana`

Once the calculation is done, the results can be analyzed in `gmx_MMPBSA_ana`. For correlation analysis, `gmx_MMPBSA_ana`
can load the system recursively from the folders with the following command-line:

``` bash
gmx_MMPBSA_ana -r
```

The following video shows how to perform correlation analysis in gmx_MMPBSA_ana.

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/QyaUTjmfYvc" frameborder="0" allowfullscreen></iframe>
</div>
    


  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [4]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool  
  [5]: ../../input_file.md#general-namelist-variables
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Correlation
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line