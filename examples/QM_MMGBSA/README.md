---
template: main.html
title: QM/MMGBSA
---

# QM/MMGBSA binding free energy calculations

!!! info
    This example can be found in the [examples/QM_MMGBSA][6] directory in the repository folder. If you didn't
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from 
    gmx_MMPBSA GitHub repository.


## Requirements
!!! danger
    The ligand mol2 file must be the Antechamber output.

In this case, `gmx_MMPBSA` requires:

| Input File required            |                             Required                              |       Type        | Description                                                                                                      |
|:-------------------------------|:-----------------------------------------------------------------:|:-----------------:|:-----------------------------------------------------------------------------------------------------------------|
| Input parameters file          |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `in`        | Input file containing all the specifications regarding the type of calculation that is going to be performed     |
| The MD Structure+mass(db) file |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |    `tpr` `pdb`    | Structure file containing the system coordinates                                                                 |
| An index file                  |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |       `ndx`       | File containing the receptor and ligand in separated groups                                                      |
| Receptor and ligand group      |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |    `integers`     | Group numbers in the index files                                                                                 |
| A trajectory file              |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     | `xtc` `pdb` `trr` | Final GROMACS MD trajectory, fitted and with no pbc.                                                             |
| Ligand parameters file         |    :octicons-check-circle-fill-16:{ .req .scale_icon_medium }     |      `mol2`       | The Antechamber output  `mol2` file of ligand parametrization                                                    |
| A topology file (not included) |  :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium }   |       `top`       | GROMACS topology file (The `* .itp` files defined in the topology must be in the same folder                     |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |       `pdb`       | Complex reference structure file (without hydrogens) with the desired assignment of chain ID and residue numbers |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

!!! tip "Remember"
    When a topology file is defined, the ligand mol2 file is not needed. The ligand mol2 file only required when  
    `gmx_MMPBSA` build the amber topology from a structure  
_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 23

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for QM/MMGBSA calculation"
Sample input file for QM/MMGBSA calculation
This input file is meant to show only that gmx_MMPBSA works. 
Although, we tried to use the input files as recommended in the
Amber manual, some parameters have been changed to perform more 
expensive calculations in a reasonable amount of time. Feel free 
to change the parameters according to what is better for your system.

&general
sys_name="QM/MMGBSA",
startframe=5,
endframe=14,
PBRadii=2,
forcefields="oldff/leaprc.ff99SB,leaprc.gaff"
/
&gb
igb=1, saltcon=0.150,
ifqnt=1, qm_theory=PM6-DH+,

# Residues to be treated with QM can be selected using different approaches. Please, make sure to include at least
# one residue from both the receptor and ligand in the qm_residues mask when using 'ifqnt'. This requirement is
# automatically fulfilled when using the within keyword https://groups.google.com/g/gmx_mmpbsa/c/GNb4q4YGCH8

# Residue selection by distance (recommended)
qm_residues="within 4"

## Explicit residue selection
#qm_residues="A/40-41,44,47,78,81-82,85,88,115,118,122,215,218-220,232 B/241"

# Residue selection with amber masks
#com_qmmask="(:44,47,85,88,218&!@N,H,CA,HA,C,O) | :241"
#rec_qmmask="(:44,47,85,88,218&!@N,H,CA,HA,C,O)"
#lig_qmmask=":1"
/
```

!!! info "Keep in mind"
    See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3]. 
    These examples are meant only to show that gmx_MMPBSA works. It is recommended to go over these variables, even 
    the ones that are not included in this input file but are available for the calculation that it's performed and
    see the values they can take (check the [input file section](../../input_file.md)). This will allow you to 
    tackle a number of potential problems or simply use fancier approximations in your calculations.

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, an MD Structure+mass(db) file (`com.tpr`), an index file (`index.ndx`),
a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the index file (`1 13`) are needed.
A ligand .mol2 file is also needed for generating the ligand topology. The `mmpbsa.in` input file will contain all 
the parameters needed for the QM/MMGBSA calculation. 10 frames are going to be used when performing QM/MMGBSA 
calculation with the igb1 (GB-HCT) model (note that `mbondi` raddi set `PBRadii=2` 
is used), **PM6-DH+** (the default dispersion- and hydrogen-bond-corrected PM6 Hamiltonian) and a salt concentration of 0.15 M.
If `qm_theory` is omitted, the same **PM6-DH+** default is used.

A plain text output file with all the statistics (default: `FINAL_RESULTS_MMPBSA.dat`) and a CSV-format 
output file containing all energy terms for every frame in every calculation will be saved. The file name in 
'-eo' flag will be forced to end in [.csv] (`FINAL_RESULTS_MMPBSA.csv` in this case). This file is only written when 
specified on the command-line.

!!! note
    Once the calculation is done, the results can be analyzed in `gmx_MMPBSA_ana` (if `-nogui` flag was not used in the command-line). 
    Please, check the [gmx_MMPBSA_ana][5] section for more information

## References for `PM6-DH+`

`PM6-DH+` is the default `qm_theory` because protein-ligand, nucleic-acid-ligand, and carbohydrate interfaces
are dominated by hydrogen bonding and dispersion - interactions that plain PM3/PM6 treat poorly. Key references:

1. **Method development:** Řezáč & Hobza, *J. Chem. Theory Comput.* **2009**, 5, 1749-1760. [doi:10.1021/ct9000922](https://doi.org/10.1021/ct9000922) — PM6-DH dispersion/H-bond corrections; tested on DNA base pairs.
2. **PM6-DH+ H-bond correction:** Korth, *J. Chem. Theory Comput.* **2010**, 6, 3808-3816. [doi:10.1021/ct100408b](https://doi.org/10.1021/ct100408b)
3. **Protein-ligand review:** Grimme & Brandenburg, *Front. Chem.* **2015**, 3, 8. [PMC4881564](https://pmc.ncbi.nlm.nih.gov/articles/PMC4881564/) — SQM-DH methods (incl. PM6-DH+) for non-covalent interactions.
4. **QM/MM-GBSA benchmark (protein-carbohydrate):** Thapa *et al.*, *J. Phys. Chem. B* **2018**, 122, 7866-7878. [doi:10.1021/acs.jpcb.8b03655](https://doi.org/10.1021/acs.jpcb.8b03655)
5. **QM/MMGBSA with gmx_MMPBSA + PM6-DH+:** *Commun. Biol.* **2025**. [doi:10.1038/s42003-025-09143-z](https://doi.org/10.1038/s42003-025-09143-z)
6. **Host-guest binding with PM6-DH+:** Muddana & Gilson, *J. Chem. Theory Comput.* **2012**, 8, 2868-2880. [doi:10.1021/ct3002738](https://doi.org/10.1021/ct3002738)

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/QM_MMGBSA
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
