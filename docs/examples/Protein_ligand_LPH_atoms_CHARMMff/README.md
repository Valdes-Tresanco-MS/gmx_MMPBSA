---
template: main.html
title: Protein-ligand LPH (Charmm)
---

!!! danger "CHARMM and MM(PB/GB)SA"
    PB model is recommended when working with CHARMMff files. Nevertheless, the combination of PB/GB models and 
    CHARMM force field hasn't been tested extensively. Please, check this [thread][1] for more information and 
    proceed with caution.

# Protein-ligand with LPH atoms BFE calculations (Single Trajectory method) -- CHARMMff files

!!! info

    This example can be found in the [docs/examples/Protein_ligand_LPH_atoms_CHARMMff][6] directory in the repository
    folder

    LPH is a positively charged virtual particle attached to halogen atoms. This strategy aims to get a better 
    representation of the halogen bond which is a highly directional, non-covalent interaction between a halogen atom 
    and another electronegative atom (See [here][8] for more info). Unfortunately, including these particles in the 
    topology will cause gmx_MMPBSA to end in an error. However, there is a way to generate the files without these 
    particles and get gmx_MMPBSA up and running.

    
!!! danger "Keep in mind"

    As the LPH particle is not considered during the calculations in gmx_MMPBSA, take the results with a grain of 
    salt, especially when working with systems where the halogen bond is determinant for the binding.


## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`             | input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb` `gro`     | Structure file containing the system coordinates|
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`          | Receptor and ligand group numbers in the index file |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `gro` `trr` | final GROMACS MD trajectory, fitted and with no pbc.|
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`            | take into account that *.itp files belonging to the topology file should be also present in the folder       |
| A Reference Structure file     | :octicons-check-circle-fill-16:{ .req_optrec .scale_icon_medium } |           `pdb`            |  Complex reference structure file with correct assignment of chain ID and residue numbers       |
              
:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][2]_

In order to generate the corresponding files (The MD Structure+mass(db), index, trajectory and the topology files) 
without the LPH particles, it's necessary to run a few commands. Bear with me!

Let's generate the index file first:

!!! important
    The main idea here is to generate a receptor group, a 
    ligand group without the LPH particles and a complex group containing both the receptor and the ligand without 
    the LPH particles. In general, index files generated with GROMACS directly will contain more detailed information 
    (_i.e._, receptor and ligand separated)

```
gmx make_ndx -f com.tpr -o index_mod_gromacs.ndx

  0 System              : 70483 atoms
  1 Protein             :  5580 atoms
  2 Protein-H           :  2817 atoms
  3 C-alpha             :   334 atoms
  4 Backbone            :  1002 atoms
  5 MainChain           :  1335 atoms
  6 MainChain+Cb        :  1654 atoms
  7 MainChain+H         :  1654 atoms
  8 SideChain           :  3926 atoms
  9 SideChain-H         :  1482 atoms
 10 Prot-Masses         :  5580 atoms
 11 non-Protein         : 64903 atoms
 12 Other               : 64903 atoms
 13 3G5                 :    32 atoms
 14 CLA                 :    62 atoms
 15 SOD                 :    63 atoms
 16 TIP3                : 64746 atoms

Splitting the ligand (group 13) by atoms
>splitat 13

Grouping both LPH particles
>47|48

Excluding both LPH particles from the ligand
>13&!49

Naming ligand as lig
>name 50 lig

Grouping rec and lig
>1|50

Cleaning
>del 17-49

save and quit
>q

This is how it should look like at the end

  0 System              : 70483 atoms
  1 Protein             :  5580 atoms
  2 Protein-H           :  2817 atoms
  3 C-alpha             :   334 atoms
  4 Backbone            :  1002 atoms
  5 MainChain           :  1335 atoms
  6 MainChain+Cb        :  1654 atoms
  7 MainChain+H         :  1654 atoms
  8 SideChain           :  3926 atoms
  9 SideChain-H         :  1482 atoms
 10 Prot-Masses         :  5580 atoms
 11 non-Protein         : 64903 atoms
 12 Other               : 64903 atoms
 13 3G5                 :    32 atoms
 14 CLA                 :    62 atoms
 15 SOD                 :    63 atoms
 16 TIP3                : 64746 atoms
 17 lig                 :    30 atoms
 18 Protein_lig         :  5610 atoms
```

!!! note
    Note that the number of atoms in the generated complex is 5610 because it doesn't include the LPH particles.

Let's generate the MD Structure+mass(db) file:

    echo 18 | gmx trjconv -s com.tpr -f traj_fit.xtc -dump 0 -o str_noLP.pdb -n index_mod_gromacs.ndx

Open `str_noLP.pdb` in your favorite visualizer and see it doesn't contain the LPH particles. Now, let's generate 
the trajectory with no LPH particles:

    echo 18 | gmx trjconv -s com.tpr -f traj_fit.xtc -o com_traj.xtc -n index_mod_gromacs.ndx

Finally, let's edit the topology file. Go inside the toppar folder and open the `HETA.itp` file. As you will see, we
deleted all the information related with LPH particles (atom numbers 31, and 32 respectively). In this case, we
deleted the information for LPH particles in `atoms` (lines 47, 48) and `pairs` (lines 124, 132, 150, 153, 154, 157, 
158, 159). Besides, delete the whole `[ virtual_sites3 ]` (lines 296-299) and `[ exclusions ]` (lines 301-318) 
fields. The original .itp (`HETA_original_with_LPH_info.itp`) is included for comparison purposes.


## Command-line
That being said, once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA -O -i mmpbsa.in -cs str_noLP.pdb -ci index_mod_gromacs.ndx -cg 1 17 -ct com_traj.xtc -cp topol.top

=== "With MPI"

        mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs str_noLP.pdb -ci index_mod_gromacs.ndx -cg 1 17 -ct com_traj.xtc -cp topol.top

where the `mmpbsa.in` input file, is a text file containing the following lines:

``` linenums="1"
Sample input file for PB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-Lig-ST",
startframe=1, endframe=11, verbose=2, solvated_trajectory=0,
/
&pb
# radiopt=0 is recommended which means using radii from the prmtop file
# for both the PB calculation and for the NP calculation

istrng=0.15, fillratio=4.0, radiopt=0, inp=1,
/
```

!!! warning "Remember"
    `radiopt = 0` is recommended which means using radii from the `prmtop` file

_See a detailed list of all the options in `gmx_MMPBSA` input file [here][3] as well as several [examples][4]_

## Considerations
In this case, a single trajectory (ST) approximation is followed, which means the receptor and ligand structures and 
trajectories will be obtained from that of the complex. To do so, a MD Structure+mass(db) file (`str_noLP.pdb`), an 
index file (`index_mod.ndx`), a trajectory file (`com_traj.xtc`), and both the receptor and ligand group numbers in the 
index file (`3 4`) are needed. The `mmpbsa.in` input file will contain all the parameters needed for the MM/PB(GB)SA 
calculation. A topology file is also needed (mandatory) in this case to generate the topology files in amber format 
with all the terms for CHARMM force field.
!!! note
    Once the calculation is done, you can analyze the results in `gmx_MMPBSA_ana` (if you didn't define `-nogui`). 
    Please see the [gmx_MMPBSA_ana][5] section for more information


  [1]: http://archive.ambermd.org/201508/0382.html 
  [2]: ../../command-line.md#gmx_mmpbsa-command-line
  [3]: ../../input_file.md#the-input-file
  [4]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples/Protein_ligand_LPH_atoms_CHARMMff
  [7]: ../../command-line.md#gmx_mmpbsa_test-command-line
  [8]: https://www.sciencedirect.com/science/article/abs/pii/S0968089616304576?via%3Dihub