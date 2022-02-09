---
template: main.html
title: Output files
---

# Output Files

## The output file

This is how a typical output file ("FINAL_RESULTS_MMPBSA.dat" by default) looks like:

```  yaml title="Output file"
| Run on Tue Feb  8 22:31:58 2022                                                     | # (1) 
|                                                                                     +       
|gmx_MMPBSA Version=v1.4.3+462.gf64aa73 based on MMPBSA.py v.16.0                     | # (2) 
|Complex Structure file                                                  com.tpr      +       
|Complex (AMBER) topology file                                        COM.prmtop      |       
|Receptor (AMBER) topology file                                       REC.prmtop      | # (3) 
|Ligand Structure file                                               ligand.mol2      |       
|Complex (AMBER) topology file                                        LIG.prmtop      |       
|Initial trajectories                                             COM_traj_0.xtc      +       
|                                                                                             
|Receptor mask                   ":1-240"                                             + # (4) 
|Ligand mask                     ":241"                                               |       
|Ligand residue name is          "RAL"                                                +       
|                                                                                             
|Calculations performed using 16 complex frames                                       + # (5) 
|C2 Entropy Std. Dev. and Conf. Interv. (95%) have been obtained by bootstrapping     |       
|with number of re-samplings = 2000                                                   |       
|                                                                                     |       
|Generalized Born ESURF calculated using 'LCPO' surface areas                         |       
|                                                                                     |       
|Using temperature = 300.00 K)                                                        |       
|All units are reported in kcal/mol.                                                  |       
|                                                                                     |       
|SD - Sample standard deviation, SEM - Sample standard error of the mean              |       
|SD(Prop.), SEM(Prop.) - SD and SEM obtained with propagation of uncertainty formula  |       
|https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae            +
       
-------------------------------------------------------------------------------       +       
-------------------------------------------------------------------------------       |        
ENTROPY RESULTS (C2 ENTROPY)                                                          |       
Model           σ(Int. Energy)    C2 Value         Std. Dev.   Conf. Interv. (95%)    | # (6) 
-------------------------------------------------------------------------------       |       
gb                  3.308           9.176           2.086         4.601-12.477        |       
-------------------------------------------------------------------------------       |       
-------------------------------------------------------------------------------       +       
                                                                                              
GENERALIZED BORN:                                                                     | # (7)        
                                                                                       
Complex:                                                                              +              
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM       |              
-------------------------------------------------------------------------------       |       
BOND                    730.92         21.73      21.73         5.43       5.43       |       
ANGLE                  2022.96         27.54      27.54         6.89       6.89       |       
DIHED                  2631.46         15.99      15.99         4.00       4.00       |       
VDWAALS               -2035.53         14.81      14.81         3.70       3.70       |       
EEL                  -16750.43         22.60      22.60         5.65       5.65       |       
1-4 VDW                 911.73         15.17      15.17         3.79       3.79       | # (8) 
1-4 EEL               10292.96         23.03      23.03         5.76       5.76       |       
EGB                   -3270.77         15.90      15.90         3.98       3.98       |       
ESURF                    96.07          0.58       0.58         0.15       0.15       |       
                                                                                      |       
GGAS                  -2195.93         54.56      33.37        13.64       8.34       |       
GSOLV                 -3174.70         15.91      15.72         3.98       3.93       |       
                                                                                      |       
TOTAL                 -5370.63         56.84      31.88        14.21       7.97       + 
                                                                                       
                                                                                              
Receptor:                                                                             +              
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM       |              
-------------------------------------------------------------------------------       |       
BOND                    719.57         21.28      21.28         5.32       5.32       |       
ANGLE                  1996.26         27.60      27.60         6.90       6.90       |       
DIHED                  2597.25         13.84      13.84         3.46       3.46       |       
VDWAALS               -1972.61         13.52      13.52         3.38       3.38       |       
EEL                  -16735.20         22.33      22.33         5.58       5.58       |       
1-4 VDW                 895.00         14.81      14.81         3.70       3.70       | # (9) 
1-4 EEL               10339.15         23.01      23.01         5.75       5.75       |       
EGB                   -3288.63         16.59      16.59         4.15       4.15       |       
ESURF                    99.77          0.61       0.61         0.15       0.15       |       
                                                                                      |       
GGAS                  -2160.58         53.26      35.11        13.31       8.78       |       
GSOLV                 -3188.86         16.60      16.39         4.15       4.10       |       
                                                                                      |       
TOTAL                 -5349.43         55.78      32.13        13.95       8.03       +
                                                                                       
                                                                                              
Ligand:                                                                               +               
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM       |               
-------------------------------------------------------------------------------       |       
BOND                     11.35          1.94       1.94         0.48       0.48       |       
ANGLE                    26.70          2.68       2.68         0.67       0.67       |       
DIHED                    34.21          3.59       3.59         0.90       0.90       |       
VDWAALS                  -4.03          1.35       1.35         0.34       0.34       |       
EEL                      15.90          0.96       0.96         0.24       0.24       |       
1-4 VDW                  16.73          1.19       1.19         0.30       0.30       | # (10)
1-4 EEL                 -46.19          0.94       0.94         0.24       0.24       |       
EGB                     -23.14          0.63       0.63         0.16       0.16       |       
ESURF                     4.52          0.02       0.02         0.00       0.00       |       
                                                                                      |       
GGAS                     54.67          5.37       3.97         1.34       0.99       |       
GSOLV                   -18.61          0.63       0.63         0.16       0.16       |       
                                                                                      |       
TOTAL                    36.05          5.41       4.16         1.35       1.04       +
                                                                                      
                                                                                              
Delta (Complex - Receptor - Ligand):                                                  + # (11)        
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM       |       
-------------------------------------------------------------------------------       |       
ΔBOND                     0.00          1.49       0.00         0.37       0.00       | # (12)
ΔANGLE                    0.00          2.62       0.00         0.65       0.00       | # (13)
ΔDIHED                   -0.00          1.44       0.00         0.36       0.00       | # (14)
ΔVDWAALS                -58.89          0.06       2.31         0.01       0.58       | # (15)
ΔEEL                    -31.13          0.69       3.04         0.17       0.76       | # (16)
Δ1-4 VDW                  0.00          0.83       0.00         0.21       0.00       | # (17)
Δ1-4 EEL                 -0.00          0.92       0.00         0.23       0.00       | # (18)
ΔEGB                     40.99          0.05       1.41         0.01       0.35       | # (19)
ΔESURF                   -8.22          0.01       0.09         0.00       0.02       | # (20)
                                                                                      |       
ΔGGAS                   -90.02          0.69       3.31         0.17       0.83       | # (21)
ΔGSOLV                   32.77          0.06       1.37         0.01       0.34       | # (22)
                                                                                      |       
ΔTOTAL                  -57.25          0.69       2.55         0.17       0.64       + # (23)
-------------------------------------------------------------------------------       
-------------------------------------------------------------------------------               
Using C2 Entropy Approximation:                                                       |        
ΔG binding =  -48.0704 +/-  2.1972                                                    | # (24)
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
```

1. Date of running
2. gmx_MMPBSA version used
3. Input files used
4. Receptor and ligand masks
5. General description of the methods used and units 
6. Entropy results in case any entropy approximation was used
7. Model used (GB in this case)
8. Energy components (complex)
9. Energy components (receptor)
10. Energy components (ligand)
11. Energy components (delta)
12. Bond potential term
13. Angle potential term
14. Dihedral potential term
15. Van der Waals contribution
16. Electrostatic contribution
17. Van der Waals 1-4 contribution
18. Electrostatic 1-4 contribution
19. Polar contribution to the solvation free energy
20. Non-polar contribution to the solvation free energy
21. = ΔBOND + ΔANGLE + ΔDIHED + ΔVDWAALS + ΔEEL
22. = ΔEGB + ΔESURF
23. = ΔGGAS + ΔGSOLV
24. Binding free energy<br>ΔG binding = ΔTOTAL - TΔS

The header of the output file will contain information about the calculation. It will also show the names of all files 
that were used in the calculation (topology files and coordinate file(s)). If the masks
were not specified, it prints its best guess so that you can verify its accuracy, along with the residue name of the
ligand (if it is only a single residue). After that, general information about methods, units, constants used is 
included. Entropy results are shown next in case any entropy approximation was used. Next, the energy and entropy 
contributions are broken up into their components as they are in `sander` and `nmode` or `cpptraj`. The contributions 
are further broken for the complex, receptor and ligand into `GGAS` and `GSOLV`. `GGAS` is the interaction energy 
and is obtained after sum the internal(bonded) components (`BOND` + `ANGLE` + `DIHED`) and the 
non-bonded (`VDWAALS` + `EEL`) components. For `GSOLV`, the polar and non-polar contributions are `EGB` (or `EPB`) 
and `ESURF` (or `ENPOLAR + EDISPER`), respectively for `GB` (or `PB`) calculations. A single trajectory protocol does 
not produce any differences between bond lengths, angles, dihedrals or 1-4 interactions between the complex and 
receptor/ligand structures. Thus, when subtracted they cancel completely. If not, these values are displayed and 
inconsistency warnings are printed. When this occurs the results are generally useless. Of course this does not 
hold for the multiple trajectory protocol as independent trajectories are used for the complex, receptor and ligand. 
Two approaches are used when calculating the standard deviation, and the standard error of the mean. The `SD` and `SEM`
values are calculated using a sample (array) of values. On the other hand, `SD(Prop.)` and `SEM(Prop.)` are 
obtained with the 
[propagation of uncertainty formula](https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae) for 
_f_ = _A_ - _B_. Check this [thread](http://archive.ambermd.org/201202/0317.html) for more details on MM/PB(GB)SA 
statistics.

## Temporary files

!!! warning
    This section does not record all the temporary files that are currently generated.

`gmx_MMPBSA` creates working files during the execution of the script beginning with the prefix `_GMXMMPBSA_`.
If `gmx_MMPBSA` does not finish successfully, several of these files may be helpful in diagnosing the problem. For that
reason, every temporary file is described below. Note that not every temporary file is generated in every simulation. At
the end of each description, the lowest value of the original “keep_files” variable that will retain this file will be
shown in parentheses. Nevertheless, in the current version, all the files are retained for plotting purposes.

`gmx_MMPBSA.log` This file contains the output coming from `gmx_MMPBSA`.

`leap.log` This file contains the output coming from tleap program.

`_GMXMMPBSA_gb.mdin` Input file that controls the GB calculation done in sander. (2)

`_GMXMMPBSA_pb.mdin` Input file that controls the PB calculation done in sander. (2)

`_GMXMMPBSA_gb_decomp_com.mdin` Input file that controls the GB decomp calculation for the complex done in sander. (2)

`_GMXMMPBSA_gb_decomp_rec.mdin` Input file that controls the GB decomp calculation for the receptor done in sander. (2)

`_GMXMMPBSA_gb_decomp_lig.mdin` Input file that controls the GB decomp calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_pb_decomp_com.mdin` Input file that controls the PB decomp calculation for the complex done in sander. (2)

`_GMXMMPBSA_pb_decomp_rec.mdin` Input file that controls the PB decomp calculation for the receptor done in sander. (2)

`_GMXMMPBSA_pb_decomp_lig.mdin` Input file that controls the PB decomp calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_com.mdin` Input file that controls the GB QM/MM calculation for the complex done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_rec.mdin` Input file that controls the GB QM/MM calculation for the receptor done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_lig.mdin` Input file that controls the GB QM/MM calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_complex.mdcrd.#` Trajectory file(s) that contains only those complex snapshots that will be processed by
MMPBSA.py. (1)

`_GMXMMPBSA_ligand.mdcrd.#` Trajectory file(s) that contains only those ligand snapshots that will be processed by
MMPBSA.py. (1)

`_GMXMMPBSA_receptor.mdcrd.#` Trajectory file(s) that contains only those receptor snapshots that will be processed by
MMPBSA.py. (1)

`_GMXMMPBSA_complex_nc.#` Same as _GMXMMPBSA_complex.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_receptor_nc.#` Same as _GMXMMPBSA_receptor.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_ligand_nc.#` Same as _GMXMMPBSA_ligand.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_dummycomplex.inpcrd` Dummy inpcrd file generated by _GMXMMPBSA_complexinpcrd.in for use with imin=5
functionality in sander. (1)

`_GMXMMPBSA_dummyreceptor.inpcrd` Same as above, but for the receptor. (1)

`_GMXMMPBSA_dummyligand.inpcrd` Same as above, but for the ligand. (1)

`_GMXMMPBSA_complex.pdb` Dummy PDB file of the complex required to set molecule up in nab programs

`_GMXMMPBSA_receptor.pdb` Dummy PDB file of the receptor required to set molecule up in nab programs

`_GMXMMPBSA_ligand.pdb` Dummy PDB file of the ligand required to set molecule up in nab programs

`_GMXMMPBSA_complex_nm.mdcrd.#` Trajectory file(s) for each thread with snapshots used for normal mode calcula- tions on
the complex. (1)

`_GMXMMPBSA_receptor_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calcula- tions on
the receptor. (1)

`_GMXMMPBSA_ligand_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calculations on the
ligand. (1)

`_GMXMMPBSA_ptrajentropy.in` Input file that calculates the entropy via the quasi-harmonic approximation. This file is
processed by ptraj. (2)

`_GMXMMPBSA_avgcomplex.pdb` PDB file containing the average positions of all complex conformations processed by

`_GMXMMPBSA_cenptraj.in.` It is used as the reference for the _GMXMMPBSA_ptrajentropy.in file above.
(1)

`_GMXMMPBSA_complex_entropy.out` File into which the entropy results from _GMXMMPBSA_ptrajentropy.in analysis on the
complex are dumped. (1)

`_GMXMMPBSA_receptor_entropy.out` Same as above, but for the receptor. (1)

`_GMXMMPBSA_ligand_entropy.out` Same as above, but for the ligand. (1)

`_GMXMMPBSA_ptraj_entropy.out` Output from running ptraj using _GMXMMPBSA_ptrajentropy.in. (1)

`_GMXMMPBSA_complex_gb.mdout.#` sander output file containing energy components of all complex snapshots done in GB. (1)

`_GMXMMPBSA_receptor_gb.mdout.#` sander output file containing energy components of all receptor snapshots done in GB. (
1)

`_GMXMMPBSA_ligand_gb.mdout.#` sander output file containing energy components of all ligand snapshots done in GB. (1)

`_GMXMMPBSA_complex_pb.mdout.#` sander output file containing energy components of all complex snapshots done in PB. (1)

`_GMXMMPBSA_receptor_pb.mdout.#` sander output file containing energy components of all receptor snapshots done in PB. (
1)

`_GMXMMPBSA_ligand_pb.mdout.#` sander output file containing energy components of all ligand snapshots done in PB. (1)

`_GMXMMPBSA_complex_rism.out.#` rism3d.snglpnt output file containing energy components of all complex snap- shots done
with 3D-RISM (1)

`_GMXMMPBSA_receptor_rism.out.#` rism3d.snglpnt output file containing energy components of all receptor snap- shots
done with 3D-RISM (1)

`_GMXMMPBSA_ligand_rism.out.#` rism3d.snglpnt output file containing energy components of all ligand snapshots done with
3D-RISM (1)

`_GMXMMPBSA_pbsanderoutput.junk.#` File containing the information dumped by sander.APBS to STDOUT. (1)

`_GMXMMPBSA_ligand_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the ligand for all
snapshots. (1)

`_GMXMMPBSA_receptor_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the receptor for
all snapshots. (1)

`_GMXMMPBSA_complex_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the com- plex for
all snapshots. (1)

`_GMXMMPBSA_mutant_...` These files are analogs of the files that only start with `_GMXMMPBSA_` described above, but
instead refer to the mutant system of alanine scanning calculations.

`_GMXMMPBSA_*out.#` These files are thread-specific files. For serial simulations, only #=0 files are created. For
parallel, #=0 through NUM_PROC - 1 are created.