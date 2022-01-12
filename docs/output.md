---
template: main.html
title: Output files
---

# Output Files

## The output file

This is how a typical output file looks like:

```
| Run on Wed Jun 16 02:53:07 2021                                                   |---> Date of running
|
|Input file:
|--------------------------------------------------------------                     | 
|Sample input file for GB calculation                                               |
|#This input file is meant to show only that gmx_MMPBSA works. Although, we         |
|#tried to use the input files as recommended in the                                |
|#Amber manual, some parameters have been changed to perform more expensive         |
|#calculations in a reasonable amount of time. Feel free to change the              |
|#parameters according to what is better for your system.                           |
|&general                                                                           |
|startframe=301, endframe=1000, verbose=2, PBRadii=3, interval=7,                   |---> Input file (*.in) 
|/                                                                                  |
|&gb                                                                                |
|igb=2, saltcon=0.150,                                                              |
|/                                                                                  |
|&pb                                                                                |
|istrng=0.15, fillratio=4.0, radiopt=0,                                             |
|/                                                                                  |
|--------------------------------------------------------------
|gmx_MMPBSA Version=v1.4.3+6.g25a02b8 based on MMPBSA.py v.16.0                     |---> gmx_MMPBSA version used
|Complex topology file:           COM.prmtop                                        |
|Receptor topology file:          REC.prmtop                                        |
|Ligand topology file:            LIG.prmtop                                        |
|Initial mdcrd(s):                COM_traj_0.xtc                                    |
|                                                                                   
|Receptor mask:                  ":1-266"                                           |  
|Ligand mask:                    ":267"                                             |---> receptor and ligand information
|Ligand residue name is "MFU"                                                       |
|                                                                                   
|Calculations performed using 100 complex frames.                                   |
|Poisson Boltzmann calculations performed using internal PBSA solver in sander.     |
|                                                                                   |---> general description of the 
|Generalized Born ESURF calculated using 'LCPO' surface areas                       |     method used and units
|                                                                                   |
|All units are reported in kcal/mole.                                               |
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

GENERALIZED BORN:                                                                   |---> results for GB calculation

Complex:
Energy Component            Average              Std. Dev.   Std. Err. of Mean      |
-------------------------------------------------------------------------------     |
BOND                       759.4785               23.3077              2.3308       |
ANGLE                     1803.1508               30.8183              3.0818       |
DIHED                     3216.1217               18.5641              1.8564       |
VDWAALS                  -2085.3102               15.8085              1.5808       |
EEL                     -18735.2906               44.2758              4.4276       |
1-4 VDW                    873.7890               11.3567              1.1357       |
1-4 EEL                  11088.2966               34.0205              3.4021       |---> Energetic components (complex)
EGB                      -2349.7950               20.1898              2.0190       |
ESURF                       74.5964                0.4837              0.0484       |
                                                                                    |
G gas                    -3079.7641               46.8525              4.6852       |
G solv                   -2275.1986               20.1875              2.0187       |
                                                                                    |
TOTAL                    -5354.9627               39.7488              3.9749       |


Receptor:
Energy Component            Average              Std. Dev.   Std. Err. of Mean      |                                   
-------------------------------------------------------------------------------     |                                   
BOND                       752.3559               22.9770              2.2977       |                                   
ANGLE                     1787.7219               30.5141              3.0514       |                                   
DIHED                     3211.0390               18.6469              1.8647       |                                   
VDWAALS                  -2066.0693               15.3943              1.5394       |                                   
EEL                     -18584.3471               43.8360              4.3836       |                                   
1-4 VDW                    873.7890               11.3567              1.1357       |                                   
1-4 EEL                  11088.2966               34.0205              3.4021       |---> Energetic components (receptor)
EGB                      -2365.3357               20.2920              2.0292       |                                   
ESURF                       75.1071                0.4790              0.0479       |                                   
                                                                                    |                                   
G gas                    -2937.2141               46.8804              4.6880       |                                   
G solv                   -2290.2286               20.2874              2.0287       |                                   
                                                                                    |                                   
TOTAL                    -5227.4426               39.5203              3.9520       |                                   


Ligand:
Energy Component            Average              Std. Dev.   Std. Err. of Mean      |                                   
-------------------------------------------------------------------------------     |                                   
BOND                         7.1226                2.2202              0.2220       |                                   
ANGLE                       15.4289                2.9490              0.2949       |                                   
DIHED                        5.0827                0.5376              0.0538       |                                   
VDWAALS                     -1.9909                0.3267              0.0327       |                                   
EEL                        -96.9811                2.5449              0.2545       |                                   
EGB                        -25.3239                1.3225              0.1322       |                                   
ESURF                        2.2611                0.0154              0.0015       |---> Energetic components (ligand)
                                                                                    |                                   
G gas                      -71.3377                5.0695              0.5070       |                                   
G solv                     -23.0628                1.3164              0.1316       |                                   
                                                                                    |                                   
TOTAL                      -94.4005                4.4250              0.4425       |                                   
                                                                                                                       
                                                                                                                       
Differences (Complex - Receptor - Ligand):                                          |---> Differences in Energetic components
Energy Component            Average              Std. Dev.   Std. Err. of Mean      |
-------------------------------------------------------------------------------     |
BOND                        -0.0000                0.0000              0.0000       |---> bond potential term
ANGLE                       -0.0000                0.0000              0.0000       |---> angle potential term
DIHED                        0.0000                0.0001              0.0000       |---> dihedral potential term
VDWAALS                    -17.2500                2.4428              0.2443       |---> van der Waals contribution
EEL                        -53.9624                4.2749              0.4275       |---> electrostatic contribution
1-4 VDW                      0.0000                0.0000              0.0000       |---> van der Waals 1-4 contribution
1-4 EEL                      0.0000                0.0000              0.0000       |---> electrostatic 1-4 contribution
EGB                         40.8646                2.0983              0.2098       |---> polar contribution to the solvation free energy
ESURF                       -2.7718                0.0462              0.0046       |---> non-polar contribution to the solvation free energy 
                                                                                    |
DELTA G gas                -71.2123                3.3059              0.3306       |---> = BOND + ANGLE + DIHED + VDWAALS + EEL
DELTA G solv                38.0928                2.1012              0.2101       |---> = EGB + ESURF
                                                                                    |
DELTA TOTAL                -33.1195                2.1204              0.2120       |---> = DELTA G gas + DELTA G solv


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
```

The header of the output file will contain information about the calculation. It will show a copy of the input file as
well as the names of all files that were used in the calculation (topology files and coordinate file(s)). If the masks
were not specified, it prints its best guess so that you can verify its accuracy, along with the residue name of the
ligand (if it is only a single residue). The energy and entropy contributions are broken up into their components as
they are in sander and `nmode` or `ptraj`. The contributions are further broken into `G gas` and `G solv`. The polar and
non-polar contributions are `EGB` (or `EPB`) and `ESURF` (or `ECAVITY` / `ENPOLAR`), respectively for `GB` (or `PB`) 
calculations. A single trajectory does not produce any differences between bond lengths, angles, or dihedrals
between the complex and receptor/ligand structures. Thus, when subtracted they cancel completely. This includes the
`BOND`, `ANGLE`, `DIHED`, and `1-4 interactions`. If inconsistencies are found, these values are displayed and inconsistency
warnings are printed. When this occurs the results are generally useless. Of course this does not hold for the multiple
trajectory protocol, and so all energy components are printed in this case. Finally, all warnings generated during the
calculation that do not result in fatal errors are printed after calculation details but before any results.


## Temporary Files

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