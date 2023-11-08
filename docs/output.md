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

This is how a typical decomposition output file (`FINAL_DECOMP_MMPBSA.dat` by default) looks like:

```  yaml title="FINAL_DECOMP_MMPBSA.dat"  
|Run on Tue Mar  9 23:48:23 2021                                                  | # (1)!
|GB non-polar solvation energies calculated with gbsa=2							  + # (2)!												
|idecomp = 2 Per-residue decomp adding 1-4 interactions to EEL and VDW.			  | 																
|Energy Decomposition Analysis (All units kcal/mol) Generalized Born solvent      +

Complex:                         | # (3)!																		
Total Energy Decomposition:      | # (4)!
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   22.06	3.40	    0.82	            -6.83	1.27	    0.31	            -32.53	0.81	    0.20	            -1.57	0.38	    0.09	            0.01	0.00	    0.00	            -18.86	3.27	    0.79	
THR  41	   19.15	2.58	    0.63	            -3.64	1.06	    0.26	            -36.16	0.96	    0.23	            -3.58	0.37	    0.09	            0.23	0.02	    0.00	            -23.99	1.97	    0.48	
ALA  44	   15.67	1.92	    0.47	            -5.42	0.72	    0.18	            -7.54	0.90	    0.22	            -2.36	0.33	    0.08	            0.00	0.00	    0.00	            0.35	2.42	    0.59	
RAL 241	   72.32	3.91	    0.95	            -16.88	2.66	    0.65	            -46.00	1.94	    0.47	            -0.74	0.89	    0.22	            0.42	0.02	    0.00	            9.12	3.79	    0.92	

Sidechain Energy Decomposition:	 | # (5)!																	
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   9.64	    2.16	    0.52	            -4.15	0.92	    0.22	            -27.51	0.45	    0.11	            4.14	0.11	    0.03	            0.01	0.00	    0.00	            -17.88	2.48	    0.60	
THR  41	   6.60	    1.19	    0.29	            -1.62	0.70	    0.17	            -26.20	0.93	    0.23	            2.38	0.25	    0.06	            0.23	0.02	    0.00	            -18.61	1.11	    0.27	
ALA  44	   2.94	    1.08	    0.26	            -2.20	0.26	    0.06	            2.19	0.13	    0.03	            -0.31	0.04	    0.01	            0.00	0.00	    0.00	            2.63	1.25	    0.30	
RAL 241	   72.32	3.91	    0.95	            -16.88	2.66	    0.65	            -46.00	1.94	    0.47	            -0.74	0.89	    0.22	            0.42	0.02	    0.00	            9.12	3.79	    0.92	

Backbone Energy Decomposition:	 | # (6)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   12.42	2.04	    0.49	            -2.67	0.66	    0.16	            -5.03	0.80	    0.19	            -5.70	0.43	    0.11	            0.00	0.00	    0.00	            -0.98	1.72	    0.42	
THR  41	   12.55	1.83	    0.44	            -2.02	0.71	    0.17	            -9.95	0.98	    0.24	            -5.96	0.34	    0.08	            0.00	0.00	    0.00	            -5.38	1.59	    0.39	
ALA  44	   12.72	1.21	    0.29	            -3.22	0.73	    0.18	            -9.73	0.87	    0.21	            -2.05	0.33	    0.08	            0.00	0.00	    0.00	            -2.28	1.44	    0.35	
RAL 241	   0.00	    0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	
																			
Receptor:                        | # (7)!																			
Total Energy Decomposition:	     | # (8)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   22.06	3.40	    0.82	            -4.37	1.27	    0.31	            -32.08	0.79	    0.19	            -2.68	0.39	    0.09	            0.23	0.02	    0.01	            -16.84	3.30	    0.80	
THR  41	   19.15	2.58	    0.63	            -1.55	1.04	    0.25	            -35.55	0.96	    0.23	            -4.84	0.33	    0.08	            0.52	0.02	    0.00	            -22.28	2.03	    0.49	
ALA  44	   15.67	1.92	    0.47	            -3.91	0.72	    0.17	            -7.61	0.86	    0.21	            -2.02	0.35	    0.08	            0.21	0.02	    0.00	            2.33	2.36	    0.57	

Sidechain Energy Decomposition:	 | # (9)!																		
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   9.64	    2.16	    0.52	            -2.70	0.89	    0.21	            -27.27	0.44	    0.11	            3.89	0.09	    0.02	            0.20	0.02	    0.01	            -16.25	2.49	    0.60	
THR  41	   6.60	    1.19	    0.29	            -0.40	0.70	    0.17	            -25.86	0.88	    0.21	            1.62	0.21	    0.05	            0.49	0.02	    0.00	            -17.55	1.10	    0.27	
ALA  44	   2.94	    1.08	    0.26	            -1.33	0.18	    0.04	            2.24	0.12	    0.03	            -0.43	0.04	    0.01	            0.20	0.01	    0.00	            3.62	1.17	    0.28	

Backbone Energy Decomposition:	 | # (10)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	   12.42	2.04	    0.49	            -1.66	0.66	    0.16	            -4.81	0.81	    0.20	            -6.57	0.42	    0.10	            0.03	0.01	    0.00	            -0.59	1.69	    0.41	
THR  41	   12.55	1.83	    0.44	            -1.15	0.72	    0.17	            -9.70	0.95	    0.23	            -6.46	0.34	    0.08	            0.02	0.01	    0.00	            -4.73	1.65	    0.40	
ALA  44	   12.72	1.21	    0.29	            -2.59	0.71	    0.17	            -9.85	0.84	    0.20	            -1.58	0.35	    0.08	            0.01	0.01	    0.00	            -1.28	1.46	    0.35	
																			
Ligand:                          | # (11)!																			
Total Energy Decomposition:	     | # (12)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
RAL   1	   72.32	3.91	    0.95	            12.61	2.08	    0.51	            -30.28	0.65	    0.16	            -23.18	0.63	    0.15	            5.65	0.03	    0.01	            37.11	4.05	    0.98	

Sidechain Energy Decomposition:	 | # (13)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
RAL   1	   72.32	3.91	    0.95	            12.61	2.08	    0.51	            -30.28	0.65	    0.16	            -23.18	0.63	    0.15	            5.65	0.03	    0.01	            37.11	4.05	    0.98	

Backbone Energy Decomposition:	 | # (14)!																			
Residue	   Internal			                        van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	       Avg.	    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
RAL   1	   0.00	    0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	
																			
DELTAS:                          | # (15)!																			
Total Energy Decomposition:	     | # (16)!																			
Residue	Location    Internal			                    van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	                Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	R LEU  40	0.00	0.00	    0.00	            -2.46	0.21	    0.05	            -0.45	0.13	    0.03	            1.11	0.18	    0.04	            -0.22	0.02	    0.01	            -2.02	0.23	    0.06
THR  41	R THR  41	0.00	0.00	    0.00	            -2.09	0.15	    0.04	            -0.60	0.12	    0.03	            1.26	0.15	    0.04	            -0.28	0.02	    0.00	            -1.72	0.17	    0.04
ALA  44	R ALA  44	0.00	0.00	    0.00	            -1.50	0.13	    0.03	            0.07	0.05	    0.01	            -0.34	0.04	    0.01	            -0.21	0.02	    0.00	            -1.98	0.17	    0.04
RAL 241	L RAL   1	0.00	0.00	    0.00	            -29.49	1.14	    0.28	            -15.72	1.60	    0.39	            22.44	0.74	    0.18	            -5.22	0.04	    0.01	            -28.00	1.36	    0.33
																			
Sidechain Energy Decomposition:	 | # (17)!																			
Residue	Location    Internal			                    van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	                Avg.	td. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	R LEU  40	0.00	0.00	    0.00	            -1.45	0.17	    0.04	            -0.23	0.10	    0.02	            0.25	0.05	    0.01	            -0.19	0.02	    0.01	            -1.63	0.15	    0.04
THR  41	R THR  41	0.00	0.00	    0.00	            -1.22	0.13	    0.03	            -0.34	0.12	    0.03	            0.76	0.15	    0.04	            -0.26	0.02	    0.00	            -1.06	0.13	    0.03
ALA  44	R ALA  44	0.00	0.00	    0.00	            -0.87	0.15	    0.04	            -0.05	0.02	    0.01	            0.13	0.02	    0.00	            -0.19	0.01	    0.00	            -0.99	0.16	    0.04
RAL 241	L RAL   1	0.00	0.00	    0.00	            -29.49	1.14	    0.28	            -15.72	1.60	    0.39	            22.44	0.74	    0.18	            -5.22	0.04	    0.01	            -28.00	1.36	    0.33
																			
Backbone Energy Decomposition:	 | # (18)!																			
Residue	Location    Internal			                    van der Waals			                Electrostatic			                Polar Solvation			                Non-Polar Solv.			                TOTAL			
	                Avg.    Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean	Avg.	Std. Dev.	Std. Err. of Mean
LEU  40	R LEU  40	0.00	0.00	    0.00	            -1.01	0.08	    0.02	            -0.22	0.10	    0.02	            0.86	0.18	    0.04	            -0.03	0.01	    0.00	            -0.40	0.11	    0.03
THR  41	R THR  41	0.00	0.00	    0.00	            -0.87	0.10	    0.02	            -0.26	0.06	    0.02	            0.50	0.04	    0.01	            -0.02	0.01	    0.00	            -0.66	0.16	    0.04
ALA  44	R ALA  44	0.00	0.00	    0.00	            -0.63	0.06	    0.01	            0.12	0.04	    0.01	            -0.47	0.04	    0.01	            -0.01	0.01	    0.00	            -1.00	0.08	    0.02
RAL 241	L RAL   1	0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00	            0.00	0.00	    0.00
```

1. Date of running
2. General description of the methods used and units
3. Data for the complex
4. Total Decomposition (TDC) by term for the complex (= Sidechain Energy Decomposition + Backbone Energy Decomposition)
5. Sidechain Energy Decomposition (SDC) by term for the complex
6. Backbone Energy Decomposition (BDC) by term for the complex
7. Data for the receptor
8. Total Decomposition (TDC) by term for the receptor (= Sidechain Energy Decomposition + Backbone Energy Decomposition)
9. Sidechain Energy Decomposition (SDC) by term for the receptor
10. Backbone Energy Decomposition (BDC) by term for the receptor
11. Data for the ligand
12. Total Decomposition (TDC) by term for the ligand (= Sidechain Energy Decomposition + Backbone Energy Decomposition)
13. Sidechain Energy Decomposition (SDC) by term for the ligand
14. Backbone Energy Decomposition (BDC) by term for the ligand
15. Delta Energies
16. Delta energy for the Total Decomposition (TDC) by term
17. Delta energy for the Sidechain Energy Decomposition (SDC) by term
18. Delta energy for the Backbone Energy Decomposition (BDC) by term

The header of the output file will contain information about the calculation and parameters specified. Next, the TDC,
SDC, and BDC data are shown for the complex, receptor and ligand, respectively. Finally, the delta energies are shown 
by terms for TDC, SDC, and BDC, respectively.

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