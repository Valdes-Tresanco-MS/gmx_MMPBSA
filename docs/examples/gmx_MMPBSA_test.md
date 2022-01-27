---
template: main.html
title: gmx_MMPBSA_test
---

## `gmx_MMPBSA_test` command-line
```
$ gmx_MMPBSA_test -h
usage: gmx_MMPBSA_test [-h] [-v] 
       [-t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18} [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18} ...]]] 
       [-f FOLDER] [-r] [-ng] [-n NUM_PROCESSORS]

This program is part of gmx_MMPBSA and will allow you to run various gmx_MMPBSA examples easily.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Test options:
  -t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18} [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18} ...]]
                        The level the test is going to be run at. Multiple systems and analysis can be run at the same 
                        time.
                              Nr. of Sys  
                        * 0      16     All -- Run all examples (Can take a long time!!!)
                        * 1      13     Minimal -- Does a minimal test with a set of systems and analyzes 
                                        that show that gmx_MMPBSA runs correctly. Only exclude 3drism, nmode
                                        protein-ligand MT because take a long time or are redundant
                        * 2       9     Fast -- Only the calculations that take a short time are run (Default)
                        [Systems]:
                             Slow Frames
                        * 3    . | 10   Protein-Ligand (Single trajectory approximation)
                        * 4    . | 10   Protein-Protein
                        * 5    . | 10   Protein-DNA
                        * 6    x |  4   Protein-Membrane
                        * 7    . | 10   Protein-Glycan
                        * 8    x |  4   Metalloprotein-Peptide
                        * 9    . | 10   Protein-DNA-RNA-IONs-Ligand
                        * 10   x |  4   Protein-Ligand (CHARMM force field)
                        * 11   x |  4   Protein-ligand complex in membrane with CHARMMff 
                        [Analysis]:
                             Slow Frames
                        * 12   . | 10   Alanine Scanning
                        * 13   . | 10   Stability calculation
                        * 14   . | 10   Decomposition Analysis
                        * 15   . | 16   Interaction Entropy approximation
                        * 16   . | 10   Protein-Ligand (Multiple trajectory approximation)
                        * 17   x |  4   Entropy calculation using Normal Mode approximation 
                        * 18   x |  4   Calculations using 3D-RISM approximation
  -f FOLDER, --folder FOLDER
                        Defines the folder to store all data
  -r, --reuse           Defines the existing test forlder will be reuse
  -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
  -n NUM_PROCESSORS, --num_processors NUM_PROCESSORS
                        Defines the number of processor cores you want to use with MPI per calculation. If the number 
                        of frames is less than the number of cpus defined, the calculation will be performed with 
                        the number of processors = number of frames.


gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
Based on MMPBSA.py (version 16.0) and AmberTools20
```

## Running gmx_MMPBSA_test
gmx_MMPBSA_test is designed to run a set of samples (all or minimal) or a specific example efficiently. 
Additionally, gmx_MMPBSA_test can run in parallel, decreasing the execution time gmx_MMPBSA_test will download the 
most recent version of the repository in the specified folder and will perform the calculations

=== "Minimal"
    
        gmx_MMPBSA_test -f /home/user/Documents -n 10
    
    Through this command-line, gmxMMPBSA_test will:
    
    * Download gmx_MMPBSA repository content in `/home/user/Documents`
    * Works with `Fast` set of examples [-t 2 is the default]
    * Perform the calculation on 9 examples sequentially, using 10 cpus each time

=== "All"
    
        gmx_MMPBSA_test -f /home/user/Documents -t 0 -n 10
    
    Through this command-line, gmxMMPBSA_test will:
    
    * Download gmx_MMPBSA repository content in `/home/user/Documents`
    * Works with `All` set of examples
    * Perform the calculation on 16 examples sequentially, using 10 cpus each time
    
=== "Specific"
    
        gmx_MMPBSA_test -f /home/user/Documents -t 3
    
    Through this command-line, gmxMMPBSA_test will:
    
    * Download gmx_MMPBSA repository content in `/home/user/Documents`
    * Execute `3` [Protein-Ligand (Single Trajectory approach)] example