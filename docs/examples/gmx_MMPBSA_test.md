---
template: main.html
title: gmx_MMPBSA_test
---

## `gmx_MMPBSA_test` command-line
<div class="termy">
    ```bash
    $ gmx_MMPBSA_test -h
    usage: gmx_MMPBSA_test [-h] [-v] 
           [-t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,101} [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,101} ...]]]
           [-f FOLDER] [-r] [-ng] [-n NUM_PROCESSORS] [-j NUM_CONCURRENT]
    
    This program is part of gmx_MMPBSA and will allow you to run various gmx_MMPBSA examples easily.
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    Test options:
      -t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,101} [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,101} ...]]
                            The level the test is going to be run at. Multiple systems and analysis can be run at the same 
                            time.
                                  Nr. of Sys  
                            * 0      23     All -- Run all examples (Can take a long time!!!)
                            * 1      12     Minimal -- Does a minimal test with a set of systems and analyzes
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
                            * 19          C2 Entropy approximation
                            * 20          LPB Calculation
                            * 21          NLPB Calculation
                            * 22          Protein-Ligand_LPH (CHARMM force field)
                            * 23          QM/MMGBSA Calculation
                            * 24          GBNSR6 Calculation
                            * 25     |  5 AMBER input files
      -f FOLDER, --folder FOLDER
                            Defines the folder to store all data
      -r, --reuse           Defines the existing test forlder will be reuse
      -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
      -n NUM_PROCESSORS, --num_processors NUM_PROCESSORS
                            Defines the number of processor cores you want to use with MPI per calculation. If the number 
                            of frames is less than the number of cpus defined, the calculation will be performed with 
                            the number of processors = number of frames.
      -j NUM_CONCURRENT, --num_concurrent NUM_CONCURRENT
                            Defines the number of examples to run concurrently. Each example can use up to
                            --num_processors MPI ranks, so the total rank count can be -j * -n.
    
    
    gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
    Based on MMPBSA.py (version 16.0) and AmberTools20
    ```
</div>

## Running gmx_MMPBSA_test
gmx_MMPBSA_test is designed to run a set of samples (all or minimal) or a specific example efficiently.
By default, examples run sequentially. To run multiple examples at the same time, use `-j/--num_concurrent`.
Each example can use up to `-n/--num_processors` MPI ranks.

!!! info "Sets in gmx_MMPBSA_test"

    === "Fast"
        
            gmx_MMPBSA_test -f /home/user/Documents -n 10
        
        Through this command-line, gmxMMPBSA_test will:
        
        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Works with `Fast` set of examples [-t 2 is the default]
        * Perform the calculation on 9 examples sequentially, using 10 cpus each time

    === "Parallel examples"

            gmx_MMPBSA_test -f /home/user/Documents -t 3 5 7 -n 4 -j 2

        Through this command-line, gmxMMPBSA_test will:

        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Execute examples `3`, `5`, and `7`
        * Run up to 2 examples at the same time
        * Use up to 4 MPI ranks per example
    
    === "Minimal"
        
            gmx_MMPBSA_test -f /home/user/Documents -n 10 -t 1
        
        Through this command-line, gmxMMPBSA_test will:
        
        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Works with `Minimal` set of examples [-t 1]
        * Perform the calculation on 12 examples sequentially, using 10 cpus each time
    
    === "All"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 0 -n 10
        
        Through this command-line, gmxMMPBSA_test will:
        
        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Works with `All` set of examples
        * Perform the calculation on 23 examples sequentially, using 10 cpus each time
        
    === "Multiple selection"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 3 5 7
        
        Through this command-line, gmxMMPBSA_test will:
        
        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Execute `3` [Protein-Ligand (Single Trajectory approach)], `5` [Protein-DNA], and `7` [Protein-Glycan]
        examples

    === "Single selection"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 3
        
        Through this command-line, gmxMMPBSA_test will:
        
        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Execute `3` [Protein-Ligand (Single Trajectory approach)] example

    === "AMBER input files"

            gmx_MMPBSA_test -f /home/user/Documents -t 25

        Through this command-line, gmxMMPBSA_test will:

        * Download gmx_MMPBSA repository content in `/home/user/Documents`
        * Execute `25` [AMBER input files] example

!!! warning "3D-RISM AmberTools runtime failures"
    Test `18` uses AmberTools 3D-RISM. If this test fails while the other examples pass and the log contains
    `Fortran runtime error: Missing comma between descriptors` from `amber_rism_interface.F90`, the failure is a
    known AmberTools/Fortran runtime compatibility issue. It has been reproduced with conda AmberTools builds linked
    against newer `libgfortran` runtimes. A known working workaround is `gmx_MMPBSA` 1.6.4 with Python 3.9 or 3.10,
    AmberTools 23, and `libgfortran5`/`libgcc-ng` 12.x, or a patched AmberTools build. See the
    [3D-RISM example](3D-RISM/README.md) for details.
