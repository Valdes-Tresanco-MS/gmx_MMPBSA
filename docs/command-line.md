---
template: main.html
title: Command-line
---


# `gmx_MMPBSA` applications command-line

## `gmx_MMPBSA` command-line
`gmx_MMPBSA` is invoked through the command line as follows:
```
$ gmx_MMPBSA -h

usage: gmx_MMPBSA [-h] [-v] [--input-file-help] [--create_input [{gb,pb,rism,ala,decomp,nmode,all}] 
                  [-O] [-prefix <file prefix>] [-i FILE] [-xvvfile XVVFILE] [-o FILE] [-do FILE] [-eo FILE]
                  [-deo FILE] [-nogui] [-s] [-cs <Structure File>] [-ci <Index File>] [-cg index index] 
                  [-ct [TRJ [TRJ ...]]] [-cp <Topology>] [-cr <PDB File>] [-rs <Structure File>] [-ri <Index File>] 
                  [-rg index] [-rt [TRJ [TRJ ...]]] [-rp <Topology>] [-lm <Structure File>] [-ls <Structure File>] 
                  [-li <Index File>] [-lg index] [-lt [TRJ [TRJ ...]]] [-lp <Topology>] [--rewrite-output] [--clean]

gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state 
free energy calculations with GROMACS files. This program is an adaptation of 
Amber's MMPBSA.py and essentially works as such. gmx_MMPBSA works with any GROMACS version.
This program will calculate binding free energies using end-state free energy methods 
on an ensemble of snapshots using a variety of implicit solvent models. This is the core 
of gmx_MMPBSA and it will do all the calculations

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --input-file-help     Print all available options in the input file. (default: False)
  --create_input        Create an new input file with selected calculation type. (default: None)
                         [{gb,pb,rism,ala,decomp,nmode,all}]

Miscellaneous Options:
  -O, --overwrite       Allow output files to be overwritten (default: False)
  -prefix <file prefix> Prefix for intermediate files. (default: _GMXMMPBSA_)

Input and Output Files:
  These options specify the input files and optional output files.

  -i FILE               MM/PBSA input file. (default: None)
  -xvvfile XVVFILE      XVV file for 3D-RISM.
                         (default: /home/mario/programs/amber20/dat/mmpbsa/spc.xvv)
  -o FILE               Output file with MM/PBSA statistics.
                         (default: FINAL_RESULTS_MMPBSA.dat)
  -do FILE              Output file for decomposition statistics summary.
                         (default: FINAL_DECOMP_MMPBSA.dat)
  -eo FILE              CSV-format output of all energy terms for every frame in
                         every calculation. File name forced to end in [.csv].
                         This file is only written when specified on the
                         command-line. (default: None)
  -deo FILE             CSV-format output of all energy terms for each printed
                         residue in decomposition calculations. File name forced
                         to end in [.csv]. This file is only written when
                         specified on the command-line. (default: None)
  -nogui                No open gmx_MMPBSA_ana after all calculations finished
                         (default: True)
  -s, --stability       Perform stability calculation. Only the complex parameters
                         are required. Only If the ligand is non-Protein (small
                         molecule) type and you not define a complex topology,
                         then ligand *.mol2 file is required. In any other case
                         receptor and ligand parameters will be ignored. See
                         description bellow (default: False)

Complex:
  Complex files and info that are needed to perform the calculation. If the
  receptor and/or the ligand info is not defined, we generate them from that of
  the complex.

  -cs <Structure File>  Structure file of the complex. If it is Protein-Ligand
                         (small molecule) complex and -cp is not defined, make
                         sure that you define -lm option. See -lm description
                         below Allowed formats: *.tpr (recommended), *.pdb,
                         *.gro (default: None)
  -ci <Index File>      Index file of the bound complex. (default: None)
  -cg index index       Groups of receptor and ligand in complex index file. The
                        notation is as follows:
                        "-cg <Receptor group> <Ligand group>", ie. -cg 1 13
                         (default: None)
  -ct [TRJ [TRJ ...]]   Complex trajectories. Make sure the trajectory is fitted
                         and pbc have been removed. Allowed formats: *.xtc 
                         (recommended), *.trr, *.pdb (specify as many as you'd 
                         like). (default: None)
  -cp <Topology>        The complex Topology file. When it is defined -lm
                         option is not needed (default: None)
  -cr <PDB File>        Complex Reference Structure file. This option is optional
                         but recommended (Use the PDB file used to generate the 
                         topology in GROMACS). If not defined, the chains ID 
                         assignment (if the structure used in -cs does not have 
                         chain IDs) will be done automatically according to the 
                         structure (can generate wrong mapping). (default: None)

Receptor:
  Receptor files and info that are needed to perform the calculation. If the
  receptor info is not defined, we generate it from that of the complex.

  -rs <Structure File>  Structure file of the unbound receptor for multiple
                         trajectory approach. Allowed formats: *.tpr (recommended),
                         *.pdb, *.gro (default: None)
  -ri <Index File>      Index file of the unbound receptor. (default: None)
  -rg index             Receptor group in receptor index file. Notation:
                         "-lg <Receptor group>", e.g. -rg 1 (default: None)
  -rt [TRJ [TRJ ...]]   Input trajectories of the unbound receptor for multiple
                         trajectory approach. Allowed formats: *.xtc (recommended),
                         *.trr, *.pdb, *.gro (specify as many as you'd like).
                         (default: None)
  -rp <Topology>        Topology file of the receptor. (default: None)

Ligand:
  Ligand files and info that are needed to perform the calculation. If the ligand
  are not defined, we generate it from that of the complex.

  -lm <Structure File>  A *.mol2 file of the unbound ligand used to parametrize
                         ligand for GROMACS using Antechamber. Must be defined
                         if Protein-Ligand (small molecule) complex was define 
                         and -cp or -lp option are not defined. No needed for 
                         Proteins, DNA, RNA, Ions, Glycans or any ligand 
                         parametrized in the Amber force fields. Must be the 
                         Antechamber output *.mol2. (default: None)
  -ls <Structure File>  Structure file of the unbound ligand. If ligand is a 
                         small molecule and -lp is not defined, make sure that you
                         define above -lm option. Allowed formats: *.tpr 
                         (recommended), *.pdb, *.gro (default: None)
  -li <Index File>      Index file of the unbound ligand. Only if tpr file was
                         define in -ls. (default: None)
  -lg index             Ligand group in ligand index file. Notation:
                         "-lg <Ligand group>", e.g. -lg 13 (default: None)
  -lt [TRJ [TRJ ...]]   Input trajectories of the unbound ligand for multiple
                         trajectory approach. Allowed formats: *.xtc
                        (recommended), *.trr, *.pdb, *.gro (specify as many as
                         you'd like). (default: None)
  -lp <Topology>        Topology file of the ligand. (default: None)

Miscellaneous Actions:
  -rewrite-output       Do not re-run any calculations, just parse the output
                         files from the previous calculation and rewrite the
                         output files. (default: False)
  --clean               Clean temporary files and quit. (default: False)

gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
Based on MMPBSA.py (version 16.0) and AmberTools20
```

### Running gmx_MMPBSA

!!! tip
    * Since version 1.4.0 we have fixed the `gmx_MMPBSA` inconsistencies when using `MPI`.
    * **We currently recommend the use of MPI since the computation time decreases considerably.**

=== "Parallel (MPI) version"
    `gmx_MMPBSA` as `MMPBSA.py` uses the `MPI` only to perform the calculations during the dynamics, the rest of 
    the process (Generation/conversion of Amber topologies, mutation, division of the trajectories, etc) occurs in 
    a single thread (See **Figure 1[^1]** for better reference). This means that it is not necessary to install any 
    program (AmberTools or GROMACS) with `MPI`, which can be used in any circumstance, and the time required to 
    process the data prior to the calculation depends on the system and will be the same for both versions (`Serial` 
    and `MPI`). 

    !!! note
        Note that `gmx_MMPBSA` processes, converts, or builds topologies from GROMACS files, so it takes slightly 
        longer than `MMPBSA.py` at the same stage of the process. However, this is not really significant.


    ???+ tip "Remember" 
        Make sure that you install the OpenMPI library
        
            sudo apt install openmpi-bin libopenmpi-dev openssh-client

    A usage example is shown below:

    === "Local"
    
            mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc
    
    === "HPC"
    
            #!/bin/sh
            #PBS -N nmode
            #PBS -o nmode.out
            #PBS -e nmode.err
            #PBS -m abe
            #PBS -M email@domain.edu
            #PBS -q brute
            #PBS -l nodes=1:surg:ppn=3
            #PBS -l pmem=1450mb or > 5gb for nmode calculation
            
            cd $PBS_O_WORKDIR
            
            mpirun -np 3 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc > progress.log

    
    !!! danger
        Unfortunately, when `gmx_MMPBSA` is run with `MPI`, GROMACS's `gmx_mpi` is incompatible. This is probably 
        because of `gmx_mpi` conflicts with mpirun. In any case, this is not a problem since `gmx` works correctly 
        and `gmx_mpi` only parallels `mdrun`, the rest of the GROMACS tools work in a single thread. See this 
        [issue](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/26) to see the output.

    !!! warning
        The nmode calculations require a considerable amount of RAM. Consider that the total amount of RAM will be:

        RAM~total~ = RAM~1_frame~ * NUM of Threads
        
        If it consumes all the RAM of the system it can cause crashes, instability or system shutdown!
        

    !!! note
        At a certain level, running RISM in parallel may actually hurt performance, since previous solutions are used 
        as an initial guess for the next frame, hastening convergence. Running in parallel loses this advantage. Also, 
        due to the overhead involved in which each thread is required to load every topology file when calculating 
        energies, parallel scaling will begin to fall off as the number of threads reaches the number of frames. 

=== "Serial version"
    This version is installed via pip as described above. `AMBERHOME` variable must be set, or it will quit with an error. 
    An example command-line call is shown below:
    
        gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc
    
    You can found test files on [GitHub][1]

  [1]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples

<figure markdown="1">
[![overview][2]][2]
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 1**. **`MPI` benchmark description from <a href="https://pubs.acs.org/doi/10.1021/ct300418h">`MMPBSA.py` paper</a>.**
  `MMPBSA.py` scaling comparison for `MM-PBSA` and `MM-GBSA` calculations on 200 frames of a 5910-atom complex. Times 
  shown are the times required for the calculation to finish. Note that `MM-GBSA` calculations are ∼5 times faster 
  than `MM-PBSA` calculations. All calculations were performed on NICS Keeneland (2 Intel Westmere 6-core CPUs per 
  node, QDR infiniband interconnect) 
  </figcaption>
</figure>

[2]: assets/images/mmpbsa_py_mpi.png

## `gmx_MMPBSA_ana` command-line
```
$ gmx_MMPBSA_ana -h

usage: run_ana.py [-h] [-v] [-f [FILES [FILES ...]]] [-r]

This program is part of gmx_MMPBSA and will show a workspace to analyze the 
gmx_MMPBSA results

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Info file:
  -f [FILES [FILES ...]], --files [FILES [FILES ...]]
                        gmx_MMPBSA info files or container folder or list of 
                        them (default: [.] Current working dir)
  -r, --recursive       Search recursively in this folder at depth = 1 
                        (default: False)

gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
Based on MMPBSA.py (version 16.0) and AmberTools20
```

### Running gmx_MMPBSA_ana
In order to analyze multiple systems in the same section and implement the correlation between them, we improved the 
file input to gmx_MMPBSA_ana. Currently, gmx_MMPBSA_ana supports info files, the folder that contains it or a list 
of them.

Additionally, you can pass as input one or more folders containing multiple systems using `recursive`. The folders 
must be the following structure:

=== "Valid folder structure"
    
    Here all folders (5 systems) will be processed by `gmx_MMPBSA_ana`
    ```
    Defined folder
      ├── System-1
      │      └── _GMXMMPBSA_info
      ├── System-2
      │      └── _GMXMMPBSA_info
      ├── System-3
      │      └── _GMXMMPBSA_info
      ├── System-4
      │      └── _GMXMMPBSA_info
      └── System-5
             └── _GMXMMPBSA_info
    ```

=== "Valid folders + files structure"

    Here all systems (4) will be processed by `gmx_MMPBSA_ana`
    ```
    System-1      
      └──_GMXMMPBSA_info
    Folder-1
      ├── System-2.1
      │      └── _GMXMMPBSA_info
      └── System-2.2
             └── _GMXMMPBSA_info
    System-3 (_GMXMMPBSA_info)
    ```
    
    !!! tip ""
        The command-line for this is:

            gmx_MMPBSA_ana -f /path/to/System-1 /path/to/System-3/_GMXMMPBSA_info /path/to/Folder-1 -r
        
        This estructure only work if recursive option was defined. See the examples below
    
    !!! note ""
        Note that:

        * `System-1` is defined as folder that contain a _GMXMMPBSA_info
        * Folder-1 contain two folder (Systems), each containing a _GMXMMPBSA_info
        * System-3 is defined as a _GMXMMPBSA_info file


=== "Wrong folder structure"
    
    Here only  4 (systems) folders will be processed by gmx_MMPBSA_ana. The systems in the `Internal folder` will be 
    ignored 
    ```
    Defined folder
      ├── System-1
      │      └── _GMXMMPBSA_info
      ├── System-2
      │      └── _GMXMMPBSA_info
      ├── Internal folder
      │      ├─X─ System-3.1
      │      │      └── _GMXMMPBSA_info
      │      └─X─ System-3.2
      │             └── _GMXMMPBSA_info
      ├── System-4
      │      └── _GMXMMPBSA_info
      └── System-5
             └── _GMXMMPBSA_info
    ```

!!! examples
        
    === "One file"
        
        Passing a _GMXMMPBSA_info file as input:
        
        * Current directory
        
                gmx_MMPBSA_ana -f _GMXMMPBSA_info
        
        * other location  :material-new-box:{: .medium .heart } Version: 1.4.0

                gmx_MMPBSA_ana -f /path/to/_GMXMMPBSA_info

    === "One folder" 
        :material-new-box:{: .medium .heart } Version: 1.4.0

        Passing a folder as input:  
        
        * Current directory
        
                gmx_MMPBSA_ana -f .
        
        * other location

                gmx_MMPBSA_ana -f /path/to/folder

        !!! warning "Remember"
            This folder must contain a valid `_GMXMMPBSA_info` file
    
    === "Files + Folders"
        :material-new-box:{: .medium .heart } Version: 1.4.0
        
            gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/_GMXMMPBSA_info-1 /path/to/folder-2
        
        !!! warning "Remember"
            * All defined folders must contain a valid `_GMXMMPBSA_info` file
            * All `_GMXMMPBSA_info` files defined must be valid
    
    === "Recursive option"
        :material-new-box:{: .medium .heart } Version: 1.4.0

        Passing a folder as input with `recursive` option:  
        
        * Current directory
        
                gmx_MMPBSA_ana -f . -r
        
        * other location

                gmx_MMPBSA_ana -f /path/to/folder --recursive

        * combine multiple folders

                gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/folder-2 /path/to/folder-3  -r

            !!! warning ""
            Folders can contain one or more systems

        * combine multiple folders and files

                gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/_GMXMMPBSA_info-1 /path/to/folder-3  -r

            !!! warning ""
            * Folders can contain one or more systems
            * Note that if you remove the option -r, each folder must contain a valid _GMXMMPBSA_info file.


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

### Running gmx_MMPBSA_test
gmx_MMPBSA_test is designed to run a set of samples (all or minimal) or a specific example efficiently. 
Additionally, gmx_MMPBSA_test can run in parallel, decreasing the execution time gmx_MMPBSA_test will download the 
most recent version of the repository in the specified folder and will perform the calculations

=== "Minimal"
    
        gmx_MMPBSA_test -f /home/user/Documents -n 10 [-t 2 is the default]
    
    Through this command-line, gmxMMPBSA_test will:
    
    * Download gmx_MMPBSA repository content in `/home/user/Documents`
    * Works with `Fast` set of examples
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

[^1]: Taken from the MMPBSA.py paper