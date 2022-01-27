---
template: main.html
title:
---

## Running gmx_MMPBSA

!!! tip
    * Since version 1.4.0 we have fixed the `gmx_MMPBSA` inconsistencies when using `MPI`.
    * **We currently recommend the use of MPI since the computation time decreases considerably.**

=== "Parallel (MPI) version"
    `gmx_MMPBSA` as `MMPBSA.py` uses the `MPI` only to perform the calculations, the rest of 
    the process (_i.e_, Generation/conversion of Amber topologies, mutation, division of the trajectories, etc) 
    occurs in a single thread (See **Figure 1[^1]** for better reference). This means that it is not necessary to 
    install any 
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
        Unfortunately, when running `gmx_MMPBSA` with `MPI`, GROMACS's `gmx_mpi` can't be used. This is probably 
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
![drawing](assets/images/mmpbsa_py_mpi.png){ width=70% style="display: block; margin: 0 auto"}
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 1**. **`MPI` benchmark description from <a href="https://pubs.acs.org/doi/10.1021/ct300418h">MMPBSA.py paper</a>.**
  `MMPBSA.py` scaling comparison for `MM-PBSA` and `MM-GBSA` calculations on 200 frames of a 5910-atom complex. Times 
  shown are the times required for the calculation to finish. Note that `MM-GBSA` calculations are ∼5 times faster 
  than `MM-PBSA` calculations. All calculations were performed on NICS Keeneland (2 Intel Westmere 6-core CPUs per 
  node, QDR infiniband interconnect) 
  </figcaption>
</figure>

[2]: assets/images/mmpbsa_py_mpi.png