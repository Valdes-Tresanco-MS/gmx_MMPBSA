---
template: main.html
title: Command-line
---


# gmx_MMPBSA command-line

## Calling gmx_MMPBSA from the command-line
`gmx_MMPBSA` is invoked through the command line as follows:
```
$ gmx_MMPBSA -h

usage: gmx_MMPBSA [-h] [-v] [--input-file-help] [-O] [-prefix <file prefix>]
                  [-i FILE] [-xvvfile XVVFILE] [-o FILE] [-do FILE] [-eo FILE]
                  [-deo FILE] [-gui] [-s] [-cs <Structure File>]
                  [-ci <Index File>] [-cg index index] [-ct [TRJ [TRJ ...]]]
                  [-cp <Topology>] [-cr <PDB File>] [-rs <Structure File>]
                  [-ri <Index File>] [-rg index] [-rt [TRJ [TRJ ...]]]
                  [-rp <Topology>] [-lm <Structure File>] [-ls <Structure File>]
                  [-li <Index File>] [-lg index] [-lt [TRJ [TRJ ...]]]
                  [-lp <Topology>] [-make-mdins] [-use-mdins] [-rewrite-output]
                  [--clean]
gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS.
This program is an adaptation of Amber's MMPBSA.py and essentially works as such. 
As gmx_MMPBSA adapts MMPBSA.py, since it has all the resources of this script and
work with any GROMACS version.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --input-file-help     Print all available options in the input file. 
                        (default: False)

Miscellaneous Options:
  -O, --overwrite       Allow output files to be overwritten (default: False)
  -prefix <file prefix>
                        Prefix for intermediate files. (default: _GMXMMPBSA_)

Input and Output Files:
  These options specify the input files and optional output files.

  -i FILE               MM/PBSA input file. (default: None)
  -xvvfile XVVFILE      XVV file for 3D-RISM. 
                        (default: /path/installed/amber20/dat/mmpbsa/spc.xvv)
  -o FILE               Output file with MM/PBSA statistics. 
                        (default: FINAL_RESULTS_MMPBSA.dat)
  -do FILE              Output file for decomposition statistics summary. 
                        (default: FINAL_DECOMP_MMPBSA.dat)
  -eo FILE              CSV-format output of all energy terms for every frame 
                        in every calculation. File name forced to end in [.csv].
                        This file is only written when specified on the 
                        command-line. (default: None)
  -deo FILE             CSV-format output of all energy terms for each printed 
                        residue in decomposition calculations. File name forced 
                        to end in [.csv]. This file is only written when 
                        specified on the command-line. (default: None)
  -nogui                Open charts application when all calculations finished 
                        (default: True)
  -stability            Perform stability calculation. Only the complex 
                        parameters are required. If ligand is non-Protein (small
                        molecule) type, then ligand *.mol2 file is required. In
                        any other case receptor and ligand parameters will be 
                        ignored. See description bellow (default: False)

Complex:
  Complex files and info that are needed to perform the calculation. If the 
  receptor and/or the ligand info is not defined, we generate them from that of
  the complex.

  -cs <Structure File>  Structure file of the complex. If it is Protein-Ligand 
                        (small molecule) complex, make sure that you define -lm
                        option. See -lm description below Allowed formats: *.tpr
                        (recommended), *.pdb, *.gro (default: None)
  -ci <Index File>      Index file of the bound complex. (default: None)
  -cg index index       Groups of receptor and ligand in complex index file. The 
                        notation is as follows: "-cg <Receptor group> <Ligand 
                        group>", ie. -cg 1 13 (default: None)
  -ct [TRJ [TRJ ...]]   Input trajectories of the complex. Make sure the 
                        trajectory is fitted and pbc have been removed. Allowed 
                        formats: *.xtc (recommended), *.trr, *.pdb (specify as 
                        many as you'd like). (default: None)
  -cp <Topology>        Topology file of the complex. (default: None)
  -cr <PDB File>        Complex Reference Structure file. This option is optional
                        but recommended (Use the PDB file used to generate the 
                        topology in GROMACS). If not defined, the chains ID 
                        assignment (if the structure used in -cs does not have 
                        chain IDs) will be done automatically according to the 
                        structure (can generate inconsistencies). (default: None)

Receptor:
  Receptor files and info that are needed to perform the calculation. If the 
  receptor info is not defined, we generate it from that of the complex.

  -rs <Structure File>  Structure file of the unbound receptor for multiple 
                        trajectory approach. Allowed formats: *.tpr (recommended), 
                        *.pdb, *.gro (default: None)
  -ri <Index File>      Index file of the unbound receptor. (default: None)
  -rg index             Receptor group in receptor index file. Notation: 
                        "-rg <Receptor group>", e.g. -rg 1 (default: None)
  -rt [TRJ [TRJ ...]]   Input trajectories of the unbound receptor for multiple 
                        trajectory approach. Allowed formats: *.xtc (recommended), 
                        *.trr, *.pdb, *.gro (specify as many as you'd like). 
                        (default: None)
  -rp <Topology>        Topology file of the receptor. (default: None)

Ligand:
  Ligand files and info that are needed to perform the calculation. If the ligand 
  are not defined, we generate it from that of the complex.

  -lm <Structure File>  A *.mol2 file of the unbound ligand used to parametrize 
                        ligand for GROMACS using Anetchamber. Must be defined if 
                        Protein-Ligand (small molecule) complex was define. No 
                        needed for Proteins, DNA, RNA, Ions, and Glycans. 
                        Antechamber output *.mol2 is recommended. (default: None)
  -ls <Structure File>  Structure file of the unbound ligand. If ligand is a small
                        molecule, make sure that you define above -lm option. 
                        Allowed formats: *.tpr (recommended), *.pdb, *.gro 
                        (default: None)
  -li <Index File>      Index file of the unbound ligand. Only if tpr file was
                        define in -ls. (default: None)
  -lg index             Ligand group in ligand index file. Notation: 
                        "-lg <Ligand group>", e.g. -lg 13 (default: None)
  -lt [TRJ [TRJ ...]]   Input trajectories of the unbound ligand for multiple 
                        trajectory approach. Allowed formats: *.xtc (recommended), 
                        *.trr, *.pdb, *.gro (specify as many as you'd like). 
                        (default: None)
  -lp <Topology>        Topology file of the ligand. (default: None)

Miscellaneous Actions:
  -make-mdins           Create the input files for each calculation and quit. This 
                        allows you to modify them and re-run using -use-mdins 
                        (default: False)
  -use-mdins            Use existing input files for each calculation. If they do
                        not exist with the appropriate names, run_cmd.py will quit 
                        in error. (default: False)
  -rewrite-output       Do not re-run any calculations, just parse the output 
                        files from the previous calculation and rewrite the output
                        files. (default: False)
  --clean               Clean temporary files and quit. (default: False)

This program will calculate binding free energies using end-state free energy 
methods on an ensemble of snapshots using a variety of implicit solvent models. 
Based on MMPBSA.py (version 16.0) and AmberTools20
```

## Running gmx_MMPBSA
### Serial version
This version is installed via pip as described above. `AMBERHOME` variable must be set, or it will quit with an error. 
An example command-line call is shown below:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc

You can found test files on [GitHub][1]

### Parallel (MPI) version
Unlike MMPBSA.py, `gmx_MMPBSA` will be installed as a separate package from the Amber installation. When installing
Amber with mpi, a MMPBSA.py version called "MMPBSA.py.MPI" will be installed as well. Since we cannot detect if 
Amber was installed one way or another, we simply decided to adapt the `gmx_MMPBSA` executable to use an argument. 
That is, `gmx_MMPBSA` is a single script that executes the serial version or the parallel version with mpi depending 
on whether the user defines the "mpi" or "MPI" argument. In principle, both the serial and parallel versions should 
work correctly when Amber was installed in parallel.

The parallel version, like MMPBSA.py.MPI requires the mpi4py module. If you did the parallel installation of Amber, it 
should be installed. In any case, it could be installed in the following way:

    amber.python -m pip install mpi4py
    
A usage example is shown below:

    mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc

or

    mpirun -np 2 gmx_MMPBSA mpi -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc

!!! note
    At a certain level, running RISM in parallel may actually hurt performance, since previous solutions are used 
    as an initial guess for the next frame, hastening convergence. Running in parallel loses this advantage. Also, 
    due to the overhead involved in which each thread is required to load every topology file when calculating 
    energies, parallel scaling will begin to fall off as the number of threads reaches the number of frames. 


  [1]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/docs/examples