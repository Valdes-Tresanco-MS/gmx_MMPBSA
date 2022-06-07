---
template: main.html
title:
---

## `gmx_MMPBSA` command-line

<div class="termy">
    ```console
    // All flags available in `gmx_MMPBSA` are shown below:
    
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
                             (default: $AMBERHOME/AmberTools/test/rism1d/tip3p-kh/tip3p.xvv.save)
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
                             *.pdb (default: None)
      -ri <Index File>      Index file of the unbound receptor. (default: None)
      -rg index             Receptor group in receptor index file. Notation:
                             "-lg <Receptor group>", e.g. -rg 1 (default: None)
      -rt [TRJ [TRJ ...]]   Input trajectories of the unbound receptor for multiple
                             trajectory approach. Allowed formats: *.xtc (recommended),
                             *.trr, *.pdb (specify as many as you'd like).
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
                             (recommended), *.pdb (default: None)
      -li <Index File>      Index file of the unbound ligand. Only if tpr file was
                             define in -ls. (default: None)
      -lg index             Ligand group in ligand index file. Notation:
                             "-lg <Ligand group>", e.g. -lg 13 (default: None)
      -lt [TRJ [TRJ ...]]   Input trajectories of the unbound ligand for multiple
                             trajectory approach. Allowed formats: *.xtc
                            (recommended), *.trr, *.pdb (specify as many as
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
</div>