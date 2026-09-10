---
template: main.html
title:
---

## `gmx_MMPBSA` command-line

<div class="termy">
    ```console
    // All flags available in `gmx_MMPBSA` are shown below:
    
    $ gmx_MMPBSA -h
    
    usage: gmx_MMPBSA [-h] [-v] [--input-file-help]
                      [--create_input [{gb,pb,pb_mem,rism,ala,decomp,nmode,gbnsr6,all} ...]]
                      [-O] [-prefix <file prefix>] [-sys_name <system name>]
                      [--progress-style {auto,rich,classic,plain,none}] [-i FILE] [-xvvfile XVVFILE] [-o FILE]
                      [-do FILE] [-eo FILE] [-deo FILE] [-nogui] [-s] [--no-error-bundle] [-cs <Structure File>]
                      [-ci <Index File>] [-cg group group]
                      [-ct [TRJ ...]] [-cp <Topology>] [-cr <PDB File>] [-rs <Structure File>] [-ri <Index File>]
                      [-rg group] [-rt [TRJ ...]] [-rp <Topology>] [-lm <Structure File>] [-ls <Structure File>]
                      [-li <Index File>] [-lg group] [-lt [TRJ ...]] [-lp <Topology>] [--rewrite-output] [--clean]
    
    gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state 
    free energy calculations with GROMACS files. This program is an adaptation of 
    Amber's MMPBSA.py and essentially works as such. gmx_MMPBSA works with any GROMACS version.
    This program will calculate binding free energies using end-state free energy methods 
    on an ensemble of snapshots using a variety of implicit solvent models. This is the core 
    of gmx_MMPBSA and it will do all the calculations
    
    options:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      --input-file-help     Print all available options in the input file. (default: False)
      --create_input [{gb,pb,pb_mem,rism,ala,decomp,nmode,gbnsr6,all} ...]
                            Create an new input file with selected calculation
                            type. (default: None)
    
    Miscellaneous Options:
      -O, --overwrite       Allow output files to be overwritten (default: False)
      -prefix <file prefix> Prefix for intermediate files. (default: _GMXMMPBSA_)
      -sys_name <system name>, --sys_name <system name>
                            System name. Overrides sys_name in the input file.
                            (default: None)
      --progress-style {auto,rich,classic,plain,none}
                            Calculation progress display. Auto uses Rich in a
                            terminal and classic otherwise. (default: auto)
    
    Input and Output Files:
      These options specify the input files and optional output files.
    
      -i FILE               MM/PBSA input file. (default: None)
      -xvvfile XVVFILE      XVV file for 3D-RISM. (default: installed bundled
                             GMXMMPBSA/data/xvv_files/tip3p.xvv, unless the
                             legacy AMBERHOME XVV file is available)
      -o FILE               Output file with MM/PBSA statistics.
                             (default: FINAL_RESULTS_MMPBSA.dat)
      -do FILE              Output file for decomposition statistics summary.
                             (default: FINAL_DECOMP_MMPBSA.dat)
      -eo FILE              CSV-format output of all energy terms for every frame in
                             every calculation. Defaults to the .csv counterpart of
                             the file specified with -o; use -eo to override its name.
                             (default: derived from -o)
      -deo FILE             CSV-format output of all energy terms for each printed
                             residue in decomposition calculations. Defaults to the
                             .csv counterpart of the file specified with -do; use -deo
                             to override its name. (default: derived from -do)
      -nogui                No open gmx_MMPBSA_ana after all calculations finished
                             (default: True)
      -s, --stability       Perform stability calculation. Only the complex parameters
                             are required. Only If the ligand is non-Protein (small
                             molecule) type and you not define a complex topology,
                             then ligand *.mol2 file is required. In any other case
                             receptor and ligand parameters will be ignored. See
                             description bellow (default: False)
      --no-error-bundle     Do not create a diagnostic zip bundle automatically when
                             gmx_MMPBSA fails. (default: False)
    
    Complex:
      Complex files and info that are needed to perform the calculation. If the
      receptor and/or the ligand info is not defined, we generate them from that of
      the complex.
    
      -cs <Structure File>  Structure file of the complex. If it is Protein-Ligand
                             (small molecule) complex and -cp is not defined, make
                             sure that you define -lm option. See -lm description
                             below. Allowed formats: *.tpr (recommended), *.pdb
                             (default: None)
      -ci <Index File>      Index file of the bound complex. (default: None)
      -cg group group       Receptor and ligand groups in the complex index file,
                            specified by zero-based group number or group name. For
                            example: -cg 1 13 or -cg Protein LIG
                             (default: None)
      -ct [TRJ ...]         Complex trajectories. Make sure the trajectory is
                            fitted and pbc have been removed. Allowed formats:
                            *.xtc (recommended), *.trr, *.pdb (specify as many as
                            you'd like). (default: None)
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
      -rg group             Receptor group in the receptor index file, specified by
                             zero-based group number or group name. For example: -rg 1
                             or -rg Protein (default: None)
      -rt [TRJ ...]         Input trajectories of the unbound receptor for
                            multiple trajectory approach. Allowed formats: *.xtc
                            (recommended), *.trr, *.pdb (specify as many as
                            you'd like). (default: None)
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
      -lg group             Ligand group in the ligand index file, specified by zero-
                             based group number or group name. For example: -lg 13 or
                             -lg LIG (default: None)
      -lt [TRJ ...]         Input trajectories of the unbound ligand for multiple
                            trajectory approach. Allowed formats: *.xtc
                            (recommended), *.trr, *.pdb (specify as many as
                            you'd like). (default: None)
      -lp <Topology>        Topology file of the ligand. (default: None)
    
    Miscellaneous Actions:
      -rewrite-output       Do not re-run any calculations, just parse the output
                             files from the previous calculation and rewrite the
                             output files. (default: False)
      --clean               Clean temporary files and quit. (default: False)
    
    gmx_MMPBSA is an effort to implement the GB/PB and others calculations in
    GROMACS. Based on MMPBSA.py (version 14.0) and AmberTools 26.0 and GROMACS
    2026.0
    ```
</div>

The help block above follows the current `gmx_MMPBSA -h` output. The bundled XVV
path is written relative to the installation here because the executable expands
it to an installation-specific absolute path at runtime. If an AMBERHOME XVV file
is available, it takes precedence; `-xvvfile` always overrides the default.

### Automatic CSV filenames

Per-frame CSV output is generated automatically. The default is the summary filename with its suffix
replaced by `.csv`. If that would name the summary itself, `.frames.csv` is used instead: `-o results.csv`
produces the text summary `results.csv` and the energy vectors `results.frames.csv`. The same rule applies
to decomposition output (`-do`/`-deo`). Explicit `-eo` and `-deo` values are preserved. Active output paths
must refer to distinct files; collisions are rejected before opening the output files.
