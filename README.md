[![pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![support](https://img.shields.io/badge/support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![python](https://img.shields.io/badge/python-v3.x-blue)]()


# Documentation
**Note: We do not intend to replace the original [MMPBSA.py](https://pubs.acs.org/doi/10.1021/ct300418h); instead, 
we have implemented and improved some functionalities, and what is most important, made this valuable tool available 
for GROMACS users. Most of the documentation below is found in the [Amber manual](https://ambermd.org/doc12/Amber20.pdf#chapter.34), 
we will point out what is new or different. Neither of these should be considered as a “black-box”, and users 
should be familiar with Amber and MM/PB(GB)SA method before at-tempting these sorts of calculations. These scripts 
automate a series of calculations, and cannot trap all the types of errors that might occur. You can review some of the 
answers to the questions that we consider most common here. If you find a bug or have any question, please consider 
opening an [issue](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues) or post in our [Google group](https://groups.google.com/g/gmx_mmpbsa)** 

# TOC
* [Requirements](#requirements)
* [Installing `gmx_MMPBSA`](#installing-gmx_mmpbsa)
    * [Update](#update)
    * [After Install](#after-install)
* [Introduction](#introduction)
    * [Literature](#literature)
* [`gmx_MMPBSA` in a nutshell](#gmx_mmpbsa-in-a-nutshell)
    * [Types of calculation you can do](#types-of-calculations-you-can-do)
    * [Comparison of `gmx_MMPBSA` vs other programs](#comparison-of-gmx_MMPBSA-vs-other-programs)
    * [Examples](#examples)
* [Calling `gmx_MMPBSA` from the command-line](#calling-gmx_mmpbsa-from-the-command-line)
* [Running `gmx_MMPBSA`](#running-gmx_mmpbsa)
    * [Serial version](#serial-version)
    * [Parallel (MPI) version](#parallel-mpi-version)
* [Input and Output](#input-and-output)
    * [The input file](#the-input-file)
        * [Sample input files](#sample-input-files)
    * [The Output File](#the-output-file)
    * [Temporary Files](#temporary-files)
* [Advanced Options](#advanced-options)

* [Python API](#python-api)
    * [Using the API](#using-the-api)
    * [Properties of mmpbsa_data](#properties-of-mmpbsa_data)
        * [Attributes](#attributes)
        * [Defined operators](#defined-operators)
    * [Example API Usage](#example-api-usage)
    * [Decomposition Data](#decomposition-data)
    

## Requirements
`gmx_MMPBSA` requires AmberTools20 to be installed in your machine with Python3 and the shell environment correctly set up for Amber. 
The AmberTools suite is free of charge and you can check [Amber Manual](https://ambermd.org/doc12/Amber20.pdf#section.2.1)
for a detailed installation guide. Of note, you can have more than one AmberTools installed in your machine. In case 
AmberTools20 is not the default Amber in your computer, just make sure to source AmberTools20 before 
installing/updating/running `gmx_MMPBSA`.
`gmx_MMPBSA` also requires GROMACS (series 4.x.x or 5.x.x or 20xx.x) to be installed in your computer and the shell 
environment correctly set up for GROMACS. `gmx_MMPBSA` has been tested with GROMACS 4.6.7, 5.1.2 and 2018.3, although it
should run smoothly with any GROMACS present in the PATH and that is compatible with the files you are using.

`gmx_MMPBSA` contains a module that allows for plotting the results (`gmx_MMPBSA_gui`). For this, it requires the
 installation of PyQt5.

    amber.python -m pip install PyQt5

## Installing `gmx_MMPBSA`
You can install `gmx_MMPBSA` from the `stable` version on PYPI:

    amber.python -m pip install gmx_MMPBSA
or the `development` version from GitHub:
    
    amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
Make sure that you have `git` installed. If not you can install it as follow:

    sudo apt install git

### Update
If you already have installed a previous `gmx_MMPBSA` version, you can update it as follows:

`stable` version (recommended):

    amber.python -m pip install gmx_MMPBSA --upgrade

`development` version from GitHub:

    amber.python -m pip intall git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA --upgrade 
Make sure that you have `git` installed.
    
**We will do our best to keep the PYPI package up to date.**

### After Install
Once the installation is completed, the following warning may appear:

    WARNING: The scripts gmx_MMPBSA and gmx_MMPBSA_gui are installed in 
    '/home/user/path_to_amber_install/amber20/miniconda/bin' which is not on PATH.
This warning is because `pip` installs the executables (`gmx_MMPBSA` and `gmx_MMPBSA_gui`) in \*/amber20/miniconda/bin.

You have two options to solve this:

Add this folder (*/amber20/miniconda/bin) to PATH:

    export PATH="/home/user/path_to_amber_install/amber20/miniconda/bin:$PATH"
*This option is more permanent and is recommended if you don't want to activate and deactivate the 
conda environment*

**Note:** Make sure to update **user** and **path_to_amber_install** in the PATH variable

or

initializing the environment of conda amber:

    amber.conda init bash

You can deactivate like this:
    
    conda deactivate

After using one of the above options, you should be able to run `gmx_MMPBSA` and `gmx_MMPBSA_gui` through the terminal


If when running `gmx_MMPBSA`, you get an error like this:

    ModuleNotFoundError: No module named 'parmed'

please see the following [issue](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/2) to see the solution

## Introduction
Molecular Mechanics / Poisson Boltzmann (or Generalized Born) Surface Area (MM/PB(GB)SA) is a post-processing
method in which representative snapshots from an ensemble of conformations are used to calculate the free energy
change between two states (typically a bound and free state of a receptor and ligand). Free energy differences are
calculated by combining the so-called gas phase energy contributions (MM term) that are independent of the chosen solvent
model as well as solvation free energy components (both polar and non-polar) calculated from an implicit solvent
model combination (PBSA or GBSA) for each species. Entropy contributions to the total free energy may be added as a further 
refinement.

The gas phase free energy contributions are calculated by sander or mmpbsa_py_energy within the AmberTools package 
according to the force field used in the MD simulation. The solvation free energy contributions may be further 
decomposed into a polar and non-polar contributions. The polar portion is calculated using the Poisson Boltzmann (PB) 
equation, the Generalized Born method, or the Reference Interaction Site Model (RISM). The PB equation is solved 
numerically by either the pbsa program included with AmberTools or by the Adaptive Poisson Boltzmann Solver (APBS) 
program (for more information, see http://www.poissonboltzmann.org/apbs). The non-polar contribution is approximated by 
the LCPO method implemented within sander or the molsurf method as implemented in cpptraj. The entropy calculations 
can be done in either a HCT Generalized Born solvation model or in the gas phase using a mmpbsa_py_nabnmode 
program written in the nab programming language, or via the quasi-harmonic approximation in ptraj as in the original 
[MMPBSA.py](https://pubs.acs.org/doi/10.1021/ct300418h). In this new module, entropy term (−TΔS) can be also estimated 
using the so-called [Interaction Entropy](https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) method, which is 
theoretically rigorous, computationally efficient, and numerically reliable for calculating entropic contribution to 
free energy in protein–ligand binding and other interaction processes.

Usually, the Single Trajectory (ST) approximation is employed when performing MM/PB(GB)SA calculations. This approximation
assumes that the configurational space explored by the systems are very similar between the bound and unbound states, 
so every snapshot for each species (_i.e._ complex, receptor, and ligand) is extracted from the same trajectory file. This
approximation improves binding free energies convergence and also reduces the computing time. However, it only should be 
applied when the molecules in the unbound state present a similar behavior to that of the bound state. On the other 
hand, in the so-called Multiple Trajectory (MT) approximation, the snapshots for each one of the species (_i.e._ complex, 
receptor, and ligand) are extracted from their own trajectory file. This approximation, theoretically more rigorous 
though, leads to higher standard deviation of the binding free energies.  

### Literature
Further information can be found in [Amber manual](https://ambermd.org/doc12/Amber20.pdf):
* [MMPBSA.py](https://ambermd.org/doc12/Amber20.pdf#chapter.34)
* [The Generalized Born/Surface Area Model](https://ambermd.org/doc12/Amber20.pdf#chapter.4)
* [PBSA](https://ambermd.org/doc12/Amber20.pdf#chapter.6)
* [Reference Interaction Site Model](https://ambermd.org/doc12/Amber20.pdf#chapter.7)
* [Generalized Born (GB) for QM/MM calculations](https://ambermd.org/doc12/Amber20.pdf#subsection.10.1.3)

and the foundational papers:
* [Srinivasan J. et al., 1998](https://pubs.acs.org/doi/abs/10.1021/ja981844+) 
* [Kollman P. A. et al., 2000](https://pubs.acs.org/doi/abs/10.1021/ar000033j) 
* [Gohlke H., Case D. A. 2004](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.10379) 

as well as some reviews:
* [Genheden S., Ryde U. 2015](https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1032936) 
* [Wang et. al., 2018](https://www.frontiersin.org/articles/10.3389/fmolb.2017.00087/full)  
* [Wang et. al., 2019](https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.9b00055) 

## `gmx_MMPBSA` in a nutshell
`gmx_MMPBSA` brings in all the [MMPBSA.py](https://pubs.acs.org/doi/10.1021/ct300418h) functionalities to GROMACS users. 
In addition, few other functionalities were implemented that eases a number of calculations (_e.g._ MM/PB(GB)SA 
with different internal dielectric constant, interaction entropy calculation). A GUI application is also incorporated 
that allows for visualizing the results and saving high-quality images.

### Types of calculations you can do
There are many different options for running `gmx_MMPBSA`. Among the types of calculations you can do are:
* **Normal binding free energies**, with either PB or GB implicit solvent models. Each can be done with either
1, 2, or 3 different trajectories, but the complex, receptor, and ligand topology files must all be defined. The
complex trajectory must always be provided. Whichever trajectories of the receptor and/or ligand that are NOT
specified will be extracted from the complex trajectory. This allows a 1-, 2-, or 3-trajectory analysis. All PB
calculations and GB models can be performed with just AmberTools via the mmpbsa_py_energy program installed with 
MMPBSA.py.
* **Stability** calculations with any calculation type.
* **Alanine scanning** with either PB or GB implicit solvent models. All trajectories will be mutated to match
the mutated topology files, and whichever calculations that would be carried out for the normal systems are
also carried out for the mutated systems. Note that only 1 mutation is allowed per simulation, and it must
be to an alanine. If mutant_only is not set to 1, differences resulting from the mutations are calculated. This
option is incompatible with intermediate NetCDF trajectories (see the netcdf = 1 option above). This has the
same program requirements as option 1 above.
* **Entropy corrections**. An entropy term can be added to the free energies calculated above using either the
quasi-harmonic, the normal mode or interaction entropy approximations. Calculations will be performed for the normal 
and mutated systems (alanine scanning) as requested. Normal mode calculations are done with the
mmpbsa_py_nabnmode program included with AmberTools.
* **Decomposition schemes**. The energy terms will be decomposed according to the decomposition scheme
outlined in the idecomp variable description. This should work with all of the above, though entropy terms
cannot be decomposed. APBS energies cannot be decomposed, either. Neither can PBSA surface area terms.
This functionality requires sander from the Amber 11 (or later) package.
* **QM/MMGBSA**. This is a binding free energy (or stability calculation) using the Generalized Born solvent
model allowing you to treat part of your system with a quantum mechanical Hamiltonian. See [“Advanced
Options”](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#Advanced-Options) for tips about optimizing this option. 
This functionality requires sander from the Amber package.
* **MM/3D-RISM**. This is a binding free energy (or stability calculation) using the 3D-RISM solvation model.
This functionality is performed with rism3d.snglpnt built with AmberTools.
* **Membrane Protein MMPBSA**. Calculate the MMPBSA binding free energy for a ligand bound to a protein
that is embedded into a membrane. Only use_sander=1 is supported.

### Comparison of `gmx_MMPBSA` vs other programs
This comparison is based on the documentation of the different programs

| Feature | [g_mmpbsa](https://github.com/RashmiKumari/g_mmpbsa) | [MMPBSA.py](https://ambermd.org/doc12/Amber20.pdf#chapter.34) <sup>1</sup> | [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/) |
|:---|:---:|:---:|:---:|
| **Normal binding free energies**| PB | PB and GB  | PB and GB |
| GB models | None | 1, 2, 3, 5, 7 and 8  | 1, 2, 3, 5, 7 and 8 |
| **Stability** |  | :heavy_check_mark: | :heavy_check_mark: |
| **Alanine scanning** | :heavy_check_mark: <sup>2</sup>| :heavy_check_mark: | :heavy_check_mark: |
| **Entropy corrections** <sup>3</sup>|   | nmode and qh | nmode, qh, and IE |
| **Decomposition schemes** | Per-Residues | Per-Residues and Per-Wise | Per-Residues and Per-Wise |
| **QM/MMGBSA** |   | :heavy_check_mark: | :heavy_check_mark: |
| **MM/3D-RISM** |   | :heavy_check_mark: | :heavy_check_mark: |
| **Membrane Protein MMPBSA** |      | :heavy_check_mark: | :heavy_check_mark: |
| **GROMACS Version** | 4.x, 5.x and 2016+ <sup>4</sup> |  --- | 4.x, 5.x and 20xx.x |
| **Approximations** | ST | ST and MT | ST and MT |
| **API** |      | :heavy_check_mark: | :heavy_check_mark: |
| **Graphical Analyzer** |    |  | :heavy_check_mark: |
| Energy to PDB | :heavy_check_mark: |     | :heavy_check_mark: |
| Energetic Terms charts | Per-Frame | Average and/or Per-Frame <sup>5</sup>  | Average and Per-frame |
| Energetic Terms charts representation | xmgrace/matplotlib/gnuplot | API and graphics library | gmx_MMPBSA_gui |
| **Externals programs** | APBS (1.2.x, 1.3.x or 1.4.x) |  AmberTools20 | AmberTools20 |
| **Parallel computation** | Depends on APBS version | :heavy_check_mark: | :heavy_check_mark: |
| **Steps for:** | | | |
| Calculation and Summary | Multiple | One | One |
| Analysis| Multiple | Multiple | One |

<sup>1</sup> MMPBSA.py is included in AMBER package.

<sup>2</sup> Without documentation.

<sup>3</sup> nmode = Normal modes approximation, qh = Quasic-Harmony approximation and IE = Interaction Entropy
approximation

<sup>4</sup> GROMACS 20xx.x is not officially supported. There is a Pull Request that offers a minimum of compatibility 
with versions higher than 2016.x but with limitations 

<sup>5</sup> The user can obtain each energetic term per frame or its average values using the API. This means that
 user must be familiar with Python to handle the API, perform custom calculations or graph such data.

### Examples...
* [Protein-DNA binding free energy calculations](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_DNA)
* [Protein-ligand binding free energy calculations (Single Trajectory method)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_ligand/ST) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.php))
* [Protein-ligand binding free energy calculations (Multiple Trajectory method)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_ligand/MT) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.php))
* [MMPBSA with membrane proteins](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_membrane)
* [Protein-protein binding free energy calculations](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_protein) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section2.php))
* [Protein-protein binding free energy calculations with MM/3D-RISM](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/3D-RISM)
* [Alanine scanning](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Alanine_scanning) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section3.php))
* [Decomposition analysis](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Decomposition_analysis) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section6.php))
* [Entropy calculations with normal modes](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Entropy_calculations/nmode) (based on this [tutorial](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section5.php))
* [Entropy calculations with Interaction Entropy](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Entropy_calculations/Interaction_Entropy)
* [Stability calculations](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Stability)
* [Protein-glycan binding free energy calculations](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files/Protein_glycan)

## Calling `gmx_MMPBSA` from the command-line
`gmx_MMPBSA` is invoked through the command line as follows:
```
usage: gmx_MMPBSA [-h] [-v] [--input-file-help] [-O] [-prefix <file prefix>] [-i FILE]
                  [-xvvfile XVVFILE] [-o FILE] [-do FILE] [-eo FILE] [-deo FILE] [-gui] [-s] 
                  [-cs <Structure File>] [-ci <Index File>] [-cg index index] [-ct [TRJ [TRJ ...]]]
                  [-rs <Structure File>] [-ri <Index File>] [-rg index] [-rt [TRJ [TRJ ...]]] 
                  [-lm <Structure File>] [-ls <Structure File>] [-li <Index File>] [-lg index] 
                  [-lt [TRJ [TRJ ...]]] [-make-mdins] [-use-mdins] [-rewrite-output]

gmx_MMPBSA is an effort to bring in Amber's MMPBSA.py functionalities and more to GROMACS users. This 
program is based on Amber's MMPBSA.py and essentially works as such. gmx_MMPBSA minimizes 
compatibility-related issues since it will process any GROMACS files compatible with the GROMACS 
in the path.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --input-file-help     Print all available options in the input file. (default: False)

Miscellaneous Options:
  -O, --overwrite       Allow output files to be overwritten (default: False)
  -prefix <file prefix>
                        Prefix for intermediate files. (default: _GMXMMPBSA_)

Input and Output Files:
  These options specify the input files and optional output files.

  -i FILE               MM/PBSA input file. (default: None)
  -xvvfile XVVFILE      XVV file for 3D-RISM. (default: $AMBERHOME/dat/mmpbsa/spc.xvv)
  -o FILE               Output file with MM/PB(GB)SA statistics. (default: FINAL_RESULTS_MMPBSA.dat)
  -do FILE              Output file for decomposition statistics. (default: FINAL_DECOMP_MMPBSA.dat)
  -eo FILE              CSV-format output of all energy terms for every frame in every calculation.
                         File name forced to end in [.csv]. This file is only written when specified
                         on the command-line. (default: None)
  -deo FILE             CSV-format output of all energy terms for each printed residue in
                         decomposition calculations. 
                         File name forced to end in [.csv]. This file is only written when specified
                         on the command-line. (default: None)
  -gui                  Open GUI plotting app when all calculations are finished (default: True)
  -s                    Perform stability calculation. Only the complex parameters are required. If
                         ligand is non-Protein (small molecule) type, then ligand *.mol2 file is 
                         required. In any other case receptor and ligand parameters will be ignored.
                         See description bellow (default: False)

Complex:
  Complex files and info that are needed to perform the calculation. If the receptor and / or the
    ligand info is not defined, we generate them from that of the complex.

  -cs <Structure File>  Structure file of the complex. If it is Protein-Ligand (small molecule)
                         complex, make sure that you define -lm option. See -lm description below
                         Allowed formats: *.tpr (recommended), *.pdb, *.gro (default: None)
  -ci <Index File>      Index file of the bound complex. (default: None)
  -cg index index       Groups of receptor and ligand in complex index file. (default: None)
                         Notation: "-cg <Receptor group> <Ligand group>", ie. -cg 1 13
  -ct [TRJ [TRJ ...]]   Input trajectories of the complex. Make sure the trajectory is fitted and
                         pbc have been removed. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                         (specify as many as you'd like). (default: None)

Receptor:
  Receptor files and info that are needed to perform the calculation. If the receptor info is not 
   defined, we generate it from that of the complex.

  -rs <Structure File>  Structure file of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.tpr (recommended), *.pdb, *.gro (default: None)
  -ri <Index File>      Index file of the unbound receptor. (default: None)
  -rg index             Receptor group in receptor index file. Notation: "-lg <Receptor group>", 
                         e.g. -rg 1 (default: None)
  -rt [TRJ [TRJ ...]]   Input trajectories of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like). (default: None)

Ligand:
  Ligand files and info that are needed to perform the calculation. If the ligand are not defined, 
   we generate it from that of the complex.

  -lm <Structure File>  A *.mol2 file of the unbound ligand used to parametrize ligand for GROMACS
                         using Anetchamber. Must be defined if Protein-Ligand (small molecule) 
                         complex was define. No needed for Proteins, DNA, RNA, Ions and Glycans.
                         Antechamber output *.mol2 is recommended. (default: None)
  -ls <Structure File>  Structure file of the unbound ligand. If ligand is a small molecule, make 
                         sure that you define above -lm option. Allowed formats: *.tpr (recommended),
                         *.pdb, *.gro (default: None)
  -li <Index File>      Index file of the unbound ligand. (default: None)
  -lg index             Ligand group in ligand index file. Notation: "-lg <Ligand group>", 
                         e.g. -lg 13 (default: None)
  -lt [TRJ [TRJ ...]]   Input trajectories of the unbound ligand for multiple trajectory approach. 
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like). (default: None)

Miscellaneous Actions:
  -make-mdins           Create the input files for each calculation and quit. This allows you to 
                         modify them and re-run using -use-mdins (default: False)
  -use-mdins            Use existing input files for each calculation. If they do not exist with the
                         appropriate names, gmx_MMPBSA will quit in error. (default: False)
  -rewrite-output       Do not re-run any calculations, just parse the output files from the 
                         previous calculation and rewrite the output files. (default: False)

This program will calculate binding free energies using end-state free energy methods on an ensemble
of snapshots using a variety of implicit solvent models
```

## Running `gmx_MMPBSA`
### Serial version
This version is installed via pip as described above. AMBERHOME variable must be set, or it will quit with an error. An example 
command-line call is shown below:

    gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc

You can found test files in [GitHub](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/test_files)

### Parallel (MPI) version
Unlike MMPBSA.py, `gmx_MMPBSA` will be installed as a separate package from the Amber installation. When installing Amber with mpi,
a MMPBSA.py version called "MMPBSA.py.MPI" will be installed as well. Since we cannot detect if Amber was installed one way or
another, we simply decided to adapt the `gmx_MMPBSA` executable to use an argument. That is, `gmx_MMPBSA` is a single script
that executes the serial version or the parallel version with mpi depending on whether the user defines the "mpi" or
"MPI" argument. In principle, both the serial and parallel versions should work correctly when Amber was installed in parallel.

The parallel version, like MMPBSA.py.MPI requires the mpi4py module. If you did the parallel installation of Amber, it 
should be installed. In any case, it could be installed in the following way:

    amber.python -m pip install mpi4py
    
A usage example is shown below:

    mpirun -np 2 gmx_MMPBSA MPI -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc

or

    mpirun -np 2 gmx_MMPBSA mpi -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc
    
One note: at a certain level, running RISM in parallel may actually hurt performance, since previous solutions are used 
as an initial guess for the next frame, hastening convergence. Running in parallel loses this advantage. Also, due to 
the overhead involved in which each thread is required to load every topology file when calculating energies, parallel 
scaling will begin to fall off as the number of threads reaches the number of frames. 

## Input and Output

### The input file
As `gmx_MMPBSA` is based on [MMPBSA.py](https://pubs.acs.org/doi/10.1021/ct300418h), it uses an input file containing 
all the specification for the MM/PB(GB)SA calculation. The input file is designed to be as syntactically similar to 
other programs in Amber as possible. The input file has the same namelist structure as both sander and pmemd. The allowed 
namelists are &general, &gb, &pb, &rism, &alanine_scanning, &nmode, and &decomp. The input variables recognized in each 
namelist are described below, but those in &general are typically variables that apply to all aspects of the calculation
or parameters required for build amber topologies from GROMACS files.
The &gb namelist is unique to Generalized Born calculations, &pb is unique to Poisson Boltzmann calculations, &rism is 
unique to 3D-RISM calculations, &alanine_scanning is unique to alanine scanning calculations, &nmode is unique to the
normal mode calculations used to approximate vibrational entropies, and &decomp is unique to the decomposition
scheme. All of the input variables are described below according to their respective namelists. Integers and floating
point variables should be typed as-is while strings should be put in either single- or double-quotes. All variables
should be set with `variable = value` and separated by commas. See several 
[examples](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA#sample-input-files) below. As you will see, several 
calculations can be performed in the same run (_i.e._ &gb and &pb, &gb and &alanine_scanning, &pb and &decomp, etc).
Variables will usually be matched to the minimum number of characters required to uniquely identify that variable 
within that namelist. Variables require at least 4 characters to be matched unless that variable name has fewer than 4 
characters (in which case the whole variable name is required). For example, “star” in &general will match “startframe”. 
However, “stare” and “sta” will match nothing.

**&general namelist variables**

`debug_printlevel` MMPBSA.py prints errors by raising exceptions, and not catching fatal errors. If debug_printlevel 
is set to 0, then detailed tracebacks (effectively the call stack showing exactly where in the program the 
error occurred) is suppressed, so only the error message is printed. If debug_printlevel is set to 1 or 
higher, all tracebacks are printed, which aids in debugging of issues. Default: 0. (Advanced Option)

`startframe` The frame from which to begin extracting snapshots from the full, concatenated trajectory comprised
of every trajectory file placed on the command-line. This is always the first frame read. (Default = 1)  
          
`endframe` The frame from which to stop extracting snapshots from the full, concatenated trajectory comprised of every 
 trajectory file supplied on the command-line. (Default = 9999999)

```diff
@@ Input variable modified. Included Interaction entropy aproximation @@
```            
`entropy` It specifies whether to perform a quasi-harmonic entropy (QH) approximation with ptraj or the 
[Interaction Entropy (IE)](https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) approximation. The allowed values are 
(default = 0): 
* 0: Don’t
* 1: perform QH
* 2: perform IE

```diff
+ New input variable added
```
`entropy_seg` Specify the representative segment (in %), starting from the `endframe`, for the calculation of 
the Interaction Entropy, _e.g._: `entropy_seg = 25` means that the last quartile of the total number of frames 
(`(endframe-startframe)/interval`) will be used to calculate the average Interaction Entropy. Default: 25 (Only if `entropy = 2`)

```diff
+ New input variable added
```
`entropy_temp` Specify the temperature to calculate the entropy term `−TΔS` (Only if `entropy` = 2). Avoid 
inconsistencies with defined internal temperature (298.15 K) when nmode is used (Default = 298.15)

```diff
+ New input variable added
```
`gmx_path` Define an additional path to search for GROMACS executables. This path takes precedence over the path 
defined in the PATH variable. In these path the following executables will be searched: `gmx`, `gmx_mpi`, `gmx_d`, 
`gmx_mpi_d` (Gromcas > 5.x.x), `make_ndx` and `trjconv` (GROMACS 4.x.x) 
   
`interval` The offset from which to choose frames from each trajectory file. For example, an interval of 2 will pull
every 2nd frame beginning at startframe and ending less than or equal to endframe. (Default = 1)

```diff
- Input variable deleted. All files are needed for plotting
``` 
~~-`keep_files` The variable that specifies which temporary files are kept. All temporary files have the prefix 
`_GMXMMPBSA_` prepended to them (unless you change the prefix on the command-line—see subsection Subsection 34.3.2 for 
details). Allowed values are 0, 1, and 2. 0: Keep no temporary files 1: Keep all generated trajectory files and mdout 
files created by sander simulations 2: Keep all temporary files. Temporary files are only deleted if MMPBSA.py 
completes successfully (Default = 1) A verbose level of 1 is sufficient to use -rewrite-output and recreate the output
 file without rerunning any simulations.~~
            
`netcdf` Specifies whether or not to use NetCDF trajectories internally rather than writing temporary ASCII trajectory 
files. For very large trajectories, this could offer significant speedups, and requires less temporary space.
However, this option is incompatible with alanine scanning. Default value is 0.
* 0: Do NOT use temporary NetCDF trajectories
* 1: Use temporary NetCDF trajectories

```diff
+ New input variable added
```
`PBRadii` PBRadii to build amber topology files (Default = 3):
* 1: bondi, recommended when igb = 7
* 2: mbondi, recommended when igb = 1
* 3: mbondi2, recommended when igb = 2 or 5
* 4: mbondi3, recommended when igb = 8

```diff
+ New input variable added
```
`protein_forcefield` Define the force field used to build Amber topology for proteins. Make sure this force field is 
the same as the one used in GROMACS (Default = "oldff/leaprc.ff99SB")
Force fields tested:
* "oldff/leaprc.ff99"
* "oldff/leaprc.ff03" 
* "oldff/leaprc.ff99SB" 
* "oldff/leaprc.ff99SBildn"
* "leaprc.protein.ff14SB"

```diff
+ New input variable added
```
`ligand_forcefield` Define the force field used to build Amber topology for small molecules or glycams. Make sure this 
force field is the same as the one used for GROMACS (Default = "leaprc.gaff").
Force fields tested:
* "leaprc.gaff"
* "leaprc.gaff2"
* "leaprc.GLYCAM_06j-1"    (Compatible with amber12SB and later)
* "leaprc.GLYCAM_06EPb"    (Compatible with amber12SB and later)
* "leaprc.GLYCAM_06h-1"    (Included in gmx_MMPBSA package. If it is selected will be copied to $AMBERHOME/dat/leap. Compatible with amber99SB and earlier)

```diff
+ New input variable added
```
`ions_parameters` Define ions parameters to build the Amber topology. (Default = 1)
* 1: frcmod.ions234lm_126_tip3p
* 2: frcmod.ions234lm_iod_tip4pew
* 3: frcmod.ions234lm_iod_spce
* 4: frcmod.ions234lm_hfe_spce
* 5: frcmod.ions234lm_126_tip4pew
* 6: frcmod.ions234lm_126_spce
* 7: frcmod.ions234lm_1264_tip4pew
* 8: frcmod.ions234lm_1264_tip3p
* 9: frcmod.ions234lm_1264_spce
* 10: frcmod.ions234lm_iod_tip3p
* 11: frcmod.ions234lm_hfe_tip4pew
* 12: frcmod.ions234lm_hfe_tip3p

```diff
+ New input variable added
```
`reuse_files` Define whether the trajectories files will be reused when the program ends in error. Note that the 
trajectories files may not be generated correctly due to internal errors or interruptions. Please use it with care. 
(Default = 0) 
* 0: Don't reuse. If there are temporary trajectory files, they will be deleted
* 1: Reuse existing trajectory file

```diff
+ New input variable added
```
`solvated_trajectory` Define if it is necessary to build a clean trajectory with no water and ions (Default = 1)
* 0: Don’t
* 1: Build clean trajectory

```diff
- Input variable deleted. ALways must be defined to get GROMACS
```
~~-`search_path` Advanced option. By default, MMPBSA.py will only search for executables in $AMBERHOME/bin .
To enable it to search for binaries in your full PATH if they can’t be found in $AMBERHOME/bin , set
search_path to 1. Default 0 (do not search through the PATH ). This is particularly useful if you are using
an older version of sander that is not in AMBERHOME .~~

`use_sander` use sander for energy calculations, even when mmpbsa_py_energy will suffice (Default = 0)
* 0: Use mmpbsa_py_energy when possible
* 1: Always use sander

`verbose` The variable that specifies how much output is printed in the output file. There are three allowed 
values (Default = 1): 
* 0: print difference terms 
* 1: print all complex, receptor, and ligand terms 
* 2: also print bonded terms if one trajectory is used 

**&gb namelist variables**

`ifqnt` Specifies whether a part of the system is treated with quantum mechanics. This functionality requires sander 
igb Generalized Born method to use (seeSection 4). Allowed values are 1, 2, 5, 7 and 8. (Default = 5) All models are 
now available with both mmpbsa_py_energy and sander.(Default = 0)
* 0: Potential function is strictly classical 
* 1: Use QM/MM

`qm_residues` Comma- or semicolon-delimited list of complex residues to treat with quantum mechanics. All
whitespace is ignored. All residues treated with quantum mechanics in the complex must be treated with
quantum mechanics in the receptor or ligand to obtain meaningful results.

```diff
+ Input variable added.
```
`intdiel` Define Internal dielectric constant without use external *.mdin file (Default = 1.0)

`qm_theory` Which semi-empirical Hamiltonian should be used for the quantum calculation. No default, this must
be specified. See its description in the QM/MM section of the manual for options.

`qmcharge_com` The charge of the quantum section for the complex. (Default = 0)

`qmcharge_lig` The charge of the quantum section of the ligand. (Default = 0)

`qmcharge_rec` The charge of the quantum section for the receptor. (Default = 0)

`qmcut` The cutoff for the qm/mm charge interactions. (Default = 9999.0)

`saltcon` Salt concentration in Molarity. (Default = 0.0)

`surfoff` Offset to correct (by addition) the value of the non-polar contribution to the solvation free energy term
(Default 0.0)

`surften` Surface tension value (Default = 0.0072). Units in kcal/mol/ Å 2

`molsurf` When set to 1, use the molsurf algorithm to calculate the surface area for the nonpolar solvation term.
When set to 0, use LCPO (Linear Combination of Pairwise Overlaps). (Default 0)

`probe` Radius of the probe molecule (supposed to be the size of a solvent molecule), in Angstroms, to use when
determining the molecular surface (only applicable when molsurf is set to 1). Default is 1.4.

`msoffset` Offset to apply to the individual atomic radii in the system when calculating the molsurf surface. See
the description of the molsurf action command in [cpptraj](https://ambermd.org/doc12/Amber20.pdf#chapter.32). Default is 0.

**&pb namelist variables**

`inp` Option to select different methods to compute non-polar solvation free energy (Default = 2).
* 0: No non-polar solvation free energy is computed
* 1: The total non-polar solvation free energy is modeled as a single term linearly proportional to the solvent accessible
surface area, as in the PARSE parameter set, that is, if INP = 1, USE_SAV must be equal to 0.
* 2: The total non-polar solvation free energy is modeled as two terms: the cavity term and the dispersion term. The 
dispersion term is computed with a surface-based integration method closely related to the PCM solvent for quantum 
chemical programs. Under this framework, the cavity term is still computed as a term linearly proportional to the 
molecular solvent-accessible-surface area (SASA) or the molecular volume enclosed by SASA.

`cavity_offset` Offset value used to correct non-polar free energy contribution (Default = -0.5692) This is not used
for APBS.

`cavity_surften` Surface tension. (Default = 0.0378 kcal/mol Angstrom 2 ). Unit conversion to kJ done automatically for 
APBS.

`exdi` External dielectric constant (Default = 80.0).

`indi` Internal dielectric constant (Default = 1.0).

`fillratio` The ratio between the longest dimension of the rectangular finite-difference grid and that of the solute
(Default = 4.0).

`scale` Resolution of the Poisson Boltzmann grid. It is equal to the reciprocal of the grid spacing. (Default = 2.0)

`istrng` Ionic strength in Molarity. It is converted to mM for PBSA and kept as M for APBS. (Default = 0.0)

`linit` Maximum number of iterations of the linear Poisson Boltzmann equation to try (Default = 1000)

`prbrad` Solvent probe radius in Angstroms. Allowed values are 1.4 and 1.6 (Default = 1.4)

`radiopt` The option to set up atomic radii according to 0: the prmtop, or 1: pre-computed values (see Amber
manual for more complete description). (Default = 1)

`sander_apbs` Option to use APBS for PB calculation instead of the built-in PBSA solver. This will work only
through the iAPBS interface built into sander.APBS. Instructions for this can be found online at the
iAPBS/APBS websites. Allowed values are 0: Don’t use APBS, or 1: Use sander.APBS. (Default = 0)

`memopt` Turn on membrane protein support (Default = 0).

`emem` Membrane dielectric constant (Default = 1.0).

`mthick` Membrane thickness (Default = 40.0).

`mctrdz` Absolute membrane center in the z-direction (Default=0.0, use protein center as the membrane center).

`poretype` Turn on the automatic membrane channel/pore finding method (Default=1).

A more thorough description of these and other options can be found [here](https://ambermd.org/doc12/Amber20.pdf#chapter.6).
Please also note that the default options have changed over time. For a detailed discussion of all related options on 
the quality of the MM/PB(GB)SA calculations, please check this [publication](https://onlinelibrary.wiley.com/doi/10.1002/jcc.24467).

**&alanine_scanning namelist variables**

`mutant_only` Option to perform specified calculations only for the mutants. Allowed values are 0: Do mutant and
original or 1: Do mutant only (Default = 0)
Note that all calculation details are controlled in the other namelists, though for alanine scanning to be performed,
the namelist must be included (blank if desired)
```diff
+Two options added to ease the alanine_scanning calculations
```
`mutant` Define whether the mutation will be perform in receptor or ligand. Allowed values are: receptor, rec, ligand or lig
 in any capitalization (Default = receptor or REC)

`mutant_res` Define the specific residue that is going to be mutated. Use the following format CHAIN:RESNUM (eg: 'A:350').
Please, make sure that your selection is correct and based on GROMACS numbering in processed files.

**&nmode namelist variables**

`dielc` Distance-dependent dielectric constant (Default = 1.0)

`drms` Convergence criteria for minimized energy gradient. (Default = 0.001)

`maxcyc` Maximum number of minimization cycles to use per snapshot in sander. (Default = 10000)

`nminterval` ∗ Offset from which to choose frames to perform nmode calculations on (Default = 1)

`nmendframe` ∗ Frame number to stop performing nmode calculations on (Default = 1000000)

`nmode_igb` Value for Generalized Born model to be used in calculations. Options are 0: Vacuum, 1: HCT GB
model (Default 1)

`nmode_istrng` Ionic strength to use in nmode calculations. Units are Molarity. Non-zero values are ignored if `nmode_igb`
 is 0 above. (Default = 0.0)

`nmstartframe` ∗ Frame number to begin performing nmode calculations on (Default = 1)

* These variables will choose a subset of the frames chosen from the variables in the &general namelist. Thus, the
“trajectory” from which snapshots will be chosen for nmode calculations will be the collection of snapshots upon
which the other calculations were performed.

**&decomp namelist variables**

`csv_format` Print the decomposition output in a Comma-Separated-Variable (CSV) file. CSV files open natively
in most spreadsheets. If set to 1, this variable will cause the data to be written out in a CSV file, and standard
error of the mean will be calculated and included for all data. If set to 0, the standard, ASCII format will be
used for the output file. Default is 1 (CSV-formatted output file)

`dec_verbose` Set the level of output to print in the decomp_output file.
* 0: DELTA energy, total contribution only
* 1: DELTA energy, total, sidechain, and backbone contributions
* 2: Complex, Receptor, Ligand, and DELTA energies, total contribution only
* 3: Complex, Receptor, Ligand, and DELTA energies, total, sidechain, and backbone contributions

Note: If the values 0 or 2 are chosen, only the Total contributions are required, so only those will be printed
to the mdout files to cut down on the size of the mdout files and the time required to parse them. Default = 0

`idecomp` Energy decomposition scheme to use:
* 1: Per-residue decomp with 1-4 terms added to internal potential terms
* 2: Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms.
* 3: Pairwise decomp with 1-4 terms added to internal potential terms
* 4: Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms

(No default. This must be specified!) This functionality requires sander.

`print_res` Select residues from the complex to print. Default is print "within 6". This variable also accepts a 
sequence of individual residues and/or ranges. The different fields must be either comma- or semicolon-delimited. 
For example: print_res = "within 6", where _within_ corresponds to the keyword and _6_ to the maximum distance 
criterion in Angstroms necessary to select the residues from both the receptor and the ligand; or 
print_res = “1, 3-10, 15, 100”, or print_res = “1; 3-10; 15; 100”. Both of these will print residues 1, 3 
through 10, 15, and 100 from the complex topology file and the corresponding residues in either the ligand 
and/or receptor topology files.

```diff
- *Please note: Using idecomp=3 or 4 (pairwise) with a very large number of printed residues and a
-  large number of frames can quickly create very, very large temporary mdout files. Large print 
-  selections also demand a large amount of memory to parse the mdout files and write 
-  decomposition output file (~500 MB for just 250 residues, since that’s 62500 pairs!) It is not
-  unusual for the output file to take a significant amount of time to print if you have a lot of
-  data. This is most applicable to pairwise decomp, since the amount of data scales as O(N 2 ).

+ Based on the above, we decided to add a new option that limits the selection of the residues 
+ that will be printed by default. We defined a selection method with the following structure:
+ print_res = "within 6" where _within_ corresponds to the keyword and _6_ to the maximum 
+ distance criterion in Angstroms necessary to select the residues from both the receptor and
+ the ligand.
```

**&rism namelist variables***

`buffer` Minimum distance between solute and edge of solvation box. Specify this with grdspc below. Mutually
exclusive with ng and solvbox. Set buffer < 0 if you wish to use ng and solvbox. (Default = 14 Å)
closure The approximation to the closure relation. Allowed choices are kh (Kovalenko-Hirata), hnc (Hypernetted-
chain), or psen (Partial Series Expansion of order-n) where “n” is a positive integer (_e.g._, “pse3”). (Default
= ‘kh’)

`closureorder` (Deprecated) The order at which the PSE-n closure is truncated if closure is specified as “pse” or
“psen” (no integers). (Default = 1)

`grdspc` Grid spacing of the solvation box. Specify this with buffer above. Mutually exclusive with ng and solvbox.
(Default = 0.5 Å)

`ng` Number of grid points to use in the x, y, and z directions. Used only if buffer < 0. Mutually exclusive with
buffer and grdspc above, and paired with solvbox below. No default, this must be set if buffer < 0. Define
like “ng=1000,1000,1000”

`polardecomp` Decompose the solvation free energy into polar and non-polar contributions. Note that this will
increase computation time by roughly 80%. 0: Don’t decompose solvation free energy. 1: Decompose
solvation free energy. (Default = 0)

`rism_verbose` Level of output in temporary RISM output files. May be helpful for debugging or following con-
vergence. Allowed values are 0 (just print the final result), 1 (additionally prints the total number of iterations
for each solution), and 2 (additionally prints the residual for each iteration and details of the MDIIS solver).
(Default = 0)

`solvbox` Length of the solvation box in the x, y, and z dimensions. Used only if buffer < 0. Mutually exclusive
with buffer and grdspc above, and paired with ng above. No default, this must be set if buffer < 0. Define
like “solvbox=20,20,20”

`solvcut` Cutoff used for solute-solvent interactions. The default is the value of buffer. Therefore, if you set buffer
< 0 and specify ng and solvbox instead, you must set solvcut to a nonzero value or the program will quit in
error. (Default = buffer )

`thermo` Which thermodynamic equation you want to use to calculate solvation properties. Options are “std”, “gf”,
or “both” (case-INsensitive). “std” uses the standard closure relation, “gf” uses the Gaussian Fluctuation
approximation, and “both” will print out separate sections for both. (Default = “std”). Note that all data are
printed out for each RISM simulation, so no choice is any more computationally demanding than another.
Also, you can change this option and use the -rewrite-output flag to obtain a different printout after-the-fact.

`tolerance` Upper bound of the precision requirement used to determine convergence of the self-consistent solution.
This has a strong effect on the cost of 3D-RISM calculations. (Default = 1e-5).

* 3D-RISM calculations are performed with the rism3d.snglpnt program built with AmberTools, written by
Tyler Luchko. It is the most expensive, yet most statistical mechanically rigorous solvation model available in
MMPBSA.py. See [Chapter 7](https://ambermd.org/doc12/Amber20.pdf#chapter.7) for a more thorough description of options 
and theory. A list of references can be found there, too. One advantage of 3D-RISM is that an arbitrary solvent can 
be chosen; you just need to change the xvvfile specified on the command line (see [34.3.2](https://ambermd.org/doc12/Amber20.pdf#subsection.34.3.2)).

#### Sample input files
```
Sample input file for GB and PB calculation
&general
startframe=5, endframe=100, interval=5,
verbose=2, protein_forcefield="oldff/leaprc.ff99SB", ligand_forcefield="leaprc.gaff",
/
&gb
igb=5, saltcon=0.150,
/
&pb
istrng=0.15, fillratio=4.0
/
--------------------------------------------------------
Sample input file for Alanine scanning
&general
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield="oldff/leaprc.ff99SB", PBRadii=4
/
&gb
igb=8, saltcon=0.150, intdiel=10
/
&alanine_scanning
#make sure to change this parameter to 'ligand' is the mutation is going to be performed 
#in the ligand
mutant='receptor'
mutant_res='B:65'
/
--------------------------------------------------------
Sample input file for entropy calculations
&general
#
startframe=5, endframe=21, verbose=2, interval=1,
#entropy variable control whether to perform a quasi-harmonic entropy (QH) approximation or the 
#Interaction Entropy (IE)(https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) approximation
protein_forcefield="oldff/leaprc.ff99SB", entropy=2, entropy_seg=25, entropy_temp=298
/
&gb
igb=2, saltcon=0.150,
/
#uncomment the next 4 lines for normal mode calculations
#&nmode
#nmstartframe=5, nmendframe=21, nminterval=2,
#maxcyc=50000, drms=0.0001,
#/
/
--------------------------------------------------------
Sample input file with decomposition analysis
#make sure to include at least one residue from both the receptor
#and ligand in the print_res mask of the &decomp section.
#http://archive.ambermd.org/201308/0075.html
&general
startframe=5, endframe=21, interval=1,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=3,
print_res="within 4"
#check _GMXMMPBSA_COM_FIXED.pdb file to select which residues are going to be printed
#in the output file
#print_res="40-41,44,47,78,81-82,85,88,115,118,122,215,218-220,232,241"
/
--------------------------------------------------------
Sample input file for QM/MMGBSA
&general
startframe=5, endframe=100, interval=5,
/
&gb
igb=5, saltcon=0.100, ifqnt=1, qmcharge_com=0,
qm_residues="100-105, 200", qm_theory="PM3"
/
--------------------------------------------------------
Sample input file for MM/3D-RISM
&general
startframe=20, endframe=100, interval=5,
/
&rism
polardecomp=1, thermo="gf"
/
--------------------------------------------------------
Sample input file for MMPBSA with membrane proteins
&general
startframe=1, endframe=100, interval=1, 
debug_printlevel=2, use_sander=1,
/
&pb
radiopt=0, indi=20.0, istrng=0.150,
fillratio=1.25, ipb=1, nfocus=1,
bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0,
npbverb=1, solvopt=2, inp=2,
memopt=1, emem=7.0, mctrdz=-10.383, mthick=36.086, poretype=1,
maxarcdot=15000
/
```

A few important notes about input files. Comments are allowed by placing a # at the beginning of the line (whites-
pace is ignored). Variable initialization may span multiple lines. In-line comments (_i.e._, putting a # for a comment
after a variable is initialized in the same line) is not allowed and will result in an input error. Variable declarations
must be comma-delimited, though all whitespace is ignored. Finally, all lines between namelists are ignored, so
comments can be added before each namelist without using #.

### The Output File
The header of the output file will contain information about the calculation. It will show a copy of the input
file as well as the names of all files that were used in the calculation (topology files and coordinate file(s)). If the
masks were not specified, it prints its best guess so that you can verify its accuracy, along with the residue name of
the ligand (if it is only a single residue).
The energy and entropy contributions are broken up into their components as they are in sander and nmode or
ptraj. The contributions are further broken into G gas and G solv . The polar and non-polar contributions are EGB (or
EPB) and ESURF (or ECAVITY / ENPOLAR), respectively for GB (or PB) calculations.
By default, bonded terms are not printed for any one-trajectory simulation. They are computed and their dif-
ferences calculated, however. They are not shown (nor included in the total) unless specifically asked for because
they should cancel completely. A single trajectory does not produce any differences between bond lengths, angles,
or dihedrals between the complex and receptor/ligand structures. Thus, when subtracted they cancel completely.
This includes the BOND, ANGLE, DIHED, and 1-4 interactions. If inconsistencies are found, these values are
displayed and inconsistency warnings are printed. When this occurs the results are generally useless. Of course
this does not hold for the multiple trajectory protocol, and so all energy components are printed in this case.
Finally, all warnings generated during the calculation that do not result in fatal errors are printed after calculation
details but before any results.

### Temporary Files
`gmx_MMPBSA` creates working files during the execution of the script beginning with the prefix `_GMXMMPBSA_`.
If `gmx_MMPBSA` does not finish successfully, several of these files may be helpful in diagnosing the problem.
For that reason, every temporary file is described below. Note that not every temporary file is generated in every 
simulation. At the end of each description, the lowest value of the original “keep_files” variable that will retain 
this file will be shown in parentheses. Nevertheless, in the current version, all the files are retained for plotting 
purposes.

`make_top.log` This file contains the output coming from all the GROMACS programs.

`leap.log` This file contains the output coming from tleap program.

`_GMXMMPBSA_gb.mdin` Input file that controls the GB calculation done in sander. (2)

`_GMXMMPBSA_pb.mdin` Input file that controls the PB calculation done in sander. (2)

`_GMXMMPBSA_gb_decomp_com.mdin` Input file that controls the GB decomp calculation for the complex done in
sander. (2)

`_GMXMMPBSA_gb_decomp_rec.mdin` Input file that controls the GB decomp calculation for the receptor done in
sander. (2)

`_GMXMMPBSA_gb_decomp_lig.mdin` Input file that controls the GB decomp calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_pb_decomp_com.mdin` Input file that controls the PB decomp calculation for the complex done in
sander. (2)

`_GMXMMPBSA_pb_decomp_rec.mdin` Input file that controls the PB decomp calculation for the receptor done in
sander. (2)

`_GMXMMPBSA_pb_decomp_lig.mdin` Input file that controls the PB decomp calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_com.mdin` Input file that controls the GB QM/MM calculation for the complex done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_rec.mdin` Input file that controls the GB QM/MM calculation for the receptor done in sander.
(2)

`_GMXMMPBSA_gb_qmmm_lig.mdin` Input file that controls the GB QM/MM calculation for the ligand done in sander.
(2)

`_GMXMMPBSA_complex.mdcrd.#` Trajectory file(s) that contains only those complex snapshots that will be processed
by MMPBSA.py. (1)

`_GMXMMPBSA_ligand.mdcrd.#` Trajectory file(s) that contains only those ligand snapshots that will be processed by
MMPBSA.py. (1)

`_GMXMMPBSA_receptor.mdcrd.#` Trajectory file(s) that contains only those receptor snapshots that will be processed
by MMPBSA.py. (1)

`_GMXMMPBSA_complex_nc.#` Same as _GMXMMPBSA_complex.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_receptor_nc.#` Same as _GMXMMPBSA_receptor.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_ligand_nc.#` Same as _GMXMMPBSA_ligand.mdcrd.#, except in the NetCDF format. (1)

`_GMXMMPBSA_dummycomplex.inpcrd` Dummy inpcrd file generated by _GMXMMPBSA_complexinpcrd.in for use with
imin=5 functionality in sander. (1)

`_GMXMMPBSA_dummyreceptor.inpcrd` Same as above, but for the receptor. (1)

`_GMXMMPBSA_dummyligand.inpcrd` Same as above, but for the ligand. (1)

`_GMXMMPBSA_complex.pdb` Dummy PDB file of the complex required to set molecule up in nab programs

`_GMXMMPBSA_receptor.pdb` Dummy PDB file of the receptor required to set molecule up in nab programs

`_GMXMMPBSA_ligand.pdb` Dummy PDB file of the ligand required to set molecule up in nab programs

`_GMXMMPBSA_complex_nm.mdcrd.#` Trajectory file(s) for each thread with snapshots used for normal mode calcula-
tions on the complex. (1)

`_GMXMMPBSA_receptor_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calcula-
tions on the receptor. (1)

`_GMXMMPBSA_ligand_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calculations
on the ligand. (1)

`_GMXMMPBSA_ptrajentropy.in` Input file that calculates the entropy via the quasi-harmonic approximation. This
file is processed by ptraj. (2)

`_GMXMMPBSA_avgcomplex.pdb` PDB file containing the average positions of all complex conformations processed by

`_GMXMMPBSA_cenptraj.in.` It is used as the reference for the _GMXMMPBSA_ptrajentropy.in file above.
(1)

`_GMXMMPBSA_complex_entropy.out` File into which the entropy results from _GMXMMPBSA_ptrajentropy.in analysis
on the complex are dumped. (1)

`_GMXMMPBSA_receptor_entropy.out` Same as above, but for the receptor. (1)

`_GMXMMPBSA_ligand_entropy.out` Same as above, but for the ligand. (1)

`_GMXMMPBSA_ptraj_entropy.out` Output from running ptraj using _GMXMMPBSA_ptrajentropy.in. (1)

`_GMXMMPBSA_complex_gb.mdout.#` sander output file containing energy components of all complex snapshots done
in GB. (1)

`_GMXMMPBSA_receptor_gb.mdout.#` sander output file containing energy components of all receptor snapshots done
in GB. (1)

`_GMXMMPBSA_ligand_gb.mdout.#` sander output file containing energy components of all ligand snapshots done in
GB. (1)

`_GMXMMPBSA_complex_pb.mdout.#` sander output file containing energy components of all complex snapshots done
in PB. (1)

`_GMXMMPBSA_receptor_pb.mdout.#` sander output file containing energy components of all receptor snapshots done
in PB. (1)

`_GMXMMPBSA_ligand_pb.mdout.#` sander output file containing energy components of all ligand snapshots done in
PB. (1)

`_GMXMMPBSA_complex_rism.out.#` rism3d.snglpnt output file containing energy components of all complex snap-
shots done with 3D-RISM (1)

`_GMXMMPBSA_receptor_rism.out.#` rism3d.snglpnt output file containing energy components of all receptor snap-
shots done with 3D-RISM (1)

`_GMXMMPBSA_ligand_rism.out.#` rism3d.snglpnt output file containing energy components of all ligand snapshots
done with 3D-RISM (1)

`_GMXMMPBSA_pbsanderoutput.junk.#` File containing the information dumped by sander.APBS to STDOUT. (1)

`_GMXMMPBSA_ligand_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the ligand
for all snapshots. (1)

`_GMXMMPBSA_receptor_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the
receptor for all snapshots. (1)

`_GMXMMPBSA_complex_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the com-
plex for all snapshots. (1)

`_GMXMMPBSA_mutant_...` These files are analogs of the files that only start with `_GMXMMPBSA_` described above, but
instead refer to the mutant system of alanine scanning calculations.

`_GMXMMPBSA_*out.#` These files are thread-specific files. For serial simulations, only #=0 files are created. For
parallel, #=0 through NUM_PROC - 1 are created.

### Advanced Options
The default values for the various parameters as well as the inclusion of some variables over others in the
general MMPBSA.py input file were chosen to cover the majority of all MM/PB(GB)SA calculations that would
be attempted while maintaining maximum simplicity. However, there are situations in which MMPBSA.py may
appear to be restrictive and ill-equipped to address. Attempts were made to maintain the simplicity described above
while easily providing users with the ability to modify most aspects of the calculation easily and without editing
the source code.

`-make-mdins` This flag will create all of the mdin and input files used by sander and nmode so that additional
control can be granted to the user beyond the variables detailed in the input file section above. The files
created are _GMXMMPBSA_gb.mdin which controls GB calculation; _GMXMMPBSA_pb.mdin which controls the
PB calculation; _GMXMMPBSA_sander_nm_min.mdin which controls the sander minimization of snapshots to
be prepared for nmode calculations; and _GMXMMPBSA_nmode.in which controls the nmode calculation. If
no input file is specified, all files above are created with default values, and _GMXMMPBSA_pb.mdin is created
for AmberTools’s pbsa. If you wish to create a file for sander.APBS, you must include an input file with
“sander_apbs=1” specified to generate the desired input file. Note that if an input file is specified, only those
mdin files pertinent to the calculation described therein will be created!

`-use-mdins` This flag will prevent MMPBSA.py from creating the input files that control the various calculations 
(_GMXMMPBSA_gb.mdin, _GMXMMPBSA_pb.mdin, _GMXMMPBSA_sander_nm_min.mdin, and _GMXMMPBSA_nmode.in). It will instead 
attempt to use existing input files (though they must have those names above!) in their place. In this way, the 
user has full control over the calculations performed, however care must be taken. The mdin files created by 
MMPBSA.py have been tested and are (generally) known to be consistent. Modifying certain variables (such as imin=5)
may prevent the script from working, so this should only be done with care. It is recommended that users start 
with the existing mdin files (generated by the -make-mdins flag above), and add and/or modify parameters from there.

`-make-mdins` and `-use-mdins` are intended to give added flexibility to user input. If the MM/PBSA input file does
not expose a variable you require, you may use the -make-mdins flag to generate the MDIN files and then quit.
Then, edit those MDIN files, changing the variables you need to, then running `gmx_MMPBSA` with -use-mdins to
use those modified files.

`QM/MMGBSA` There are a lot of options for QM/MM calculations in sander, but not all of those options were
made available via options in the MMPBSA.py input file. In order to take advantage of these other options,
you’ll have to make use of the -make-mdins and -use-mdins flags as detailed above and change the resulting
_GMXMMPBSA_gb_qmmm_com/rec/lig.mdin files to fit your desired calculation. Additionally, MMPBSA.py
suffers all shortcomings of sander, one of those being that PB and QM/MM are incompatible. Therefore,
only QM/MMGBSA is a valid option right now.

### Python API
The aim of the MMPBSA.py API is to provide you with direct access to the raw data produced during a 
MMPBSA.py calculation. By default, MMPBSA.py calculates an average, standard deviation, and standard error
of the mean for all the generated data sets, but it does not support custom analyses. The API reads an
_MMPBSA_info file, from which it will determine what kind of calculation was performed, then automatically
parse the output files and load the data into arrays.

It currently does NOT load decomposition data into available data structures. The topology files you used in 
the `gmx_MMPBSA` calculation must also be available in the location specified in the _GMXMMPBSA_info file.

#### Using the API
We have derived a new API to reorganize the data so that it is arranged more hierarchically. This makes
easier to transform the data into graphs in the `gmx_MMPBSA_gui`. **The original and the current API 
only differ in the name of the callable function, the disposition of the data in Per-wise decomposition analysis
and in the new 'delta' key. If you want to use the original, see the [Amber manual](http://ambermd.org/doc12/Amber20.pdf#section.34.4)**

The function `load_gmxmmpbsa_info` takes the name of a gmx_MMPBSA info file (typically _GMXMMPBSA_info)
and returns a populated `mmpbsa_data` instance with all the parsed data. An example code snippet that
creates a `mmpbsa_data` instance from the information in _GMXMMPBSA_info is shown below.

```python
from GMXMMPBSA import API as gmxMMPBSAapi
data = gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA_info")
```

#### Properties of `mmpbsa_data`

The `mmpbsa_data` class is a nested dictionary structure (`mmpbsa_data` is actually derived from `dict`). The
various attributes of `mmpbsa_data` are described below followed by the defined operators.

##### Attributes
If the numpy package is installed and available, all data arrays will be `numpy.ndarray` instances. 
Otherwise, all data arrays will be `array.array` instances with the ’d’ data type specifier (for a double
precision float). The data is organized in an `mmpbsa_data` instance in the following manner:

```python
mmpbsa_data_instance[calc_key][system_component][energy_term]
```

In this example, `calc_key` is a `dict` key that is paired to another `dict` (`mmpbsa_data_instance` is the
first-level `dict`, in this case) (Table 2). The keys of these second-level dict instances
(`system_component`) pair to another `dict` (Table 3).

**Table 2. List and description of `calc_key` dict keys that may be present in instances of the `mmpbsa_data`
class.**

| Dictionary Key (calc_key) | Calculation Type |
|:---:|:---|
| `gb`| Generalized Born Results |
| `pb` | Poisson-Boltzmann Results |
| `rism gf` | Gaussian Fluctuation 3D-RISM Results |
| `rism std` | Standard 3D-RISM Results |
| `nmode` | Normal Mode Analysis Results |
| `qh` | Quasi-harmonic Approximation Results |

**Table 3. List and description of system_component keys that may be present in instances of the mmpbsa_data
class.**

| Dictionary Key (system_component) | Description |
|:---:|:---|
| `complex` | Data sets for the complex. (Stability & Binding)|
| `receptor` | Data sets for the receptor. (Binding only)|
| `ligand` | Data sets for the ligand. (Binding only)|
| `delta`| Data sets for the delta. (Binding only) |

The keys of these inner-most (third-level) `dict` instances are paired with the data arrays for that energy
term (Table 4). The various dictionary keys are listed below for each level. If alanine scanning was performed, the
`mmpbsa_data_instance` also has a `"mutant"` attribute that contains the same dictionary structure as
`mmpbsa_data` does for the normal system. If not, the mutant attribute is None. The only difference is
that the data is accessed as follows:

```python
mmpbsa_data_instance.mutant[calc_key][system_component][energy_term]
```
Note, all keys are case-sensitive, and if a space appears in the key, it must be present in your program.
Also, if polar/non-polar decomposition is not performed for `3D-RISM`, then the `POLAR SOLV` and `APOLAR
SOLV` keys are replaced with the single key `ERISM`

**Table 4. List and description of energy_term keys that may be present in instances of the mmpbsa_data class.
The allowed values of energy_term depend on the value of calc_key above in Table 2. The
`energy_term` keys are listed for each `calc_key` enumerated above, accompanied by a description.
The RISM keys are the same for both `rism gf` and `rism std` although the value of `POLAR
SOLV` and `APOLAR SOLV` will differ depending on the method chosen. Those keys marked with *
are specific to the CHARMM force field used through chamber. Those arrays are all 0 for normal
Amber topology files.**

| Description | `gb` | `pb` | `RISM` |
|:---|:---:|:---:|:---:|
| Bond energy | `BOND` | `BOND` |  `BOND` |
| Angle energy | `ANGLE` | `ANGLE` |  `ANGLE` |
| Dihedral Energy | `DIHED` | `DIHED` |  `DIHED` |
| Urey-Bradley* | `UB` | `UB` |  — |
| Improper Dihedrals* | `IMP` | `IMP` |  — |
| Correction Map* | `CMAP` | `CMAP` |  — |
| 1-4 van der Waals energy | `1-4 VDW` | `1-4 VDW` |  `1-4 VDW` |
| 1-4 Electrostatic energy | `1-4 EEL` | `1-4 EEL` |  `1-4 EEL` |
| van der Waals energy | `VDWAALS` | `VDWAALS` |  `VDWAALS` |
| Electrostatic energy | `EEL` | `EEL` |  `EEL` |
| Polar solvation energy | `EGB` | `EPB` |  `POLAR SOLV` |
| Non-polar solvation energy | `ESURF` | `ENPOLAR` |  `APOLAR SOLV` |
| Total solvation free energy | `G solv` | `G solv` |  `G solv` |
| Total gas phase free energy | `G gas` | `G gas` |  `G gas` |
| Total energy | `TOTAL` | `TOTAL` |  `TOTAL` |

**Table 5. Same as Table 4 for the entropy data.**

| Description | `nmode`| `qh`|
|:---|:---:|:---:|
| Translational entropy | `Translational` | `Translational` |
| Rotational entropy | `Rotational` | `Rotational` |
| Vibrational entropy | `Vibrational` | `Vibrational` |
| Total entropy | `Total` | `Total` |


#### Defined operators
In-place addition: It extends all of the arrays that are common to both `mmpbsa_data` instances. This is
useful if, for instance, you run two `gmx_MMPBSA` calculations, and you use -prefix <new_prefix> for the
second simulation. Assuming that <new_prefix> is `_GMXMMPBSA2_` for the second `gmx_MMPBSA` calculation,
the following pseudo-code will generate an mmpbsa_data instance with all of the data in concatenated
arrays.
The pseudo-code assumes GMXMMPBSA.API was imported as demonstrated below.
```python
data = gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA_info")
data += gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA2_info")
```

### Example API Usage
In many cases, the autocorrelation function of the energy can aid in the analysis of MM/PBSA data, since it
provides a way of determining the statistical independence of your data points. For example, 1000
correlated snapshots provide less information, and therefore less statistical certainty, than 1000
uncorrelated snapshots. The standard error of the mean calculation performed by `gmx_MMPBSA` assumes a
completely uncorrelated set of snapshots, which means that it is a lower bound of the true standard error
of the mean, and a plot of the autocorrelation function may help determine the actual value.

The example program below will calculate the autocorrelation function of the total energy (complex only for
both the normal and alanine mutant systems) from a GB calculation and plot the resulting code using matplotlib.

```python
import os
import sys
# append AMBERHOME/bin to sys.path
sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))
# Now import the MMPBSA API
from GMXMMPBSA import API as gmxMMPBSAapi
import matplotlib.pyplot as plt
import numpy as np

data = gmxMMPBSAapi.load_gmxmmpbsa_info('_GMXMMPBSA_info')
total = data['gb']['complex']['TOTAL'].copy()
data = gmxMMPBSAapi.load_gmxmmpbsa_info('_GMXMMPBSA_info')
total_mut = data.mutant['gb']['complex']['TOTAL'].copy()
# Create a second copy of the data set. The np.correlate function does not
# normalize the correlation function, so we modify total and total2 to get
# that effect
total -= total.mean()
total /= total.std()
total2 = total.copy() / len(total)
acor = np.correlate(total, total2, 'full')
total_mut -= total_mut.mean()
total_mut /= total_mut.std()
total2_mut = total_mut.copy() / len(total_mut)
acor_mut = np.correlate(total_mut, total2_mut, 'full')
# Now generate the 'lag' axis
xdata = np.arange(0, len(total))
# The acor data set is symmetric about the origin, so only accept the
# positive lag times. Graph the result
plt.plot(xdata, acor[len(acor)//2:], xdata, acor_mut[len(acor)//2:])
plt.show()
```

### Decomposition Data
When performing decomposition analysis, the various decomp data is stored in a separate tree of dicts 
referenced with the `decomp` key. The key sequence is similar to the sequence for the `normal` data
described above, where `decomp` is followed by the solvent model (`GB` or `PB`), followed by the species
(`complex`, `receptor`, or `ligand`) (additionally, we include `delta` key), followed by the decomposition
components (total, backbone, or sidechain), followed by the residue number (and residue pair for pairwise 
decomposition), finally followed by the contribution (`internal`, `van der Waals`, `electrostatics`, etc.)
The available keys are shown in Figure 1 (and each key is described afterwards).

**Decomp Key Descriptions**
* `gb` All Generalized Born results
* `pb` All Poisson-Boltzmann results
* `complex` All results from the complex trajectory
* `receptor` All results from the receptor trajectory
* `ligand` All results from the ligand trajectory
* `delta` All results from delta total decomposition [ complex TDC - (receptor TDC + ligand TDC) ]  
* `TDC` All results from the total decomposition
* `SDC` All results from the sidechain decomposition
* `BDC` All results from the backbone decomposition
* `#` All data from residue number `#` in per-residue and per-wise decomposition (same residue numbering
  scheme as in each respective topology file)
    * `##` All interaction energies between residues `##` and their respective pair `#` in per-wise 
  decomposition (same residue numbering scheme as in each respective topology file)
* `int` Internal energy contributions (see the idecomp variable description above)
* `vdw` van der Waals energy contributions
* `eel` Electrostatic energy contributions
* `pol` Polar solvation free energy contributions
* `sas` Non-polar solvation free energy contributions
* `tot` Total free energy contributions (sum of previous 5).

Figure 1: Tree of `dict` keys following the `decomp` key in a `mmpbsa_data` instance.
![](./doc_files/decomp_dict_keys.png)


