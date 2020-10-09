# GMX-MMPBSA

## Install
### Requirements

## Run GMX-MMPBSA
amber.python ../../GMX_MMPBSA.py -cs COM_ion_em.tpr -ci index.ndx -cg 1 13 -ct traj.xtc -i mmpbsa_igb2.in -lm ligand.mol2

## Test

## Documentation
# GMX-MMPBSA.py
**Note: Here we intend to clarify the new implementations and the differences with the previous one. We do not intend 
to replace the original MMPBSA, we have only extended and made this valuable tool easier for Gromacs users. Most of the 
documentation is found in the Amber manual (AmberTools20), we will point out what is new or different.
Note: You can review some of the answers to the questions that we consider most common here. If you find a bug or have 
any doubts, consider opening an issue**

Neither of these should be considered as a “black-box”, and users should be familiar with Amber before at-
tempting these sorts of calculations. These scripts automate a series of calculations, and cannot trap all the types
of errors that might occur. You should be sure that you know how to carry out an MM-PBSA calculation “by
hand” (i.e., without using the scripts); if you don’t understand in detail what is going on, you will have no good
reason to trust the results. Also, if something goes awry (and this is not all that uncommon), you will need to run
and examine the individual steps to carry out useful debugging.

## Introduction
This section describes the use of the python script GMX-MMPBSA.py [676] to perform Molecular Mechanics / Pois-
son Boltzmann (or Generalized Born) Surface Area (MM/PB(GB)SA) calculations for Gromacs. This is a post-processing
method in which representative snapshots from an ensemble of conformations are used to calculate the free energy
change between two states (typically a bound and free state of a receptor and ligand). Free energy differences are
calculated by combining the so-called gas phase energy contributions that are independent of the chosen solvent
model as well as solvation free energy components (both polar and non-polar) calculated from an implicit solvent
model for each species. Entropy contributions to the total free energy may be added as a further refinement. The
entropy calculations can be done in either a HCT Generalized Born solvation model [185, 196] or in the gas phase
using a mmpbsa_py_nabnmode program written in the nab programming language, or via the quasi-harmonic
approximation in ptraj.
The gas phase free energy contributions are calculated by sander within the Amber program suite or
mmpbsa_py_energy within the AmberTools package according to the force field with which the topology files were
created. The solvation free energy contributions may be further decomposed into an electrostatic and hydrophobic
contribution. The electrostatic portion is calculated using the Poisson Boltzmann (PB) equation, the Generalized
Born method, or the Reference Interaction Site Model (RISM). The PB equation is solved numerically by either the
pbsa program included with AmberTools or by the Adaptive Poisson Boltzmann Solver (APBS) program through
the iAPBS interface[446] with Amber (for more information, see http://www.poissonboltzmann.org/apbs). The
hydrophobic contribution is approximated by the LCPO method [170] implemented within sander or the molsurf
method as implemented in cpptraj.
MM/PB(GB)SA typically employs the approximation that the configurational space explored by the systems are
very similar between the bound and unbound states, so every snapshot for each species is extracted from the same
trajectory file, although GMX-MMPBSA.py will accept separate trajectory files for each species. Furthermore, explicit
solvent and ions are stripped from the trajectory file(s) to hasten convergence by preventing solvent-solvent inter-
actions from dominating the energy terms. A more detailed explanation of the theory can be found in Srinivasan,
et. al.[677] You may also wish to refer to reviews summarizing many of the applications of this model,[678, 679]
as well as to papers describing some of its applications.[680–684]
Many popular types of MM/PBSA calculations can be performed using just AmberTools, while some of the
more advanced functionality requires the sander (now included in AmberTools) program from Amber.

## Preparing for an MM/PB(GB)SA calculation
MM/PB(GB)SA is often a very useful tool for obtaining relative free energies of binding when comparing
ligands. Perhaps its biggest advantage is that it is very computationally inexpensive compared to other free energy
calculations, such as TI or FEP. Following the advice given below before any MD simulations are run will make
running MMPBSA.py successfully much easier.

### Building Topology Files
As converting a topology from Gromacs to Amber is somewhat unpredictable and can be inconsistent, so we decided to 
build the topology from a structure. To build the topology, GMX-MMPBSA requires the following (see options):
* Structure file (tpr)
* Index file (ndx)
* Receptor and ligand group number in the index file

With these requirements we can build the topology and extract the dry trajectory that are necessary to be able to carry
 out the calculation with MMPBSA.

The most common procedure is to obtain the receptor and ligand topologies from the complex, so we have directly 
included a functionality similar to ante-MMPBSA.py (see Amber manual). Additionally we define explicitly the stability 
calculation mode, so GMX-MMPBSA only requires the complex structure(tpr) or if you prefer the complex, receptor and 
ligand structures (tpr and mol2 if ligand is a small molecule). In case of only defining the complex structure like 
Protein-Protein type, we will build the receptor and ligand topologies from it. If the complex is like Protein-ligand 
type (small molecule), the mol2 file used to parameterize the ligand must be defined.

You will always need a topology for the entire complex, but you can define also receptor and ligand. Moreover, they 
must be compatible with one another (i.e., each must have the same charges for the same atoms and the same force field 
must be used for all three of the required structures (tpr)).

Unlike MMPBSA, it can perform an alanine scan without having to provide a topology of the mutant system. Simply define 
the amino acid you want to mutate and we build the topology.

## Running MMPBSA.py
### The input file
The input file was designed to be as syntactically similar to other programs in Amber as possible. The input
file has the same namelist structure as both sander and pmemd. The allowed namelists are &general, &gb, &pb,
&rism, &alanine_scanning, &nmode, and &decomp. The input variables recognized in each namelist are described
below, but those in &general are typically variables that apply to all aspects of the calculation. The &gb namelist
is unique to Generalized Born calculations, &pb is unique to Poisson Boltzmann calculations, &rism is unique to
3D-RISM calculations, &alanine_scanning is unique to alanine scanning calculations, &nmode is unique to the
normal mode calculations used to approximate vibrational entropies, and &decomp is unique to the decomposition
scheme. All of the input variables are described below according to their respective namelists. Integers and floating
point variables should be typed as-is while strings should be put in either single- or double-quotes. All variables
should be set with “variable = value” and separated by commas. See the examples below. Variables will usually
be matched to the minimum number of characters required to uniquely identify that variable within that namelist.
Variables require at least 4 characters to be matched unless that variable name has fewer than 4 characters (in which
case the whole variable name is required). For example, “star” in &general will match “startframe”. However,
“stare” and “sta” will match nothing.


**&general namelist variables**

**debug_printlevel** MMPBSA.py prints errors by raising exceptions, and not catching fatal errors. If debug_printlevel 
            is set to 0, then detailed tracebacks (effectively the call stack showing exactly where in the program the 
            error occurred) is suppressed, so only the error message is printed. If debug_printlevel is set to 1 or 
            higher, all tracebacks are printed, which aids in debugging of issues. Default: 0. (Advanced Option)
            
`endframe` The frame from which to stop extracting snapshots from the full, concatenated trajectory comprised of every 
            trajectory file supplied on the command-line. (Default = 9999999)
```diff
@@ Input variable modified. Included Interaction entropy aproximation @@
```            
`entropy` It specifies whether to make a quasi-harmonic entropy (QH) approximation with ptraj or the interaction 
            entropy (IE) approximation. The allowed values are: 0 Don’t. 1: make QH. 2 make IE (default = 0)
```diff
+ New input variable added
```
`entropy_seg` Specify the representative segment (%), starting from the end of the trajectory, for the calculation of 
            the interaction entropy, ie: 25 is equivalent to the frames contained in the last quartile of the 
            trajectory. Default: 25 (Only if `entropy` = 2)
```diff
+ New input variable added
```
`entropy_temp` Specify the temperature to calculate the entropy (Only if `entropy` = 2). Avoid inconsistencies with 
            defined internal temperature (298.15 K) when nmode is used. 
      
`interval` The offset from which to choose frames from each trajectory file. For example, an interval of 2 will pull
            every 2nd frame beginning at startframe and ending less than or equal to endframe. (Default = 1)

```diff
- Input variable deleted. All files are needed for analysis with charts
``` 
~~-`keep_files` The variable that specifies which temporary files are kept. All temporary files have the prefix 
`_MMPBSA_` prepended to them (unless you change the prefix on the command-line—see subsection Subsection 34.3.2 for 
details). Allowed values are 0, 1, and 2. 0: Keep no temporary files 1: Keep all generated trajectory files and mdout 
files created by sander simulations 2: Keep all temporary files. Temporary files are only deleted if MMPBSA.py 
completes successfully (Default = 1) A verbose level of 1 is sufficient to use -rewrite-output and recreate the output
 file without rerunning any simulations.~~

`ligand_mask` The mask that specifies the ligand residues within the complex prmtop (NOT the solvated prmtop if there 
            is one). The default guess is generally sufficient and will only fail as stated above. You should use the 
            default mask assignment if possible because it provides a good error catch. This follows the same 
            description as the receptor_mask above.
            
`netcdf` Specifies whether or not to use NetCDF trajectories internally rather than writing temporary ASCII trajectory 
            files. NOTE: NetCDF trajectories can be used as input for MMPBSA.py regardless of what this variable is set 
            to, but NetCDF trajectories are faster to write and read. For very large trajectories, this could offer 
            significant speedups, and requires less temporary space. However, this option is incompatible with alanine
            scanning. Default value is 0.<br>
            0: Do NOT use temporary NetCDF trajectories<br>
            1: Use temporary NetCDF trajectories
        
`receptor_mask` The mask that specifies the receptor residues within the complex prmtop (NOT the solvated prm-
top if there is one). The default guess is generally sufficient and will only fail if the ligand residues are not
found in succession within the complex prmtop. You should use the default mask assignment if possible
because it provides a good error catch. It uses the “Amber mask” syntax described elsewhere in this manual.
This will be replaced with the default receptor_mask if ligand_mask (below) is not also set.

```diff
- Input variable deleted. ALways must be defined to get gromacs
```
~~-`search_path` Advanced option. By default, MMPBSA.py will only search for executables in $AMBERHOME/bin .
To enable it to search for binaries in your full PATH if they can’t be found in $AMBERHOME/bin , set
search_path to 1. Default 0 (do not search through the PATH ). This is particularly useful if you are using
an older version of sander that is not in AMBERHOME .~~


`startframe` The frame from which to begin extracting snapshots from the full, concatenated trajectory comprised
of every trajectory file placed on the command-line. This is always the first frame read. (Default = 1)

`strip_mask` The variable that specifies which atoms are stripped from the trajectory file if a solvated_prmtop is
provided on the command-line. See 34.3.2. (Default = “:WAT:CL:CIO:CS:IB:K:LI:MG:NA:RB”)

`use_sander` Forces MMPBSA.py to use sander for energy calculations, even when mmpbsa_py_energy will suf-
fice (Default 0)
0 - Use mmpbsa_py_energy when possible
1 - Always use sander

`full_traj` This variable is for calculations performed in parallel to control whether complete trajectories are made of
the complex, receptor, and ligand. In parallel calculations, a different trajectory is made for each processor to
analyze only the selected frames for that processor. A value of 0 will only create the intermediate trajectories
analyzed by each processor, while a value of 1 will additionally combine those trajectories to make a single
trajectory of all frames analyzed across all processors for the complex, receptor, and ligand. (Default = 0)

`verbose` The variable that specifies how much output is printed in the output file. There are three allowed values:
0, 1, and 2. A value of 0 will simply print difference terms, 1 will print all complex, receptor, and ligand
terms, and 2 will also print bonded terms if one trajectory is used. (Default = 1)
```diff
+ Input variable added. Define Internal dielectric constant without use external mdin file
```
`intdiel` Define a new intenal dielectric constant (Default=1.0)

**&gb namelist variables**

`ifqnt` Specifies whether a part of the system is treated with quantum mechanics. 1: Use QM/MM, 0: Potential
function is strictly classical (Default = 0). This functionality requires sander
igb Generalized Born method to use (seeSection 4). Allowed values are 1, 2, 5, 7 and 8. (Default = 5) All models
are now available with both mmpbsa_py_energy and sander

`qm_residues` Comma- or semicolon-delimited list of complex residues to treat with quantum mechanics. All
whitespace is ignored. All residues treated with quantum mechanics in the complex must be treated with
quantum mechanics in the receptor or ligand to obtain meaningful results. If the default masks are used,
then MMPBSA.py will figure out which residues should be treated with QM in the receptor and ligand.
Otherwise, skeleton mdin files will be created and you will have to manually enter qmmask in the ligand and
receptor topology files. There is no default, this must be specified.

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

the description of the molsurf action command in cpptraj. Default is 0.

**&pb namelist variables**

`inp` Nonpolar optimization method (Default = 2).

`cavity_offset` Offset value used to correct non-polar free energy contribution (Default = -0.5692) This is not used
for APBS.

`cavity_surften` Surface tension. (Default = 0.0378 kcal/mol Angstrom 2 ). Unit conversion to kJ done automati-
cally for APBS.

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
through the iAPBS interface[446] built into sander.APBS. Instructions for this can be found online at the
iAPBS/APBS websites. Allowed values are 0: Don’t use APBS, or 1: Use sander.APBS. (Default = 0)

`memopt` Turn on membrane protein support (Default = 0).

`emem` Membrane dielectric constant (Default = 1.0).

`mthick` Membrane thickness (Default = 40.0).

`mctrdz` Absolute membrane center in the z-direction (Default=0.0, use protein center as the membrane center).

`poretype` Turn on the automatic membrane channel/pore finding method (Default=1).
A more thorough description of these and other options can be found in Chapter 6. Please also note that the default
options have changed over time. For a detailed discussion of all related options on the quality of the MM/PBSA
calculations, please refer to our recent publication [234].

**&alanine_scanning namelist variables**

`mutant_only` Option to perform specified calculations only for the mutants. Allowed values are 0: Do mutant and
original or 1: Do mutant only (Default = 0)
Note that all calculation details are controlled in the other namelists, though for alanine scanning to be performed,
the namelist must be included (blank if desired)
```diff
+Added two option to make more aesy the alanine_scanning performace
```
`mutant` Define if the moutation will be make in receptor or ligand. Allowed values are: receptor, rec, ligand or lig
 in any capitalization (Default = receptor or REC)

`mutant_res` Define the specific residue that do you want to mutate in format CHAIN:RESNUM (eg: 'A:350'). Please, make 
sure that your selection is based on consistent structure file.

**&nmode namelist variables**

`dielc` Distance-dependent dielectric constant (Default = 1.0)

`drms` Convergence criteria for minimized energy gradient. (Default = 0.001)

`maxcyc` Maximum number of minimization cycles to use per snapshot in sander. (Default = 10000)

`nminterval` ∗ Offset from which to choose frames to perform nmode calculations on (Default = 1)

`nmendframe` ∗ Frame number to stop performing nmode calculations on (Default = 1000000)

`nmode_igb` Value for Generalized Born model to be used in calculations. Options are 0: Vacuum, 1: HCT GB
model [185, 196] (Default 1)

`nmode_istrng` Ionic strength to use in nmode calculations. Units are Molarity. Non-zero values are ignored if

`nmode_igb` is 0 above. (Default = 0.0)

`nmstartframe` ∗ Frame number to begin performing nmode calculations on (Default = 1)
* These variables will choose a subset of the frames chosen from the variables in the &general namelist. Thus, the
“trajectory” from which snapshots will be chosen for nmode calculations will be the collection of snapshots upon
which the other calculations were performed.

**&decomp namelist variables**

`csv_format` Print the decomposition output in a Comma-Separated-Variable (CSV) file. CSV files open natively
in most spreadsheets. If set to 1, this variable will cause the data to be written out in a CSV file, and standard
error of the mean will be calculated and included for all data. If set to 0, the standard, ASCII format will be
used for the output file. Default is 1 (CSV-formatted output file)

`dec_verbose` Set the level of output to print in the decmop_output file.
0 - DELTA energy, total contribution only
1 - DELTA energy, total, sidechain, and backbone contributions
2 - Complex, Receptor, Ligand, and DELTA energies, total contribution only
3 - Complex, Receptor, Ligand, and DELTA energies, total, sidechain, and backbone contributions
Note: If the values 0 or 2 are chosen, only the Total contributions are required, so only those will be printed
to the mdout files to cut down on the size of the mdout files and the time required to parse them. However,
this means that -rewrite-output cannot be used to change the default verbosity to print out sidechain and/or
backbone energies, but it can be used to reduce the amount of information printed to the final output. The
parser will extract as much information from the mdout files as it can, but will complain and quit if it cannot
find everything it’s being asked for.
Default = 0

`idecomp` Energy decomposition scheme to use:
1. Per-residue decomp with 1-4 terms added to internal potential terms
2. Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms.
3. Pairwise decomp with 1-4 terms added to internal potential terms
4. Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms

(No default. This must be specified!) This functionality requires sander.

`print_res` Select residues from the complex prmtop to print. The receptor/ligand residues will be automatically
figured out if the default mask assignments are used. If you specify your own masks, you will need to
modify the mdin files created by MMPBSA.py and rerun MMPBSA.py with the -use-mdins flag. Note that
the DELTAs will not be computed in this case. This variable accepts a sequence of individual residues and/or
ranges. The different fields must be either comma- or semicolon-delimited. For example: print_res = “1,
3-10, 15, 100”, or print_res = “1; 3-10; 15; 100”. Both of these will print residues 1, 3 through 10, 15, and
100 from the complex prmtop and the corresponding residues in either the ligand and/or receptor prmtops.
```diff
-(Default: print all residues)*
- *Please note: Using idecomp=3 or 4 (pairwise) with a very large number of printed residues and a large number
-of frames can quickly create very, very large temporary mdout files. Large print selections also demand a large
-amount of memory to parse the mdout files and write decomposition output file (~500 MB for just 250 residues,
-since that’s 62500 pairs!) It is not unusual for the output file to take a significant amount of time to print if you
-have a lot of data. This is most applicable to pairwise decomp, since the amount of data scales as O(N 2 ).

+ Based on the above, we decided to limit the selection of the residues that will be displayed by default. In turn, 
+ we facilitate this selection since now it is not necessary to explicitly define the list of amino acids. We define a 
+ selection method with the following structure:
+ "within 6"
+ where within corresponds to the keyword and 6 to the maximum distance criterion in angstroms necessary to select any
+ residue from  receptor and ligand
```
(Default: print_res="within 6")

**&rism namelist variables***

`buffer` Minimum distance between solute and edge of solvation box. Specify this with grdspc below. Mutually
exclusive with ng and solvbox. Set buffer < 0 if you wish to use ng and solvbox. (Default = 14 Å)
closure The approximation to the closure relation. Allowed choices are kh (Kovalenko-Hirata), hnc (Hypernetted-
chain), or psen (Partial Series Expansion of order-n) where “n” is a positive integer (e.g., “pse3”). (Default
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
MMPBSA.py. See Chapter 7 for a more thorough description of options and theory. A list of references can be
found there, too. One advantage of 3D-RISM is that an arbitrary solvent can be chosen; you just need to change
the xvvfile specified on the command line (see 34.3.2).

#### Sample input files
```
Sample input file for GB and PB calculation
&general
startframe=5, endframe=100, interval=5,
verbose=2
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
verbose=2,
/
&gb
igb=2, saltcon=0.10
/
&alanine_scanning
mutant='receptor'
mutant_res='A:350'
/
--------------------------------------------------------
Sample input file with nmode analysis
&general
startframe=5, endframe=100, interval=5,
verbose=2
/
&gb
igb=5, saltcon=0.150,
/
&nmode
nmstartframe=2, nmendframe=20, nminterval=2,
maxcyc=50000, drms=0.0001,
/
--------------------------------------------------------
Sample input file with decomposition analysis
&general
startframe=5, endframe=100, interval=5,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=3,
print_res="20, 40-80, 200"
/
--------------------------------------------------------
Sample input file for QM/MMGBSA
&general
startframe=5, endframe=100, interval=5,
/
&gb
igb=5, saltcon=0.100, ifqnt=1, qmcharge=0,
qm_residues="100-105, 200", qm_theory="PM3"
/
--------------------------------------------------------
Sample input file for MM/3D-RISM
&general
startframe=5, endframe=100, interval=5,
/
&rism
polardecomp=1, thermo="gf"
/
--------------------------------------------------------
Sample input file for MMPBSA with membrane proteins
&general
use_sander=1,
startframe=1, endframe=100, interval=1, debug_printlevel=2
/
&pb
radiopt=0, indi=20.0, istrng=0.150,
fillratio=1.25, ipb=1, nfocus=1,
bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0,
npbverb=1, solvopt=2, inp=1,
memopt=1, emem=7.0, mctrdz=-10.383, mthick=36.086, poretype=1,
maxarcdot=15000
/
```

A few important notes about input files. Comments are allowed by placing a # at the beginning of the line (whites-
pace is ignored). Variable initialization may span multiple lines. In-line comments (i.e., putting a # for a comment
after a variable is initialized in the same line) is not allowed and will result in an input error. Variable declarations
must be comma-delimited, though all whitespace is ignored. Finally, all lines between namelists are ignored, so
comments may be put before each namelist without using #.

### Calling MMPBSA.py from the command-line
MMPBSA.py is invoked through the command line as follows:
``` diff
Usage: MMPBSA.py [Options]
Options:
usage: GMX_MMPBSA.py [-h] [-v] [--input-file-help] [-O]
                     [-prefix <file prefix>] [-i FILE] [-xvvfile XVVFILE]
                     [-o FILE] [-do FILE] [-eo FILE] [-deo FILE] [-s]
                     [-pff <Forcefield>] [-lff <Forcefield>] [-st]
                     [-cs <Structure File>] [-ci <Index File>]
                     [-cg index index] [-ct [TRJ [TRJ ...]]]
                     [-rs <Structure File>] [-ri <Index File>] [-rg index]
                     [-rt [TRJ [TRJ ...]]] [-lm <Structure File>]
                     [-ls <Structure File>] [-li <Index File>] [-lg index]
                     [-lt [TRJ [TRJ ...]]] [-make-mdins] [-use-mdins]
                     [-rewrite-output] [--clean]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --input-file-help     Print all available options in the input file.

Miscellaneous Options:
  -O, --overwrite       Allow output files to be overwritten
  -prefix <file prefix>
                        Prefix for intermediate files.

Input and Output Files:
  These options specify the input files and optional 
  output files.

  -i FILE               MM/PBSA input file.
  -xvvfile XVVFILE      XVV file for 3D-RISM.
  -o FILE               Output file with MM/PBSA statistics.
  -do FILE              Output file for decomposition statistics summary.
  -eo FILE              CSV-format output of all energy terms for every frame
                        in every calculation. File name forced to end in
                        [.csv]. This file is only written when specified on
                        the command-line.
  -deo FILE             CSV-format output of all energy terms for each printed
                        residue in decomposition calculations. File name
                        forced to end in [.csv]. This file is only written
                        when specified on the command-line.

Options:
  These options specify explicit calculation type and forcefield to prepare the 
  Amber topologies

  -s                    Perform stability calculation. Only the complex
                        parameter are required. Only, if ligand is non-Protein
                        (small molecule) type will required the ligand
                        parameters. IN any other case receptor and ligand
                        parameters will be ignored
  -pff <Forcefield>     Used forcefield to make the protein MD in Gromacs.
                        Allowed: amber14sb, amber99sb-ildn, amber99sb,
                        amber03, amber99, amber96, amber94. (Default:
                        amber14sb)
  -lff <Forcefield>     Used forcefield to make the ligand MD in Gromacs.
                        Allowed: gaff, gaff2. (Default: gaff)
  -st                   Define if complex, receptor and ligand trajectories is
                        solvated. We assume that the path contains ions and
                        water

Complex:
  
  Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
  groups names in the index file and forcefield to make Amber topology. 
  If is necessary we split complex to generate Receptor and Ligand topology like 
  ante-MMPBSA

  -cs <Structure File>  Structure file of the bound complex. If it is Protein-
                        Ligand (small molecule) complex, make sure that you
                        define -lm option
  -ci <Index File>      Index file of the bound complex.
  -cg index index       Groups of the receptor and ligand in complex index
                        file. The notation is as follows: "-cg <Receptor
                        group> <Ligand group>"
  -ct [TRJ [TRJ ...]]   Input trajectories of the complex. Allowed formats:
                        *.xtc (recommended), *.trr, *.pdb (specify as many as
                        you'd like).

Receptor:
  
  Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
  groups names in the index file and forcefield to make Amber topology. 
  If is necessary we split complex to generate Receptor and Ligand topology like 
  ante-MMPBSA

  -rs <Structure File>  Structure file of the unbound receptor. If omitted and
                        stability is False, the structure from complex will be
                        used.
  -ri <Index File>      Index file of the unbound receptor.
  -rg index             Group of the receptor in receptor index file.
  -rt [TRJ [TRJ ...]]   Input trajectories of the unbound receptor. Allowed
                        formats: *.xtc (recommended), *.trr, *.pdb. (specify
                        as many as you'd like).

Ligand:
  
  Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
  groups names in the index file and forcefield to make Amber topology. 
  If is necessary we split complex to generate Receptor and Ligand topology like 
  ante-MMPBSA

  -lm <Structure File>  Structure file of the unbound ligand used to
                        parametrize ligand for Gromacs. Most be defined if
                        Protein-Ligand (small molecule) complex was define.
  -ls <Structure File>  Structure file of the unbound ligand. If omitted, is
                        protein-like and stability is False the structure from
                        complex will be used. If ligand is a small molecule,
                        make sure that you definde above -lm option
  -li <Index File>      Index file of the unbound ligand. Only if tpr file was
                        define in -ls.
  -lg index             Group of the ligand in ligand index file. Only if tpr
                        file was define in -ls.
  -lt [TRJ [TRJ ...]]   Input trajectories of the unbound ligand. Allowed
                        formats: *.xtc (recommended), *.trr, *.pdb (specify as
                        many as you'd like).

Miscellaneous Actions:
  -make-mdins           Create the input files for each calculation and quit.
                        This allows you to modify them and re-run using -use-
                        mdins
  -use-mdins            Use existing input files for each calculation. If they
                        do not exist with the appropriate names, GMX_MMPBSA.py
                        will quit in error.
  -rewrite-output       Do not re-run any calculations, just parse the output
                        files from the previous calculation and rewrite the
                        output files.
  --clean               Clean temporary files and quit.
```

`-make-mdins` and `-use-mdins` are intended to give added flexibility to user input. If the MM/PBSA input file does
not expose a variable you require, you may use the -make-mdins flag to generate the MDIN files and then quit.
Then, edit those MDIN files, changing the variables you need to, then running MMPBSA.py with -use-mdins to
use those modified files.

`--clean` will remove all temporary files created by MMPBSA.py in a previous calculation.

`--version` will display the program version and exit.

### Running MMPBSA.py
#### Serial version
This version is installed with Amber during the serial install of AmberTools. AMBERHOME must be set, or it will
quit on error. If any changes are made to the modules, MMPBSA.py must be remade so the updated modules are
found by MMPBSA.py. An example command-line call is shown below:

`MMPBSA.py -O -i mmpbsa.in -cp com.top -rp rec.top -lp lig.top -y traj.crd`

The tests, found in ${AMBERHOME}/test/mmpbsa_py provide good examples for running MMPBSA.py calcula-
tions.
#### Parallel (MPI) version

This version is installed with Amber during the parallel install. The python package mpi4py is included with
the MMPBSA.py source code and must be successfully installed in order to run the MPI version of MMPBSA.py.
It is run in the same way that the serial version is above, except MPI directions must be given on the command
line as well. Note, if mpi4py does not install correctly, you must install it yourself in order to use
MMPBSA.py.MPI. One note: at a certain level, running RISM in parallel may actually hurt performance, since
previous solutions are used as an initial guess for the next frame, hastening convergence. Running in parallel loses
this advantage. Also, due to the overhead involved in which each thread is required to load every topology file
when calculating energies, parallel scaling will begin to fall off as the number of threads reaches the number of
frames. A usage example is shown below:

`mpirun -np 2 MMPBSA.py.MPI -O -i mmpbsa.in -cp com.top -rp rec.top -lp lig.top -y traj.crd`

### Types of calculations you can do
There are many different options for running MMPBSA.py. Among the types of calculations you can do are:
* Normal binding free energies, with either PB or GB implicit solvent models. Each can be done with either
1, 2, or 3 different trajectories, but the complex, receptor, and ligand topology files must all be defined. The
complex mdcrd must always be provided. Whichever trajectories of the receptor and/or ligand that are NOT
specified will be extracted from the complex trajectory. This allows a 1-, 2-, or 3-trajectory analysis. All PB
calculations and GB models can be performed with just AmberTools via the mmpbsa_py_energy program
installed with MMPBSA.py.
* Stability calculations with any calculation type. If you only specify the complex prmtop (and leave receptor
and ligand prmtop options blank), then a “stability” calculation will be performed, and you will get statistics
based on only a single system. Any additional receptor or ligand information given will be ignored, but note
that if receptor and/or ligand topologies are given, it will no longer be considered a stability calculation. The
previous statement refers principally to mutated receptor/ligand files or extra ligand/receptor trajectory files.
* Alanine scanning with either PB or GB implicit solvent models. All trajectories will be mutated to match
the mutated topology files, and whichever calculations that would be carried out for the normal systems are
also carried out for the mutated systems. Note that only 1 mutation is allowed per simulation, and it must
be to an alanine. If mutant_only is not set to 1, differences resulting from the mutations are calculated. This
option is incompatible with intermediate NetCDF trajectories (see the netcdf = 1 option above). This has the
same program requirements as option 1 above.
* Entropy corrections. An entropy term can be added to the free energies calculated above using either the
quasi-harmonic approximation or the normal mode approximation. Calculations will be done for the nor-
mal and mutated systems (alanine scanning) as requested. Normal mode calculations are done with the
mmpbsa_py_nabnmode program included with AmberTools.
* Decomposition schemes. The energy terms will be decomposed according to the decomposition scheme
outlined in the idecomp variable description. This should work with all of the above, though entropy terms
cannot be decomposed. APBS energies cannot be decomposed, either. Neither can PBSA surface area terms.
This functionality requires sander from the Amber 11 (or later) package.
* QM/MMGBSA. This is a binding free energy (or stability calculation) using the Generalized Born solvent
model allowing you to treat part of your system with a quantum mechanical Hamiltonian. See “Advanced
Options” for tips about optimizing this option. This functionality requires sander from the Amber package.
* MM/3D-RISM. This is a binding free energy (or stability calculation) using the 3D-RISM solvation model.
This functionality is performed with rism3d.snglpnt built with AmberTools.
* Membrane Protein MMPBSA. Calculate the MMPBSA binding free energy for a ligand bound to a protein
that is embedded into a membrane. Only use_sander=1 is supported.

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
MMPBSA.py creates working files during the execution of the script beginning with the prefix _MMPBSA_.
The variable “keep_files” controls how many of these files are kept after the script finishes successfully. If the
script quits in error, all files will be kept. You can clean all temporary files from a directory by running MMPBSA
–clean described above.
If MMPBSA.py does not finish successfully, several of these files may be helpful in diagnosing the problem.
For that reason, every temporary file is described below. Note that not every temporary file is generated in every simulation. At the end of each description, the lowest value of “keep_files” that will retain this file will be shown
in parentheses.

`_MMPBSA_gb.mdin` Input file that controls the GB calculation done in sander. (2)

`_MMPBSA_pb.mdin` Input file that controls the PB calculation done in sander. (2)

`_MMPBSA_gb_decomp_com.mdin` Input file that controls the GB decomp calculation for the complex done in
sander. (2)

`_MMPBSA_gb_decomp_rec.mdin` Input file that controls the GB decomp calculation for the receptor done in
sander. (2)

`_MMPBSA_gb_decomp_lig.mdin` Input file that controls the GB decomp calculation for the ligand done in sander.
(2)

`_MMPBSA_pb_decomp_com.mdin` Input file that controls the PB decomp calculation for the complex done in
sander. (2)

`_MMPBSA_pb_decomp_rec.mdin` Input file that controls the PB decomp calculation for the receptor done in
sander. (2)

`_MMPBSA_pb_decomp_lig.mdin` Input file that controls the PB decomp calculation for the ligand done in sander.
(2)

`_MMPBSA_gb_qmmm_com.mdin` Input file that controls the GB QM/MM calculation for the complex done in sander.
(2)

`_MMPBSA_gb_qmmm_rec.mdin` Input file that controls the GB QM/MM calculation for the receptor done in sander.
(2)

`_MMPBSA_gb_qmmm_lig.mdin` Input file that controls the GB QM/MM calculation for the ligand done in sander.
(2)

`_MMPBSA_complex.mdcrd.#` Trajectory file(s) that contains only those complex snapshots that will be processed
by MMPBSA.py. (1)

`_MMPBSA_ligand.mdcrd.#` Trajectory file(s) that contains only those ligand snapshots that will be processed by
MMPBSA.py. (1)

`_MMPBSA_receptor.mdcrd.#` Trajectory file(s) that contains only those receptor snapshots that will be processed
by MMPBSA.py. (1)

`_MMPBSA_complex_nc.#` Same as _MMPBSA_complex.mdcrd.#, except in the NetCDF format. (1)

`_MMPBSA_receptor_nc.#` Same as _MMPBSA_receptor.mdcrd.#, except in the NetCDF format. (1)

`_MMPBSA_ligand_nc.#` Same as _MMPBSA_ligand.mdcrd.#, except in the NetCDF format. (1)

`_MMPBSA_dummycomplex.inpcrd` Dummy inpcrd file generated by _MMPBSA_complexinpcrd.in for use with
imin=5 functionality in sander. (1)

`_MMPBSA_dummyreceptor.inpcrd` Same as above, but for the receptor. (1)

`_MMPBSA_dummyligand.inpcrd` Same as above, but for the ligand. (1)

`_MMPBSA_complex.pdb` Dummy PDB file of the complex required to set molecule up in nab programs

`_MMPBSA_receptor.pdb` Dummy PDB file of the receptor required to set molecule up in nab programs

`_MMPBSA_ligand.pdb` Dummy PDB file of the ligand required to set molecule up in nab programs

`_MMPBSA_complex_nm.mdcrd.#` Trajectory file(s) for each thread with snapshots used for normal mode calcula-
tions on the complex. (1)

`_MMPBSA_receptor_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calcula-
tions on the receptor. (1)

`_MMPBSA_ligand_nm.mdcrd.#` Trajectory file for each thread with snapshots used for normal mode calculations
on the ligand. (1)

`_MMPBSA_ptrajentropy.in` Input file that calculates the entropy via the quasi-harmonic approximation. This
file is processed by ptraj. (2)

`_MMPBSA_avgcomplex.pdb` PDB file containing the average positions of all complex conformations processed by

`_MMPBSA_cenptraj.in.` It is used as the reference for the _MMPBSA_ptrajentropy.in file above.
(1)

`_MMPBSA_complex_entropy.out` File into which the entropy results from _MMPBSA_ptrajentropy.in analysis
on the complex are dumped. (1)

`_MMPBSA_receptor_entropy.out` Same as above, but for the receptor. (1)

`_MMPBSA_ligand_entropy.out` Same as above, but for the ligand. (1)

`_MMPBSA_ptraj_entropy.out` Output from running ptraj using _MMPBSA_ptrajentropy.in. (1)

`_MMPBSA_complex_gb.mdout.#` sander output file containing energy components of all complex snapshots done
in GB. (1)

`_MMPBSA_receptor_gb.mdout.#` sander output file containing energy components of all receptor snapshots done
in GB. (1)

`_MMPBSA_ligand_gb.mdout.#` sander output file containing energy components of all ligand snapshots done in
GB. (1)

`_MMPBSA_complex_pb.mdout.#` sander output file containing energy components of all complex snapshots done
in PB. (1)

`_MMPBSA_receptor_pb.mdout.#` sander output file containing energy components of all receptor snapshots done
in PB. (1)

`_MMPBSA_ligand_pb.mdout.#` sander output file containing energy components of all ligand snapshots done in
PB. (1)

`_MMPBSA_complex_rism.out.#` rism3d.snglpnt output file containing energy components of all complex snap-
shots done with 3D-RISM (1)

`_MMPBSA_receptor_rism.out.#` rism3d.snglpnt output file containing energy components of all receptor snap-
shots done with 3D-RISM (1)

`_MMPBSA_ligand_rism.out.#` rism3d.snglpnt output file containing energy components of all ligand snapshots
done with 3D-RISM (1)

`_MMPBSA_pbsanderoutput.junk.#` File containing the information dumped by sander.APBS to STDOUT. (1)

`_MMPBSA_ligand_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the ligand
for all snapshots. (1)

`_MMPBSA_receptor_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the
receptor for all snapshots. (1)

`_MMPBSA_complex_nm.out.#` Output file from mmpbsa_py_nabnmode that contains the entropy data for the com-
plex for all snapshots. (1)

`_MMPBSA_mutant_...` These files are analogs of the files that only start with _MMPBSA_ described above, but
instead refer to the mutant system of alanine scanning calculations.

`_MMPBSA_*out.#` These files are thread-specific files. For serial simulations, only #=0 files are created. For
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
created are _MMPBSA_gb.mdin which controls GB calculation; _MMPBSA_pb.mdin which controls the
PB calculation; _MMPBSA_sander_nm_min.mdin which controls the sander minimization of snapshots to
be prepared for nmode calculations; and _MMPBSA_nmode.in which controls the nmode calculation. If
no input file is specified, all files above are created with default values, and _MMPBSA_pb.mdin is created
for AmberTools’s pbsa. If you wish to create a file for sander.APBS, you must include an input file with
“sander_apbs=1” specified to generate the desired input file. Note that if an input file is specified, only those
mdin files pertinent to the calculation described therein will be created!

`-use-mdins` This flag will prevent MMPBSA.py from creating the input files that control the vari-
ous calculations (_MMPBSA_gb.mdin, _MMPBSA_pb.mdin, _MMPBSA_sander_nm_min.mdin, and
_MMPBSA_nmode.in). It will instead attempt to use existing input files (though they must have those
names above!) in their place. In this way, the user has full control over the calculations performed, however
care must be taken. The mdin files created by MMPBSA.py have been tested and are (generally) known to
be consistent. Modifying certain variables (such as imin=5) may prevent the script from working, so this
should only be done with care. It is recommended that users start with the existing mdin files (generated by
the -make-mdins flag above), and add and/or modify parameters from there.

`strip_mask` This input variable allows users to control which atoms are stripped from the trajectory files associated
with solvated_prmtop. In general, counterions and water molecules are stripped, and the complex is centered
and imaged (so that if iwrap caused the ligand to “jump” to the other side of the periodic box, it is replaced
inside the active site). If there is a specific metal ion that you wish to include in the calculation, you can
prevent ptraj from stripping this ion by NOT specifying it in strip_mask. Note that strip_mask does nothing
if no solvated_prmtop is provided.

`QM/MMGBSA` There are a lot of options for QM/MM calculations in sander, but not all of those options were
made available via options in the MMPBSA.py input file. In order to take advantage of these other options,
you’ll have to make use of the -make-mdins and -use-mdins flags as detailed above and change the resulting
_MMPBSA_gb_qmmm_com/rec/lig.mdin files to fit your desired calculation. Additionally, MMPBSA.py
suffers all shortcomings of sander, one of those being that PB and QM/MM are incompatible. Therefore,
only QM/MMGBSA is a valid option right now.