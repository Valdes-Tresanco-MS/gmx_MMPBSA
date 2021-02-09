---
template: main.html
title: The input file
---

# The input file

## Description

As `gmx_MMPBSA` is based on [MMPBSA.py][1], it uses an input file containing all the specification for the MM/PB(GB)
SA calculation. The input file is designed to be as syntactically similar to other programs in Amber as possible. 
The input file has the same namelist structure as both sander and pmemd. The allowed namelists are `&general`, `&gb`, 
`&pb`, `&rism`, `&alanine_scanning`, `&nmode`, and `&decomp`. The input variables recognized in each namelist are 
described below, but those in `&general` are typically variables that apply to all aspects of the calculation or 
parameters required for build amber topologies from GROMACS files. The &gb namelist is unique to Generalized Born 
calculations, `&pb` is unique to Poisson Boltzmann calculations, `&rism` is unique to 3D-RISM calculations, 
`&alanine_scanning` is unique to alanine scanning calculations, `&nmode` is unique to the normal mode calculations 
used to approximate vibrational entropies, and `&decomp` is unique to the decomposition scheme. All of the input 
variables are described below according to their respective namelists. Integers and floating point variables should 
be typed as-is while strings should be put in either single- or double-quotes. All variables should be set with 
`variable = value` and separated by commas. See several [examples][2] below. As you will see, several calculations 
can be performed in the same run (_i.e._ `&gb` and `&pb`, `&gb` and `&alanine_scanning`, `&pb` and `&decomp`, etc). 
Variables will usually be matched to the minimum number of characters required to uniquely identify that variable 
within that namelist. Variables require at least 4 characters to be matched unless that variable name has fewer than 4 
characters (in which case the whole variable name is required). For example, "star" in &general will match 
`startframe`. However, "stare" and "sta" will match nothing.

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [2]: #sample-input-files

### **`&general` namelist variables**

`assign_chainID` Defines the chains ID assignment mode. _It is ignored when defining a reference structure
(recommended)_. If `assign_chainID = 1`, gmx_MMPBSA check if the structure has no chains ID and it is assigned according
to the structure[^1]. If `assign_chainID = 2`, `gmx_MMPBSA` re-assign the chains ID, exist or not, according to the
structure[^1] (can generate inconsistencies). If a `*.gro` file was used for complex structure
(`-cs` flag) and not reference structure was provided, `gmx_MMPBSA` assume `assign_chainID = 1`. (Default = 0)

  [^1]: _The chain ID is assigned according to two criteria: **terminal amino acids** and **residue numbering**. If
        both criteria or residue numbering changes are present, we assign a new chain ID. If there are terminal 
        amino acids but the numbering of the residue continues, we do not change the ID of the chain._


`debug_printlevel` MMPBSA.py prints errors by raising exceptions, and not catching fatal errors. If debug_printlevel is
set to 0, then detailed tracebacks (effectively the call stack showing exactly where in the program the error occurred)
is suppressed, so only the error message is printed. If debug_printlevel is set to 1 or higher, all tracebacks are
printed, which aids in debugging of issues. (Default = 0) (Advanced Option)

!!! tip
    _Now `gmx_MMPBSA` shows the command-line used to build AMBER topologies when `debug_printlevel = 1 or higher`._ 

`startframe` (Default = 1)
:   The frame from which to begin extracting snapshots from the full, concatenated trajectory comprised of
    every trajectory file placed on the command-line. This is always the first frame read.

`endframe` (Default = 9999999)
:   The frame from which to stop extracting snapshots from the full, concatenated trajectory comprised of every
    trajectory file supplied on the command-line. 

`entropy` (default = 0) 
:    It specifies whether to perform a quasi-harmonic entropy (QH) approximation with ptraj or the
     [Interaction Entropy (IE)][3] approximation. The allowed values are:
     
     * 0: Don’t
     * 1: perform QH
     * 2: perform IE

  [3]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682

!!! tip "Improvement" 
    _Included Interaction entropy aproximation._
 

`entropy_seg` (Default = 25)
:    Specify the representative segment (in %), starting from the `endframe`, for the calculation of the
     Interaction Entropy, _e.g._: `entropy_seg = 25` means that the last quartile of the total number of frames
     (`(endframe-startframe)/interval`) will be used to calculate the average Interaction Entropy. (Only
      if `entropy = 2`)

`entropy_temp` (Default = 298.15)
:    Specify the temperature to calculate the entropy term `−TΔS` (Only if `entropy` = 2). Avoid inconsistencies 
      with defined internal temperature (298.15 K) when nmode is used.

`gmx_path` 
:     Define an additional path to search for GROMACS executables. This path takes precedence over the path defined
      in the PATH variable. In these path the following executables will be searched: `gmx`, `gmx_mpi`, `gmx_d`,
      `gmx_mpi_d` (Gromcas > 5.x.x), `make_ndx`, `editconf` and `trjconv` (GROMACS 4.x.x)

`interval` (Default = 1)
:     The offset from which to choose frames from each trajectory file. For example, an interval of 2 will pull
      every 2nd frame beginning at startframe and ending less than or equal to endframe. 

`netcdf` (Default = 0)
:     Specifies whether or not to use NetCDF trajectories internally rather than writing temporary ASCII trajectory
      files. For very large trajectories, this could offer significant speedups, and requires less temporary space. 
      However, this option is incompatible with alanine scanning.

      * 0: Do NOT use temporary NetCDF trajectories
      * 1: Use temporary NetCDF trajectories

`overwrite_data` (Default = 0)
:   Defines whether the gmxMMPBSA data will be overwritten

    * 0: don't
    * 1: overwrite gxmMMPBSA data if exist

    !!! tip Keep in mind
        We recommend activating this option with each new release because there may be changes and/or corrections in 
        the gmxMMPBSA data files.

`PBRadii` (Default = 3)
:   PBRadii to build amber topology files:

    * 1: bondi, recommended when igb = 7
    * 2: mbondi, recommended when igb = 1
    * 3: mbondi2, recommended when igb = 2 or 5
    * 4: mbondi3, recommended when igb = 8

`protein_forcefield` (Default = "oldff/leaprc.ff99SB")
:   Define the force field used to build Amber topology for proteins. Make sure this force field is the same as the 
    one used in GROMACS. Force fields tested:

    * "oldff/leaprc.ff99"
    * "oldff/leaprc.ff03"
    * "oldff/leaprc.ff99SB"
    * "oldff/leaprc.ff99SBildn"
    * "leaprc.protein.ff14SB"

    !!! important "Keep in mind"
        * You don't need to define it when you use a topology. Please refer to the section "How gmx_MMPBSA works"
        * _This notation format is the one used in tleap._


`ligand_forcefield` (Default = "leaprc.gaff")
:   Define the force field used to build Amber topology for small molecules or glycams. Make sure this force field 
    is the same as the one used for GROMACS . Force fields tested:

    * "leaprc.gaff"
    * "leaprc.gaff2"
    * "leaprc.GLYCAM_06j-1"    (Compatible with amber12SB and later)
    * "leaprc.GLYCAM_06EPb"    (Compatible with amber12SB and later)
    * "gmxMMPBSA/leaprc.GLYCAM_06h-1"    `*`(Included in gmx_MMPBSA package. Compatible with amber99SB and earlier)
    * "gmxMMPBSA/leaprc.zaa99SB"    `*`Parameters for Zwitterionic amino acids. (Included in gmx_MMPBSA package. 
                                       Compatible with amber 99SB)
  
    !!! important "Keep in mind"
        * You don't need to define it when you use a topology or the ligand is protein-like type. Please refer to the 
        section "How gmx_MMPBSA works"
        * _This notation format is the one used in tleap._

        `*` We create a new folder (named _gmxMMPBSA_) in each one of the Amber's parameter folders 
            ($AMBERHOME/dat/leap/[cmd, prep, lib, parm]/gmxMMPBSA). This way, we keep `gmx_MMPBSA` data separated from 
            Amber's.


`ions_parameters` (Default = 1)
:   Define ions parameters to build the Amber topology. 

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

    !!! important "Keep in mind"
        * You don't need to define it when you use a topology. Please refer to the section "How gmx_MMPBSA works"    
        * This notation is simpler since these parameter files are generally the same for all systems        


`reuse_files` (Default = 0)
:   Define whether the trajectories files will be reused when the program ends in error. 

    * 0: Don't reuse. If there are temporary trajectory files, they will be deleted
    * 1: Reuse existing trajectory file

    !!! danger
        Note that the trajectories files may not be generated correctly due to internal errors or interruptions. 
        Please use it with care.

`solvated_trajectory` (Default = 1)
:   Define if it is necessary to build a clean trajectory with no water and ions
    
    * 0: Don’t
    * 1: Build clean trajectory


`use_sander` (Default = 0)
:   use sander for energy calculations, even when `mmpbsa_py_energy` will suffice 
    
    * 0: Use `mmpbsa_py_energy` when possible
    * 1: Always use sander

    !!! note
        Sander is always used when building the Amber topology from a Gromacs topology. This, because the conversion 
        can generate parameters that are not recognized by `mmpbsa_py_energy` 

`verbose` (Default = 1)
:   The variable that specifies how much output is printed in the output file. There are three allowed values:

    * 0: print difference terms
    * 1: print all complex, receptor, and ligand terms
    * 2: also print bonded terms if one trajectory is used

### **`&gb` namelist variables**

`ifqnt` (Default = 0)
:   Specifies whether a part of the system is treated with quantum mechanics. 
    
    * 1: Use QM/MM
    * 0: Potential function is strictly classical 

    !!! note
        This functionality requires `sander`

`igb` (Default = 5)
:   Generalized Born method to use (seeSection 4).  Allowed values are 1, 2, 5, 7 and 8

    !!! note
        All models are now available with both `mmpbsa_py_energy` and `sander`

`qm_residues`
:   Comma- or semicolon-delimited list of complex residues to treat with quantum mechanics. All whitespace is 
    ignored. All residues treated with quantum mechanics in the complex must be treated with quantum mechanics in the
    receptor or ligand to obtain meaningful results.

`intdiel` (Default = 1.0)
:   Define Internal dielectric constant without use external *.mdin file

`qm_theory` 
:   Which semi-empirical Hamiltonian should be used for the quantum calculation. See its description in the QM/MM 
    section of the manual for options.

    !!! danger
         No default, this must be specified.

`qmcharge_com` (Default = 0)
:   The charge of the quantum section for the complex. 

`qmcharge_lig` (Default = 0)
:   The charge of the quantum section of the ligand.

`qmcharge_rec` (Default = 0)
:   The charge of the quantum section for the receptor.

`qmcut` (Default = 9999.0)
:   The cutoff for the qm/mm charge interactions. 

`saltcon` (Default = 0.0)
:   Salt concentration in Molarity. 

`surfoff` (Default 0.0)
:   Offset to correct (by addition) the value of the non-polar contribution to the solvation free energy term

`surften` (Default = 0.0072)
:   Surface tension value. Units in kcal/mol/Å^2^

`molsurf` (Default 0)
:   Define the algorithm to calculate the surface area for the nonpolar solvation termWhen 
    
    * 0: LCPO (Linear Combination of Pairwise Overlaps)
    * 1: molsurf algorithm

`probe` (Default = 1.4)
:   Radius of the probe molecule (supposed to be the size of a solvent molecule), in Angstroms, to use when
    determining the molecular surface. 
    
    !!! note
        only applicable when `molsurf` is set to 1

`msoffset` (Default = 0) 
:   Offset to apply to the individual atomic radii in the system when calculating the `molsurf` surface. See the
    description of the `molsurf` action command in [cpptraj][4]. 
    
  [4]: https://ambermd.org/doc12/Amber20.pdf#chapter.32

### **`&pb` namelist variables**

`inp` (Default = 2) 
:   Option to select different methods to compute non-polar solvation free energy.

    * 0: No non-polar solvation free energy is computed
    * 1: The total non-polar solvation free energy is modeled as a single term linearly proportional to the solvent
      accessible surface area, as in the PARSE parameter set, that is, if INP = 1, USE_SAV must be equal to 0.
    * 2: The total non-polar solvation free energy is modeled as two terms: the cavity term and the dispersion term. The
      dispersion term is computed with a surface-based integration method closely related to the PCM solvent for quantum
      chemical programs. Under this framework, the cavity term is still computed as a term linearly proportional to the
      molecular solvent-accessible-surface area (SASA) or the molecular volume enclosed by SASA.

`cavity_offset` (Default = -0.5692)
:   Offset value used to correct non-polar free energy contribution. 

    !!! note
        This is not used for `APBS`.

`cavity_surften` (Default = 0.0378  [ kcal/mol Å^2^ ] )
:   Surface tension. Unit conversion to kJ done automatically for `APBS`.

`exdi` (Default = 80.0)
:   External dielectric constant.

`indi` (Default = 1.0)
:   Internal dielectric constant.

`fillratio` (Default = 4.0)
:   The ratio between the longest dimension of the rectangular finite-difference grid and that of the solute.

`scale` (Default = 2.0)
:   Resolution of the Poisson Boltzmann grid. It is equal to the reciprocal of the grid spacing. 

`istrng` (Default = 0.0)
:   Ionic strength in Molarity. It is converted to mM for `PBSA` and kept as M for `APBS`.

`linit` (Default = 1000) 
:   Maximum number of iterations of the linear Poisson Boltzmann equation to try 

`prbrad` (Default = 1.4)
:   Solvent probe radius in Angstroms. Allowed values are 1.4 and 1.6 

`radiopt` (Default = 1)
:   The option to set up atomic radii according to: 

    * 0: the prmtop, or 
    * 1: pre-computed values 
    
    !!! warning
        * radiopt=0 is recommended which means using radii from the prmtop file for both the PB calculation and for 
        the NP. Check this [thread](http://archive.ambermd.org/201303/0548.html) and See Amber manual for more 
        complete description.

`sander_apbs` (Default = 0)
:   Option to use `APBS` for `PB` calculation instead of the built-in `PBSA` solver. This will work only through the
    `iAPBS` interface built into `sander.APBS`. Instructions for this can be found online at the iAPBS/APBS websites.
    
    * 0: Don’t use `APBS`
    * 1: Use `sander.APBS` 

`memopt` (Default = 0)
:   Turn on membrane protein support.

`emem` (Default = 1.0)
:   Membrane dielectric constant.

`mthick` (Default = 40.0)
:   Membrane thickness.

`mctrdz` (Default=0.0, use protein center as the membrane center)
:   Absolute membrane center in the z-direction.

`poretype` (Default=1)
:   Turn on the automatic membrane channel/pore finding method.

!!! note
    A more thorough description of these and other options can be found [here][5]. Please also note that the default 
    options have changed over time. For a detailed discussion of all related options on the quality of the 
    MM/PB(GB)SA calculations, please check this [publication][6].

  [5]: https://ambermd.org/doc12/Amber20.pdf#chapter.6
  [6]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.24467

### **`&alanine_scanning` namelist variables**

`mutant_only`  (Default = 0)
:   Option to perform specified calculations only for the mutants. 

    * 0: Do mutant and original
    * 1: Do mutant only
    
    !!! note
        Note that all calculation details are controlled in the other namelists, though for alanine scanning to be 
        performed, the namelist must be included (blank if desired)

`mutant` (Default = "ALA") 
:   Defines the residue by which it is going to mutate.
    
    !!! note
        Allowed values are: "ALA" or "A" for Alanine scanning and "GLY" or "G" for Glycine scanning.

`mutant_res`
:   Define the specific residue that is going to be mutated. Use the following format CHAIN:RESNUM (eg: 'A: 350'). 

    !!! important
        * We recommend using the reference structure (-cr) to ensure the perfect match between the selected residue in 
        the defined structure or topology 
        * This option allow `gmx_MMPBSA` to do the mutation. This way the user does not have to provide the mutant 
        topology

### **`&nmode` namelist variables**

`dielc` (Default = 1.0)
:   Distance-dependent dielectric constant 

`drms` (Default = 0.001)
:   Convergence criteria for minimized energy gradient.

`maxcyc` (Default = 10000)
:   Maximum number of minimization cycles to use per snapshot in sander.

`nminterval`[^2] (Default = 1)
:   Offset from which to choose frames to perform `nmode` calculations on

`nmendframe`[^2] (Default = 1000000)
:   Frame number to stop performing `nmode` calculations on 

`nmode_igb` (Default = 1)
:   Value for Generalized Born model to be used in calculations. Options are:
    
    * 0: Vacuum
    * 1: HCT GB model 

`nmode_istrng` (Default = 0.0)
:   Ionic strength to use in `nmode` calculations. Units are Molarity. Non-zero values are ignored if `nmode_igb`
    is 0 above. 

`nmstartframe`[^2]
:   Frame number to begin performing `nmode` calculations on 

  [^2]: _These variables will choose a subset of the frames chosen from the variables in the `&general` namelist. Thus,
        the "trajectory" from which snapshots will be chosen for `nmode` calculations will be the collection of 
        snapshots upon which the other calculations were performed._

### **`&decomp` namelist variables**

`csv_format`  (Default = 1 [CSV-formatted output file])
:   Print the decomposition output in a Comma-Separated-Variable (CSV) file. CSV files open natively in most
spreadsheets. 
    * If set to 1, this variable will cause the data to be written out in a CSV file, and standard error of the mean 
    will be calculated and included for all data. 
    * If set to 0, the standard, ASCII format will be used for the output file.

`dec_verbose` (Default = 0)
:   Set the level of output to print in the decomp_output file.

    * 0: DELTA energy, total contribution only
    * 1: DELTA energy, total, sidechain, and backbone contributions
    * 2: Complex, Receptor, Ligand, and DELTA energies, total contribution only
    * 3: Complex, Receptor, Ligand, and DELTA energies, total, sidechain, and backbone contributions

    !!! note
        If the values 0 or 2 are chosen, only the Total contributions are required, so only those will be printed to the
        mdout files to cut down on the size of the mdout files and the time required to parse them.


`idecomp`
:   Energy decomposition scheme to use:
    
    * 1: Per-residue decomp with 1-4 terms added to internal potential terms
    * 2: Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms.
    * 3: Pairwise decomp with 1-4 terms added to internal potential terms
    * 4: Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms

    !!! warning
        * No default. This must be specified!.
        * This functionality requires sander.


`print_res` (Default = "within 6")
:   Select residues from the complex to print. This variable also accepts a sequence of individual residues and/or 
    ranges. The different fields must be either comma- or semicolon-delimited. For example: print_res = "within 6", 
    where _within_ corresponds to the keyword and _6_ to the maximum distance criterion in Angstroms necessary to 
    select the residues from both the receptor and the ligand; or print_res = "1,3-10, 15, 100", or 
    print_res = "1; 3-10; 15; 100". Both of these will print residues 1, 3 through 10, 15, and 100 from the complex 
    topology file and the corresponding residues in either the ligand and/or receptor topology files.

    ???+ warning
        Using idecomp=3 or 4 (pairwise) with a very large number of printed residues and a large number of frames 
        can quickly create very, very large temporary mdout files. Large print selections also demand a large amount 
        of memory to parse the mdout files and write decomposition output file (~500 MB for just 250 residues, since 
        that’s 62500 pairs!) It is not unusual for the output file to take a significant amount of time to print if 
        you have a lot of data. This is most applicable to pairwise decomp, since the amount of data scales as O(N 2 ).

    ???+ check
        Based on the above, we decided to add a new option that limits the selection of the residues that will be 
        printed by default. The new option (within 6) considerably reduces the number of residues to be printed 
        reducing the computational cost. Additionally, it selects the interaction interface residues (according to 
        the selected value) automatically, without the user needing to define them explicitly.


### **`&rism` namelist variables**

???+ warning
    `3D-RISM` calculations are performed with the `rism3d.snglpnt` program built with AmberTools, written by Tyler 
    Luchko. It is the most expensive, yet most statistical mechanically rigorous solvation model available in 
    `MMPBSA.py`. See [Chapter 7][7] for a more thorough description of options and theory. 
    A list of references can be found there, too. One advantage of `3D-RISM` is that an arbitrary solvent can be chosen; 
    you just need to change the `xvvfile` specified on the command line (see [34.3.2][8]).

  [7]: https://ambermd.org/doc12/Amber20.pdf#chapter.7
  [8]: https://ambermd.org/doc12/Amber20.pdf#subsection.34.3.2

`buffer` (Default = 14 Å)
:   Minimum distance between solute and edge of solvation box. Specify this with `grdspc` below. Mutually exclusive
    with ng and solvbox. Set buffer < 0 if you wish to use `ng` and `solvbox`. 

`closure` (Default = "kh")
:   The approximation to the closure relation. Allowed choices are `kh` (Kovalenko-Hirata), `hnc` (Hypernetted- chain),
    or `psen` (Partial Series Expansion of order-n) where "n" is a positive integer (_e.g._, "pse3").

???+ warning "Deprecated"
    `closureorder` (Default = 1)
    :    The order at which the PSE-n closure is truncated if closure is specified as "pse" or "psen" (no integers).

`grdspc`(Default = 0.5 Å)
:   Grid spacing of the solvation box. Specify this with buffer above. Mutually exclusive with `ng` and `solvbox`.


`ng` 
:   Number of grid points to use in the x, y, and z directions. Used only if buffer < 0. Mutually exclusive with `buffer`
    and `grdspc` above, and paired with `solvbox` below. 

    !!! warning 
        No default, this must be set if buffer < 0. Define like "ng=1000,1000,1000"

`polardecomp` (Default = 0)
:   Decompose the solvation free energy into polar and non-polar contributions. Note that this will increase 
    computation time by roughly 80%. 
    
    * 0: Don’t decompose solvation free energy. 
    * 1: Decompose solvation free energy. 

`rism_verbose` (Default = 0)
:   Level of output in temporary RISM output files. May be helpful for debugging or following convergence. 

    * 0: just print the final result
    * 1: additionally prints the total number of iterations for each solution
    * 2: additionally prints the residual for each iteration and details of the MDIIS solver


`solvbox`
:   Length of the solvation box in the x, y, and z dimensions. Used only if buffer < 0. Mutually exclusive with
    `buffer` and `grdspc` above, and paired with `ng` above. 

    !!! warining 
        No default, this must be set if buffer < 0. Define like "solvbox=20,20,20"

`solvcut`  (Default = buffer )
:   Cutoff used for solute-solvent interactions. The default is the value of buffer. Therefore, if you set `buffer` < 
    0 and specify `ng` and `solvbox` instead, you must set `solvcut` to a nonzero value, or the program will quit in 
    error.

`thermo` (Default = "std")
:   Which thermodynamic equation you want to use to calculate solvation properties. Options are "std", "gf", or 
    "both" (case-INsensitive). "std" uses the standard closure relation, "gf" uses the Gaussian Fluctuation 
    approximation,and "both" will print out separate sections for both. . 
    
    !!! note
        Note that all data are printed out for each RISM simulation, so no choice is any more computationally demanding 
        than another. Also, you can change this option and use the -rewrite-output flag to obtain a different 
        printout after-the-fact.

`tolerance` (Default = 1e-5)
:   Upper bound of the precision requirement used to determine convergence of the self-consistent solution. This has 
    a strong effect on the cost of 3D-RISM calculations. 

## Sample input files

!!! tip
    You can refer to the [examples](../Examples/) to understand the input file in a practical way.

### GB and PB

``` linenums="1"
Sample input file for GB and PB calculation building the Amber topologies
from structures. Please refer to the section "How gmx_MMPBSA works"

&general
startframe=5, endframe=100, interval=5, verbose=2, 
protein_forcefield="oldff/leaprc.ff99SB", ligand_forcefield="leaprc.gaff"
/

&gb
igb=5, saltcon=0.150,
/

&pb
istrng=0.15, fillratio=4.0
/
```

### Alanine scanning

``` linenums="1"
Sample input file for Alanine scanning

&general
startframe=5, endframe=21, verbose=2, interval=1,
protein_forcefield="oldff/leaprc.ff99SB", PBRadii=4
/

&gb
igb=8, saltcon=0.150, intdiel=10
/

&alanine_scanning
mutant='ALA'
mutant_res='B:12'
/
```

### Entropy

``` linenums="1"
Sample input file for entropy calculations

&general
startframe=5, endframe=21, verbose=2, interval=1,
# `entropy` variable control whether to perform a quasi-harmonic entropy (QH)
# approximation or the Interaction Entropy (IE)
# (https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) approximation
protein_forcefield="oldff/leaprc.ff99SB", entropy=2, entropy_seg=25,
entropy_temp=298
/

&gb
igb=2, saltcon=0.150,
/

uncomment the next 4 lines for normal mode calculations
#&nmode
#nmstartframe=5, nmendframe=21, nminterval=2,
#maxcyc=50000, drms=0.0001,
#/
```

### Decomposition analysis

```
Sample input file with decomposition analysis
Make sure to include at least one residue from both the receptor
and ligand in the print_res mask of the &decomp section.
http://archive.ambermd.org/201308/0075.html

&general
startframe=5, endframe=21, interval=1,
/

&gb
igb=5, saltcon=0.150,
/

&decomp
idecomp=2, dec_verbose=3,
# This will print all residues that are less than 4 angstroms between
# the receptor and the ligand
print_res="within 4"
/
```

### QM/MMGBSA

```
Sample input file for QM/MMGBSA

&general
startframe=5, endframe=100, interval=5,
/

&gb
igb=5, saltcon=0.100, ifqnt=1, qmcharge_com=0,
qm_residues="100-105, 200", qm_theory="PM3"
/
```

### MM/3D-RISM

```
Sample input file for MM/3D-RISM

&general
startframe=20, endframe=100, interval=5,
/

&rism
polardecomp=1, thermo="gf"
/
```

### MMPBSA with membrane proteins

```
Sample input file for MMPBSA with membrane proteins

&general
startframe=1, endframe=100, interval=1, 
debug_printlevel=2, use_sander=1,
/

&pb
radiopt=0, indi=20.0, istrng=0.150, fillratio=1.25, ipb=1, 
nfocus=1, bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0, npbverb=1,
solvopt=2, inp=2, memopt=1, emem=7.0, mctrdz=-10.383, 
mthick=36.086, poretype=1, maxarcdot=15000
/
```

!!! info
    Comments are allowed by placing a # at the beginning of the line (whites-space are ignored). Variable 
    initialization may span multiple lines. In-line comments (_i.e._, putting a # for a comment after a variable is 
    initialized in the same line) is not allowed and will result in an input error. Variable declarations must be 
    comma-delimited, though all whitespace is ignored. Finally, all lines between namelists are ignored, so comments can
    be added before each namelist without using #.
