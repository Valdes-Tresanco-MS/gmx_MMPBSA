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

`assign_chainID` (Default = 0) 
:   Defines the chains ID assignment mode. _It is ignored when defining a reference structure
(recommended)_. If `assign_chainID = 1`, gmx_MMPBSA check if the structure has no chains ID and it is assigned according
to the structure[^1]. If `assign_chainID = 2`, `gmx_MMPBSA` assign the chains ID, exist or not, according to the
structure[^1] (can generate inconsistencies). If a `*.gro` file was used for complex structure
(`-cs` flag) and not reference structure was provided, `gmx_MMPBSA` assume `assign_chainID = 1`. 

    _New in v1.2.0_

  [^1]: _The chain ID is assigned according to two criteria: **terminal amino acids** and **residue numbering**. If
        both criteria or residue numbering changes are present, we assign a new chain ID. If there are terminal 
        amino acids but the numbering of the residue continues, we do not change the ID of the chain._


`debug_printlevel`
:   gmx_MMPBSA prints errors by raising exceptions, and not catching fatal errors. If `debug_printlevel` is
set to 0, then detailed tracebacks (effectively the call stack showing exactly where in the program the error occurred)
is suppressed, so only the error message is printed. If `debug_printlevel` is set to 1 or higher, all tracebacks are
printed, which aids in debugging of issues. (Default = 0) (Advanced Option)

    _Changed in v1.2.0: Now `gmx_MMPBSA` shows the command-line used to build AMBER topologies when 
    `debug_printlevel = 1` or higher_

`startframe` (Default = 1)
:   The frame from which to begin extracting snapshots from the full, concatenated trajectory comprised of
    every trajectory file placed on the command-line. This is always the first frame read.

`endframe` (Default = 9999999)
:   The frame from which to stop extracting snapshots from the full, concatenated trajectory comprised of every
    trajectory file supplied on the command-line. 

`qh_entropy` (Default = 0) 
:    It specifies whether to perform a quasi-harmonic entropy (QH) approximation with `ptraj` or not.
     
     * 0: Don’t
     * 1: perform QH

    _New in v1.4.2: Equivalent to (Deprecate) `entropy = 1`_

???+ warning "Deprecated in v1.4.2: It will be removed in next version (v1.5.0). Use `qh_entropy` or `interaction_entropy` instead"

    `entropy` (default = 0) 
    :    ~~It specifies whether to perform a quasi-harmonic entropy (QH) approximation with ptraj or the
         [Interaction Entropy (IE)][3] approximation. The allowed values are:~~
         
         * ~~0: Don’t~~
         * ~~1: perform QH~~
         * ~~2: perform IE~~

        _Deprecated since v1.4.2: Divided in `qh_entropy` and `interaction_entropy`_
    
        _Changed in v1.0.0: Include Interaction Entropy approximation_    

`interaction_entropy` (default = 0) 
:    It specifies whether to use the [Interaction Entropy (IE)][3] approximation.
     
     * 0: Don’t
     * 1: perform IE

    !!! tip "Keep in mind"
        - The standard deviation of the interaction energy (σIE) should always be reported.
        - The Interaction Entropy method should be avoided if σIE > ~ 3.6 kcal/mol because it is impossible to 
        converge the exponential average.
        - It is advisable to study how the Interaction Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, as reported in the original [IE publication][3], seems to be 3–40 times too 
        dense. A sampling frequency of 0.1 ps would be more appropriate.
        - The Interaction Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation. 

        Please, consult this [paper][10] for further details.

    _New in v1.4.2: Equivalent to (Deprecate) `entropy = 2`_

  [3]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682
  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374


`ie_segment` (Default = 25)
:    Representative segment (in %), starting from the last frame, for the calculation of the
     Interaction Entropy, _e.g._: `ie_segment = 25` means that the last quartile of the total number of frames
     (`(endframe-startframe)/interval`) will be used to calculate the average Interaction Entropy.

    _New in v1.4.2_

???+ warning "Deprecated in v1.4.2: It will be removed in next version (v1.5.0). Use `ie_segment` instead"

    `entropy_seg` (Default = 25)
    :    ~~Representative segment (in %), starting from the `endframe`, for the calculation of the
         Interaction Entropy, _e.g._: `entropy_seg = 25` means that the last quartile of the total number of frames
         (`(endframe-startframe)/interval`) will be used to calculate the average Interaction Entropy. (Only
          if `entropy = 2`)~~
    
        _Deprecated since v1.4.2: Replaced by `ie_segment`_
        
        _New in v1.0.0_

???+ warning "Deprecated in v1.4.0: It will be removed in next version. Use `temperature` instead"
    
    `entropy_temp` (Default = 298.15)
    :    ~~Specify the temperature to calculate the entropy term `−TΔS` (Only if `entropy` = 2). Avoid inconsistencies 
          with defined internal temperature (298.15 K) when `nmode` is used.~~

`c2_entropy` (default = 0) 
:    It specifies whether to use the [C2 Entropy][11] approximation.
     
     * 0: Don’t
     * 1: perform C2

    !!! tip "Keep in mind"
        - The standard deviation of the interaction energy (σIE) should always be reported.
        - The C2 Entropy method should be avoided if σIE > ~ 6 kcal/mol because it gives unrealistically large entropies.
        - It is advisable to study how the C2 Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, seems to be 3–40 times too dense. A sampling frequency of 0.1 ps would be more 
        appropriate.
        - The C2 Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation. 

        Please, consult this [paper][10] for further details.

    _New in v1.5.0_

  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374
  [11]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00418


`c2_segment` (Default = 25)
:    Representative segment (in %), starting from the last frame, for the calculation of the
     C2 Entropy, _e.g._: `ie_segment = 25` means that the last quartile of the total number of frames
     (`(endframe-startframe)/interval`) will be used to calculate the C2 Entropy.

    _New in v1.5.0_

`exp_ki` (Default = 0.0)
:   Specify the experimental Ki in nM for correlations analysis. If not defined or exp_ki = 0 then this system will be 
omitted in the correlation analysis

    _New in v1.4.0_

`gmx_path` 
:     Define an additional path to search for GROMACS executables. This path takes precedence over the path defined
      in the PATH variable. In these path the following executables will be searched: `gmx`, `gmx_mpi`, `gmx_d`,
      `gmx_mpi_d` (GROMACS > 5.x.x), `make_ndx`, `editconf` and `trjconv` (GROMACS 4.x.x)

    _New in v1.1.1_
    
`interval` (Default = 1)
:     The offset from which to choose frames from each trajectory file. For example, an interval of 2 will pull
      every 2nd frame beginning at startframe and ending less than or equal to endframe. 

`netcdf` (Default = 0)
:     Specifies whether or not to use NetCDF trajectories internally rather than writing temporary ASCII trajectory
      files. For very large trajectories, this could offer significant speedups, and requires less temporary space. 
      However, this option is incompatible with alanine scanning.

      * 0: Do NOT use temporary NetCDF trajectories
      * 1: Use temporary NetCDF trajectories

??? danger "`overwrite_data` Removed" 
    `overwrite_data` (Default = 0)
    :   Defines whether the gmx_MMPBSA data will be overwritten. `gmx_MMPBSA` detects if the gmxMMPBSA data files exist 
        before copying them. This option allows the user to control the copy process and prevents the system from being 
        overloaded by copying files unnecessarily.
    
        * 0: don't
        * 1: overwrite gxmMMPBSA data if exist
    
        !!! tip Keep in mind
            We recommend activating this option with each new release because there may be changes and/or corrections in 
            the gmxMMPBSA data files.
        
        _Removed in v1.4.3_

        _New in v1.3.1_
        

`PBRadii` (Default = 3)
:   PBRadii to build amber topology files:

    * 1: bondi, recommended when igb = 7
    * 2: mbondi, recommended when igb = 1
    * 3: mbondi2, recommended when igb = 2 or 5
    * 4: mbondi3, recommended when igb = 8

`forcefields` (Default = "oldff/leaprc.ff99SB,leaprc.gaff")
:   Comma-separated list of force fields used to build Amber topologies. This variable is more flexible than the 
existing ones (`protein_forcefield` and `ligand_forcefield`). The goal of this variable is to provide convenient 
support for complex systems like this one: [5O8F](https://www.rcsb.org/3d-view/5o8f). It supports all force fields 
tested in `protein_forcefield` and `ligand_forcefield` variables:
    
    !!! inline end important "Keep in mind"
        * You don't need to define it when you use a topology. Please refer to the section 
          ["How gmx_MMPBSA works"](howworks.md#how-gmx_mmpbsa-works)
        * _This notation format is the one used in tleap._

    Protein and Nucleic Acids

    * "oldff/leaprc.ff99"
    * "oldff/leaprc.ff03"
    * "oldff/leaprc.ff99SB"
    * "oldff/leaprc.ff99SBildn"
    * "leaprc.protein.ff14SB"
    
    
    !!! inline end important "Keep in mind"
        * You don't need to define it when you use a topology or the ligand is protein-like type. Please refer to the 
        section ["How gmx_MMPBSA works"](howworks.md#how-gmx_mmpbsa-works)
        * _This notation format is the one used in tleap._

        `*` We added the gmxMMPBSA data to the tleap path. This way, we keep `gmx_MMPBSA` data separated from 
            Amber's.

    Organic molecules, Glycan and Zwitterionic amino acids

    * "leaprc.gaff"
    * "leaprc.gaff2"
    * "leaprc.GLYCAM_06j-1"    (_Compatible with amber12SB and later_)
    * "leaprc.GLYCAM_06EPb"    (_Compatible with amber12SB and later_)
    * "gmxMMPBSA/leaprc.GLYCAM_06h-1"    `*`(_Compatible with amber99SB and earlier_)
    * "gmxMMPBSA/leaprc.zaa99SB"    `*`(_Parameters for Zwitterionic amino acids. Compatible with amber 99SB_)
  
    !!! warning
        Be careful defining this variable since you can define two forces fields with a similar purpose which can 
        generate inconsistencies.    

    _New in v1.4.1_

    _Modiefied in v1.4.3. Internal change_


???+ warning "Deprecated since v1.4.1: It will be removed in v1.5.0. Use `forcefields` instead"
   
    `protein_forcefield` (Default = "oldff/leaprc.ff99SB")
    :   Define the force field used to build Amber topology for proteins. Make sure this force field is the same as the 
        one used in GROMACS. Force fields tested:
    
        * "oldff/leaprc.ff99"
        * "oldff/leaprc.ff03"
        * "oldff/leaprc.ff99SB"
        * "oldff/leaprc.ff99SBildn"
        * "leaprc.protein.ff14SB"
    
        !!! important "Keep in mind"
            * You don't need to define it when you use a topology. Please refer to the section 
              ["How gmx_MMPBSA works"](howworks.md#how-gmx_mmpbsa-works)
            * _This notation format is the one used in tleap._
    
    
???+ warning "Deprecated since v1.4.1: It will be removed in v1.5.0. Use `forcefields` instead"
    
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
            section ["How gmx_MMPBSA works"](howworks.md#how-gmx_mmpbsa-works)
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
        * You don't need to define it when you use a topology. Please refer to the section 
          ["How gmx_MMPBSA works"](howworks.md#how-gmx_mmpbsa-works)   
        * This notation is simpler since these parameter files are generally the same for all systems

`sys_name` (Default = None) (Optional)
:   Define the System Name. This is useful when trying to analyze several systems at the same time or calculating 
the correlation between the predicted and the experimental energies. If the name is not defined, one will be 
assigned when loading it in gmx_MMPBSA_ana according to the order in this is done.

    !!! tip 
        The definition of the system name is entirely optional, however it can provide a better clarity during 
        the results analysis. All files associated with this system will be saved using its name.

    _New in v1.4.0_   

`solvated_trajectory` (Default = 1)
:   Define if it is necessary to build a clean trajectory with no water and ions
    
    * 0: Don’t
    * 1: Build clean trajectory

    _New in v1.3.0_

`temperature` (Default = 298.15)  
: Specify the temperature to calculate experimental ΔG from Ki and the Interaction Entropy (Only if `entropy` = 2). 
Avoid inconsistencies with defined internal temperature (298.15 K) when `nmode` or `qh` are used.  

    !!! warning "Keep in mind"
        Note that the `nmode` as `qh` model to calculate the entropy term are parameterized and will only work at 298.15 K
   
    _New in v1.4.0: Replace `entropy_temp`_

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

`intdiel` (Default = 1.0)
:   Define Internal dielectric constant without use external *.mdin file

`igb` (Default = 5)
:   Generalized Born method to use (seeSection 4).  Allowed values are 1, 2, 5, 7 and 8

    !!! note
        All models are now available with both `mmpbsa_py_energy` and `sander`

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

`ifqnt` (Default = 0)
:   Specifies whether a part of the system is treated with quantum mechanics.
    
    * 0: Potential function is strictly classical
    * 1: Use QM/MM

    !!! note
        This functionality requires `sander`

`qm_residues`
:   Complex residues to treat with quantum mechanics. All residues treated with quantum mechanics in the complex 
    must be treated with quantum mechanics in the receptor or ligand to obtain meaningful results. This notation is 
    the same used for `print_res` variable 

    Notation: [ `CHAIN`/(`RESNUM` or `RESNUM-RESNUM`) ]
    :    Treat with quantum mechanics residues individual or ranges. This notation also supports insertion codes, in 
    which case you must define them individually

    !!! example
        `qm_residues="A/1,3-10,15,100"` This treat with quantum mechanic Chain A residues 1, 3 through 10, 15, and 
        100 from the complex topology file and the corresponding residues in either the ligand and/or receptor 
        topology files.

        Suppost that we can have the following sequence: A:LEU:5, A:GLY:6:A, A:THR:6:B, A:SER:6:C A:ASP:6D, A:ILE:7
        
        === "Right notation"
            
            **Ranges selection**
            :   `qm_residues="A/5-7` Will treat with quantum mechanic all mentioned residues because all residues with 
            insertion code are contained in the range
            
            **Individual selection**
            :   `qm_residues="A/5,6:B,6:C,7` Will treat with quantum mechanic all mentioned residues except the 
            residues A:6:A and A:6:D
            
            **Multiple chain selection**
            :   `qm_residues="A/5-10,100 B/34,56` Will treat with quantum mechanic residues 3 through 10, 100 from 
            chain A and residues 34 and 56 from Chain B.

        === "Wrong notation"
            `qm_residues="A/5-6B,6D-7` Will end in error.

`qm_theory` 
:   Which semi-empirical Hamiltonian should be used for the quantum calculation. Options are PM3, AM1, MNDO, PDDG-PM3, 
    PM3PDDG, PDDG-MNDO, PDDGMNDO, PM3-CARB1, PM3CARB1, DFTB, SCC-DFTB, RM1, PM6, PM3-ZnB, PM3-MAIS, PM3ZNB, MNDO/D, 
    MNDOD. The dispersion correction can be switched on for AM1 and PM6 by choosing AM1-D* and PM6-D, respectively. The 
    dispersion and hydrogen bond correction will be applied for AM1-DH+ and PM6-DH+.

    !!! danger
         No default, this must be specified if QM/MM calculations are going to be performed.

`qmcharge_com` (Default = 0)
:   The charge of the quantum section for the complex.

`qmcharge_lig` (Default = 0)
:   The charge of the quantum section of the ligand.

`qmcharge_rec` (Default = 0)
:   The charge of the quantum section for the receptor.

`qmcut` (Default = 9999.0)
:   The cutoff for the qm/mm charge interactions.
    
  [4]: https://ambermd.org/doc12/Amber20.pdf#chapter.32

### **`&pb` namelist variables**

!!! note
    PB calculations can be performed either with **_mmpbsa_py_energy_** or **_sander_**. Of note, **_mmpbsa_py_energy_**
    works with a lower number of variables compared to **_sander_**. Nevertheless, **_mmpbsa_py_energy_** will suffice 
    for most routine applications of PB. 

    **_mmpbsa_py_energy_** variables: `inp`, `smoothopt`, `radiopt`, `npbopt`, `solvopt`, `maxitn (linit)`, `nfocus`, 
    `fscale`, `epsin (indi)`, `epsout (exdi)`, `istrng`, `dprob (prbrad)`, `iprob`, `accept`, `fillratio`, 
    `space (1/scale)`, `bcopt`, `eneopt`, `cutnb`, `sprob`, `cavity_surften`, `cavity_offset`

    On the other hand, **_sander_** offers access to all [pbsa][5] functionalities.

    **_sander_** variables: `ntb`, `cut`, `nsnb`, `imin`, `maxcyc`, `ipb`, `inp`, `ioutfm`, `ntx`, `epsin (indi)`, 
    `epsout (exdi)`, `istrng`, `radiopt`, `sprob`, `dprob (prbrad)`, `space (1/scale)`, `maxitn (linit)`, 
    `cavity_surften`, `cavity_offset`, `fillratio`, `epsmem (emen)`, `membraneopt (memopt)`, `sasopt`, `mthick`, 
    `maxarcdot`, `solvopt`, `nfocus`, `bcopt`, `eneopt`, `frcopt`, `cutfd`, `cutnb`, `mctrdz`, `poretype`, `npbverb`, 
    `npbopt`, `pbtemp0`, `iprob`, `arcres`, `mprob`, `accept`, `nbuffer`, `npbgrid`, `scalec`, `nsnba`, `phiout`, 
    `phiform`, `decompopt`, `use_rmin`, `vprob`, `rhow_effect`, `use_sav`, `maxsph`

    Hereafter, a selected group of variables is presented, which should suffice for most PB calculations.  The default 
    values for these parameters are appropriate for most calculations on solvated molecular systems. Also note that 
    the default options may have changed over time. A more thorough description of all the options can be found 
    [here][5]. For a detailed discussion of all related options on the quality of the MM/PB(GB)SA calculations, 
    please check this [publication][6].

`npbopt` (Default = 0) 
:   Option to select the linear, or the full nonlinear PB equation.

    * 0: Linear PB equation is solved
    * 1: Nonlinear PB equation is solved

    !!! note
        While the linear PB equation will suffice for most calculations, the nonlinear PB equation is recommended 
        for highly charged systems. It is supported in both **_mmpbsa_py_energy_** and **_sander_**. Parameters such 
        as `eneopt` or `cutnb` should be adjusted accordingly when using the nonlinear PB equation. Check the following 
        threads on how to proceed in [mmbsa_py_energy](http://archive.ambermd.org/201203/0191.html) or 
        [sander](http://archive.ambermd.org/201610/0114.html) when using the nonlinear PB equation. Last but not 
        least, take into account that using nonlinear PB equation can significantly increase the calculation time 
        required for PB calculation.

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
        This is not used for `APBS`

`cavity_surften` (Default = 0.0378  [ kcal/mol Å^2^ ] )
:   Surface tension. Unit conversion to kJ done automatically for `APBS`.

`exdi` (Default = 80.0)
:   External dielectric constant. This corresponds to `epsout` in [pbsa][5].

`indi` (Default = 1.0)
:   Internal dielectric constant. This corresponds to `epsin` in [pbsa][5].

`fillratio` (Default = 4.0)
:   The ratio between the longest dimension of the rectangular finite-difference grid and that of the solute.

`scale` (Default = 2.0)
:   Resolution of the Poisson Boltzmann grid. It is equal to the reciprocal of the grid spacing (`space` in [pbsa][5]). 

`istrng` (Default = 0.0)
:   Ionic strength in Molarity. It is converted to mM for `PBSA` and kept as M for `APBS`.

`linit` (Default = 1000) 
:   Maximum number of iterations of the linear Poisson Boltzmann equation to try. This corresponds to `maxitn` 
    in [pbsa][5].

`prbrad` (Default = 1.4)
:   Solvent probe radius in Angstroms. Allowed values are 1.4 and 1.6. This corresponds to `dprob` in [pbsa][5].

`radiopt` (Default = 1)
:   The option to set up atomic radii according to:

    * 0: the prmtop, or 
    * 1: pre-computed values 
    
    !!! warning
        `radiopt=0` is recommended which means using radii from the prmtop file for both the PB calculation and for 
        the NP. Check this [thread](http://archive.ambermd.org/201303/0548.html) and See Amber manual for more 
        complete description.

`sander_apbs` (Default = 0)
:   Option to use `APBS` for `PB` calculation instead of the built-in `PBSA` solver. This will work only through the
    `iAPBS` interface built into `sander.APBS`. Instructions for this can be found online at the iAPBS/APBS websites.
    
    * 0: Don’t use `APBS`
    * 1: Use `sander.APBS` 

`memopt` (Default = 0)
:   Turn on membrane protein support. This corresponds to `membraneopt` in [pbsa][5].

`emem` (Default = 1.0)
:   Membrane dielectric constant. This corresponds to `epsmem` in [pbsa][5].

`mthick` (Default = 40.0)
:   Membrane thickness.

`mctrdz` (Default=0.0, use protein center as the membrane center)
:   Absolute membrane center in the z-direction.

`poretype` (Default=1)
:   Turn on the automatic membrane channel/pore finding method.

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
:   Defines the residue by which it is going to mutate. Allowed values are: `"ALA"` or `"A"` for Alanine scanning and 
    `"GLY"` or `"G"` for Glycine scanning.

    _Changed in v1.3.0: Change mol (receptor or ligand) by mutant aminoacid (ALA or GLY)_

`mutant_res` (Default = None. Most be defined)
:   Define the specific residue that is going to be mutated. Use the following format CHAIN:RESNUM (eg: 'A:350') or 
CHAIN:RESNUM:INSERTION_CODE if applicable (eg: "A:27:B"). 

    !!! important
        * Only one residue for mutation is supported!
        * We recommend using the reference structure (-cr) to ensure the perfect match between the selected residue in 
        the defined structure or topology 
        * This option allow `gmx_MMPBSA` to do the mutation. This way the user does not have to provide the mutant 
        topology
    
    _Changed in v1.4.0: Allow mutation in antibodies since it support insertion code notation_

`cas_intdiel` (Default = 0)
:   The dielectric constant (`intdiel`(GB)/`indi`(PB)) will be modified depending on the nature of the residue to be 
mutated. 
    
    * 0: Don’t
    * 1: Adaptative `intdiel` assignation

    !!! important
        * Works with the GB and PB calculations
        * It is ignored when `intdiel`(GB)/`indi`(PB) has been explicitly defined, that is, it is ignored if 
        `intdiel != 1.0`/`indi != 1.0` (default values)
        * Dielectric constant values has been assigned according to [Yan et al., 2017][9]

  [9]: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00734
    
    _New in v1.4.2_

`intdiel_nonpolar` (Default = 1)
:   Define the `intdiel`(GB)/`indi`(PB) value for nonpolar residues (`PHE`, `TRP`, `VAL`, `ILE`, `LEU`, `MET`, `PRO`,
`CYX`, `ALA`, `GLY`, `PRO`)
    
    _New in v1.4.2_

`intdiel_polar` (Default = 3)
:   Define the `intdiel`(GB)/`indi`(PB) value for polar residues (`TYR`, `SER`, `THR`, `CYM`, `CYS`, `HIE`, `HID`, 
`ASN`, `GLN`, `ASH`, `GLH`, `LYN`)
    
    _New in v1.4.2_

`intdiel_positive` (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for positive charged residues (`LYS`, `ARG`, `HIP`)
    
    _New in v1.4.2_

`intdiel_negative` (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for negative charged residues (`GLU`, `ASP`)
    
    _New in v1.4.2_

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

`print_res` (Default = "within 6")
:   Select residues from the complex to print. The default selection should be sufficient in most cases, however we 
have added several additional notations
    
    === "By Distance"
        Notation: [ `within` `distance` ]
        :   `within` corresponds to the keyword and `distance` to the maximum distance criterion in 
            Angstroms necessary to select the residues from both the receptor and the ligand

        !!! example
            `print_res="within 6"` Will print all residues within 6 Angstroms between receptor and 
            ligand including both.

    === "Amino acid selection"
        Notation: [ `CHAIN`/(`RESNUM` or `RESNUM-RESNUM`) ]
        :   Print residues individual or ranges. This notation also supports insertion codes, in which case you must 
            define them individually

        !!! example
            `print_res="A/1,3-10,15,100 B/25"` This will print Chain A residues 1, 3 through 10, 15, and 100 along with 
            chain B residue 25 from the complex topology file and the corresponding residues in either the ligand and/or 
            receptor topology files.

            !!! danger
                make sure to include at least one residue from both the receptor and ligand in the `print_res` mask of the 
                `&decomp` section. Check http://archive.ambermd.org/201308/0075.html

            Suppost that we can have the following sequence where chain A is the receptor and B is the ligand: 
            A:LEU:5, A:GLY:6:A, A:THR:6:B, A:SER:6:C A:ASP:6D, A:ILE:7 , B:25
            
            === "Supported notation"
                
                **Ranges selection**
                :   `print_res="A/5-7 B/25` Will print all mentioned residues because all residues with insertion code are 
                    contained in the range
                
                **Individual selection**
                :   `print_res="A/5,6:B,6:C,7 B/25` Will print all mentioned residues except the residues A:6:A and A:6:D

            === "Wrong notation"
                `print_res="A/5-6B,6D-7` Will end in error.

    === "All"

        Notation: `all`
        :   will print all residues. This option is often not recommended since most residues contribution is zero and 
            it is just going to be a waste of time and computational resources.

        !!! danger
            Using idecomp=3 or 4 (pairwise) with a very large number of printed residues and a large number of frames 
            can quickly create very, very large temporary mdout files. Large print selections also demand a large amount 
            of memory to parse the mdout files and write decomposition output file (~500 MB for just 250 residues, since 
            that’s 62500 pairs!) It is not unusual for the output file to take a significant amount of time to print if 
            you have a lot of data. This is most applicable to pairwise decomp, since the amount of data scales as  
            O(N^2^).
 
    !!! important
        We recommend using the reference structure (-cr) to ensure the perfect match between the selected residue in 
        the defined structure or topology 
        
    
    _Changed in v1.4.0: Improve residue selection_

### **`&rism` namelist variables**

???+ warning
    `3D-RISM` calculations are performed with the `rism3d.snglpnt` program built with AmberTools, written by Tyler 
    Luchko. It is the most expensive, yet most statistical mechanically rigorous solvation model. See [RISM chapter][7] 
    for a thorough description of options and theory. A list of references can be found there, too.
    
    We have included more variables in 3D-RISM calculations than the ones available in the MMPBSA.py original code. That 
    way, users can be more in control and tackle various issues (e.g., convergence problems).

    **3D-RISM variables and their default values:** 

    `closureorder= 1`, `asympcorr= 1`, `buffer= 14.0`, `solvcut= None (If -1 or no value is specified then the 
    buffer distance is used)`, 
    `grdspc= 0.5,0.5,0.5`, `ng= -1,-1,-1`, `solvbox= -1,-1,-1`, `tolerance= 0.00001 (1.0e-5)`, `mdiis_del= 0.7`, 
    `mdiis_restart= 10.0`, `mdiis_nvec= 5`, `maxstep= 10000`, `npropagate= 5`, `centering= 1`, `polarDecomp= 0`, 
    `entropicDecomp= 0`, `gf= 0`, `pc+= 0`, `uccoeff= 0.0,0.0,0.0,0.0`, `rism_verbose= 0`, `treeDCF= 1`, `treeTCF= 1`, 
    `treeCoulomb= 0`, `treeDCFOrder= 2`, `treeTCFOrder= 2`, `treeCoulombOrder= 2`, `treeDCFN0= 500`, `treeTCFN0= 500`, 
    `treeCoulombN0= 500`, `treeDCFMAC= 0.1`, `treeTCFMAC= 0.1`, `treeCoulombMAC= 0.1`, `asympKSpaceTolerance= -1.0`, 
    `ljTolerance= -1.0`

    One advantage of `3D-RISM` is that an arbitrary solvent can be chosen; you just need to change the `xvvfile` 
    specified on the command line (see [34.3.2][8]).

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
        No default, this must be set if buffer < 0. Define like `ng=1000,1000,1000`

`solvbox`
:   Length of the solvation box in the x, y, and z dimensions. Used only if buffer < 0. Mutually exclusive with
    `buffer` and `grdspc` above, and paired with `ng` above. 

    !!! warning 
        No default, this must be set if buffer < 0. Define like `solvbox=20,20,20`

`polardecomp` (Default = 0)
:   Decompose the solvation free energy into polar and non-polar contributions. Note that this will increase 
    computation time by roughly 80%. 
    
    * 0: Don’t decompose solvation free energy. 
    * 1: Decompose solvation free energy. 

`rism_verbose` (Default = 0)
:   Level of output in temporary RISM output files. May be helpful for debugging or following convergence. 

    * 0: just print the final result
    * 1: additionally prints the total number of iterations for each solution
    * 2: additionally prints the residual for each iteration and details of the MDIIS solver (useful for debugging 
    and convergence analyses)

`solvcut`  (Default = buffer )
:   Cutoff used for solute-solvent interactions. The default is the value of buffer. Therefore, if you set `buffer` < 
    0 and specify `ng` and `solvbox` instead, you must set `solvcut` to a nonzero value, or the program will quit in 
    error.

`thermo` (Default = "std")
:   Which thermodynamic equation you want to use to calculate solvation properties. Options are "std", "gf", or 
    "both" (case-INsensitive). "std" uses the standard closure relation, "gf" uses the Gaussian Fluctuation 
    approximation, and "both" will print out separate sections for both.
    
    !!! note
        Note that all data are printed out for each RISM simulation, so no choice is any more computationally demanding 
        than another.

`tolerance` (Default = 1e-5)
:   Upper bound of the precision requirement used to determine convergence of the self-consistent solution. This has 
    a strong effect on the cost of 3D-RISM calculations (smaller value for tolerance -> more computation).

## Sample input files

!!! tip
    You can refer to the [examples](examples/3D-RISM/README.md) to understand the input file in a practical way.

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
mutant='ALA', mutant_res='B:12'
/
```

### Entropy

``` linenums="1"
Sample input file for entropy calculations

&general
startframe=5, endframe=21, verbose=2, interval=1,
# `entropy` variable control whether to perform a quasi-harmonic entropy (QH)
# or the Interaction Entropy (IE)
# (https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) approximation
protein_forcefield="oldff/leaprc.ff99SB", entropy=2, entropy_seg=25,
temperature=298
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
http://archive.ambermd.org/201308/0075.html. This is automally
guaranteed if use "within" keyword.

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
qm_residues="B/240-251", qm_theory="PM3"
/
```

### MM/3D-RISM

```
Sample input file for 3D-RISM

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
