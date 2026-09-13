---
template: main.html
title: The input file
---

# The input file

## Description

The `gmx_MMPBSA` input file contains the settings for each calculation. Its syntax is similar to that used by other
Amber programs, but gmx_MMPBSA also accepts a format resembling GROMACS `*.mdp` files (see below). The input file is
divided into sections called namelists, in which the variables for each calculation are defined. The allowed
namelists are:

- [`&general`](input_file.md#general-namelist-variables): contains variables that apply to all aspects of the 
  calculation or parameters required for building AMBER topologies from GROMACS files.
- [`&gb`](input_file.md#gb-namelist-variables): variables specific to Generalized Born (GB) calculations.
- [`&gbnsr6`](input_file.md#gbnsr6-namelist-variables): variables specific to GBNSR6 calculations.
- [`&pb`](input_file.md#pb-namelist-variables): variables specific to Poisson–Boltzmann (PB) calculations.
- [`&rism`](input_file.md#rism-namelist-variables): variables specific to 3D-RISM calculations.
- [`&alanine_scanning`](input_file.md#alanine_scanning-namelist-variables): variables specific to alanine-scanning
  calculations.
- [`&decomp`](input_file.md#decomp-namelist-variables): variables specific to decomposition calculations.
- [`&nmode`](input_file.md#nmode-namelist-variables): variables specific to normal-mode (NMODE) calculations used to
  approximate vibrational entropies.

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h

## Generation of input files with gmx_MMPBSA
The input file can be created using `gmx_MMPBSA` by selecting the calculations you want to perform.

``` title="Command-line"
gmx_MMPBSA --create_input args

where `args` can be:  gb, gbnsr6, pb, pb_mem, rism, ala, decomp, nmode, all

`pb_mem` creates a membrane-oriented `&pb` template with automatic membrane
center/thickness detection from phosphorus atoms in `-ct` and membrane-specific
PB defaults.
```

Example:
=== "GB calculation"
        
        gmx_MMPBSA --create_input gb
    
=== "PB calculation"
    
        gmx_MMPBSA --create_input pb

=== "Membrane PB calculation"

        gmx_MMPBSA --create_input pb_mem

=== "GB, PB and Decomposition calculations"
    
        gmx_MMPBSA --create_input gb pb decomp

=== "All calculations"

        gmx_MMPBSA --create_input
     
    or 
        
        gmx_MMPBSA --create_input all
        
!!! Danger 
    Note that several variables must be explicitly defined in the input file

!!! note "Historical version markers"
    Version markers such as “Implemented in v1.5.0” below record when an option was introduced or changed. They are
    historical provenance, not claims that the option is new in the 1.7.0 release.

_Introduced in v1.5.0 (historical)_

!!! tip "Viewing wide examples"
    Generated input examples and option tables are intentionally wide. On a narrow screen, scroll horizontally to
    view the complete lines and values.

## Format
The input variables are described below by namelist. Some descriptions are adapted from the original sources. Enter
integer and floating-point values without quotes, and enclose strings in either single or double quotes. Assign each
variable with `variable = value`. Variables on the same line must be separated by commas; commas are optional when each
variable is placed on a separate line. See the [sample input files](#sample-input-files) below. Several calculations can
be performed in one run (_e.g._, `&gb` with `&pb`, `&gb` with `&alanine_scanning`, or `&pb` with `&decomp`). You can
also generate an input file with the `--create_input` option. gmx_MMPBSA accepts both the traditional Amber/Fortran
namelist style and the GROMACS-like multiline style, and the two styles can be combined:

=== "Multiline format style "
    ``` title="Multiline format style input file example"
            
    # General namelist variables
    &general
      sys_name             = ""                      # System name
      startframe           = 1                       # First frame to analyze
      endframe             = 9999999                 # Last frame to analyze
      ...
      interval              = 1                      # The offset from which to choose frames from each trajectory file
    /
    
    # Generalized-Born namelist variables
    &gb
      igb                  = 8                       # GB model to use
      ...
      probe                = 1.4                     # Solvent probe radius (Å) for surface-area calculation
    /
    ```

=== "Old format style"
    ``` title="Old format style Input file example"
            
    # General namelist variables
    &general
      sys_name = "", startframe = 1, endframe = 9999999
      ...
      interval = 1
    /
    
    # Generalized-Born namelist variables
    &gb
      igb = 8,
      ...
      probe = 1.4
    /
    ```

## Namelists

### **`&general` namelist variables**

#### **Basic input options**

`sys_name` (Default = None)
:   Define the System Name. This is useful when trying to analyze several systems at the same time or calculating 
the correlation between the predicted and the experimental energies. If the name is not defined, one will be assigned 
when loading the system in `gmx_MMPBSA_ana` on a first-come, first-served basis.

    !!! tip 
        Defining the system name is optional, but a descriptive name can make result analysis clearer. All files
        associated with the system will be saved using this name.

    _Implemented in v1.4.0_  

`startframe` (Default = 1)
:   The frame from which to begin extracting snapshots from the full, concatenated trajectory composed of
every trajectory file placed on the command-line, in command-line order. This is always the first frame read.

`endframe` (Default = 9999999)
:   The frame at which to stop extracting snapshots from the full, concatenated trajectory composed of every
trajectory file supplied on the command-line.

`interval` (Default = 1)
:   The offset from which to choose frames from the combined trajectory stream. For example, an interval of 2 will
pull every 2nd frame beginning at `startframe` and ending less than or equal to `endframe`; the interval is not
restarted when the next trajectory file is read.

    !!! warning "Multiple complex trajectories"
        Multiple files supplied with `-ct` are concatenated in the order given and analyzed as one trajectory.
        `startframe`, `endframe`, and `interval` are applied to this combined frame stream. The resulting energy
        statistics are pooled over the concatenated frames; they are not calculated as independent-replica
        statistics. If the input files are replicas, this concatenated interpretation should be reported explicitly.

`explicit_waters` (Default = 0)
:   Number of explicit water molecules to keep in the working complex trajectory and assign to the receptor. A value of
0 disables this feature and preserves the standard dry MM/PB(GB)SA workflow.

    !!! important "Scope"
        Explicit receptor waters are currently supported for single-trajectory GB, GBNSR6, PB, RISM, and normal-mode
        entropy calculations. Quasi-harmonic entropy is not available for new calculations, and multi-trajectory
        receptor/ligand inputs are not supported with this mode.
        Extra-point water models such as OPC or TIP4P can fail in `sander` because of their virtual-site atoms. By
        default, `gmx_MMPBSA` stops when these atoms are found. Set `explicit_waters_extra_points="strip"` only if you
        intentionally want to remove the virtual sites and use the result as an approximate relative comparison.

    !!! note "Multiple trajectories with explicit waters"
        When multiple `-ct` files are supplied, the concatenated frame selection is performed before
        `cpptraj closest` selects the nearest explicit waters. The reference/interface mask remains static, while
        the selected water identities may change from frame to frame.

`explicit_waters_mask` (Default = "")
:   Reference selection used to choose the closest explicit waters when `explicit_waters > 0`. Accepted values are:

    * An Amber residue mask, for example `":4,6,7,8,9"`
    * A decomposition-style distance selection, for example `"within 4"`
    * `"dASA"` to identify interface residues with cpptraj dASA before selecting closest waters

    The selected interface/reference residues are static for the calculation. The water identities can change from frame
    to frame because `cpptraj closest` is applied during trajectory processing.

`explicit_waters_group` (Default = "")
:   Solvent group name in the complex index file. When empty, `gmx_MMPBSA` looks for common solvent group names and
water model names such as `SOLV`, `SOL`, `Water`, `WAT`, `TP3`, `TIP3P`, `SPC`, and `OPC`. Set this option only when
the solvent group in the index file uses a custom name.

`explicit_waters_dasa_cutoff` (Default = 0.5)
:   dASA cutoff used when `explicit_waters_mask="dASA"`. This option is ignored for Amber mask and
    `within <distance>` selections.

`explicit_waters_as` (Default = "receptor")
:   Molecule that receives the explicit waters. The only supported value is `"receptor"` in the current implementation.

`explicit_waters_extra_points` (Default = "error")
:   How to handle virtual-site/extra-point atoms in selected explicit waters. Accepted values are:

    * `"error"`: stop when extra-point atoms are found.
    * `"strip"`: remove extra-point atoms from the explicit-water topologies and trajectories, with a warning.

    The `"strip"` option changes the electrostatics of OPC/TIP4P-style waters and should be used only for controlled
    relative comparisons where this approximation is acceptable.

#### **Topology and parameter provenance**

The calculation uses the supplied topology rather than rebuilding parameters from structures. Bonded terms,
nonbonded terms, charges, ligand parameters, and ion parameters must already be present in the topology used for the
MD system. For GROMACS calculations, provide the complex topology with `-cp` and the unbound topologies with `-rp`
and `-lp` when using the multiple-trajectory protocol. Native AMBER calculations use the corresponding AMBER
topologies.

`forcefields`, `protein_forcefield`, `ligand_forcefield`, and `ions_parameters`
:   Legacy topology-preparation settings. They are retained only for compatibility with older input files and are not
    used to parameterize current topology-based calculations. Do not add them to new input files. Prepare the
    topology externally with the force field and ion model intended for the MD simulation.

    Historical structure-to-topology workflows using `tleap`, `parmchk2`, and `-lm` are not part of the current
    GROMACS calculation path. A small-molecule ligand must already be represented in the supplied GROMACS topology.

`PBRadii` (Default = 4)
:   Continuum-radius set used by the implicit-solvent calculation. The value may be given as the numeric code or the
    corresponding named set (for example, `PBRadii=mbondi3`). It does not select a bonded or nonbonded force field.

    This is the continuum-radius set used by GB, PB, and GBNSR6 topology preparation. Continuum radii are distinct
    from Lennard-Jones radii and are part of the scoring-model parameterization.

    The default value is `4` (`mbondi3`) to match the default GB-Neck2 model (`igb = 8`), so it does not need to be
    specified for standard GB calculations using the default GB model.

    * 1: bondi, recommended when igb = 7
    * 2: mbondi, recommended when igb = 1
    * 3: mbondi2, recommended when igb = 2 or 5
    * 4: mbondi3, recommended when igb = 8

    !!! note "Native AMBER topologies"
        When using `amber_MMPBSA` with a native AMBER `prmtop`, the per-atom `RADII` and `SCREEN` values already
        stored in the topology are preserved. In this case, `PBRadii` does not rebuild or replace the normal working
        topologies; select the radius set during `tleap` preparation. `amber_MMPBSA` warns when the recognized
        topology `RADIUS_SET` does not conventionally match the selected `igb`, but it does not change or reject the
        combination automatically. For GROMACS inputs, `gmx_MMPBSA` applies the selected `PBRadii` to the generated
        AMBER topology and gives the same advisory warning when it differs from the conventional choice for `igb`.

    * 5: mbondi_pb2

        !!! note "_mbondi_pb2_ radii set"
            This radii set was added in _gmx_MMPBSA v1.5.0_. It is based on the _mbondi_ radii set and contains a
            new optimized set of halogen PB radii for halogenated compounds (without extra point (EP) of charge) 
            parametrized with General Amber Force Field (GAFF):

            Values from Table 3 in [§3.1 Halogen Radii Optimization Without EP][300]:

            * Cl: 1.76
            * Br: 1.97
            * I: 2.09
    
            !!! warning 
                Make sure that you installed the updated ParmEd
                ```
                python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236
                ```

            This radii set should be used with the following PBSA setup:
    
            ```
            Sample input file for PB calculation with halogenated compounds
            
            &general
            sys_name="PB_Halogens",
            PBRadii=5,
            /
            &pb
            radiopt=0, istrng=0.150, inp=1,
            /
            ```

    * 6: mbondi_pb3

        !!! note "_mbondi_pb3_ radii set"
            This radii set was added in _gmx_MMPBSA v1.5.x_. It is based on the _mbondi_ radii set and contains a
            new optimized set of halogen PB radii for halogenated compounds (without extra point (EP) of charge) 
            parametrized with General Amber Force Field (GAFF):

            Values from Table 3 in [§3.1 Halogen Radii Optimization Without EP][300]:

            * Cl: 2.20
            * Br: 2.04
            * I: 2.19
                    
            !!! warning 
                Make sure that you installed the updated ParmEd
                ```
                python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236
                ```

            This radii set should be used with the following PBSA setup:
    
            ```
            Sample input file for PB calculation with halogenated compounds
            
            &general
            sys_name="PB_Halogens",
            PBRadii=6,
            /
            &pb
            radiopt=0, istrng=0.150, inp=2,
            /
            ```
            
  [300]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    * 7: charmm_radii (compatible only with &pb)
        
        !!! note "_charmm_radii_ radii set"
            This radii set was added in _gmx_MMPBSA v1.5.0_. **Use it only with systems prepared with CHARMM
            force fields**. This atomic radii set for Poisson-Boltzmann calculations has been derived from average 
            solvent electrostatic charge distribution with explicit solvent. The accuracy has been tested with free 
            energy perturbation with explicit solvent. Most of the values were taken from a _*radii.str_ file used in 
            PBEQ Solver in [charmm-gui](https://www.charmm-gui.org/?doc=input/pbeqsolver).

            * Radii for protein atoms in 20 standard amino acids from 
            [Nina, Belogv, and Roux](https://pubs.acs.org/doi/10.1021/jp970736r)
            * Radii for nucleic acid atoms (RNA and DNA) from 
            [Banavali and Roux](https://pubs.acs.org/doi/abs/10.1021/jp025852v)
            * Halogens and other atoms from [Fortuna and Costa](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
            
            !!! warning 
                Make sure that you installed the updated ParmEd
                ```
                python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236
                ```

    _Updated in v1.5.0: New PB radii sets have been added_

`radii_audit` (Default = 0)
:   Set to `1` to write `GMXMMPBSA_radii_<component>.csv` files with per-atom continuum-radius provenance. Normal runs
    always write the compact `GMXMMPBSA_radii.json` summary and include the same summary in `_GMXMMPBSA_info`.
    The audit flags metals, dummy/extra-point atoms, unknown elements, nonstandard residues, and assignments whose
    exact rule cannot be established. It records the final topology arrays; it does not change them.

`source_force_field` (Default = "auto")
:   Optional source force-field family override used for provenance and advisories. Allowed values are `auto`,
    `amber`, `charmm`, `opls`, `gromos`, and `other`. Automatic classification is conservative; use an override when
    the topology does not identify the source family clearly. This setting never selects a radius set or changes a
    topology.

`temperature` (Default = 298.15)  
:   Specify the temperature (in K) used in the calculations.
   
    _Implemented in v1.4.0: Replace `entropy_temp`_

    _Reworked in v1.5.0: Temperature value used not only for entropy but all calculations

#### **Entropy options**

`qh_entropy` (Default = 0)
:    Legacy quasi-harmonic entropy (QH) calculations are no longer supported for new runs. Keep this value at `0`.
     Historical result files containing QH data remain readable by the analyzer in this release only. All QH
     compatibility support will be removed after this release.

    !!! important "Keep in mind"
        QH output is retained only for backward-compatible reading of existing results in this final compatibility
        release. Use NMODE, interaction entropy, or C2 entropy for new calculations. QH will not be available after
        this release.

    _Implemented in v1.4.2: Equivalent to (Removed) `entropy = 1`_

`interaction_entropy` (default = 0)
:    It specifies whether to use the [Interaction Entropy (IE)][3] approximation.
     
     * 0: Do not perform IE calculation
     * 1: Perform IE calculation

    !!! note "Keep in mind"
        - The Interaction Entropy can be calculated independently of the solvent model used.
        - A sample Interaction Entropy input file is shown [here](input_file.md#interaction-entropy)
        - A tutorial on the use of Interaction Entropy is 
        available [here](examples/Entropy_calculations/Interaction_Entropy/README.md)
        - The standard deviation of the interaction energy (`σIE`) should always be reported when using the Interaction 
        Entropy method.
        - The Interaction Entropy method should be avoided if `σIE > ~ 3.6 kcal/mol` because it is impossible to 
        converge the exponential average.
        - It is advisable to study how the Interaction Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, as reported in the original [IE publication][3], seems to be 3-40 times too 
        dense. A sampling frequency of 0.1 ps would be more appropriate.
        - The Interaction Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation. 

        See this [paper][10] for further details.

    _Implemented in v1.4.2: Equivalent to (Removed) `entropy = 2`_

    _Updated in v1.5.0: Now reports the `σIE`. Chart improved in `gmx_MMPBSA_ana`_

  [3]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682
  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374


`ie_segment` (Default = 25)
:    Convergence-diagnostic segment (in %), starting from the last frame. The Interaction Entropy reported as the
primary result is calculated from the full selected ensemble. For example, `ie_segment = 25` retains the last
quartile of the cumulative IE curve to report its tail mean and fluctuation as a diagnostic; it does not replace the
full-ensemble IE estimate. Nonoverlapping block results at several values of `N` are also written to the IE output.

    _Implemented in v1.4.2_

`c2_entropy` (default = 0) 
:    It specifies whether to use the [C2 Entropy][11] approximation.
     
     * 0: Do not perform C2 calculation
     * 1: Perform C2 calculation

    !!! note "Keep in mind"
        - The C2 Entropy can be calculated independently of the solvent model used.
        - A tutorial on the use of C2 Entropy is 
        available [here](examples/Entropy_calculations/C2_Entropy/README.md)
        - The standard deviation of the interaction energy (`σIE`) should always be reported.
        - C2 Entropy magnitudes may become unrealistic if `σIE > ~ 6.0 kcal/mol` (approximately 25 kJ/mol), even
        though C2 converges more readily than IE.
        - It is advisable to study how the C2 Entropy depends on N by block averaging (which also provide an 
        estimate of the precision of the calculated entropies).
        - A sampling frequency of 10 fs, seems to be 3-40 times too dense. A sampling frequency of 0.1 ps would be more 
        appropriate.
        - The C2 Entropy results may vary depending on the system flexibility or whether constraints were used 
        or not in the MD simulation.

        See this [paper][10] for further details.

    _Implemented in v1.5.0_

  [10]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00374
  [11]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00418

#### **Miscellaneous options**

`assign_chainID` (Default = 0) 
:   Defines the chains ID assignment mode. _It is ignored when defining a reference structure
(recommended)_. If `assign_chainID = 1`, `gmx_MMPBSA` check if the structure has no chains ID, and it is assigned 
according to the structure[^1]. If `assign_chainID = 2`, `gmx_MMPBSA` assign the chains ID, exist or not, 
according to the structure[^1] (can generate inconsistencies).

    _Implemented in v1.2.0_

    _Updated in v1.5.0: Internal changes_

  [^1]: _The chain ID is assigned according to two criteria: **terminal amino acids** and **residue numbering**. If
        both criteria or residue numbering changes are present, we assign a new chain ID. If there are terminal 
        amino acids, but the numbering of the residue continues, we do not change the ID of the chain._

`exp_ki` (Default = 0.0)
:   Specify the experimental Ki (in nM) for correlations analysis. If not defined or exp_ki = 0 then this system 
will be omitted in the correlation analysis

    _Implemented in v1.4.0_

`full_traj` (Default = 0)
:   Print trajectories

    * 0: Print only MPI-rank trajectories in *.mdcrd format
    * 1: Print a full trajectory and the MPI-rank trajectories in *.mdcrd format

`gmx_path` 
:   Define a path to search for GROMACS executables. This path takes precedence over the path defined
in the PATH variable. In this path the following executables will be searched: `gmx`, `gmx_mpi`, `gmx_d`, or
`gmx_mpi_d` (GROMACS > 5.x.x), and `make_ndx`, `editconf` and `trjconv` (GROMACS 4.x.x)

    !!! note "Keep in mind"
        Use this variable when the GROMACS installation used to run the simulation differs from the installation used
        for the analysis. Set it to the GROMACS `bin` directory in which the executables should be located.
        An example of the use of this variable is given below:

            &general
            sys_name="my_system",
            verbose=2
            gmx_path="/home/programs/gromacs/bin"
            /
            &gb
            igb=8, saltcon=0.150
            /
           
            # replace this "/home/programs/gromacs/bin" with the path to the GROMACS you want to use.

    _Implemented in v1.1.1_

`keep_files` (Default = 2)
:   Specifies which files are kept.

    * 0: Keep only the normalized input metadata and compact result (_GMXMMPBSA_info_ and _COMPACT_MMXSA_RESULTS.mmxsa_)
    * 1: Keep all temporary files (_prefix_*)
    * 2: Keep all temporary files (_prefix_*) and binary file

    !!! note "Keep in mind"
        The binary file contains all the information necessary to analyze the data with gmx_MMPBSA_ana. Its 
        use is only recommended in special cases where the original 
        raw data are not present. We recommend using `keep_files = 1` or `keep_files = 2` to
        avoid errors when using `gmx_MMPBSA_ana`.

`netcdf` (Default = 0)
:   Specifies whether to use NetCDF trajectories internally rather than writing temporary ASCII trajectory
files. For very large trajectories, this could offer significant speedups, and requires less temporary space. 
However, this option is incompatible with alanine scanning.

    * 0: Do NOT use temporary NetCDF trajectories
    * 1: Use temporary NetCDF trajectories

`solvated_trajectory` (Default = 1)
:   Define if it is necessary to generate a clean trajectory with no water and ions.
    
    * 0: Do not generate clean trajectory
    * 1: Generate clean trajectory

    _Implemented in v1.3.0_

    _Updated in v1.5.0. Bugs fixed_

`use_sander`
:   Use sander for energy calculations, even when `mmpbsa_py_energy` will suffice.

    _Removed in v1.5.0: Now `sander` is used in all calculations_

`verbose` (Default = 1)
:   Specifies how much output is printed in the output file. The accepted levels are 0, 1, and 2.

    * 0: Print only difference terms
    * 1: Print all complex, receptor, ligand, and difference terms
    * 2: Print the full component terms, including the terms that cancel in a single-trajectory calculation

    _Updated in v1.5.0_

### **`&gb` namelist variables**

!!! note "Keep in mind"
    * A default GB input file can be created as follows:

        ```
        gmx_MMPBSA --create_input gb
        ```
    
    * A sample GB input file is shown [here](input_file.md#gb)
    * A tutorial on binding free energy calculation with GB model is available 
    [here](examples/Protein_ligand/ST/README.md)

#### **Basic input options**

`igb` (Default = 8)
:   Generalized Born method to use (see [§4](https://ambermd.org/doc12/Amber21.pdf#chapter.4) for more info).

    * 1: The [Hawkins, Cramer, Truhlar][191] pairwise GB model (GB-HCT)
    * 2: Modified GB model 1 developed by [A. Onufriev, D. Bashford and D.A. Case][188] (GB-OBC1)
    * 5: Modified GB model 2 developed by [A. Onufriev, D. Bashford and D.A. Case][188] (GB-OBC2)
    * 7: GBn model described by [Mongan, Simmerling, McCammon, Case and Onufriev][206] (GB-Neck)
    * 8: Same GB functional form as the GBn model (igb=7), but with different parameters. Developed by [Nguyen, Pérez, 
         Bermeo, and Simmerling][200] (GB-Neck2)

  [191]: https://pubs.acs.org/doi/10.1021/jp961710n
  [188]: https://onlinelibrary.wiley.com/doi/10.1002/prot.20033
  [206]: https://pubs.acs.org/doi/10.1021/ct600085e
  [200]: https://pubs.acs.org/doi/10.1021/acs.jctc.5b00271

`alpb` (Default = 0)
:   Use [Analytical Linearized Poisson-Boltzmann (ALPB)][209] approximation to handle electrostatic interactions 
within the implicit solvent model (see [§4.2](https://ambermd.org/doc12/Amber21.pdf#section.4.2)):

    <img src="../assets/images/alpb_Eqn.svg" align="center"/>

    where β = ε<sub>in</sub>/ε<sub>ex</sub> is the ratio of the internal and external dielectrics, α=0.571412, and A 
    is the so-called effective electrostatic size of the molecule (see `arad_method` below). The ALPB requires one 
    of the analytical GB models to be set, that is igb = 1, 2, 5, or 7, for computing the effective Born radii. It uses 
    the same sets of radii as required by the particular GB model.

    * 0: Do not use ALPB
    * 1: Use ALPB

  [209]: https://aip.scitation.org/doi/10.1063/1.1857811

`arad_method` (Default = 1)
:   Method used to estimate the effective electrostatic size/radius (`A` in ALPB equation) of the molecule 
(See [Sigalov, Fenley, and Onufriev](https://aip.scitation.org/doi/10.1063/1.2177251)).

    * 1: Use structural invariants
    * 2: Use elementary functions
    * 3: Use elliptic integral (numerical)

`intdiel` (Default = 1.0)
:   Define Internal dielectric constant.

`extdiel` (Default = 78.5)
:   Define External dielectric constant.

    _Implemented in v1.5.0: This variable `= 80.0` in previous versions_

`saltcon` (Default = 0.0)
:   Salt concentration in Molarity (M).

`surften` (Default = 0.0072)
:   Surface tension value. Units are kcal/mol/Å².

`surfoff` (Default = 0.0)
:   Offset to correct (by addition) the value of the non-polar contribution to the solvation free energy term.

`molsurf` (Default = 0)
:   Define the algorithm to calculate the surface area for the non-polar solvation term.
    
    * 0: Use LCPO (Linear Combination of Pairwise Overlaps)
    * 1: Use molsurf algorithm

`msoffset` (Default = 0) 
:   Offset to apply to the individual atomic radii in the system when calculating the `molsurf` surface. See the
description of the `molsurf` action command in [cpptraj][4].

`probe` (Default = 1.4)
:   Radius in Å of the probe molecule (supposed to be the size of a solvent molecule), to use when determining the 
molecular surface.
    
    !!! note
        Only applicable when `molsurf` is set to 1

#### **QM options**

`ifqnt` (Default = 0)
:   Specifies whether a part of the system is treated with quantum mechanics.
    
    * 0: Potential function is strictly classical
    * 1: Use QM/MM

    !!! note "Keep in mind"
        * Calculations where part of the system is treated with quantum mechanics can be performed only with GB
        * QM/MMGBSA can be combined with the explicit receptor-water workflow in single-trajectory mode
        * A sample QM/MMGBSA input file is shown [here](input_file.md#qmmmgbsa)
        * A tutorial on binding free energy calculation with QM/MMGBSA is available 
        [here](examples/QM_MMGBSA/README.md)

`qm_theory` 
:   Which semi-empirical Hamiltonian should be used for the quantum calculation. The values accepted by the current
validator are `PM3`, `AM1`, `RM1`, `MNDO`, `PM3-PDDG`, `PM3-PDDG_08`, `MNDO-PDDG`, `PM3-CARB1`, `PM3-ZNB`,
`PM3-MAIS`, `AM1-D*`, `AM1-DH+`, `MNDO/D`, `AM1/D`, `PM6`, `PM6-D`, `PM6-DH+`, `DFTB`, `DFTB2`, and `DFTB3`.
The dispersion correction can be switched on for `AM1` and `PM6` by choosing `AM1-D*` and `PM6-D`, respectively.
The dispersion and hydrogen bond correction will be applied for `AM1-DH+` and `PM6-DH+`.

    These names are case-sensitive and must be used exactly as shown. Historical aliases such as `PDDG-PM3`,
    `SCC-DFTB`, and `PM3-ZnB` are not accepted by the current QM/MMGBSA input validator. `DFTB`, `DFTB2`,
    and `DFTB3` also require the corresponding Amber parameter data to be available. Individual Hamiltonians can
    have element-specific parameter limitations; for example, PM3-MAIS does not provide parameters for every element.

    The following is a limited smoke-test matrix, not a guarantee of universal support. The current development
    checkout was tested with a one-frame QM/MMGBSA run using the prepared test system. The following values completed
    successfully in the configured SANDER environment: `PM3`, `AM1`, `RM1`, `MNDO`,
    `PM3-PDDG`, `PM3-PDDG_08`, `MNDO-PDDG`, `PM3-CARB1`, `PM3-ZNB`, `AM1-D*`, `AM1-DH+`, `MNDO/D`, `AM1/D`,
    `PM6`, `PM6-D`, `PM6-DH+`, `DFTB`, `DFTB2`, and `DFTB3`. The DFTB variants required their Amber parameter data
    to be available. `PM3-MAIS` is accepted by the validator but failed for the tested QM region because SANDER
    reported missing PM3-MAIS parameters for nitrogen; it should therefore be treated as element-limited rather
    than universally runnable. `gmx_MMPBSA` validates the method name and writes the QM/MM input; SANDER performs
    the quantum calculation. Actual availability therefore depends on the SANDER build, Amber parameter data,
    QM-region elements and charge, and SCF convergence. A successful validator check or one successful smoke test
    does not guarantee that a Hamiltonian will run for every user-defined QM region.

    The default is `PM6-DH+` which includes dispersion and hydrogen-bond corrections missing from plain PM3/PM6,and have been used in the past to study protein-ligand interactions. 
    (see [Řezáč & Hobza, *JCTC* **2009**, 5, 1749](https://doi.org/10.1021/ct9000922); 
    [Korth, *JCTC* **2010**, 6, 3808](https://doi.org/10.1021/ct100408b)). 
    For QM/MMGBSA applications, see also 
    [Thapa *et al.*, *J. Phys. Chem. B* **2018**, 122, 7866](https://doi.org/10.1021/acs.jpcb.8b03655) and 
    [*Commun. Biol.* **2025**](https://doi.org/10.1038/s42003-025-09143-z).

    If `qm_theory` is omitted while `ifqnt` = 1, it defaults to `PM6-DH+`.

    gmx_MMPBSA checks the SANDER output before accepting a QM/MM energy calculation. Fatal diagnostics such as SCF
    nonconvergence, missing Hamiltonian parameters, missing dispersion-correction parameters, or missing DFTB
    Slater-Koster files stop the calculation and identify the mdout file, diagnostic class, and suggested remedy.
    The numerical-derivative message for d orbitals is reported as a performance warning because SANDER can still
    complete the calculation with numerical derivatives. A final binding result is not considered valid when a
    QM/MM energy component contains one of the fatal diagnostics.

`qm_residues`
:   Complex residues to treat with quantum mechanics. All residues treated with quantum mechanics in the complex 
must be treated with quantum mechanics in the receptor or ligand to obtain meaningful results. This notation is 
the same used for `print_res` variable in `&decomp` namelist.

    !!! danger
         No `qm_residues` default, this must be specified if `ifqnt` = 1 except when using user defined masks.

    !!! example "Selection schemes"

        === "By Distance (recommended)"
            Notation: [ `within` `distance` ]
            :   `within` corresponds to the keyword and `distance` to the maximum distance criterion in Å necessary to 
                select the residues from both the receptor and the ligand. In case the cutoff used is so small that 
                the number of qm_residues = 0, the cutoff value will be increased by 0.1 until the number of 
                qm_residues > 0.
    
            !!! example
                `qm_residues="within 5"` treats receptor and ligand residues within 5 Å of one another with quantum
                mechanics.

        === "Amino acid selection"
            Notation: [ `CHAIN`/(`RESNUM`, `RESNUM` + `INSERTION_CODE`, or `RESNUM-RESNUM`) ]
            :    Treat individual residues or residue ranges with quantum mechanics. An insertion code is appended directly
            to the residue number, for example `A/27B`; the colon form `A/27:B` is invalid. Numeric ranges select only
            residues with a blank insertion code, so insertion-coded residues must be listed individually.

            `qm_residues="A/1,3-10,15,100"` treats chain A residues 1, 3 through 10, 15, and
            100 from the complex topology file and the corresponding residues in either the ligand and/or receptor 
            topology files.
    
            Let's suppose that we can have the following sequence: - A:LEU:5 - A:GLY:6:A - A:THR:6:B - A:SER:6:C - 
            A:ASP:6:D - A:ILE:7
    
            with the format `CHAIN`/`RESNUMBER` followed immediately by `INSERTION_CODE`
            
            === "Right notation"
                
                **Ranges selection**
                :   `qm_residues="A/5-7"` selects only the blank-code residues A/5 and A/7. Numeric ranges do not include
                    insertion-coded residues such as A/6A through A/6D.
                
                **Individual selection**
                :   `qm_residues="A/5,6A,6B,6C,6D,7"` treats all listed residues with quantum mechanics, including
                    every insertion-coded residue in the sequence.
                
                **Multiple chain selection**
                :   `qm_residues="A/5-10,100 B/34,56"` treats residues 5 through 10 and 100 from
                chain A, and residues 34 and 56 from Chain B.
    
            === "Wrong notation"
                `qm_residues="A/27:B"` is invalid; use `qm_residues="A/27B"` for an insertion code.

                `qm_residues="A/5-6B,6D-7"` produces an error.

`com_qmmask` (Default = '')
:   Amber mask specifying the quantum atoms in the complex. If defined, residues selected with `qm_residues` variable are 
ignored. When using user defined masks, `com_qmmask`, `rec_qmmask`, and `lig_qmmask` must be defined.
    
    !!! danger
        When using user-defined masks, automatic assignment of `qmcharge_com` is disabled and the default or
        user-defined `qmcharge_com` value is used.

`rec_qmmask` (Default = '')
:   Amber mask specifying the quantum atoms in the receptor. When using user defined masks, `com_qmmask`, `rec_qmmask`, 
and `lig_qmmask` must be defined.
    
    !!! danger
        When using user-defined masks, automatic assignment of `qmcharge_rec` is disabled and the default or
        user-defined `qmcharge_rec` value is used.

`lig_qmmask` (Default = '')
:   Amber mask specifying the quantum atoms in the ligand. When using user defined masks, `com_qmmask`, `rec_qmmask`,
and `lig_qmmask` must be defined.
    
    !!! danger
        When using user-defined masks, automatic assignment of `qmcharge_lig` is disabled and the default or
        user-defined `qmcharge_lig` value is used.

`qmcharge_com` (Default = 0)
:   Charge of the quantum region in the complex. `qmcharge_com` is assigned automatically from `qm_residues`. When
using user-defined masks, automatic assignment is disabled and the default or user-defined `qmcharge_com` value is used.

`qmcharge_rec` (Default = 0)
:   Charge of the quantum region in the receptor. `qmcharge_rec` is assigned automatically from `qm_residues`. When
using user-defined masks, automatic assignment is disabled and the default or user-defined `qmcharge_rec` value is used.

`qmcharge_lig` (Default = 0)
:   Charge of the quantum region in the ligand. `qmcharge_lig` is assigned automatically from `qm_residues`. When
using user-defined masks, automatic assignment is disabled and the default or user-defined `qmcharge_lig` value is used.

`qmcut` (Default = 9999.0)
:   The QM/MM charge-interaction cutoff in Å.

`scfconv` (Default = 1.0e-8)
:   Controls the convergence criteria for the SCF calculation, in kcal/mol. The tighter the 
convergence the longer the calculation will take. Values tighter than 1.0e-11 are not recommended as these can lead 
to oscillations in the SCF, due to limitations in machine precision, that can lead to convergence failures.

`itrmax` (Default = 1000)
:   Maximum number of SCF iterations allowed for each QM/MM step. If SANDER reports an unconverged SCF step,
the calculation stops rather than including the unconverged energy in the binding-energy result. Increase `itrmax`
or adjust `scfconv` only after checking the QM/MM output for stable convergence.

`ndiis_attempts` (Default = None)
:   Optional maximum number of DIIS attempts used by SANDER during each QM/MM SCF calculation. `None` leaves the
SANDER default unchanged; explicit values from 0 to 1000 are accepted. This option is not available for DFTB, DFTB2,
or DFTB3. Set it only when the QM/MM output shows repeatable SCF convergence difficulty; increasing it can increase
runtime. For example, `ndiis_attempts=700` recovered convergence for the Fig3 test system. It does not override the
fatal-diagnostic check: unconverged frames are still rejected.

`writepdb` (Default = 1)
:   Write the atoms selected for the QM region to `qmmm_region.pdb` during the first step. This file makes it easier to
verify which atoms were included in the QM region.

    * 0: Do not write a PDB file of the selected QM region
    * 1: Write a PDB file of the selected QM region

`peptide_corr` (Default = 0)
:   Apply MM correction to peptide linkages. This correction is of the form: 

[comment]: <> (    <img src="https://latex.codecogs.com/svg.)

[comment]: <> (    image?E_{scf}&space;=&space;E_{scf}&space;&plus;&space;h_{type}&#40;i_{type}&#41;)

[comment]: <> (    *sin^{2}\phi" title="https://latex.codecogs.com/svg.image?E_{scf} = E_{scf} + h_{type}&#40;i_{type}&#41;*sin^{2}\phi )

[comment]: <> (    align="center""/>)

    <img src="../assets/images/peptide_correction.svg" align="center"/>

    where _ϕ_ is the dihedral angle of the H-N-C-O linkage and h<sub>type</sub> is a constant dependent on the 
    Hamiltonian used. Recommended, except for DFTB, DFTB2, and DFTB3.

    * 0: Do not apply a MM correction to peptide linkages
    * 1: Apply a MM correction to peptide linkages

`verbosity` (Default = 0)
:   Controls the verbosity of QM/MM related output. Values of 2 or higher will produce a lot of output.

    * 0: only minimal information is printed - Initial QM geometry and link atom positions as
    well as the SCF energy at every ntpr steps.
    * 1: Print SCF energy at every step to many more significant figures than usual. Also print the
    number of SCF cycles needed on each step.
    * 2: As 1 and also print info about memory reallocations, number of pairs per QM atom, QM core -
    QM core energy, QM core - MM atom energy, and total energy.
    * 3: As 2 and also print SCF convergence information at every step.
    * 4: As 3 and also print forces on the QM atoms due to the SCF calculation and the coordinates of
    the link atoms at every step.
    * 5: As 4 and also print all of the info in kJ/mol as well as kcal/mol.
    
  [4]: https://ambermd.org/doc12/Amber21.pdf#subsection.34.11.49

### **`&gbnsr6` namelist variables**

!!! note "Keep in mind"
    * GBNSR6 is an implementation of the Generalized Born (GB) model in which the effective Born radii are computed 
    numerically, via the so-called "R6" integration ([ref.][222]) over molecular surface of the solute. In contrast to 
    most GB practical models, GBNSR6 model is parameter free in the same sense as the numerical PB framework is. Thus, 
    accuracy of GBNSR6 relative to the PB standard is virtually unaffected by the choice of input atomic radii. Check
    Chapter [§5](https://ambermd.org/doc12/Amber21.pdf#chapter.5) in Amber manual for a more thorough description 
    of the GBNSR6 model and its parameters.
    * `gmx_MMPBSA` prepares temporary GBNSR6 topology copies for the energy calculation and keeps the original
    COM/REC/LIG topology files for output parsing. This compatibility step is automatic and does not require extra
    input options.
    * A default GBNSR6 input file can be created as follows:

        ```
        gmx_MMPBSA --create_input gbnsr6
        ```
    
    * A sample GBNSR6 input file is shown [here](input_file.md#gbnsr6)
    * A tutorial on binding free energy calculation with GBNSR6 model is available 
    [here](examples/GBNSR6/README.md)

    !!! warning "GBNSR6 + decomposition is hybrid"
        Each frame runs **sander (MM)** and **`gbnsr6`**, then merges the outputs.

        * **System totals:** `EEL`, `1-4 EEL`, and `EGB` come from GBNSR6.
        * **Decomposition columns:** `internal`, `vdw`, **`eel`**, and `sas` come from sander; only **`pol`**
          is replaced with GBNSR6 residue/pair energies from `DGij`. Amber GBNSR6 does not provide a
          sander-style gas-phase Coulomb decomposition, so residue `eel` is not a per-residue split of the
          GBNSR6 total `EEL`.
        * `epsin` is copied to sander `intdiel` (and related MM settings) so the sander Coulomb terms use the
          same solute dielectric as GBNSR6.

        Do not treat decomp `eel + pol` as a pure GBNSR6 electrostatic partition when comparing to totals.

  [222]: https://pubs.acs.org/doi/abs/10.1021/ct200786m

#### **Basic input options**

`epsin` (Default = 1.0)
:   Dielectric constant of the solute region. When decomposition is enabled, this value is also passed to the
sander MM input as `intdiel` so gas-phase Coulomb (`eel`) uses the same dielectric as GBNSR6.

`epsout` (Default = 78.5)
:   Implicit solvent dielectric constant for the solvent.

`istrng` (Default = 0.0)
:   Ionic strength in M for the GBNSR6 equation.
                           
`dprob` (Default = 1.4)
:   Radius of the solvent probe.

`cavity_surften` (Default = 0.005)
:   Surface tension parameter for nonpolar solvation calculation.

#### **Options to select numerical procedures**

`space` (Default = 0.5)
:   Sets the grid spacing that determines the resolution of the solute molecular surface. Note that memory footprint of 
this grid-based implementation of GBNSR6 may become large for large structures, e.g. the nucleosome (about 25,000 
atoms) will take close to 2 GB of RAM when the default grid spacing is used. For very large structures, one may 
consider increasing the value of space, which will reduce the memory footprint and execution time; however, the 
accuracy will also decrease.

`arcres` (Default = 0.2)
:   Arc resolution used for numerical integration over molecular surface.

`b` (Default = 0.028)
:   Specifies the value of uniform offset to the (inverse) effective radii, the default value 0.028 gives 
better agreement with the PB model, regardless of the structure size. For best agreement with the explicit solvent 
(TIP3P) solvation energies, optimal value of B depends on the structure size: for small molecules (number of atoms 
less than 50), B=0 is recommended. With -chagb option, B is calculated automatically based on the solute size.

`alpb` (Default = 1)
:   Specifies if ALPB correction is to be used.

    * 0: Use canonical GB
    * 1: Use ALPB approximation (default)

#### **Options for CHAGB model**

`chagb` (Default = 0)
:   Define if CHAGB is used.

    * 0: Do not use CHAGB.
    * 1: Use CHAGB.

`rs` (Default = 0.52)
:   Dielectric boundary shift compared to the molecular surface.

`radiopt` (Default = 0)
:   Set of intrinsic atomic radii to be used.

    * 0: uses hardcoded intrisic radii optimized for small drug like molecules, and single amino acid
    dipeptides ([ref.][215])
    * 1: intrinsic radii are read from the topology file. Note that the dielectric surface defined using
    these radii is then shifted outwards by Rs relative to the molecular surface. The option is not
    recommended unless you are planning to re-optimize the input radii set for your problem.

  [215]: https://pubs.acs.org/doi/full/10.1021/ct4010917

`roh` (Default = 1)
:   Sets the value of R<sup>z</sup><sub>OH</sub> for CHAGB model, the default is 0.586Å. This parameter defines which 
explicit water model is being mimicked with respect to its propensity to cause charge hydration asymmetry. A perfectly 
tetrahedral water, which cannot cause charge-hydration asymmetry, would have R<sup>z</sup><sub>OH</sub> = 0. The
options for `roh` are:

    * 1: R<sup>z</sup><sub>OH</sub> = 0.586Å corresponds to TIP3P and SPC/E. 
    * 2: R<sup>z</sup><sub>OH</sub> = 0.699Å for OPC.
    * 3: R<sup>z</sup><sub>OH</sub> = 0.734Å for TIP4P 
    * 4: R<sup>z</sup><sub>OH</sub> = 0.183Å for TIP5P/E. 

`tau` (Default = 1.47)
:   Value of τ in the CHAGB model. This dimensionless parameter controls the effective range of the neighboring 
charges (_j_) affecting the CHA of atom (_i_), see ([ref.][215]) for details.

### **`&pb` namelist variables**

!!! note "Keep in mind"
    * `gmx_MMPBSA` uses **_sander_** to perform PB calculations. **_sander_** offers access to all [pbsa][5] 
    functionalities. The default values for the variables described below are appropriate for most calculations 
    on solvated molecular systems. Also note that the default options may have changed over time. A more thorough 
    description of all the options can be found [here][5]. For a detailed discussion of all related options on 
    the quality of the MM/PB(GB)SA calculations, please check this [publication][6].
    * A default PB input file can be created as follows:

        ```
        gmx_MMPBSA --create_input pb
        ```
    
    * A sample PB input file is shown [here](input_file.md#pb)
    * A tutorial on binding free energy calculation with PB model is available 
    [here](examples/Linear_PB_solver/README.md)

#### **Basic input options**

`ipb` (Default = 2)
:   Option to set up a dielectric model for all numerical PB procedures. `ipb = 1` corresponds to a classical geometric 
method, while a level-set based algebraic method is used when `ipb > 2`.

    * 0: No electrostatic solvation free energy is computed.
    * 1: The dielectric interface between solvent and solute is built with a geometric approach. ([ref.][217])
    * 2: The dielectric interface is implemented with the level set function. Use of a level set function
    simplifies the calculation of the intersection points of the molecular surface and grid edges and
    leads to more stable numerical calculations. ([ref.][239])
    * 4: The dielectric interface is also implemented with the level set function. However, the linear
    equations on the grid points nearby the dielectric boundary are constructed using the IIM. In this
    option, The dielectric constant do not need to be smoothed, that is, `smoothopt` is useless.
    Only the linear PB equation is supported, that is, `npbopt = 0`. Starting from the Amber 2018
    release, `solvopt` is no longer relevant as only one stable solver is supported. ([ref.][233])
    * 6: The dielectric interface is implemented analytically with the revised density function approach
    (`sasopt = 2`). The linear equations on the irregular points are constructed using the IIM and
    fully utilizing the analytical surface. Otherwise, it is exactly the same as `ipb = 4`. ([ref.][240])
    * 7: The dielectric interface is implemented analytically with the revised density function approach
    (`sasopt = 2`). The linear equations on the irregular points are constructed using the Χ-factor
    harmonic average method. ([ref.][241])
    * 8: The dielectric interface is implemented analytically with the revised density function approach
    (`sasopt = 2`). The linear equations on the irregular points are constructed using the secondorder harmonic 
    average method. ([ref.][241])

  [217]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.10120
  [239]: https://pubs.acs.org/doi/10.1021/ct300341d
  [233]: https://www.sciencedirect.com/science/article/abs/pii/S0009261408016539?via%3Dihub
  [240]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.25783
  [241]: https://pubs.acs.org/doi/10.1021/acs.jctc.9b00602

`inp` (Default = 1)
:   Option to select different methods to compute non-polar solvation free energy.

    * 1: The total non-polar solvation free energy is modeled as a single term linearly proportional to the
    solvent accessible surface area ([ref.][227]). When using `inp = 1`, `gmx_MMPBSA` aligns the INPUT
    nonpolar parameters with Amber PBSA before writing mdins and parsing decomposition:

        * `sprob` is reset to 1.4
        * `cavity_surften` is reset to 0.005
        * `cavity_offset` is reset to 0.000
        * `radiopt` is reset to 0 (use prmtop / `PBRadii` radii)

    This keeps `_GMXMMPBSA_info` and decomp SAS scaling consistent with the γ Amber actually used for
    ENPOLAR/ECAVITY (Amber also applies the cavity/sprob reset at runtime; `radiopt` is only warned about
    unless gmx_MMPBSA aligns it).

    * 2: The total non-polar solvation free energy is modeled as two terms: the cavity term and the
    dispersion term. The dispersion term is computed with a surface-based integration method
    ([ref.][227]) closely related to the PCM solvent for quantum chemical programs. ([ref.][229]) Under this
    framework, the cavity term is still computed as a term linearly proportional to the molecular
    solvent-accessible-surface area (SASA) or the molecular volume enclosed by SASA.

    !!! info "Keep in mind"
        Sometimes, high values for the solvation energy are obtained using `inp=2`. Check 
        this [section](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/Q%26A/calculations/#possible-solutions_2) to 
        see a workaround.

  [227]: https://pubs.acs.org/doi/abs/10.1021/jp073399n
  [229]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.540100504

`sander_apbs` (Default = 0)
:   Deprecated legacy compatibility option. New calculations reject `sander_apbs = 1`; keep it at `0` and use the
    built-in `PBSA` solver. Historical APBS result files remain readable where the analyzer still needs this value.

#### **Options to define the physical constants**

`indi` (Default = 1.0)
:   Internal dielectric constant. This corresponds to `epsin` in [pbsa][5].

`exdi` (Default = 78.5)
:   External dielectric constant. This corresponds to `epsout` in [pbsa][5].

`emem` (Default = 4.0)
:   Sets the membrane dielectric constant. Only used if `memopt` > 0, does nothing otherwise. Value
used should be between `indi` and `exdi` or there may be errors. This corresponds to `epsmem` in [pbsa][5].

`smoothopt` (Default = 1)
:   Instructs PB how to set up dielectric values for finite-difference grid edges that are located across the
solute/solvent dielectric boundary.

    * 0: The dielectric constants of the boundary grid edges are always set to the equal-weight harmonic
    average of `indi` and `exdi`.
    * 1: A weighted harmonic average of `indi` and `exdi` is used for boundary grid edges. The
    weights for `indi` and `exdi` are fractions of the boundary grid edges that are inside or
    outside the solute surface. ([ref.][243])
    * 2: The dielectric constants of the boundary grid edges are set to either `indi` or `exdi` depending on whether 
    the midpoints of the grid edges are inside or outside the solute surface.

  [243]: https://pubs.acs.org/doi/10.1021/cr00101a005

    _Implemented in v1.5.0_

`istrng` (Default = 0.0)
:   Ionic strength in Molarity (M). It is converted to mM for `PBSA`.

`radiopt` (Default = 0)
:   The option to set up atomic radii.

    * 0: Use radii from the prmtop file for both the PB calculation and for the non-polar calculation (see `inp`).
      This is required when `inp = 1` (Amber PBSA warns otherwise and still runs Tan–Luo radii if `radiopt = 1`).
      Topology conversion (`-cp` / amber conversion) also forces `radiopt = 0` so the PB solver uses the
      `PBRadii` / ChRad radii written into the generated prmtops.
    * 1: Use atom-type/charge-based radii by Tan and Luo ([ref.][244]) for the PB calculation. Note that the
    radii are optimized for Amber atom types as in standard residues from the Amber database and should work fine for
    `standard` complexes such as protein-protein, protein-DNA. On the other hand, if a molecule in your system was 
    built by antechamber, _i.e._, if GAFF atom types are used, or any other extrenal software, radii from the prmtop 
    file should be used (`radiopt = 0`). Check this [thread](http://archive.ambermd.org/201303/0548.html) for more info.

    When `inp = 1`, `gmx_MMPBSA` aligns `radiopt` to `0` during input processing if it was left at `1`.

  [244]: https://pubs.acs.org/doi/abs/10.1021/jp063479b

`prbrad` (Default = 1.4)
:   Solvent probe radius (in Å). Allowed values are 1.4 and 1.6. This corresponds to `dprob` in [pbsa][5].

`iprob` (Default = 2.0)
:   Mobile ion probe radius (in Å) for ion accessible surface used to define the Stern layer.

    _Implemented in v1.5.0_

`sasopt` (Default = 0)
:   Option to determine which kind of molecular surfaces to be used in the Poisson-Boltzmann implicit solvent model.

    * 0: Use the solvent excluded surface as implemented by ([ref.][239])
    * 1: Use the solvent accessible surface. Apparently, this reduces to the van der Waals surface when
    the `prbrad` is set to zero.
    * 2: Use the smooth surface defined by a revised density function. ([ref.][245]) This must be combined with
    `ipb > 2.

  [245]: https://pubs.acs.org/doi/10.1021/ct900318u

`arcres` (Default = 0.25)
:   The `arcres` keyword gives the resolution (in Å) of dots used to represent solvent accessible arcs. More
generally, `arcres` should be set to max(0.125 Å, 0.5h) (h is the grid spacing). ([ref.][239])

    _Implemented in v1.5.0_

#### **Options for implicit membranes**

`memopt` (Default = 0)
:   Option to turn the implicit membrane on and off. The membrane is implemented as a slab like region with a uniform 
or heterogeneous dielectric constant depth profile. Details of the implicit membrane setup can be 
found [here](https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00382).

    * 0: No implicit membrane used.
    * 1: Use a uniform membrane dielectric constant in a slab-like implicit membrane. ([ref.][246])
    * 2: Use a heterogeneous membrane dielectric constant in a slab-like implicit membrane. The dielectric constant 
    varies with depth from a value of 1 in the membrane center to 80 at the membrane
    periphery. The dielectric constant depth profile was implemented using the PCHIP fitting. ([ref.][247])
    * 3: Use a heterogeneous membrane dielectric constant in a slab-like implicit membrane. The dielectric constant 
    varies with depth from a value of 1 in the membrane center to 80 at the membrane periphery. The dielectric constant 
    depth profile was implemented using the Spline fitting. ([ref.][247])

    !!! note "Keep in mind"
        * Calculations for implicit membranes can be performed only with PB
        * A sample input file is shown [here](input_file.md#mmpbsa-with-membrane-proteins)
        * A tutorial on binding free energy calculation for membrane proteins is available 
        [here](examples/Protein_membrane/README.md)
        * Check this thread for more info on [Parameters for Implicit 
        Membranes](http://archive.ambermd.org/202006/0088.html)

  [246]: https://www.sciencedirect.com/science/article/abs/pii/S0009261412012808?via%3Dihub
  [247]: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00363

`mprob` (Default = 2.70)
:   Membrane probe radius (in Å). This is used to specify the highly different lipid molecule accessibility versus 
that of the water. ([ref.][248])

  [248]: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00382

    _Implemented in v1.5.0_

[<img src="../assets/prot_memb.png" height="200" width="258" align="right"/>]()

`mthick` (Default = `automatic`)
:   Membrane thickness (in Å), or `automatic` to calculate it from the selected membrane atom coordinates in `-ct`.
    The automatic value is the mean z coordinate of the upper leaflet minus the mean z coordinate of the lower
    leaflet. A numeric value can be supplied independently of `mctrdz`.

`mctrdz` (Default = `automatic`)
:   Membrane center (in Å) in the z direction, or `automatic` to use the mean z coordinate of the selected
    membrane atoms in `-ct`. A numeric value can be supplied independently of `mthick`.

`membrane_atoms` (Default = `P`)
:   Semicolon-separated atom names used for automatic membrane parameters, for example `membrane_atoms="P;N"`.
    These are atom names, not cpptraj masks, and must be present in the full structure/trajectory supplied through
    `-cs`/`-ct`. For GROMACS inputs, automatic detection reads the original unstripped `-ct` trajectory, so the
    membrane does not need to be included in the receptor group given through `-cg`. The automatic calculation uses
    only `-ct`, including the selected `startframe`, `endframe`, and `interval`.

    Explicit lipids are retained only when they are intentionally included in the selected `-cg` receptor group;
    their residue names should not be changed.

    When automatic calculation is used, `GMXMMPBSA_membrane_parameters.csv` and
    `GMXMMPBSA_membrane_parameters.png` are retained regardless of `keep_files`. The CSV contains per-frame
    diagnostics and the resolved values; the PNG shows the selected atom z-coordinate distributions, the two
    leaflets, the resolved center and slab boundaries, and frame-wise center/thickness stability.
    The trajectory must already be oriented with the membrane normal along z and should be made continuous across
    periodic boundaries before running the calculation; automatic detection does not reorient or unwrap `-ct`.

`poretype` (Default = 1)
:   Turn on and off the automatic depth-first search method to identify the pore. ([ref.][248])

    * 0: Do not turn on the pore searching algorithm.
    * 1: Turn on the pore searching algorithm.

#### **Options to select numerical procedures**

`npbopt` (Default = 0) {#npbopt}
:   Option to select the linear, or the full nonlinear PB equation.

    * 0: Linear PB equation (LPBE) is solved
    * 1: Nonlinear PB equation (NLPBE) is solved

    !!! note
        While the linear PB equation (see [tutorial](examples/Linear_PB_solver/README.md)) will suffice for most
        calculations, the nonlinear PB equation (see [tutorial](examples/NonLinear_PB_solver/README.md)) is recommended
        for highly charged systems. Take into account that NLPBE can significantly increase PB wall time.

        Amber PBSA **forces `eneopt = 1` (P3M)** when `npbopt = 1` and requires a nonzero `cutnb`. That is Amber's
        bookkeeping for nonlinear PB, not a gmx_MMPBSA preference: `EPB` is reported as zero and reaction-field plus
        Coulombic energy are combined in `EEL`, so ΔGGAS / ΔGSOLV are not a meaningful gas/solv split. Use ΔTOTAL.
        See `eneopt` below and Amber threads ([T1](http://archive.ambermd.org/201203/0191.html),
        [T2](http://archive.ambermd.org/201610/0114.html)).

    _Implemented in v1.5.0_

`solvopt` (Default = 1)
:   Option to select iterative solvers.

    * 1 Modified ICCG or Periodic (PICCG) if `bcopt = 10`.
    * 2 Geometric multigrid. A four-level v-cycle implementation is applied by default.
    * 3 Conjugate gradient (Periodic version available under `bcopt = 10`). This option requires a large
    `linit` to converge.
    * 4 SOR. This option requires a large `linit` to converge.
    * 5 Adaptive SOR. This is only compatible with `npbopt = 1`. This option requires a large `linit` 
    converge. ([ref.][219])
    * 6 Damped SOR. This is only compatible with `npbopt = 1`. This option requires a large `linit` to 
    converge. ([ref.][219])

  [219]: https://pubs.acs.org/doi/10.1021/ct900381r

`accept` (Default = 0.001)
:   Sets the iteration convergence criterion (relative to the initial residue).

    _Implemented in v1.5.0_

`linit` (Default = 1000) 
:   Sets the maximum number of iterations for the finite difference solvers. Note that `linit` has to be set to a 
much larger value, _e.g._ 10000, for the less efficient solvers, such as conjugate gradient and SOR, to converge. 
This corresponds to `maxitn` in [pbsa][5].

`fillratio` (Default = 4.0) 
:   The ratio between the longest dimension of the rectangular finite-difference grid and that of the solute. For 
macromolecules is fine to use 4, or a smaller value like 2. A default value of 4 is large enough to be used for a 
small solute, such as a ligand molecule. Using a smaller value for `fillratio` may cause part of the small solute 
to lie outside the finite-difference grid, causing the finite-difference solvers to fail. 

`scale` (Default = 2.0)
:   Resolution of the Poisson Boltzmann grid. It is equal to the reciprocal of the grid spacing (`space` in [pbsa][5]).

`nbuffer` (Default = 0)
:   Sets how far away (in grid units) the boundary of the finite difference grid is away from the solute
surface; _i.e._, automatically set to be at least a solvent probe or ion probe (diameter) away from the solute surface.

    _Implemented in v1.5.0_

`nfocus` (Default = 2)
:   Set how many successive FD calculations will be used to perform an electrostatic focussing calculation on a 
molecule. When `nfocus` = 1, no focusing is used. It is recommended that `nfocus = 1` when the multigrid solver is used.

`fscale` (Default = 8)
:   Set the ratio between the coarse and fine grid spacings in an electrostatic focussing calculation.

    _Implemented in v1.5.0_

`npbgrid` (Default = 1)
:   Sets how often the finite-difference grid is regenerated.

    _Implemented in v1.5.0_

#### **Options to compute energy and forces**

`bcopt` (Default = 5)
:   Boundary condition options.

    * 1: Boundary grid potentials are set as zero, _i.e._ conductor. Total electrostatic potentials and energy
    are computed.
    * 5: Computation of boundary grid potentials using all grid charges. Total electrostatic potentials
    and energy are computed.
    * 6: Computation of boundary grid potentials using all grid charges. Reaction field potentials and
    energy are computed with the charge singularity free formalism. ([ref.][236])
    * 10: Periodic boundary condition is used. Total electrostatic potentials and energy are computed.
    Can be used with `solvopt = 1, 2, 3, or 4` and `ipb = 1 or 2`. It should only be used on charge-neutral 
    systems. If the system net charge is detected to be nonzero, it will be neutralized by
    applying a small neutralizing charge on each grid (_i.e._ a uniform plasma) before solving.

  [236]: https://aip.scitation.org/doi/abs/10.1063/1.3099708

`eneopt` (Default = 2) {#eneopt}
:   Option to compute total electrostatic energy and forces.

    * 1: Compute total electrostatic energy and forces with the particle-particle particle-mesh (P3M)
    procedure outlined in Lu and Luo. ([ref.][223]) In doing so, energy term EPB in the output file is set
    to zero, while EEL includes both the reaction field energy and the Coulombic energy. The van
    der Waals energy is computed along with the particle-particle portion of the Coulombic energy.
    The electrostatic forces and dielectric boundary forces can also be computed. ([ref.][223]) This option
    requires a nonzero `cutnb` and `bcopt = 5` for soluble proteins / `bcopt = 10` for membrane proteins.
    * 2: Use dielectric boundary surface charges to compute the reaction field energy. Both
    the Coulombic energy and the van der Waals energy are computed via summation of pairwise
    atomic interactions. Energy term EPB in the output file is the reaction field energy. EEL is the
    Coulombic energy.
    * 3: Similar to the first option above, a P3M procedure is applied for both solvation and Coulombic
    energy and forces for larger systems.
    * 4: Similar to the third option above, a P3M procedure for the full nonlinear PB equation is applied
    for both solvation and Coulombic energy and forces for larger systems. A more robust and
    clean set of routines were used for the P3M and dielectric surface force calculations.

    !!! warning "Reporting with `eneopt = 1` (and similar P3M modes 3/4)"
        `EPB ≈ 0` and a larger `EEL` are expected Amber output, not a missing polar term. gmx_MMPBSA still
        places `EEL` under GGAS and `EPB` under GSOLV, so ΔGGAS / ΔGSOLV are **not** the usual gas/solv
        partition; ΔTOTAL remains valid. The same applies when Amber resets `eneopt` to `1` for
        `npbopt = 1` (NLPBE). Linear PB with `eneopt = 1` (for example some membrane setups) has the same
        reporting implication without NLPB.

  [223]: https://aip.scitation.org/doi/10.1063/1.1622376

`frcopt` (Default = 0)
:   Option to compute and output electrostatic forces to a file named force.dat in the working directory.

    * 0: Do not compute or output atomic and total electrostatic forces.
    * 1: Reaction field forces are computed by trilinear interpolation. Dielectric boundary forces are
    computed using the electric field on dielectric boundary. The forces are output in the unit of
    kcal/mol·Å.
    * 2: Use dielectric boundary surface polarized charges to compute the reaction field forces and dielectric 
    boundary forces ([ref.][237]) The forces are output in the unit of kcal/mol·Å.
    * 3: Reaction field forces are computed using dielectric boundary polarized charge. Dielectric boundary forces 
    are computed using the electric field on dielectric boundary. ([ref.][249]) The forces are output in kcal/mol·Å.

  [237]: https://www.sciencedirect.com/science/article/abs/pii/S0009261411010487?via%3Dihub
  [249]: https://pubs.rsc.org/en/content/articlelanding/2012/cp/c2cp43237d

`scalec` (Default = 0)
:   Option to compute reaction field energy and forces.

    * 0: Do not scale dielectric boundary surface charges before computing reaction field energy and
    forces.
    * 1: Scale dielectric boundary surface charges using Gauss’s law before computing reaction field
    energy and forces.

    _Implemented in v1.5.0_

`cutfd` (Default = 5.0)
:   Atom-based cutoff distance to remove short-range finite-difference interactions, and to add pairwise
charge-based interactions. This is used for both energy and force calculations. See Eqn (20) in 
Lu and Luo. ([ref.][223])

`cutnb` (Default = 0.0)
:   Atom-based cutoff distance for van der Waals interactions, and pairwise Coulombic interactions when `eneopt` = 2.
When `cutnb` is set to the default value of 0, no cutoff will be used for van der Waals and Coulombic interactions, 
_i.e._, all pairwise interactions will be included. When `eneopt = 1`, this is the cutoff distance used for van der 
Waals interactions only. The particle-particle portion of the Coulombic interactions is computed with the cutoff of 
`cutfd`._

`nsnba` (Default = 1)
:   Sets how often (steps) atom-based pairlist is generated.

    _Implemented in v1.5.0_

#### **Options to select a non-polar solvation treatment**

`decompopt` (Default = 2)
:   Option to select different decomposition schemes when `inp = 2`. See ([ref.][227]) for a detailed discussion
of the different schemes. The _σ_ decomposition scheme is the best of the three schemes studied. ([ref.][227]) As 
discussed in ([ref.][227]), `decompopt = 1` is not a very accurate approach even if it is more straightforward to 
understand the decomposition.

    * 1: Use the 6/12 decomposition scheme
    * 2: Use the _σ_ decomposition scheme
    * 3: Use the WCA decomposition scheme

    _Implemented in v1.5.0_

`use_rmin` (Default = 1)
:   The option to set up van der Waals radii. The default is to use van der Waals _rmin_ to improve the agreement with
TIP3P. ([ref.][227])

    * 0: Use atomic van der Waals _σ_ values.
    * 1: Use atomic van der Waals _rmin_ values.

    _Implemented in v1.5.0_

`sprob` (Default = 0.557)
:   Solvent probe radius (in Å) for solvent accessible surface area (SASA) used to compute the dispersion term,
default to 0.557 Å in the _σ_ decomposition scheme as optimized in ([ref.][227]) with respect to the
TIP3P solvent and the PME treatment. Recommended values for other decomposition schemes can
be found in Table 4 of ([ref.][227]). If `use_sav = 0` (see below), `sprob` can be used to compute SASA
for the cavity term as well. Unfortunately, the recommended value is different from that used in the
dispersion term calculation as documented in ([ref.][227]). Thus, two separate calculations are
needed when `use_sav = 0`, one for the dispersion term and one for the cavity term. Therefore,
please carefully read ([ref.][227]) before proceeding with the option of `use_sav = 0`. Note that
`sprob` was used for ALL three terms of solvation free energies, _i.e._, electrostatic, attractive, and
repulsive terms in previous releases in Amber. However, it was found in the more recent study ([ref.][227])
that it was impossible to use the same probe radii for all three terms after each term was calibrated
and validated with respect to the TIP3P solvent. ([ref.][227])

    _Implemented in v1.5.0_

`vprob` (Default = 1.300)
:   Solvent probe radius (in Å) for molecular volume (the volume enclosed by SASA) used to compute non-polar cavity 
solvation free energy, default to 1.300 Å, the value optimized in ([ref.][227]) with respect to the TIP3P solvent. 
Recommended values for other decomposition schemes can be found in Tables 1-3 of ([ref.][227]).

    _Implemented in v1.5.0_

`rhow_effect` (Default = 1.129)
:   Effective water density used in the non-polar dispersion term calculation, default to 1.129 for `decompopt = 2`, the 
_σ_ scheme. This was optimized in ([ref.][227]) with respect to the TIP3P solvent in PME. Optimized values for other 
decomposition schemes can be found in Table 4 of ([ref.][227]).

    _Implemented in v1.5.0_

`use_sav` (Default = 1)
:   The option to use molecular volume (the volume enclosed by SASA) or to use molecular surface (SASA) for cavity term 
calculation. Recent study shows that the molecular volume approach transfers better from small training molecules to 
biomacromolecules.

    * 0: Use SASA to estimate cavity free energy
    * 1: Use the molecular volume enclosed by SASA

    _Implemented in v1.5.0_

`cavity_surften` (Default = 0.0378)
:   The regression coefficient for the linear relation between the total non-polar solvation free energy (`inp` = 1), or
the cavity free energy (`inp = 2`) and SASA/volume enclosed by SASA. The parser default is for `inp = 2` and set to the
best of three tested schemes as reported in ([ref.][227]), _i.e._ `decompopt = 2`, `use_rmin = 1`, and `use_sav = 1`.
When `inp = 1` (the PB default), `gmx_MMPBSA` resets this value to `0.005` during input processing. See
recommended values in Tables 1-3 for other schemes.

`cavity_offset` (Default = -0.5692)
:   The regression offset for the linear relation between the total non-polar solvation free energy (`inp`= 1), or
the cavity free energy (`inp = 2`) and SASA/volume enclosed by SASA. The parser default is for `inp` = 2 and set to
the best of three tested schemes as reported in ([ref.][227]), _i.e._ `decompopt = 2`, `use_rmin = 1`, and `use_sav = 1`.
When `inp = 1`, `gmx_MMPBSA` resets this value to `0.0` during input processing.
See recommended values in Tables 1-3 for other schemes.

`maxsph` (Default = 400)
:   Approximate number of dots to represent the maximum atomic solvent accessible surface. These dots are first checked 
against covalently bonded atoms to see whether any of the dots are buried. The exposed dots from the first step are 
then checked against a non-bonded pair list with a cutoff distance of 9 Å to see whether any of the exposed dots 
from the first step are buried. The exposed dots of each atom after the second step then represent the solvent 
accessible portion of the atom and are used to compute the SASA of the atom. The molecular SASA is simply a 
summation of the atomic SASA’s. A molecular SASA is used for both PB dielectric map assignment and for NP calculations.

  [5]: https://ambermd.org/doc12/Amber21.pdf#chapter.6
  [6]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.24467

    _Implemented in v1.5.0_

`maxarcdot` (Default = 1500)
:   Number of dots used to store arc dots per atom.

#### **Options for output**

`npbverb` (Default = 0)
:   Verbose mode.

    * 0: Off
    * 1: On

    _Implemented in v1.5.0_

### **`&rism` namelist variables**

!!! note "Keep in mind"
    * A default 3drism input file can be created as follows:

        ```
        gmx_MMPBSA --create_input rism
        ```
    
    * `3D-RISM` calculations are launched through the AmberTools `sander` backend by the
    `EnergyCalculation` wrapper. For each processed frame, `gmx_MMPBSA` generates a Sander RISM input
    file and collects rank-specific `*_rism.mdout.<rank>` output. The wrapper distributes frame work
    across MPI ranks; this is separate from any solver-level parallelism provided by the AmberTools
    backend. `rism3d.snglpnt` is retained below only as the historical pre-v1.5.2 backend reference.
    3D-RISM is the most expensive, yet most statistical mechanically rigorous solvation model. See
        * [Introduction to RISM](https://ambermd.org/doc12/Amber21.pdf#section.7.1) for a thorough description RISM 
        theory.
        * [General workflow for using 3D-RISM](https://ambermd.org/doc12/Amber21.pdf#section.7.3)
        * Practical considerations on:
            * [Computational Requirements and Parallel Scaling of RISM](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.1)
            * [Numerical Accuracy of RISM](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3)
            * [Convergence issues](https://ambermd.org/doc12/Amber21.pdf#subsection.7.3.1)
    * A sample 3drism input file is shown [here](input_file.md#mm3d-rism)
    * A tutorial on binding free energy calculation with 3D-RISM is available [here](examples/3D-RISM/README.md)
    * We have included more variables in 3D-RISM calculations than the ones available in the MMPBSA.py original code. 
    That way, users can be more in control and tackle various issues (_e.g._, convergence issues).
    * One advantage of `3D-RISM` is that an arbitrary solvent can be chosen; change the `-xvvfile` command-line option
    to point to the desired XVV file. By default, the command first uses
    `$AMBERHOME/AmberTools/test/rism1d/tip3p-kh/tip3p.xvv.save` when that file exists; otherwise it uses the bundled
    `GMXMMPBSA/data/xvv_files/tip3p.xvv` file resolved from the installed package. You can provide any other
    precomputed `.xvv` file explicitly with `-xvvfile`.

  [7]: https://ambermd.org/doc12/Amber21.pdf#chapter.7
  [8]: https://ambermd.org/doc12/Amber21.pdf#subsection.36.3.2

#### **Closure approximations**

`closure` (Default = "kh")
:   Comma separate list of closure approximations. If more than one closure is provided, the 3D-RISM solver will use 
the closures in order to obtain a solution for the last closure in the list when no previous solutions are available.
The solution for the last closure in the list is used for all output. The use of several closures combined with 
different tolerances can be useful to overcome convergence issues (see [§7.3.1](https://ambermd.org/doc12/Amber21.
pdf#subsection.7.3.1))

    * "kh": Kovalenko-Hirata
    * "hnc": Hyper-netted chain equation
    * "psen": Partial Series Expansion of order-n where “n” is a positive integer (_e.g._, "pse3")

    !!! example "Examples"

            === "v1.5.2"
                === "One closure"
                         closure="pse3"
                === "Several closures"
                         closure="kh","pse3"
            === "< v1.5.2"
                === "One closure"
                         closure="pse3"
                === "Several closures"
                         closure="kh,pse3"

#### **Solvation free energy corrections**

=== "v1.5.2"
    
    !!! info
        The `thermo` variable has been removed in v1.5.2. Now the standard closure relation is always reported. Use 
        `gfcorrection` and/or `pcpluscorrection` variables to compute additional excess chemical potential functionals.

    `gfcorrection` (Default = 0)
    :    Compute the Gaussian fluctuation excess chemical potential functional. 
    See [§7.1.2](https://ambermd.org/doc12/Amber21.pdf#subsection.7.1.2)

        * 0: Off
        * 1: On

    `pcpluscorrection` (Default = 0)
    :    Compute the PC+/3D-RISM excess chemical potential functional.
    See [§7.2.4](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.4)

        * 0: Off
        * 1: On

=== "< v1.5.2"
    `thermo` (Default = "std")
    :   Which thermodynamic equation you want to use to calculate solvation properties. Options are:
    
        * "std": uses the standard closure relation
        * "gf": Compute the Gaussian fluctuation excess chemical potential functional
        * "both": print out separate sections for all
        
        !!! note
            Note that all data are printed out for each RISM simulation, so no choice is any more computationally 
            demanding than another.

#### **Long-range asymptotics**

!!! info
    Long-range asymptotics are used to analytically account for solvent distribution beyond the solvent box. 
    Long-range asymptotics are always used when calculating a solution but can be omitted for
    the subsequent thermodynamic calculations, though it is not recommended.

`noasympcorr` (Default = 1)
:   Disable long-range asymptotic corrections for thermodynamic output only. Long-range asymptotics are still used to
    calculate the 3D-RISM solution. This follows the AmberTools `rism3d.snglpnt` definition of the `--noasympcorr`
    switch.

    * 0: Use long-range asymptotic corrections for thermodynamic output
    * 1: Disable long-range asymptotic corrections for thermodynamic output

    _Implemented in v1.5.0_

`treeDCF` (Default = 1)
:   Use direct sum, or the treecode approximation to calculate the direct correlation function long-range asymptotic 
correction.

    * 0: Use direct sum
    * 1: Use treecode approximation

    _Implemented in v1.5.0_

`treeTCF` (Default = 1)
:   Use direct sum, or the treecode approximation to calculate the total correlation function long-range asymptotic 
correction.

    * 0: Use direct sum
    * 1: Use treecode approximation

    _Implemented in v1.5.0_

`treeCoulomb` (Default = 0)
:   Use direct sum, or the treecode approximation to calculate the Coulomb potential energy.

    * 0: Use direct sum
    * 1: Use treecode approximation

    _Implemented in v1.5.0_

`treeDCFMAC` (Default = 0.1)
:   Treecode multipole acceptance criterion for the direct correlation function long-range asymptotic correction.

    _Implemented in v1.5.0_

`treeTCFMAC` (Default = 0.1)
:   Treecode multipole acceptance criterion for the total correlation function long-range asymptotic correction.

    _Implemented in v1.5.0_

`treeCoulombMAC` (Default = 0.1)
:   Treecode multipole acceptance criterion for the Coulomb potential energy.

    _Implemented in v1.5.0_

`treeDCFOrder` (Default = 2)
:   Treecode Taylor series order for the direct correlation function long-range asymptotic correction.

    _Implemented in v1.5.0_

`treeTCFOrder` (Default = 2)
:   Treecode Taylor series order for the total correlation function long-range asymptotic correction. Note that the 
Taylor expansion used does not converge exactly to the TCF long-range asymptotic correction, so a very high order 
will not necessarily increase accuracy.

    _Implemented in v1.5.0_

`treeCoulombOrder` (Default = 2)
:   Treecode Taylor series order for the Coulomb potential energy.

    _Implemented in v1.5.0_

`treeDCFN0` (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the direct correlation function 
long-range asymptotic correction. This sets the depth of the hierarchical octtree.

    _Implemented in v1.5.0_

`treeTCFN0` (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the total correlation function 
long-range asymptotic correction. This sets the depth of the hierarchical octtree.

    _Implemented in v1.5.0_

`treeCoulombN0` (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the Coulomb potential energy. This 
sets the depth of the hierarchical octtree.

    _Implemented in v1.5.0_

#### **Solvation box**

!!! info
    The non-periodic solvation box super-cell can be defined as variable or fixed in size. When a
    variable box size is used, the box size will be adjusted to maintain a minimum buffer distance between the atoms
    of the solute and the box boundary. This has the advantage of maintaining the smallest possible box size while
    adapting to change of solute shape and orientation. Alternatively, the box size can be specified at run-time. This
    box size will be used for the duration of the sander calculation. Solvent box dimensions have a strong effect on 
    the numerical precision of 3D-RISM. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for 
    recommendation on selecting an appropriate box size and resolution.

##### **Variable box size**

!!! info "Keep in mind"
    It is recommended to avoid specifying a large, prime number of processes (≥ 7) when using a variable solvation 
    box size.

`buffer` (Default = 14)
:   Minimum distance (in Å) between solute and edge of solvation box. Specify this with `grdspc` below. Mutually 
exclusive with `ng` and `solvbox`. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on 
how this affects numerical accuracy and how this interacts with `ljTolerance`, and `tolerance`

    * when < 0: Use fixed box size (see `ng` and `solvbox` below)
    * when >= 0: Use `buffer` distance

`grdspc`(Default = 0.5,0.5,0.5)
:   Grid spacing (in Å) of the solvation box. Specify this with `buffer` above. Mutually exclusive with `ng` and 
`solvbox`.

##### **Fixed box size**

`ng` (Default = -1,-1,-1)
:   Comma separated number of grid points to use in the x, y, and z directions. Used only if buffer < 0. Mutually 
exclusive with `buffer` and `grdspc` above, and paired with `solvbox` below.

    !!! warning 
        No default, this must be set if buffer < 0. As a general requirement, the number of grids points in each 
        dimension must be divisible by two, and the number of grid points in the z-axis must be divisible by the
        number of MPI processes.

        As an example: define like `ng=1000,1000,1000`, where all numbers are divisible by two 
        and you can use 1, 2, 4, 5, 8, 10... processors, all divisors of 1000 (value in the z-axis).

        Take into account that at a certain level, running RISM in 
        parallel may actually hurt performance, since previous solutions are used 
        as an initial guess for the next frame, hastening convergence. Running in parallel loses this advantage. Also, 
        due to the overhead involved in which each MPI rank is required to load every topology file when calculating
        energies, parallel scaling will begin to fall off as the number of MPI ranks reaches the number of frames.

`solvbox` (Default = -1,-1,-1)
:    Sets the size in Å of the fixed size solvation box. Used only if `buffer` < 0. Mutually exclusive with `buffer` 
and `grdspc` above, and paired with `ng` above. 

    !!! warning 
        No default, this must be set if buffer < 0. Define like `solvbox=20,20,20`

`solvcut` (Default = -1)
:   Cutoff used for solute-solvent interactions. The parser default `-1` delegates the cutoff to the RISM engine,
which resolves it from `buffer` when a variable-size box is used. If `buffer < 0`, specify a positive `solvcut` together
with `ng` and `solvbox`; otherwise the program stops before the calculation.

#### **Solution convergence**

`tolerance` (Default = 0.00001)
:   A comma-separated list of maximum residual values for solution convergence. This has a strong effect on the 
cost of 3D-RISM calculations (smaller value for tolerance -> more computation). When used in combination with a list 
of closures it is possible to define different tolerances for each of the closures. This can be useful for difficult 
to converge calculations (see [§7.4.1](https://ambermd.org/doc12/Amber21.pdf#page=120&zoom=100,96,798)). For the sake of 
efficiency, it is best to use as high a tolerance as possible for all but the last closure. 
See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects numerical 
accuracy and how this interacts with `ljTolerance`, `buffer`, and `solvbox`. Three formats of list are possible:

    * one tolerance: All closures but the last use a tolerance of 1. The last tolerance in the list is used
    by the last closure. In practice this is the most efficient.
    * two tolerances: All closures but the last use the first tolerance in the list. The last tolerance in the
    list is used by the last closure.
    * n tolerances: Tolerances from the list are assigned to the closure list in order.

    !!! example "Examples"
        === "v1.5.2"
            === "One closure/One tolerance"
                    closure="pse3", tolerance=0.00001
                
                A tolerance of 0.00001 will be used for closure "pse3".
            === "Several closures/One tolerance"
                     closure="kh","pse3", tolerance=0.00001
    
                A tolerance of 1 will be used for closure "kh", while 0.00001 will be used for closure "pse3".
                Equivalent to `closure="kh", "pse3", tolerance=1,0.00001`
            === "Several closures/Two tolerances"
                     closure="kh","pse2","pse3", tolerance=0.01,0.00001
    
                A tolerance of 0.01 will be used for closures "kh" and "pse2", while 0.00001 will be used for closure
                "pse3". Equivalent to `closure="kh","pse2","pse3", tolerance=0.01,0.01,0.00001`
            === "Several closures/Several tolerances"
                     closure="kh","pse2","pse3", tolerance=0.1,0.01,0.00001
    
                A tolerance of 0.1 will be used for closure "kh", 0.01 will be used for closure "pse2", while 0.00001
                will be used for closure "pse3".

        === "<v1.5.2"
            === "One closure/One tolerance"
                    closure="pse3", tolerance=0.00001
                
                A tolerance of 0.00001 will be used for closure "pse3".
            === "Several closures/One tolerance"
                     closure="kh,pse3", tolerance=0.00001
    
                A tolerance of 1 will be used for closure "kh", while 0.00001 will be used for closure "pse3".
                Equivalent to `closure="kh, pse3", tolerance=1,0.00001`
            === "Several closures/Two tolerances"
                     closure="kh,pse2,pse3", tolerance=0.01,0.00001
    
                A tolerance of 0.01 will be used for closures "kh" and "pse2", while 0.00001 will be used for closure
                "pse3". Equivalent to `closure="kh,pse2,pse3", tolerance=0.01,0.01,0.00001`
            === "Several closures/Several tolerances"
                     closure="kh,pse2,pse3", tolerance=0.1,0.01,0.00001
    
                A tolerance of 0.1 will be used for closure "kh", 0.01 will be used for closure "pse2", while 0.00001
                will be used for closure "pse3".

`ljTolerance` (Default = -1)
:   Lennard-Jones accuracy (Optional.) Determines the Lennard-Jones cutoff distance based on the desired accuracy of 
the calculation. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects 
numerical accuracy and how this interacts with `tolerance`, `buffer`, and `solvbox`.

    _Implemented in v1.5.0_

`asympKSpaceTolerance` (Default = -1)
:   Tolerance reciprocal space long range asymptotics accuracy (Optional.) Determines the reciprocal space long 
range asymptotic cutoff distance based on the desired accuracy of the calculation. 
See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects numerical 
accuracy. Possible values are:

    * when < 0: asympKSpaceTolerance = tolerance/10
    * when = 0: no cutoff
    * when > 0: given value determines the maximum error in the reciprocal-space long range asymptotics calculations

    _Implemented in v1.5.0_

`mdiis_del` (Default = 0.7)
:   MDIIS step size.

    _Implemented in v1.5.0_

`mdiis_nvec` (Default = 5)
:   Number of previous iterations MDIIS uses to predict a new solution.

    _Implemented in v1.5.0_

`mdiis_restart` (Default = 10)
:   If the current residual is mdiis_restart times larger than the smallest residual in memory, then the MDIIS 
procedure is restarted using the lowest residual solution stored in memory. Increasing this number can sometimes 
help convergence.

    _Implemented in v1.5.0_

`maxstep` (Default = 10000)
:   Maximum number of iterations allowed to converge on a solution.

    _Implemented in v1.5.0_

`npropagate` (Default = 5)
:   Number of previous solutions propagated forward to create an initial guess for this solute atom configuration.

    * =0: Do not use any previous solutions
    * = 1..5: Values greater than 0 but less than 4 or 5 will use less system memory but may introduce artifacts to 
    the solution (_e.g._, energy drift).

    _Implemented in v1.5.0_

#### **Output**

`polardecomp` (Default = 0)
:   Decomposes solvation free energy into polar and non-polar components. Note that this typically requires 80% more 
computation time.

    * 0: Do not decompose solvation free energy into polar and non-polar components. 
    * 1: Decompose solvation free energy into polar and non-polar components.

`entropicdecomp` (Default = 0)
:   Decomposes solvation free energy into energy and entropy components. Also performs temperature derivatives of other 
calculated quantities. Note that this typically requires 80% more computation time and requires a `.xvv` file version 
1.000 or higher (available within `GMXMMPBSA` data folder). 
See [§7.1.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.1.3) and 
[§7.3](https://ambermd.org/doc12/Amber21.pdf#section.7.3)

    * 0: No entropic decomposition
    * 1: Entropic decomposition

`rism_verbose` (Default = 0)
:   Level of output in temporary RISM output files. May be helpful for debugging or following convergence. 

    * 0: just print the final result
    * 1: additionally prints the total number of iterations for each solution
    * 2: additionally prints the residual for each iteration and details of the MDIIS solver (useful for debugging 
    and convergence analyses)

    !!! warning "AmberTools runtime compatibility"
        If a 3D-RISM job stops before the calculation starts with
        `Fortran runtime error: Missing comma between descriptors` from `amber_rism_interface.F90`, this is a known
        AmberTools/Fortran runtime compatibility issue. Changing `rism_verbose` or other `&rism` input options does
        not resolve this failure. A known working workaround is `gmx_MMPBSA` 1.6.4 with Python 3.9 or 3.10,
        AmberTools 23, and `libgfortran5`/`libgcc-ng` 12.x, or a patched AmberTools build.

### **`&alanine_scanning` namelist variables**

!!! note "Keep in mind"
    * A default alanine scanning input file can be created as follows:

        ```
        gmx_MMPBSA --create_input ala
        ```
    
    * A sample alanine scanning input file is shown [here](input_file.md#alanine-scanning)
    * A tutorial on alanine scanning is available [here](examples/Alanine_scanning/README.md)

`mutant_res` (Default = None. Must be defined)
:   Select one or more residues using `CHAIN/RESNUM`. Append an insertion code directly to the residue number
    (`A/27B`); the colon form (`A/27:B`) is invalid. Numeric ranges select only blank-code residues, so list
    insertion-coded residues individually. For example, `A/13,25` mutates two residues together in one composite mutant.
    Selections across chains are allowed only when all selected residues belong to the same component.

    !!! important
        * All selected residues must belong to the receptor or all to the ligand; mixing components is rejected.
        * Every selected residue uses the same `mutant` target (`ALA` or `GLY`). A composite run reports one
          combined mutation effect, not separate single-residue scans.
        * Composite mutations require `cas_intdiel=0`; set the desired dielectric explicitly in the GB/PB namelist.
        * Use a reference structure (`-cr`) when specific chain IDs, residue numbers, or insertion codes are needed.
        * The program builds mutant topologies and trajectories from the selected residues.

`mutant` (Default = "ALA") 
:   Defines the residue that it is going to be mutated for. Allowed values are: 

    * `"ALA"` or `"A"`: Alanine scanning
    * `"GLY"` or `"G"`: Glycine scanning

    _Changed in v1.3.0: Change mol (receptor or ligand) by mutant aminoacid (ALA or GLY)_

`mutant_only`  (Default = 0)
:   Option to perform specified calculations only for the mutants. 

    * 0: Perform calcultion on mutant and original
    * 1: Perform calcultion on mutant only
    
    !!! note
        Note that all calculation details are controlled in the other namelists, though for alanine scanning to be 
        performed, the namelist must be included (blank if desired)

`cas_intdiel` (Default = 0)
:   The dielectric constant (`intdiel`(GB)/`indi`(PB)) will be modified depending on the nature of the residue to be 
mutated. 
    
    * 0: Do not use adaptative `intdiel` assignation
    * 1: Use adaptative `intdiel` assignation

    !!! important
        * Works with the GB and PB calculations
        * It is ignored when `intdiel`(GB)/`indi`(PB) has been explicitly defined, that is, it is ignored if 
        `intdiel != 1.0`/`indi != 1.0` (default values)
        * Dielectric constant values has been assigned according to [Yan et al., 2017][9]

  [9]: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00734
    
    _Implemented in v1.4.2_

`intdiel_nonpolar` (Default = 1)
:   Define the `intdiel`(GB)/`indi`(PB) value for non-polar residues (`PHE`, `TRP`, `VAL`, `ILE`, `LEU`, `MET`, `PRO`,
`CYX`, `ALA`, `GLY`, `PRO`)
    
    _Implemented in v1.4.2_

`intdiel_polar` (Default = 3)
:   Define the `intdiel`(GB)/`indi`(PB) value for polar residues (`TYR`, `SER`, `THR`, `CYM`, `CYS`, `HIE`, `HID`, 
`ASN`, `GLN`, `ASH`, `GLH`, `LYN`)
    
    _Implemented in v1.4.2_

`intdiel_positive` (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for positive charged residues (`LYS`, `ARG`, `HIP`)
    
    _Implemented in v1.4.2_

`intdiel_negative` (Default = 5)
:   Define the `intdiel`(GB)/`indi`(PB) value for negative charged residues (`GLU`, `ASP`)
    
    _Implemented in v1.4.2_

### **`&decomp` namelist variables**

!!! note "Keep in mind"
    * A default decomp input file can be created as follows:

        ```
        gmx_MMPBSA --create_input decomp
        ```
    
    * A sample decomp input file is shown [here](input_file.md#decomposition-analysis)
    * A tutorial on binding free energy decomposition is available [here](examples/Decomposition_analysis/README.md)
    * With **`&gbnsr6`**, decomposition is a hybrid sander + GBNSR6 merge: only the polar (`pol`) column comes
      from GBNSR6. See the [GBNSR6 decomposition warning](input_file.md#gbnsr6-namelist-variables).

`idecomp` (Default = 2)
:   Energy decomposition scheme to use:
    
    * 1: Per-residue decomp with 1-4 terms added to internal potential terms
    * 2: Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms
    * 3: Pairwise decomp with 1-4 terms added to internal potential terms
    * 4: Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms

`dec_verbose` (Default = 1)
:   Set the level of output to print in the decomp_output file.

    * 0: DELTA energy, total contribution only
    * 1: DELTA energy, total, sidechain, and backbone contributions
    * 2: Complex, Receptor, Ligand, and DELTA energies, total contribution only
    * 3: Complex, Receptor, Ligand, and DELTA energies, total, sidechain, and backbone contributions

    !!! note
        If the values 0 or 2 are chosen, only the Total contributions are required, so only those will be printed to the
        mdout files to cut down on the size of the mdout files and the time required to parse them.

`print_res` (Default = "within 6")
:   Select residues whose information is going to be printed in the output file. The default selection should be 
sufficient in most cases, however we have added several additional notations
    
    !!! example "Selection schemes"

        === "By Distance (recommended)"
            Notation: [ `within` `distance` ]
            :   `within` corresponds to the keyword and `distance` to the maximum distance criterion in Å necessary to 
                select the residues from both the receptor and the ligand. In case the cutoff used is so small that 
                the number of decomp residues to print < 2, the cutoff value will be increased by 0.1 until the 
                number of decomp residues to print >= 2.
    
            !!! example
                `print_res="within 6"` Will print all residues within 6 Å between receptor and 
                ligand including both.
    
        === "Amino acid selection"
            Notation: [ `CHAIN`/(`RESNUM`, `RESNUM` + `INSERTION_CODE`, or `RESNUM-RESNUM`) ]
            :   Print individual residues or ranges. Append an insertion code directly to the residue number, for example
            `A/27B`; the colon form `A/27:B` is invalid. Numeric ranges select only blank-code residues, so insertion-coded
            residues must be defined individually.
    
            !!! example
                `print_res="A/1,3-10,15,100 B/25"` This will print Chain A residues 1, 3 through 10, 15, and 100 along with 
                chain B residue 25 from the complex topology file and the corresponding residues in either the ligand and/or 
                receptor topology files.
    
                !!! danger
                    make sure to include at least one residue from both the receptor and ligand in the `print_res` mask of 
                    the `&decomp` section. Check http://archive.ambermd.org/201308/0075.html
    
                Suppose that the following sequence contains chain A as the receptor and chain B as the ligand:
                A:LEU:5, A:GLY:6:A, A:THR:6:B, A:SER:6:C, A:ASP:6:D, A:ILE:7, B:25
                
                === "Supported notation"
                    
                    **Ranges selection**
                    :   `print_res="A/5-7 B/25"` prints only the blank-code residues A/5 and A/7 from chain A, plus
                        B/25. It does not include insertion-coded residues A/6A through A/6D.
                    
                    **Individual selection**
                    :   `print_res="A/5,6A,6B,6C,6D,7 B/25"` prints all listed residues, including the insertion-coded
                        residues from chain A.
    
                === "Wrong notation"
                    `print_res="A/5-6B,6D-7"` produces an error.
    
        === "All"
    
            Notation: `all`
            :   Prints all residues. This option is generally not recommended because most residue contributions are zero,
                so it can waste time and computational resources.
    
            !!! danger
                Using `idecomp=3 or 4` (pairwise) with a very large number of printed residues and a large number of frames 
                can quickly create very, very large temporary mdout files. Large print selections also demand a large amount 
                of memory to parse the mdout files and write decomposition output file (~500 MB for just 250 residues, since 
                that’s 62500 pairs!) It is not unusual for the output file to take a significant amount of time to print if 
                you have a lot of data. This is most applicable to pairwise decomp, since the amount of data scales as  
                O(N^2^).
     
    !!! important
        We recommend using the reference structure (-cr) to ensure the perfect match between the selected residue in 
        the defined structure or topology 
        
    
    _Changed in v1.4.0: Improve residue selection_

`csv_format`  (Default = 1)
:   Print the decomposition output in a Comma-Separated-Values (CSV) file. CSV files open natively in most
spreadsheets. 

    * 0: Data to be written out in the standard ASCII format.
    * 1: Data to be written out in a CSV file, and standard error of the mean will be calculated and included for all 
    data.

### **`&nmode` namelist variables**

!!! note "Keep in mind"
    * A default nmode input file can be created as follows:

        ```
        gmx_MMPBSA --create_input nmode
        ```
    
    * A sample nmode input file is shown [here](input_file.md#entropy-with-nmode)
    * A tutorial on normal mode analysis is available [here](examples/Entropy_calculations/nmode/README.md)

#### **Basic input options**

`nmstartframe`[^2] (Default = 1)
:   Frame number to begin performing `nmode` calculations on.

  [^2]: _These variables will choose a subset of the frames chosen from the variables in the `&general` namelist. Thus,
        the "trajectory" from which snapshots will be chosen for `nmode` calculations will be the collection of 
        snapshots upon which the other calculations were performed._

`nmendframe`[^2] (Default = 1000000)
:   Frame number to stop performing `nmode` calculations on 

`nminterval`[^2] (Default = 1)
:   Offset from which to choose frames to perform `nmode` calculations on

#### **Parameter options**

`nmode_igb` (Default = 1)
:   Value for Generalized Born model to be used in calculations. Options are:
    
    * 0: Vacuum
    * 1: HCT GB model 

`nmode_istrng` (Default = 0.0)
:   Ionic strength to use in `nmode` calculations. Units are Molarity (M). Non-zero values are ignored if `nmode_igb`
is 0 above.

`dielc` (Default = 1.0)
:   Distance-dependent dielectric constant 

`drms` (Default = 0.001)
:   Convergence criteria for minimized energy gradient.

`maxcyc` (Default = 10000)
:   Maximum number of minimization cycles to use per snapshot in sander.

    Frames that fail the minimized-energy-gradient criterion (`drms` / `maxcyc`) remain `NaN` and are omitted from
    NMODE averages and uncertainties. Increase `drms` or `maxcyc` if more frames should converge.

## Sample input files

!!! tip
    You can refer to the [examples](examples/README.md) to understand the input file in a practical way.

### GB

``` linenums="1"
Sample input file for GB calculation building the Amber topologies
from structures. Please refer to the section "How gmx_MMPBSA works"

&general
startframe=5, endframe=100, interval=5, verbose=2
/

&gb
igb=8, saltcon=0.150,
/
```

### GBNSR6

``` linenums="1"
Sample input file for GBNSR6 calculation building the Amber topologies
from structures. Please refer to the section "How gmx_MMPBSA works"

&general
startframe=5, endframe=100, interval=5, verbose=2
/

&gbnsr6
epsin=1.0, istrng=0.150,
/
```

### QM/MMGBSA

```
Sample input file for QM/MMGBSA

&general
startframe=5, endframe=100, interval=5,
/

&gb
igb=8, saltcon=0.100, ifqnt=1,
qm_residues="A/240-251 B/297", qm_theory="PM6-DH+"
/
```

### PB

``` linenums="1"
Sample input file for PB calculation building the Amber topologies
from structures. Please refer to the section "How gmx_MMPBSA works"

&general
startframe=5, endframe=100, interval=5
/

&pb
istrng=0.15, fillratio=4.0
/
```

### MMPBSA with membrane proteins

```
Sample input file for MMPBSA with membrane proteins

&general
startframe=1, endframe=100, interval=1,
/

&pb
memopt=1, emem=7.0, indi=4.0,
mctrdz=automatic, mthick=automatic, membrane_atoms="P", poretype=1,
radiopt=0, istrng=0.150, fillratio=1.25, inp=2,
sasopt=0, solvopt=2, ipb=1, bcopt=10, nfocus=1, linit=1000,
eneopt=1, cutfd=7.0, cutnb=99.0,
maxarcdot=15000,
npbverb=1,
/
```

### MM/3D-RISM

```
Sample input file for 3D-RISM

&general
startframe=20, endframe=100, interval=5,
/

&rism
polardecomp=1, gfcorrection=1
/
```

### Alanine scanning

``` linenums="1"
Sample input file for Alanine scanning

&general
startframe=5, endframe=21, interval=1, PBRadii=4
/

&gb
igb=8, saltcon=0.150, intdiel=10
/

&alanine_scanning
mutant='ALA', mutant_res='B:12'
/
```

### Decomposition analysis

```
Sample input file for decomposition analysis
Make sure to include at least one residue from both the receptor
and ligand in the print_res mask of the &decomp section.
http://archive.ambermd.org/201308/0075.html. This is automally
guaranteed when using "within" keyword.

&general
startframe=5, endframe=21, interval=1,
/

&gb
igb=8, saltcon=0.150,
/

&decomp
idecomp=2, dec_verbose=3,
# This will print all residues that are less than 4 Å between
# the receptor and the ligand
print_res="within 4"
/
```

### Entropy with NMODE

``` linenums="1"
Sample input file for entropy calculations

&general
startframe=5, endframe=21, interval=1,
temperature=298
/

&gb
igb=2, saltcon=0.150,
/

&nmode
nmstartframe=10, nmendframe=21, nminterval=2,
maxcyc=50000, drms=0.0001,
/
```

### Interaction Entropy

``` linenums="1"
Sample input file for entropy calculations

&general
startframe=5, endframe=21, interval=1,
# Interaction Entropy (IE)
# (https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682) approximation
interaction_entropy=1, ie_segment=25,
temperature=298
/

&gb
igb=2, saltcon=0.150,
/
```

!!! info "Accepted input syntax"
    The parser accepts one namelist per block. Start a block with `&general`, `&gb`, `&pb`, `&rism`, `&gbnsr6`,
    `&alanine_scanning`, `&decomp`, or `&nmode` on its own line and close it with `/` or `&end` on its own line.
    Variable assignments use `name = value`; numeric values are unquoted and strings may use single or double quotes.
    Separate assignments with commas, or put one assignment on each line without commas. Comma-separated values are
    supported for list variables. Lines beginning with `#` or `!` are ignored, and `#` starts an inline comment inside a
    namelist. Text outside namelist blocks is ignored. Unknown variables, duplicate namelists, and unterminated blocks
    are rejected.
