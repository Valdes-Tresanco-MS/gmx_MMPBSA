---
template: main.html
title: Compatibility
---

# Upgrading
Upgrade to the latest version with:

    amber.python -m pip install --upgrade gmx_MMPBSA

Inspect the currently installed version with:

    amber.python -m pip show gmx_MMPBSA  or  gmx_MMPBSA -v    

## Upgrading from 1.4.3 to 1.5.0
**The new gmx_MMPBSA version 1.5.0 includes major changes in calculation and processing modules, making it 
incompatible with previous versions.**

New variables have been included in the input file to perform new calculations or simply provide the user with full 
control over the ones already available in gmx_MMPBSA. A number of variables has been also updated, reworked or 
removed. See the list below for more details:

### Variables
* :material-new-box:{: .heart } New variables
    * [&general namelist](input_file.md#general-namelist-variables)
        - `c2_entropy`
    * [&gb namelist](input_file.md#gb-namelist-variables)
        - `extdiel`
    * [&pb namelist](input_file.md#pb-namelist-variables)
        - `smoothopt` , `iprob` , `arcres` , `mprob` , `npbopt` , `accept` , `nbuffer` , `fscale` , `npbgrid` , 
        `scalec` , `nsnba` , `decompopt` , `use_rmin` , `sprob` , `vprob` , `rhow_effect` , `use_sav` , `maxsph` , 
        `npbverb` 
    * [&rism namelist](input_file.md#rism-namelist-variables)
        - `noasympcorr` , `ljTolerance` , `asympKSpaceTolerance` , `treeDCF` , `treeTCF` , `treeCoulomb` , `treeDCFMAC` , 
        `treeTCFMAC` , `treeCoulombMAC` , `treeDCFOrder` , `treeTCFOrder` , `treeCoulombOrder` , `treeDCFN0` , 
        `treeTCFN0` , `treeCoulombN0` , `mdiis_del` , `mdiis_nvec` , `mdiis_restart` , `maxstep` , `npropagate` , 
        `entropicDecomp`

* :material-update:{: .heart } Updated variables
    * [&general namelist](input_file.md#general-namelist-variables)
        * `PBRadii` , `interaction_entropy` , `assign_chainID` , `solvated_trajectory` , `verbose`

* :material-tools:{: .heart } Reworked variables
    * [&general namelist](input_file.md#general-namelist-variables)
        * `temperature`
    
* :material-trash-can:{: .heart } Removed variables
    * [&general namelist](input_file.md#general-namelist-variables)
        * `protein_forcefield` , `ligand_forcefield` , `use_sander`
    
### Calculations
* :material-new-box:{: .heart } New calculations
    * In the new version, it is possible to perform new calculations such as 
    [C2 Entropy](examples/Entropy_calculations/C2_Entropy/README.md) and 
    [Binding free energy calculation with non-linear PB solver](examples/NonLinear_PB_solver/README.md)
      
* :material-update:{: .heart } Updated calculations
    * The addition of new variables (see above) provide more control over &pb and &rism calculations
    
### gmx_MMPBSA_ana
gmx_MMPBSA_ana has been completely reworked, and it doesn't support files from previous versions. New functions for 
customizing/exporting graphs, change number of frames and showing/exporting data have been added.

---
## Upgrading from 1.3.x to 1.4.x
The differences between both versions are small, you can see them below
### Variables
* New variables in input file
    - `qh_entropy` replace `entropy = 1` 
      
        ( :material-new-box:{: .heart } _Since v1.4.2_)
      
    - `interaction_entrpy` replace `entropy = 2`
        
        ( :material-new-box:{: .heart } _Since v1.4.2_)
      
    - `ie_segment` replace `entropy_seg`  
    
        ( :material-new-box:{: .heart } _Since v1.4.2_)
      
    - `forcefields` replace `protein_forcefield` and `ligand_forcefield`
        
        ( :material-new-box:{: .heart } _Since v1.4.1_)
      
    - `temperature` replace `entropy_temp`
        
        ( :material-new-box:{: .heart } _Since v1.4.1_)
      
    - `sys_name`
        
        ( :material-new-box:{: .heart } _Since v1.4.0_)
      
    - `exp_ki` 
        
        ( :material-new-box:{: .heart } _Since v1.4.0_)
      
* Internal variables (the user only interact with it in the info file)
    - `complex_fixed` (_Since v1.4.0_)
* Modified variables
    - `print_res` (_Since v1.4.0_)
* Deprecated variables
    - `entropy` use `qh_entropy` or `interaction_entropy` instead
        
        ( :octicons-archive-24: _Deprecated: v1.4.2_ · :octicons-trash-24: _Removed: v1.5.0_ )
    
    - `entropy_seg` use `ie_segment` instead
      
        ( :octicons-archive-24: _Deprecated: v1.4.2_ · :octicons-trash-24: _Removed: v1.5.0_ )
    
    - `entropy_temp` use `temperature` instaed
      
        ( :octicons-archive-24: _Deprecated: v1.4.1_ · :octicons-trash-24: _Removed: v1.5.0_ )
      
    - `protein_forcefield` use `forcefields` instead
      
        ( :octicons-archive-24: _Deprecated: v1.4.1_ · :octicons-trash-24: _Removed: v1.5.0_ )
      
    - `ligand_forcefield` use `forcefields` instead
      
        ( :octicons-archive-24: _Deprecated: v1.4.1_ · :octicons-trash-24: _Removed: v1.5.0_ )

!!! tip
    Check the changes in [`&general namelist variables`](input_file.md#general-namelist-variables) section

### Command-line
* `gmx_MMPBSA_ana` changes the `-p` option by `-f` with more flexibility. Please check the [gmx_MMPBSA_ana 
  command-line](gmx_MMPBSA_command-line.md#gmx_mmpbsa_ana-command-line). (_Changed in v1.4.0_)

### Results and info file
We have ensured backwards compatibility with `gmx_MMPBSA`, however there are some changes you can make

=== "Old results"

    Since the calculations are done, we only have two options.
    
    ^^Define the variables in `gmx_MMPBSA_ana`^^
    :  As we described above, these variables can be defined or modified in the `gmx_MMPBSA_ana` start dialog (See the 
    [`gmx_MMPBSA_ana` documentation](analyzer.md))
    
    ^^Modify the `*info` file (usually `_GMXMMPBSA_info`)^^
    :  Added the new variables as we describe below. _We have modified the description in the `_GMXMMPBSA_info` file a 
    bit to aid editing._
        
        !!! example "Add variables to `_GMXMMPBSA_info` file"
            
            \# You can alter the variables below
            
            INPUT['debug_printlevel'] = 0
        
            INPUT['verbose'] = 2
        
            INPUT['csv_format'] = 1
        
            INPUT['dec_verbose'] = 0
            
            {++INPUT['temperature'] = 298.15++}
            
            {++INPUT['exp_ki'] = 0.0++}
            
            {++INPUT['sys_name'] = 'My System'++}
            
            INPUT['entropy_seg'] = 25
    
    
        !!! warning
            `complex_fixed` is an internal variable that cannot be defined in the input file. This variable stores the PDB
            file name of the fixed complex. This structure corresponds to the complex with the proper chain identifiers and 
            amino acid numbers. If it does not exist, a warning will be displayed, and the complex structure extracted 
            from the structure file defined with the -cs option will be used.
            
            Note that this can lead to inconsistencies, for example: if the file defined in the -cs option is in GRO 
            format, it will not have the string IDs.
    
            You can generate a structure for this variable as follows:
    
            * generate a copy of the `_GMXMMPBSA_COM.pdb` file
            * Open it with your preferred editor (we recommend one with column editing capabilities, such as Kate, Geany,
                or Sublime Text)
            * Check the chains ID, if they do not exist add them
            * Do not change the numbering of amino acids
            * Save the document as `_GMXMMPBSA_COM_FIXED.pdb`
            * In the `_GMXMMPBSA_info` file add the following line
    
                    FILES.complex_fixed = '_GMXMMPBSA_COM_FIXED.pdb'
        
            * save the document
    
        !!! tip "Removing `entropy_temp` variable"
            After defining the `temperature` variable, you can remove the `entropy_temp` variable. This avoids getting 
            the related warning.
    
=== "New calculation"
    
    As we describe in the [`&general namelist variables`](input_file.md#general-namelist-variables) section in the input 
    file, these three variables are optional since they can be defined or modified in the `gmx_MMPBSA_ana` start dialog. 
    However, modifying a set of systems can be cumbersome. **We recommend defining them in the input file**
    
    ???+ example
        
        ```
        &general
        sys_name="Protein-Ligand",
        temperature=310
        exp_ki=10
        /
        ```     
    
    !!! note
        `exp_ki` Only needed when performing correlation analysis

