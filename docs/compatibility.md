---
template: main.html
title: Compatibility
---

# Upgrading
Upgrade to the latest version with:

    amber.python -m pip install --upgrade gmx_MMPBSA

Inspect the currently installed version with:

    amber.python -m pip show gmx_MMPBSA

or

    gmx_MMPBSA -v    


## Upgrading from 1.3.x to 1.4.x

### What's new?

### Changes

??? example "`temperature` instead `entropy_temp`"
    sdfsdfsdf

??? example "`sys_name`"
    sdfsdfsdf

??? example "`temperature` instead `entropy_temp`"
    sdfsdfsdf

??? example "`temperature` instead `entropy_temp`"
    sdfsdfsdf
    

=== "v1.4.x"

    Temperature to calculate Experimental Energy from Ki and Iteraction Entropy 
        
                

    * Sed sagittis eleifend rutrum
    * Donec vitae suscipit est
    * Nulla tempor lobortis orci

=== "v1.3.x"

    1. Sed sagittis eleifend rutrum
    2. Donec vitae suscipit est
    3. Nulla tempor lobortis orci


### Changes to `_GMXMMPBA_info` file

:octicons-archive-24: Deprecated: 5.5 · :octicons-trash-24: Removed: 6.0



`temperature` variable
:   Temperature to calculate Experimental Energy from Ki and Interaction Entropy
    
    === "v1.4.x"
        _Temperature to calculate Experimental Energy from Ki and Iteraction Entropy_
        
        `INPUT['temperature']` = 298.15
    
    === "v1.3.x"
        _Temperature to calculate Iteraction Entropy_
        
        `INPUT['entropy_temp']` = 298.15

### Additions to `_GMXMMPBA_info` file

`sys_name` variable
:   Define the System name. Useful when analyze multiple systems.

    

`exp_ki` variable
:   Specify experimental Ki in nM. Required to do correlation analysis.

    !!! tip
        This variable is optional since it can be defined in the start dialog for each system

### Changes to `mmpbsa.in` input file
**`sys_name`**  :material-new-box:{: .medium .heart } New: 1.4.x
:   Define the System name. Useful when analyze multiple systems.

    !!! warning
        If not defined, `gmx_MMPBSA_ana` will assign one according to the order in which it was loaded.
        
    !!! danger
        Currently, `gmx_MMPBSA_ana` cannot deal with duplicate names. Take this into account to avoid unforced errors


`print_res` variable   
:   Select residues from the complex to print.

    === "v1.4.x"
        This variable accepts three notations.

        !!! example inline end
            `print_res="within 6"` Will print all residues within 6 Angstroms between receptor and ligand
    
        **Distance criterion** 
        :   (`within` `distance`) --> `within` corresponds to the keyword and `6` to the maximum distance criterion in 
            Angstroms necessary to select the residues from both the receptor and the ligand;
        
        !!! example inline end
            `print_res="within 6 not lig"` Will print all residues within 6 Angstroms between receptor and 
            ligand omitting the ligand ones.
        
        **distance criterion with exclusion** 
        :   (`within` `distance` `not` (`receptor` or `ligand`)) --> here we print out all residues within 6 
            Angstroms between receptor and ligand omitting the selected component ones

        !!! example inline end
            `print_res="A/1,3-10,15,100"` This will print Chain A residues 1, 3 through 10, 15, and 100 from the 
            complex topology file and the corresponding residues in either the ligand and/or receptor topology files.

        **Individual residues or ranges**
        :   (`CHAIN`/`RESNUM1`,`RESNUM2`-`RESNUM10`,`RESNUM11`:`ICODE`) --> Select residues indivual or ranges. This 
            notation also supports insertion codes, in which case you must define them individually
        
        !!! warning
            If an amino acid with an insertion code is found in a defined range, it will end in error


    === "v1.3.x"
        This variable also accepts a sequence of individual residues and/or ranges. The different fields must be 
        either comma- or semicolon-delimited. For example: print_res = "within 6", where within corresponds to the 
        keyword and 6 to the maximum distance criterion in Angstroms necessary to select the residues from both the
        receptor and the ligand; or print_res = "1,3-10,15,100", or print_res = "1;3-10;15;100". Both of these will 
        print residues 1, 3 through 10, 15, and 100 from the complex topology file and the corresponding residues in 
        either the ligand and/or receptor topology files.

        !!! warning inline end "Don't work" 
            _These numbers correspond to the residue number according to Amber that works with the renumbered 
            residues starting at 1_        
        
        Allowed values:

        * `print_res`="within 6"
        * `print_res`="1,3-10,15,100" or
        * `print_res` = "1;3-10;15;100"

        