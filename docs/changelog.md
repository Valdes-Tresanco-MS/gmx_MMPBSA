---
template: main.html
title: Changelog
---

# Changelog
## [gxm_MMPBSA v1.4.1 (08/04/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.1)
### Additions
- Added a new class `Residue` to handle residues selection in Gromacs format with Amber index
- Verification of water molecules in user defined groups 
- Gromacs timer

### Fixes
- Gromacs topology conversion
- `qm_residues` notation
- Default path in `gmx_MMPBSA_test`
- The Entropy representation in `gmx_MMPBSA_ana`
- Bug when the structure has insertion code
- Improved ΔG Binding plot representation

### Changes
- Now `forcefields` variable unified `protein_forcefield` and `ligand_forcefield`
- `protein_forcefield` and `ligand_forcefield` pass to deprecate state. They will be removed in the next version.
- Improved documentation
  - Examples
  - Command-line
    - MPI
    - Examples
  - Links and references
  - Updated to material 7.1.0
    - Dark mode
    - Material "Back to Top" button
    - Grammatical corrections
  - Installation section
  - Figures caption
- The Ambiguous name for Entropy term in output files

## [gxm_MMPBSA v1.4.0 (22/03/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.0)
This release focuses almost entirely on `gmx_MMPBSA_ana` with minor issues fixed in `gmx_MMPBSA`
### Additions
- New start window to select options  
    - Option to make correlation (Pearson and Spearman coefficients)
    - Option to hide decomposition data
    - Option to not compute charts with a non-significant contribution
    - Option to not include terms with a non-significant contribution in bar charts
    - Selection of the components to display in addition to Delta (i.e. complex, receptor, and ligand)
    - Toggle the chart toolbar for a cleaner visualization
    - An informative table with selected systems data
        - Option to exclude any system
        - Option to change:
            - The system name
            - The experimental Ki for correlation
            - The temperature to calculate the Experimental Energy and the Interaction Entropy
    - Data reader with progress bar and multiprocessing
- Multiple systems to analyze in the same session  
- Correlation dock
    - Multiple models at the same time
    - Graphs and correlation data for each calculated energy term (ΔH, ΔH+IE, ΔH+NMODE and ΔH+QH)  
    - Table with the experimental energy of the systems, and the data of the selected model
- New arguments flags for gmx_MMPBSA_ana (See the [gmx_MMPBSA_ana documentation](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/analyzer/))
    - replaced `-p` by `-f`
    - `-f` accept a folder, single info file or a list of them
    - New flag `-r`. This flag allows to load all the systems inside a selected folder
- New graphical options
    - A new set of graphics (heatmap)
        - Per-frame when analyzing Per-residue and pair in Per-wise
        - Relation matrix for Per-wise
    - Interactive visualization of PDB files with per-residue energies with `PyMOL` (up to 5 instances).
    - Regression plot for correlation
- Plot features
    - Added Standard deviation to bar plots  
    - Added rolling average to line plots
    - Added indicators for the selected interval and average value in IE chart
    - Added crosshair cursor for better analysis on charts
- Multiprocessing application for testing (`gmx_MMPBSA_test`) 
- Embed Youtube videos for `gmx_MMPBSA_ana`

### Fixes
- Now `gmx_MMPBSA_ana` shows stability results as expected
- Errors in the documentation
- MPI

### Changes
- Converted analyzer.py into a sub-module for more flexibility, organization and portability
- Residues notation for mutation: CHAIN:RESNAME:RESNUMBER:ICODE instead of Amber residue index
- Improve the selection method in decomposition calculation
- Replaced variable `entropy_temp` (deprecated) by `temperature`
- IE in API
- Color Palette used in graphs
- Use seaborn and matplotlib for charts
- Use Pandas Dataframe and numpy to store data
- Changed the data structure to implement dynamic selection of frames in future versions
- Improved data export: now any item can be exported as CSV file
- Improved the documentation
    - Improve examples documentation
    - Added changelog button at home
    - Separated changelog in a new header
    - Added tags to mark the history of changes of variables and functionalities
  

## [gxm_MMPBSA v1.3.3 (09/03/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.3)
### Fixes
* fixed Boltzmann constant for IE
* fix mutation in ligand
* fixed analyzer error when interval > 1
* fixed residue selection within
* fixed ChainID assignation when reference structure is defined
* fixed the selection to print when decomposition

## [gxm_MMPBSA v1.3.2 (01/03/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.2) 
### Additions
* **Now, gmx_MMPBSA is in Zenodo 
  [![DOI](https://zenodo.org/badge/295050575.svg)](https://zenodo.org/badge/latestdoi/295050575). You can refer to 
  us in this way in what we publish the article**
* Added Interaction Entropy to gmx_MMPBSA output file
* Added a new class to save IE in a csv file
* Added "Go to Top" button to documentation HTML.
### Fixes
* Error when ligand and/or receptor are discontinuous (Testing it)
* Error when ligand and/or receptor are discontinuous and numbered non-consecutively
* Non-critical errors and inconsistencies in documentation
### Changes
* `receptor_mask` and `ligand_mask` have been removed from input file variables. Now we extract the amber mask directly 
  based on the GROMACS index file
* The receptor and ligand mapping in the complex was improved. Now we use a method based on the GROMACS index file
* The method `Map` of the `system_MMPBSA` class has been restructured. Now always processes amber masks
* Changing the IE calculation function to a class


## [gmx_MMPBSA-v1.3.1](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.1)
### Additions
* New variable (`overwrite_data`) to overwrite gmx_MMPBSA data. 
* More informative message when sander fail. Useful for PB calculation

### Fixes
* Protein-ligand with charmm force field example
* Stability calculation
* gmx path error
* leaprc.GLYCAM_06h-1 file
* Protein-glycan example

### Changes
* Documentation banner

## [gmx_MMPBSA-v1.3.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.0)

### Additions
* Documentation at Github pages
* Charmm force field support
* Amber topology generation from GROMACS topology. Work for Charmm and Amber force fields
  * New flag for topologies (`-cp`, `-rp` and `-lp`) added
* Now `gmx_MMPBSA` supports discontinuous receptor and ligand.
* Glycine scanning
* Autocompletion script for both `gmx_MMPBSA` and `gmx_MMPBSA_ana`
* Versioneer to control the semantic version.
* Argument type checker for the command-line 

### Fixes
* Alanine scanning tutorial 
* GROMACS executable path 
* The `-gui` option has been changed by `-nogui` and fixed when it is defined 
* Improvement on documentation

### Changes
* Documentation theme. Now we use Material
* Alanine scanning variable. Now `mutant` correspond to mutant amino acid (ALA and GLY)
* The `gmx_MMPBSA_gui` was changed by `gmx_MMPBSA_ana` 
* Improvement on the topologies construction process
* Order in which the trajectories are cleaned. Now, the topology is built and finally, the trajectories are cleaned

## [gmx_MMPBSA v1.2.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.2.0)
### Additions

* New ligand force field (Zwitterionic amino acids)
* A new flag (-cr) added for defining a reference structure (guarantee a better consistency in generated PDB files)
* API documentation

### Fixes
* More comprehensive output log file
* Best handling of structure files

### Changes
`gmx` `editconf` is used to generate PDB files instead of `gmx` `trjconv` (#14)
`gmx_MMPBSA` data is copied in AMBER as an independent folder
*.gro files can be used as a MD Structure+mass(db) file
Updated tutorial list in README (Protein_DNA_RNA_Ion_ligand BFE calculations)

## [gmx_MMPBSA v1.1.1](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.1.1)
### Additions
* New tutorial added (see Protein_DNA_RNA_Ion_ligand tutorial)
 
### Fixes
* Support various metallo-complexes formats

### Changes
* Keep all the temporary files in the folder for debugging purposes

## [gxm_MMPBSA v.1.1.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v.1.1.0)
### Additions
* Now supports carbohydrates as ligand. See this tutorial 
* Now supports metalloprotein-ligand complexes. See this tutorial.
* We have added data folder to gmx_MMPBSA module. This folder contains the GLYCAM_06h-1 force field files 
  (Compatible with amber99sb and early, see at http://glycam.org) which is not in AmberTools.

### Fixes
* Minor bugs

### Changes
* We changed the notation of the force fields, now the user can define any force field (We have only tested Amber and 
GLYCAM force fields) available in AmberTools. Charmm is not yet supported. See this section

## [gxm_MMPBSA first release](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.0.0)
**This version includes:**

* Compatibility with different Gromacs versions
* Support for all types of calculations available in MMPBSA.py
* Graphical user interface for results analysis (gmx_MMPBSA_ana)
* API modified to get more information
* Some new facilities and types of calculations