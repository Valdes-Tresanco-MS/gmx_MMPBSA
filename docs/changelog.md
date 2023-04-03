---
template: main.html
title: Changelog
---
# Changelog

## gmx_MMPBSA v1.6.1 (04/04/2023)

### Fixes

- PB decomp ends in error (https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/352)
- Better parsing for errors

## gmx_MMPBSA v1.6.0 (02/19/2023)

### Additions
- GBNSR6 model implementation
    - Enthalpy
    - Per-residue and Per-wise decomposition
    - Documentation and examples
- Added git dependency to env.yml file
- Groups can now be defined as either the number or the name (#157)
- Added error tracking for groups that are not in the index file

### Fixes
- Incompatible `compilers` version
- Fixed bug when CAS and PB/RISM calculations are performed
- Fixed duplicated options and args checking
- Documentation improved
- Fixed numpy deprecated np.float (#316)
- Fixed reproducible command line in debugging in log file (#335)
- GROMACS 2023 compatibility (#327, #335)
- Fixed IE and C2 Entropy (#325)
- PyMOL doesn't open within gmx_MMPBSA environment (#331)
- Inaccessible _temp_top.top when the topology file is in different directory (#299)

### Changes
- A copy old info file is created when rewrite-output (for backward compatibility)
- Setting inp = 1 by default (#329)

---

## gmx_MMPBSA v1.5.7 (09/10/2022)

### Additions
- Progress bar implementation
- verbose options for API and info logging for IE and C2 entropy calculations
- Added logging info for MPI calculations

### Fixes

- QH entropy calculation
- Decomposition bug (#269)
- Alanine scanning in Terminal residues
- Regression to output files subwindow in gmx_MMPBSA_ana
- Bug when using the compact results version
- Seaborn compatibility
- `print_res` in `&decomp` namelist does not accept "all" (#292)
- error when output files are open and any plot property is changed
- error if at least one residue of both the receptor and the ligand are not included in the selection for decomposition


### Changes
- now use `mpi4py` when `mpi_size > 1` otherwise use `fake_MPI`. Now "MPI" flag is not required 

---

## gmx_MMPBSA v1.5.6 (07/06/2022)

### Fix

#### gmx_MMPBSA

- Inconsistent decomposition output (#246, #248)
- IE and C2 entropies calculation fail when rewrite-output (#253)
- Delta Delta Entropies' values are incorrect (#257)
- Loading binary file fails when alanine scanning is performed (#258)

#### gmx_MMPBSA_ana

- Mutant-Normal data is not computed or displayed (#254)
- Inconsistent energy values in PyMol visualization (#255)

### Changes

#### gmx_MMPBSA

- Re-implementation of compact binary output file (#249)

---

## gmx_MMPBSA v1.5.5 (06/10/2022)

### Additions

#### gmx_MMPBSA

- Improved the API (#209)
    - Added a new method to get the Enthalpy 
    - Added a new method to get the Entropy
    - Added a new method to get the Binding
    - Added a new method to get the Decomposition energy
    - Added a new method to get the gmx_MMPBSA_ana data

#### gmx_MMPBSA_ana
- Added multiprocessing and multithreading (#209)
- Wait indicator
- Added Error line cap size
- Systems, subsystems and component selection options
- Added threshold to remove energetic terms or residues
- Correlation (#137)
    - Normal correlation using ΔG
    - Mutant correlation using ΔΔG
    - Energy values table
    - Interactive systems selection and Experimental Ki editor
    - New regression plot with bivariate and univariate graphs
    - Chart options panel
        - General options
        - Regression options
            - Scatter
            - Distribution

#### Documentation
- News section
- Animated command-line code block with terminal
- Added versioning
- Added installation using *.yml configuration file (#221)
- Added a new tutorial for correlation analysis (#194)

### Fix

#### gmx_MMPBSA_ana
- Close gmx_MMPBSA_ana properly when the selection dialog is rejected
- Compatibility with previous versions (`< v1.5.2`)
- Added several tooltips
- Improved general performance. See the [benchmark](analyzer.md#fast-like-a-rocket) (#209, #230, #243)
    - Improved the reading of output files and optimized data storage. (#228)
    - Eliminated the recalculation of IE and C2 entropies before opening the GUI. Now read from outputs files
    - Optimized data storage and access to subsets in panda's Dataframes
    - Now file processing and graphics generation does not freeze the GUI.
    - Removed redundant steps and data
    - Removed line graphs for components in the per-wise decomposition schema
    - Data access, processing, and storage are done in the API.
    - Added the option to temporarily store data on the hard disk instead of memory
    - Removed several pop-up windows
    - Decreased RAM consumption

- Removed H5 support 
- Invalid periodicity in the prmtop (#211)
- Buffer saturation (#213)


### Changes

#### gmx_MMPBSA
- Improved IE and C2 calculation

#### gmx_MMPBSA_ana

- Change the IE representation
- Removed summary widget
- Now the bar plot data contain the Average, SD, and SEM instead of all frames
- Chart Options
    - Error line color
    - Error line representation (SD or SEM)

- Restructured the init (system selection) dialog
    - Improved energy and decomposition options selection
    - Improved the correlation options
    - Added a new section for performance configuration
    - Added RAM consumption estimator
    - Improved the systems' tree representation and edition options


#### gmx_MMPBSA_test
- The examples folder moved up to the root directory

---

## gmx_MMPBSA v1.5.2 (03/23/2022)

### Additions

#### gmx_MMPBSA

- ALPB (Analytical Linearized Poisson-Boltzmann) approximation (#170)
- PC+ correction in 3D-RISM (#187)
- Warning when using -deo option without &decomp namelist in the input file (#191)

#### gmx_MMPBSA_ana

- Tooltips for different options (#151)

### Fix

####  gmx_MMPBSA_ana

- Chart options (#185)

### Changes

#### gmx_MMPBSA

- 3D-RISM calculated with sander instead of rism3d.snglpnt (#153)

#### Documentation

- Update 3D-RISM doc (#190)

## gmx_MMPBSA v1.5.1 (03/10/2022)

### Additions

#### `gmx_MMPBSA`

- Add new variables for QM/MM calculations (#171)
- Improvement of "print_res" function for &decomp and QM (#150)
- Support for OPLS force field (#160)
- New precalculated .xvv files for 3D-RISM

#### `gmx_MMPBSA_ana`

- Support for PyQt6

#### `Documentation`

- gmx_MMPBSA_ana documentation (v1.5.0) (#152)
- [Tutorial for psf_dcd files](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/#support-for-psf_dcd-files)
- [Tutorial for OPLS ff files](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/#oplsff-support)

### Fixes

#### `gmx_MMPBSA`

- Add ERROR when ligand.mol2 defined and no gaff or gaff2 in forcefields variable (#175)
- Fixed incorrect path search when the gmx_path is defined
- Fixed error when no ions are included in the groups defined

#### `gmx_MMPBSA_ana`

- Remove empty terms issue (#179)
- Fixed chart title for decomposition

## gmx_MMPBSA v1.5.0.3 (02/26/2022)
### Fixes
#### `gmx_MMPBSA`
- Improved nmode and sd in the output (#169)
- Improved autocompletion script
- Updated 3D-RISM, NMODE, Alascan, Decomp, and Linear/Non-Linear PB solver tutorials

## gmx_MMPBSA v1.5.0.1 (02/22/2022)
### Fixes
#### `gmx_MMPBSA`
- Allow index group name for group selection (#156)

#### `gmx_MMPBSA_ana`
- gmx_MMPBSA_ana does not open at the end of the calculation (#158)

## gmx_MMPBSA v1.5.0 (02/22/2022)

### Additions
#### `gmx_MMPBSA`
- Support for all files generated with CHARMM-GUI for Amber force fields (ff14SB and ff19SB)
- enabled NonLinear PB solver in `sander` (#63)
- enabled all `pbsa` options in `sander` (#64)
- enabled all `3D-RISM` variables (#68)
- Improve IE output data (#91)
- New C2 Entropy method added (#73)
- Option to create an [input file](input_file.md#generation-of-input-files-with-gmx_mmpbsa) (`--create_input`)
- New format file (h5) to store all the data result (`Experimental`)
- New method to get the `decomp` data without the API classes
- Report energy differences for every term in alanine scanning
- Automatic calculation of charge for `com` , `rec`, and `lig` groups in QM/MMGBSA
- Improve residue selection in QM/MMGBSA
- New modified PB Radii sets for GAFF (#140) and CHARMM (#141) force fields
- conda package (#55)
- Checking structure consistency method (#44)
#### `gmx_MMPBSA_ana`
- Chart properties selector. The following properties can be changed:
  - Font size:
    - Axes (Ticks, Labels)
    - Title, Subtitle and legend
    - Color bar
  - Pallete, color or theme:
    - Line plot
    - Bar plot
    - Heatmap plot
    - PyMOL visualization
  - Rotation, padding and number of ticks
  - Figure size, dpi and format
  - Chart type specific properties
  - Highlight or split components
- Frame range and interval selection
- Chart data can be visualized in table format and copy directly to excel
- Summary table for complex, receptor, ligand, etc.
- Now the outputs files can be visualized as document in their own sub-window
- New button to launch all the graphics and functions associated to item
- Highlighted Bar and Heatmap plot
- User setting file per-system
- Frames to Time conversion
#### Documentation
- Code block names and annotations

### Fixes
#### `gmx_MMPBSA`
- Protonated residues (#42)
- Inconsistent energy between tleap and ParmEd topology (#147)
- Error with alphanumeric residue numbers (#134)
- `print_res` selection scheme
- Error when `startframe=0` (#66)
- `parmchk2` always uses GAFF as force field (#45)
- IE is calculated for PB and GB independently (#72)
- Alanine scanning with CHARMM ff (#88)
- MT approach (#78)
- `get_num_terms` function run forever if TDC term not found. (#98)
- Check if the groups defined for receptor and ligand are the same (#86)
- Add SD and SEM calculated with [propagation formula](https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae) (#105)
- Structure consistency (#80, #79)
- gmx_MMPBSA launches an error when there is an `OverflowError` on IE calculation (#57)
- Improve `gmx_MMPBSA.log` file (#108)
- Inconsistency with multiple trajectories (#120)
- Alanine Scanning ERROR on THR to ALA mutation (#139)
#### `gmx_MMPBSA_ana`
- Tick labels in line plots (#65)
- Improved PyMOL 3D visualization (#85)
- Improved the system options in the init dialog
#### `gmx_MMPBSA_test`
- Error when `-f` option it not defined (#46)

### Changes
#### `gmx_MMPBSA`
- Recalculate the PB energy with --rewrite-output changing the value of `inp` (#144)
- Removed [deprecated variables](compatibility.md#variables)
- Input file format. Although it kept the structure of the previous version, the current one is more GROMACS alike
- `EnergyVector` changed to `ndarray` subclass
- Regen expression for `mutant_res`
- Now the COM, REC and LIG trajectories must have the same length when using MT approach
- Improve verbose logging in gmx_MMPBSA.log file (#108)
- Removed `*.gro` file support in `-cs` , `-rs` and `-ls` flags
- Added `trjconv` to avoid the PBC in the tpr file (#43)
#### `gmx_MMPBSA_ana`
- New set of chart buttons
- IE plot
#### `gmx_MMPBSA_test`
- Improved parallel processing
- Change command-line
#### Documentation
- Updated packages dependency
- Input and Output file pages

---

## [gxm_MMPBSA v1.4.3 (05/26/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.3)
### Additions
- Added two new tutorials [Protein_ligand_LPH_atoms_CHARMMff](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand_LPH_atoms_CHARMMff/) and [QM/MMGBSA calculations](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/examples/QM_MMGBSA/)
- Now the program reports the `p-value` associated with the correlation coefficient when performing the correlation analysis
- Google Analytics is used as a third-party tracking service to improve documentation. Check our [Private Policy](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/private_policy/) for more details

### Fixes
- Minor fixes in the documentation
- Improved parsing of `forcefields` variable
- Fixed bug when `gmx_MMPBSA_ana` runs without arguments
- Fixed compatibility with older files (`< v1.4.0`)
- Fixed error when `debug_printlevel > 1` in tleap command

### Changes
- Now the command line used is added to the log file 
- The `gmx_MMPBSA data` folder is exported directly rather than copied in the `$AMBERHOME` data folder.

---

## [gxm_MMPBSA v1.4.2 (05/01/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.2)
### Additions
- Added Covid-19 and other complex systems as examples in the documentation
- Added Q&A section to the documentation
- Implemented an adaptive `intdiel`(GB)/`indi`(PB) for Alanine scanning (Check the [input file](input_file.md#alanine_scanning-namelist-variables) section) 

### Fixes
- Fixed bug when `startframe = 0` (#33)
- Fixed bug when blank lines exist in [molecule] section in topology file
- Fixed pipe command-line for Gromacs execution in macOS 
- Fixed compatibility issues with v1.3.x
- Improved and fixed the documentation
- Improved output file information related to ΔG binding
- Improved calculation with different entropy approximations simultaneously

### Changes
- Changed `protein_forcefields` and `ligand_forcefield` by `forcefields` variable in all examples
- Now QH and IE can be calculated at the same time
- `entropy` variable was separated in `qh_entropy` and `interaction_entropy`. _The `entropy` variable is deprecated 
  and will be removed in the next major version (v1.5.0)._
- `entropy_seg` was replaced by `ie_segment`. _The `entropy_seg` variable is deprecated and will be removed in the 
  next major version (v1.5.0)._
  
---

## [gxm_MMPBSA v1.4.1 (04/08/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.1)
**Additions**

- New class `Residue` added to handle residues selection in Gromacs format with Amber index
- Verification of the presence of water molecules in receptor and ligand groups 
- Gromacs timer added

**Fixes**

- Gromacs topology conversion
- `qm_residues` notation
- Default path in `gmx_MMPBSA_test`
- The Entropy representation in `gmx_MMPBSA_ana`
- Bug when the structure has insertion code
- Improved ΔG Binding plot representation

**Changes**

- Now `forcefields` variable unifies `protein_forcefield` and `ligand_forcefield`. These variables `protein_forcefield`
  and `ligand_forcefield` are deprecated and will be removed in the next major version (v1.5.0).
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

---

## [gxm_MMPBSA v1.4.0 (03/22/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.4.0)
This release focuses almost entirely on `gmx_MMPBSA_ana` with minor issues fixed in `gmx_MMPBSA`

**Additions**

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
- New arguments flags for gmx_MMPBSA_ana (See the [gmx_MMPBSA_ana documentation](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/analyzer/))
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
- Embed YouTube videos for `gmx_MMPBSA_ana`

**Fixes**

- Now `gmx_MMPBSA_ana` shows stability results as expected
- Errors in the documentation
- MPI

**Changes**

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
  
---

## [gxm_MMPBSA v1.3.3 (03/09/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.3)
**Fixes**

* fixed Boltzmann constant for IE
* fix mutation in ligand
* fixed analyzer error when interval > 1
* fixed residue selection within
* fixed ChainID assignation when reference structure is defined
* fixed the selection to print when decomposition

---

## [gxm_MMPBSA v1.3.2 (03/01/2021)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.2) 
**Additions**

* **Now, gmx_MMPBSA is in Zenodo 
  [![DOI](https://zenodo.org/badge/295050575.svg)](https://zenodo.org/badge/latestdoi/295050575). You can refer to 
  us in this way in what we publish the article**
* Added Interaction Entropy to gmx_MMPBSA output file
* Added a new class to save IE in a csv file
* Added "Go to Top" button to documentation HTML.

**Fixes**
  
* Error when ligand and/or receptor are discontinuous (Testing it)
* Error when ligand and/or receptor are discontinuous and numbered non-consecutively
* Non-critical errors and inconsistencies in documentation

**Changes**
  
* `receptor_mask` and `ligand_mask` have been removed from input file variables. Now we extract the amber mask directly 
  based on the GROMACS index file
* The receptor and ligand mapping in the complex was improved. Now we use a method based on the GROMACS index file
* The method `Map` of the `system_MMPBSA` class has been restructured. Now always processes amber masks
* Changing the IE calculation function to a class

---

## [gmx_MMPBSA-v1.3.1](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.1)
**Additions**

* New variable (`overwrite_data`) to overwrite gmx_MMPBSA data. 
* More informative message when sander fail. Useful for PB calculation

**Fixes**

* Protein-ligand with charmm force field example
* Stability calculation
* gmx path error
* leaprc.GLYCAM_06h-1 file
* Protein-glycan example

**Changes**

* Documentation banner

---

## [gmx_MMPBSA-v1.3.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.3.0)

**Additions**

* Documentation at GitHub pages
* CHARMM force field support
* Amber topology generation from GROMACS topology. Work for CHARMM and Amber force fields
  * New flag for topologies (`-cp`, `-rp` and `-lp`) added
* Now `gmx_MMPBSA` supports discontinuous receptor and ligand.
* Glycine scanning
* Autocompletion script for both `gmx_MMPBSA` and `gmx_MMPBSA_ana`
* Versioneer to control the semantic version.
* Argument type checker for the command-line 

**Fixes**

* Alanine scanning tutorial 
* GROMACS executable path 
* The `-gui` option has been changed by `-nogui` and fixed when it is defined 
* Improvement on documentation

**Changes**

* Documentation theme. Now we use Material
* Alanine scanning variable. Now `mutant` correspond to mutant amino acid (ALA and GLY)
* The `gmx_MMPBSA_gui` was changed by `gmx_MMPBSA_ana` 
* Improvement on the topologies' construction process
* Order in which the trajectories are cleaned. Now, the topology is built and finally, the trajectories are cleaned

---

## [gmx_MMPBSA v1.2.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.2.0)
**Additions**

* New ligand force field (Zwitterionic amino acids)
* A new flag (-cr) added for defining a reference structure (guarantee a better consistency in generated PDB files)
* API documentation

**Fixes**

* More comprehensive output log file
* Best handling of structure files

**Changes**

`gmx` `editconf` is used to generate PDB files instead of `gmx` `trjconv` (#14)
`gmx_MMPBSA` data is copied in AMBER as an independent folder
*.gro files can be used as a MD Structure+mass(db) file
Updated tutorial list in README (Protein_DNA_RNA_Ion_ligand BFE calculations)

---

## [gmx_MMPBSA v1.1.1](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.1.1)
**Additions**

* New tutorial added (see Protein_DNA_RNA_Ion_ligand tutorial)
 
**Fixes**

* Support various metallo-complexes formats

**Changes**

* Keep all the temporary files in the folder for debugging purposes

---

## [gxm_MMPBSA v.1.1.0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v.1.1.0)
**Additions**

* Now supports carbohydrates as ligand. See this tutorial 
* Now supports metalloprotein-ligand complexes. See this tutorial.
* We have added data folder to gmx_MMPBSA module. This folder contains the GLYCAM_06h-1 force field files 
  (Compatible with amber99sb and early, see at http://glycam.org) which is not in AmberTools.

**Fixes**

* Minor bugs

**Changes**

* We changed the notation of the force fields, now the user can define any force field (We have only tested Amber and 
GLYCAM force fields) available in AmberTools. Charmm is not yet supported. See this section
  
---

## [gxm_MMPBSA first release](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/v1.0.0)
**This version includes:**

* Compatibility with different Gromacs versions
* Support for all types of calculations available in MMPBSA.py
* Graphical user interface for results analysis (gmx_MMPBSA_ana)
* API modified to get more information
* Some new facilities and types of calculations

---