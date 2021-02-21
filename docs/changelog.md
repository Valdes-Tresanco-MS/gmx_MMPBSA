---
template: main.html
title: Changelog
---

# Changelog
## Upcoming release (minor)
This version is almost entirely focused on gmx_MMPBSA_ana
### Additions
- New set of graphics (heatmap)
- Interactive visualization
- Multiple systems to analyze in the same session
- Calculations and plotting of the correlation between the systems

### Fixes
- Graphics improvements

## Upcoming release (path)
### Fixes
* Error when ligand and/or receptor are discontinuous (Testing it)

### Changes
* `receptor_mask` and `ligand_mask` have been removed from input file variables. Now we extract the amber mask directly 
  based on the GROMACS index file
* The receptor and ligand mapping in the complex was improved. Now we use a method based on the GROMACS index file
* The method Map of the system_MMPBSA class has been restructured. Now always processes amber masks


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

* New ligand force field (Zwitterionic aminoacids)
* A new flag (-cr) added for defining a reference structure (guarantee a better consintency in generated PDB files)
* API documentation

### Fixes
* More comprehensive output log file
* Best handling of structure files

### Changes
gmx editconf is used to generate PDB files instead of gmx trjconv (#14)
gmx_MMPBSA data is copied in AMBER as an independent folder
*gro files can be used as a MD Structure+mass(db) file
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