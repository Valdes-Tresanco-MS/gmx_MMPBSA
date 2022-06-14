---
template: main.html
title: gmx_MMPBSA vs other programs
---

# Comparison of `gmx_MMPBSA` vs other programs
This comparison is based on the documentation of the different programs


## Calculation features
| Feature                          |                 [g_mmpbsa][1]                  |             [GMXPBSA 2.1][2]              |              MMPBSA.py [^1]               |              [gmx_MMPBSA][3]              |
|:---------------------------------|:----------------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|
| **Normal binding free energies** |                       PB                       |                    PB                     |              PB [^0] and GB               |              PB [^0] and GB               |
| * GB models                      |                                                |                                           |             1, 2, 5, 7 and 8              |          1, 2, 5, 7, 8 and NSR6           |
| **Stability**                    |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Alanine scanning**             | :material-check-bold:{.scale_icon_medium} [^2] | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Entropy corrections** [^3]     |                                                |                                           |               NMODE and QH                |           NMODE, QH, IE, and C2           |
| **Decomposition schemes**        |                  Per-Residues                  |                                           |         Per-Residues and Per-Wise         |         Per-Residues and Per-Wise         |
| **QM/MMGBSA**                    |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **MM/3D-RISM**                   |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Support Membrane Protein**     |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Approximations**               |                       ST                       |                 ST and MT                 |                 ST and MT                 |                 ST and MT                 |

## Analysis features
| Feature                         |                [g_mmpbsa][1]                 | [GMXPBSA 2.1][2] |              MMPBSA.py [^1]               |              [gmx_MMPBSA][3]              |
|:--------------------------------|:--------------------------------------------:|:----------------:|:-----------------------------------------:|:-----------------------------------------:|
| **API**                         |                                              |                  | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Analyzer Tool**               |                                              |                  |                   [^4]                    | :material-check-bold:{.scale_icon_medium} |
| * Multiple systems at same time |                                              |                  |                                           | :material-check-bold:{.scale_icon_medium} |
| * Correlation between systems   |  :material-check-bold:{.scale_icon_medium}   |                  |                                           | :material-check-bold:{.scale_icon_medium} |
| * Per-residue energies to PDB   |  :material-check-bold:{.scale_icon_medium}   |                  |                                           | :material-check-bold:{.scale_icon_medium} |
| * Interactive visualization     |                                              |                  |                                           | :material-check-bold:{.scale_icon_medium} |
| ** _3D Molecular Visualization_ |                                              |                  |                                           |                   PyMOL                   |
| ** _Interactive Charts_         |                 static image                 |                  |                                           | :material-check-bold:{.scale_icon_medium} |
| * Plotting tool                 |                internal tools                |                  |       API and graphics library [^5]       |              gmx_MMPBSA_ana               |
| * Energetic Terms charts        | ΔG~polar~, ΔG~nonpolar~, ΔE~MM~ and ΔG~bind~ |                  |                                           |                    All                    |
| * Export data to CSV file       |                                              |                  | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| ** _Energy Summary_             |                                              |                  | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| ** _Individual Energetic Terms_ |                                              |                  |                                           | :material-check-bold:{.scale_icon_medium} |

## Technical features
| Feature                   |        [g_mmpbsa][1]         |               [GMXPBSA 2.1][2]               |              MMPBSA.py [^1]               |              [gmx_MMPBSA][3]              |
|:--------------------------|:----------------------------:|:--------------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|
| **GROMACS Version**       |   4.x, 5.x and 2016+ [^6]    |           4.x, 5.x and 20xx.x [^7]           |                    ---                    |            4.x, 5.x and 20xx.x            |
| **Dependencies**          | APBS (1.2.x, 1.3.x or 1.4.x) |              APBS (1.x.x) [^8]               |              AmberTools20/21              |           AmberTools20/21 [^9]            |
| **Parallel computation**  |       Depends on APBS        | Locally using APBS or in HPC divided in jobs | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Steps for:**            |                              |                                              |                                           |                                           |
| * Calculation and Summary |           Multiple           |                   Multiple                   |                    One                    |                    One                    |
| * Analysis                |           Multiple           |                   Multiple                   |                 Multiple                  |                    One                    |



  [^1]: [MMPBSA.py][4] is included in AMBER package
  [^2]: Without documentation
  [^3]: NMODE = Normal modes approximation, QH = Quasic-Harmony approximation, IE = Interaction Entropy
approximation, and C2 = C2 Entropy
  [^4]: We plan to extend gmx_MMPBSA compatibility to MMPBSA.py's results
  [^5]: Currently there is a repository ([AmberUtils][5]) for analysing the results
  [^6]: GROMACS 20xx.x is not officially supported. There is a Pull Request that offers a minimum compatibility 
with versions higher than 2016.x one, but still with limitations
  [^7]: It is not clear whether it does support GROMACS versions 20xx.x or not, but we assume that it does since 
it is script-based
  [^8]: It is not clear whether it does support APBS versions 3.x.x
  [^9]: gmx_MMPBSA is compatible with AmberTools20, AmberTools21, and AmberTools22
  [^0]: gmx_MMPBSA supports Linear and Non-Linear PB equations. [MMPBSA.py][4] on the other hand, requires the user to 
modify manually the *.mdin input files 
  

  [1]: https://github.com/RashmiKumari/g_mmpbsa
  [2]: https://github.com/aspitaleri/gmxpbsa
  [3]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
  [4]: https://ambermd.org/doc12/Amber21.pdf#chapter.36
  [5]: https://github.com/williamdlees/AmberUtils

  