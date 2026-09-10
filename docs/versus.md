---
template: main.html
title: gmx_MMPBSA vs other programs
---

# Comparison of `gmx_MMPBSA` with other programs
This comparison is based on each program's documentation.


## Calculation features
| Feature                          |                 [g_mmpbsa][1]                  |             [GMXPBSA 2.1][2]              |              MMPBSA.py [^1]               |              [gmx_MMPBSA][3]              |
|:---------------------------------|:----------------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|
| **Normal binding free energies** |                       PB                       |                    PB                     |              PB [^0] and GB               |              PB [^0] and GB               |
| * GB models                      |                                                |                                           |             1, 2, 5, 7 and 8              |          1, 2, 5, 7, 8 and NSR6           |
| **Stability**                    |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Alanine scanning**             | :material-check-bold:{.scale_icon_medium} [^2] | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Entropy corrections** [^3]     |                                                |                                           |               NMODE and QH                |       NMODE, IE, and C2; legacy QH reader |
| **Decomposition schemes**        |                  Per-Residues                  |                                           |         Per-Residues and Per-Wise         |         Per-Residues and Per-Wise         |
| **QM/MMGBSA**                    |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **MM/3D-RISM**                   |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Membrane-protein support**     |                                                |                                           | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Approximations**               |                       ST                       |                 ST and MT                 |                 ST and MT                 |                 ST and MT                 |

## Analysis features
| Feature                         |                [g_mmpbsa][1]                 | [GMXPBSA 2.1][2] |              MMPBSA.py [^1]               |              [gmx_MMPBSA][3]              |
|:--------------------------------|:--------------------------------------------:|:----------------:|:-----------------------------------------:|:-----------------------------------------:|
| **API**                         |                                              |                  | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Analyzer Tool**               |                                              |                  |                   [^4]                    | :material-check-bold:{.scale_icon_medium} |
| * Multiple systems at once      |                                              |                  |                                           | :material-check-bold:{.scale_icon_medium} |
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
| **GROMACS Version**       |   4.x, 5.x and 2016+ [^6]    |           4.x, 5.x and 20xx.x [^7]           |                    ---                    |              `>=2022,<2027` [^10]        |
| **Dependencies**          | APBS (1.2.x, 1.3.x or 1.4.x) |              APBS (1.x.x) [^8]               |              AmberTools                   |             AmberTools [^9]               |
| **Parallel computation**  |       Depends on APBS        | Locally using APBS or in HPC divided in jobs | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |
| **Steps for:**            |                              |                                              |                                           |                                           |
| * Calculation and Summary |           Multiple           |                   Multiple                   |                    One                    |                    One                    |
| * Analysis                |           Multiple           |                   Multiple                   |                 Multiple                  |                    One                    |



  [^1]: [MMPBSA.py][4] is included in the AMBER package
  [^2]: Without documentation
  [^3]: NMODE = normal-mode approximation, QH = quasi-harmonic approximation, IE = interaction entropy
approximation, and C2 = C2 Entropy. In 1.7.0, new QH calculations are disabled;
historical QH result files remain readable only during the compatibility window.
  [^4]: We plan to extend gmx_MMPBSA compatibility to MMPBSA.py's results
  [^5]: The [AmberUtils][5] repository provides tools for analyzing the results
  [^6]: GROMACS 20xx.x is not officially supported. A pull request provides limited compatibility with versions
later than 2016.x
  [^7]: Support for GROMACS 20xx.x is not documented; the table assumes compatibility because the tool is script-based
  [^8]: Support for APBS 3.x.x is not documented
  [^9]: The recommended conda dependency boundary is AmberTools `>=24.8,<27`; older compatible AmberTools
versions may also work when their Python and compiled dependency stack is consistent.
  [^10]: This is the 1.7.0 tested environment boundary, not a claim that every GROMACS release is supported.
Conversion paths and force-field/model restrictions still apply; see [compatibility and upgrades](compatibility.md).
  [^0]: gmx_MMPBSA supports linear and nonlinear PB equations. [MMPBSA.py][4], by contrast, requires the user to
modify the `*.mdin` input files manually


  [1]: https://github.com/RashmiKumari/g_mmpbsa
  [2]: https://github.com/aspitaleri/gmxpbsa
  [3]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
  [4]: https://ambermd.org/doc12/Amber21.pdf#chapter.36
  [5]: https://github.com/williamdlees/AmberUtils
