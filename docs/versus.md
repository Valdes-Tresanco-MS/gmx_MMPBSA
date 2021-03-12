---
template: main.html
title: gmx_MMPBSA vs other programs
---

# Comparison of `gmx_MMPBSA` vs other programs
This comparison is based on the documentation of the different programs


## Calculation features
| Feature                               |        [g_mmpbsa][1]         |       [GMXPBSA 2.1][2]       |                     MMPBSA.py[^1]                     |      [gmx_MMPBSA][3]      |
|:--------------------------------------|:----------------------------:|:----------------------------:|:-----------------------------------------------------:|:-------------------------:|
| **Normal binding free energies**      |              PB              |              PB              |                       PB and GB                       |         PB and GB         |
| * GB models                           |                              |                              |                  1, 2, 3, 5, 7 and 8                  |    1, 2, 3, 5, 7 and 8    |
| **Stability**                         |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Alanine scanning**                  |   :heavy_check_mark: [^2]    |      :heavy_check_mark:      |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Entropy corrections** [^3]          |                              |                              |                     nmode and qh                      |     nmode, qh, and IE     |
| **Decomposition schemes**             |         Per-Residues         |                              |               Per-Residues and Per-Wise               | Per-Residues and Per-Wise |
| **QM/MMGBSA**                         |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **MM/3D-RISM**                        |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Membrane Protein MMPBSA**           |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Approximations**                    |              ST              |              ST              |                       ST and MT                       |         ST and MT         |

## Analysis features
| Feature                               |        [g_mmpbsa][1]         |       [GMXPBSA 2.1][2]       |                     MMPBSA.py[^1]                     |      [gmx_MMPBSA][3]      |
|:--------------------------------------|:----------------------------:|:----------------------------:|:-----------------------------------------------------:|:-------------------------:|
| **API**                               |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Analyzer Tool**                     |                              |                              |                         [^4]                              |    :heavy_check_mark:     |
| * Multiple systems at same time       |                              |                              |                                                       |    :heavy_check_mark:     |
| * Correlation between systems         |      :heavy_check_mark:      |                              |                                                       |    :heavy_check_mark:     |
| * Energy to PDB                       |      :heavy_check_mark:      |                              |                                                       |    :heavy_check_mark:     |
| * Interactive visualization           |                              |                              |                                                       |    :heavy_check_mark:     |
|   ** _PDB_                            |                              |                              |                                                       |           PyMOL           |
|   ** _Charts_                         |        static image          |                              |                                                       |    :heavy_check_mark:     |
| * Energetic Terms charts              |       internal tools         |                              |               API and graphics library [^5]           |      gmx_MMPBSA_ana       |
| * Export to CSV file                  |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
|   ** _Energy Summary_                 |                              |                              |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
|   ** _Individual Energetic Terms_     |                              |                              |                                                       |    :heavy_check_mark:     |

## Technical features
| Feature                               |        [g_mmpbsa][1]         |       [GMXPBSA 2.1][2]       |                     MMPBSA.py[^1]                     |      [gmx_MMPBSA][3]      |
|:--------------------------------------|:----------------------------:|:----------------------------:|:-----------------------------------------------------:|:-------------------------:|
| **GROMACS Version**                   |   4.x, 5.x and 2016+ [^6]    |   4.x, 5.x and 20xx.x [^7]   |                          ---                          |    4.x, 5.x and 20xx.x    |
| **Externals programs**                | APBS (1.2.x, 1.3.x or 1.4.x) |      APBS (1.x.x) [^8]       |                     AmberTools20                      |       AmberTools20        |
| **Parallel computation**              |   Depends on APBS version    |      :heavy_check_mark:      |                  :heavy_check_mark:                   |    :heavy_check_mark:     |
| **Steps for:**                        |                              |                              |                                                       |                           |
| * Calculation and Summary             |           Multiple           |           Multiple           |                          One                          |            One            |
| * Analysis                            |           Multiple           |           Multiple           |                       Multiple                        |            One            |



  [^1]: [MMPBSA.py][4] is included in AMBER package.
  [^2]: Without documentation.
  [^3]: nmode = Normal modes approximation, qh = Quasic-Harmony approximation and IE = Interaction Entropy
approximation
  [^4]: We plan that gmx_MMPBSA will have compatibility with the results of MMPBSA.py
  [^5]: Currently there is a repository ([AmberUtils][5]) quite robust for the analysis of the results obtained.
  [^6]: GROMACS 20xx.x is not officially supported. There is a Pull Request that offers a minimum of compatibility 
with versions higher than 2016.x but with limitations
  [^7]: It is not clear if it supports the GROMACS versions 20xx, but we assume that it does because 
it is script-based.
  [^8]: It is not clear if it supports the APBS versions 3.x.x , but we assume that it does because it is 
script-based.
  

  [1]: https://github.com/RashmiKumari/g_mmpbsa
  [2]: https://github.com/aspitaleri/gmxpbsa
  [3]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
  [4]: https://ambermd.org/doc12/Amber20.pdf#chapter.34
  [5]: https://github.com/williamdlees/AmberUtils

  