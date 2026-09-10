---
template: main.html
title: Introduction
---

# Introduction

The MM/PB(GB)SA method can be used to calculate the binding free energies of noncovalently bound complexes.

<figure markdown="1">
![Thermodynamic cycle for binding free energy calculations](assets/images/cycle.png){ width=50% style="display: block; margin: 0 auto"}
  <figcaption markdown="1" style="margin-top:0;">
**Figure 1.** Thermodynamic cycle for binding free energy calculations
  </figcaption>
</figure>

[16]: assets/images/cycle.png

The binding free energy of a complex can be estimated as follows:

<p align="center">
    ∆𝐺<sub>𝑏𝑖𝑛𝑑</sub> = 〈𝐺<sub>𝐶𝑂𝑀</sub>〉−〈𝐺<sub>𝑅𝐸𝐶</sub>〉−〈𝐺<sub>𝐿𝐼𝐺</sub>〉
</p>
<p align="center">
    (1)
</p>

where each term on the right-hand side is given by:

<p align="center">
〈𝐺<sub>𝑥</sub>〉 = 〈𝐸<sub>𝑀𝑀</sub>〉 + 〈𝐺<sub>𝑠𝑜𝑙</sub>〉 − 〈𝑇𝑆〉
</p>
<p align="center">
    (2)
</p>

In turn, ∆𝐺<sub>𝑏𝑖𝑛𝑑</sub> can also be represented as:

<p align="center">
∆𝐺<sub>𝑏𝑖𝑛𝑑</sub> = ∆𝐻 − 𝑇∆𝑆
</p>
<p align="center">
    (3)
</p>

where ∆𝐻 corresponds to the enthalpy of binding and −𝑇∆𝑆 to the conformational entropy after ligand binding. When the
entropic term is omitted, the computed value is an enthalpy-like effective binding estimate. It can be useful for a
defined relative-comparison protocol, but omitting entropy is an approximation whose adequacy depends on the systems,
sampling, and scientific question; it is not generally sufficient by itself for relative affinity claims.

The enthalpy, ∆𝐻, can be decomposed into different terms:

<p align="center">
∆𝐻 = ∆𝐸<sub>𝑀𝑀</sub> + ∆𝐺<sub>𝑠𝑜𝑙</sub>
</p>
<p align="center">
    (4)
</p>

where:

<p align="center">
∆𝐸<sub>𝑀𝑀</sub> = ∆𝐸<sub>𝑏𝑜𝑛𝑑</sub> + ∆𝐸<sub>𝑎𝑛𝑔𝑙𝑒</sub> + ∆𝐸<sub>𝑑𝑖ℎ𝑒𝑑</sub> + ∆𝐸<sub>𝑣𝑑𝑊</sub> + ∆𝐸<sub>𝑒𝑙𝑒</sub> + ∆𝐸<sub>1-4 VDW</sub> + ∆𝐸<sub>1-4 EEL</sub>
</p>
<p align="center">
    (5)
</p>

The gas-phase contributions are calculated by `sander` within AmberTools according to the force field and method. The
ordinary terms above are supplemented by `UB`, `IMP`, and `CMAP` for applicable CHARMM calculations and by `ESCF` for
QM/MMGBSA. These terms are included when the output parser reports them. In ST, component differences for matching
topologies can cancel; that cancellation does not redefine the component totals or apply automatically to MT.

The ∆𝐺<sub>𝑠𝑜𝑙</sub> is given by:

<p align="center">
∆𝐺<sub>𝑠𝑜𝑙</sub> = ∆𝐺<sub>𝑝𝑜𝑙</sub> + ∆𝐺<sub>𝑛𝑜𝑛−𝑝𝑜𝑙</sub> = ∆𝐺<sub>𝑃𝐵/𝐺𝐵</sub> + ∆𝐺<sub>𝑛𝑜𝑛−𝑝𝑜𝑙</sub>
</p>
<p align="center">
    (6)
</p>

where:

<p align="center">
∆𝐺<sub>𝑛𝑜𝑛−𝑝𝑜𝑙𝑎𝑟</sub> = 𝑁𝑃<sub>𝑇𝐸𝑁𝑆𝐼𝑂𝑁</sub> ∗ ∆𝑆𝐴𝑆𝐴 + 𝑁𝑃<sub>𝑂𝐹𝐹𝑆𝐸𝑇</sub>
</p>
<p align="center">
    (7)
</p>

or,

<p align="center">
∆𝐺<sub>𝑛𝑜𝑛−𝑝𝑜𝑙</sub> = ∆𝐺<sub>𝑑𝑖𝑠𝑝</sub> + ∆𝐺<sub>𝑐𝑎𝑣𝑖𝑡𝑦</sub> = ∆𝐺<sub>𝑑𝑖𝑠𝑝</sub> + (𝐶𝐴𝑉𝐼𝑇𝑌<sub>𝑇𝐸𝑁𝑆𝐼𝑂𝑁</sub> ∗ 
∆𝑆𝐴𝑆𝐴 + 𝐶𝐴𝑉𝐼𝑇𝑌<sub>𝑂𝐹𝐹𝑆𝐸𝑇</sub>)
</p>
<p align="center">
    (8)
</p>

In the above equations, ∆𝐸<sub>𝑀𝑀</sub> corresponds to the molecular mechanical energy changes in the
gas phase. ∆𝐸<sub>𝑀𝑀</sub> includes ∆𝐸<sub>𝑏𝑜𝑛𝑑𝑒𝑑</sub>, also known as internal energy, and 
∆𝐸<sub>𝑛𝑜𝑛𝑏𝑜𝑛𝑑𝑒𝑑</sub>, corresponding to the van der Waals and electrostatic contributions. The solvation energy is 
determined differently depending on the method employed. In the 3D-RISM model, both the polar and nonpolar components
of the solvation energy are calculated. However, the PB and GB models estimate only the polar component of the 
solvation energy. The nonpolar component is usually assumed to be proportional to the molecule's total solvent-accessible
surface area (SASA), with a proportionality constant derived from experimental solvation energies of small nonpolar
molecules (Eq. 7). Alternatively, a modern approach that separates nonpolar solvation free energies into cavity and
dispersion terms can be used. In this approach, SASA is used to correlate the cavity term only, while a 
surface-integration method is employed to compute the dispersion term (Eq. 8).

Furthermore, the entropic component can be estimated with normal-mode analysis (NMODE). NMODE is Hessian-based: the
energy is minimized and a mass-weighted Hessian is diagonalized around the minimized structure to obtain vibrational
modes and frequencies. It is therefore distinct from quasi-harmonic (QH) analysis, which estimates fluctuations from
a coordinate covariance matrix over a sampled trajectory. NMODE can be computationally expensive, although truncated
systems can reduce the cost. New QH calculations are not supported in 1.7.0; historical
QH results remain readable for compatibility only.

Interaction Entropy (IE) estimates an entropic contribution from the fluctuation of the interaction energy along an MD
trajectory and has low additional post-processing cost. Its numerical behavior depends strongly on the distribution,
fluctuations, and convergence of the sampled interaction energies; it is not universally superior to NMODE and should
be checked with block or cumulative convergence diagnostics before interpretation. Multiple-trajectory IE/C2 use is
experimental in this release because the bound and unbound trajectories are independently sampled. See the
[GROMACS normal-mode reference](https://manual.gromacs.org/current/reference-manual/algorithms/normal-mode-analysis.html)
and [Ekberg and Ryde (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8389774/) for methodological context.

Typically, MM/PB(GB)SA calculations use one of two approaches: the single-trajectory protocol (STP) or the
multiple-trajectory protocol (MTP). In STP, both the receptor and ligand trajectories are extracted
from the complex trajectory. This approach is valid when the bound and unbound states of the receptor and ligand
are similar. It is computationally less expensive than the MTP approach since only a simulation of the complex is 
required. Additionally, the potential internal terms (_e.g._, bonds, angles, and dihedrals) cancel exactly in STP 
since these terms are the same in both bound and unbound states. On the other hand, the MTP is a more realistic 
approach because it considers separate trajectories for the complex, receptor, and ligand. However, subtracting
energies from independently sampled conformations can introduce substantial uncertainty. In practice, the system
must be studied carefully to select the appropriate approach.


## Literature
Further information can be found in [Amber manual][3]:

* [MMPBSA.py][4]
* [The Generalized Born/Surface Area Model][5]
* [PBSA][6]
* [Reference Interaction Site Model][7]
* [Generalized Born (GB) for QM/MM calculations][8]

and the foundational papers:

* [Srinivasan J. et al., 1998][9] 
* [Kollman P. A. et al., 2000][10] 
* [Gohlke H., Case D. A. 2004][11] 

as well as some reviews and expert opinions:

* [Genheden S., Ryde U. 2015][12] 
* [Wang et. al., 2018][13]  
* [Wang et. al., 2019][14]
* [Tuccinardi, 2021][16]

  [1]: https://pubs.acs.org/doi/10.1021/ct300418h
  [2]: https://pubs.acs.org/doi/abs/10.1021/jacs.6b02682

  [3]: https://ambermd.org/doc12/Amber21.pdf
  [4]: https://ambermd.org/doc12/Amber21.pdf#chapter.36
  [5]: https://ambermd.org/doc12/Amber21.pdf#chapter.4
  [6]: https://ambermd.org/doc12/Amber21.pdf#chapter.6
  [7]: https://ambermd.org/doc12/Amber21.pdf#chapter.7
  [8]: https://ambermd.org/doc12/Amber21.pdf#subsection.11.1.3
  [9]: https://pubs.acs.org/doi/abs/10.1021/ja981844+
  [10]: https://pubs.acs.org/doi/abs/10.1021/ar000033j
  [11]: https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.10379
  [12]: https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1032936
  [13]: https://www.frontiersin.org/articles/10.3389/fmolb.2017.00087/full
  [14]: https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.9b00055
  [15]: https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00418
  [16]: https://www.tandfonline.com/doi/pdf/10.1080/17460441.2021.1942836
