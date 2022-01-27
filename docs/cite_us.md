---
template: main.html
title: Cite us
---

<a href="https://www.scimagojr.com/journalsearch.php?q=5100155074&amp;tip=sid&amp;exact=no" title="SCImago Journal 
&amp; Country Rank"><img border="0" align="right" width=175 src="https://www.scimagojr.com/journal_img.php?id=5100155074" 
alt="SCImago Journal &amp; Country Rank"  /></a>

**If you found gmx_MMPBSA useful for your research, please cite:**

Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. gmx_MMPBSA: A New Tool to Perform 
End-State Free Energy Calculations with GROMACS. _Journal of Chemical Theory and Computation_, 2021 17 (10), 6281-6291. 
[https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645][1]. Download | [*.bib](gmx_MMPBSA_citation.bib)
| [*.ris](gmx_MMPBSA_citation.ris)

**Please also consider citing MMPBSA.py's paper:**

Bill R. Miller, T. Dwight McGee, Jason M. Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg. MMPBSA.
py: An Efficient Program for End-State Free Energy Calculations.  _Journal of Chemical Theory and Computation_, 
2012 8 (9), 3314-3321. [https://pubs.acs.org/doi/10.1021/ct300418h][2]. Download | [*.bib](MMPBSA_py_citation.bib)
| [*.ris](MMPBSA_py_citation.ris) | [*.xml](MMPBSA_py_citation.xml)

[1]: https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645
[2]: https://pubs.acs.org/doi/10.1021/ct300418h

## Example

!!! note "Important"
    This does not constitute by any means the only way to cite `gmx_MMPBSA` and programs/methods implemented in it. It 
    is just meant to serve as a guidance.

Here there is an example on how to cite `gmx_MMPBSA` and programs/methods implemented in it:

**MM/GBSA calculations**

PBC conditions were removed from GROMACS output trajectory before running the calculations with gmx_MMPBSA.[^1^][1]^,^[^2^][2]
Energetically relevant residues within 5 Å at the interface were predicted using the per-residue effective free energy 
decomposition (prEFED) protocol.[^3^][3] The AMBER99SB force field[^5^][5] was used to calculate the internal
term (ΔE~int~) as well as van der Waals (ΔE~vdW~) and electrostatic (ΔE~ele~) energies. The GB-Neck2 model 
(igb = 8)[^6^][6] was used to estimate the polar component of the solvation energy (ΔG~GB~) while the non-polar solvation
free energy (∆𝐺~𝑆𝐴~) was obtained by the equation:

<p align="center">
    ∆𝐺<sub>𝑆𝐴</sub> = 𝛾 · ∆𝑆𝐴𝑆𝐴 + 𝛽
</p>

where ∆𝑆𝐴𝑆𝐴 represents the solvent-accessible surface area variation of the solute molecule  upon complex formation, 
and 𝛾 and 𝛽 are empiric constants whose values for GB models are 0.0072 kcal·Å^-2^·mol^-1^ and 0, 
respectively.[^7^][7]^,^[^8^][8] The entropic term was calculated by the Interaction Entropy method.[^9^][9] The input 
file for gmx_MMPBSA decomposition calculation is shown below:

```
============================
Sample input file with decomposition analysis
&general
startframe=1750, endframe=2400, interval=1, PBRadii=4,
forcefields="oldff/leaprc.ff99SB"
/
&gb
igb=8, saltcon=0.150, intdiel=5,
/
&decomp
idecomp=2, dec_verbose=3,
print_res="within 5"
============================
```

Computational alanine scanning[^10^][10] was performed for five residues (TP62, EP68, EP70, HP76, and EP83) with a 
specific internal dielectric constant as suggested by Yan et al.[^11^][11] and according to the chemical-physical 
properties of the mutated amino acid (_i.e._, _e~i~_ = 5 for charged residues; _e~i~_ = 3 for polar residues; 
and _e~i~_ = 1 for hydrophobic residues). An example of the input file for gmx_MMPBSA alanine scanning 
(EP68 residue) calculation is shown below:

```
============================
Sample input file for alanine scanning analysis
&general
startframe=1750, endframe=2400, interval=1, PBRadii=4,
forcefields="oldff/leaprc.ff99SB", interaction_entropy=1, ie_segment=25, temperature=298
/
&gb
igb=8, saltcon=0.150,
/
&alanine_scanning
mutant='ALA', mutant_res='C:68', cas_intdiel=1
/
============================
```

[3]: https://onlinelibrary.wiley.com/doi/10.1002/jcc.20290
[5]: https://onlinelibrary.wiley.com/doi/10.1002/prot.21123
[6]: https://pubs.acs.org/doi/10.1021/ct3010485
[7]: https://pubs.acs.org/doi/10.1021/j100058a043
[8]: https://pubs.acs.org/doi/10.1021/jp073399n
[9]: https://pubs.acs.org/doi/10.1021/jacs.6b02682
[10]: https://www.sciencedirect.com/science/article/pii/S0022283603006107?via%3Dihub
[11]: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00734

## Authors

* Mario Sergio Valdés Tresanco, PhD Student. _University of Medellin, Colombia_
* Mario Ernesto Valdés Tresanco, PhD Student. _University of Calgary, Canada._
* Pedro Alberto Valiente, PhD. _University of Toronto, Canada_
* Ernesto Moreno Frías, PhD. _University of Medellin, Colombia_
