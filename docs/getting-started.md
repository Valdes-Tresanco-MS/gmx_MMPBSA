---
template: main.html
title: Getting started
---

[![pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![support](https://img.shields.io/badge/support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![python](https://img.shields.io/badge/python-v3.x-blue)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4569307.svg)](http://doi.org/10.5281/zenodo.4569307)

[<img src="../assets/logo.svg" height="120" width="240" align="right"/>]()

# Getting started

gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations with GROMACS 
files.

## Citing us

At the moment we only have Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4569307.svg)](https://doi.org/10.5281/zenodo.4569307)

!!! note "Cite as"
    Mario S. Valdés Tresanco, Mario E. Valdes-Tresanco, Pedro A. Valiente, & Ernesto Moreno Frías.  
    gmx_MMPBSA (Version v1.4.2). Zenodo. http://doi.org/10.5281/zenodo.4569307

**We will have the article available soon**

Please also consider citing MMPBSA.py's paper:

**MMPBSA.py: An Efficient Program for End-State Free Energy Calculations**. Bill R. Miller, T. Dwight McGee, Jason M.
Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg. _Journal of Chemical Theory and Computation_, 2012 8 
(9), 3314-3321. [**DOI:** 10.1021/ct300418h][1]

## Authors

* Mario Sergio Valdés Tresanco, PhD Student. _University of Medellin, Colombia_
* Mario Ernesto Valdés Tresanco, PhD Student. _University of Calgary, Canada._
* Pedro Alberto Valiente, PhD. _University of Toronto, Canada_
* Ernesto Moreno Frías, PhD. _University of Medellin, Colombia_

## About the gmx_MMPBSA implementation
Most of the documentation below is found in the [Amber manual][2], we will point out what is new or different. 
Neither of these should be considered as a “black-box”, and users should be familiar with Amber and MM/PB(GB)SA 
method before at-tempting these sorts of calculations. These scripts automate a series of calculations, and cannot 
trap all the types of errors that might occur. 

!!! important
    We do not intend to replace the original [MMPBSA.py][3]; instead, we have implemented and improved some 
    functionalities, and what is most important, made this valuable tool available for GROMACS users. 

  [1]: https://doi.org/10.1021/ct300418h
  [2]: https://ambermd.org/doc12/Amber20.pdf#chapter.34
  [3]: https://pubs.acs.org/doi/10.1021/ct300418h
