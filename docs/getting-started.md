---
template: main.html
title: Getting started
---

[![pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![support](https://img.shields.io/badge/support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![python](https://img.shields.io/badge/python-v3.x-blue)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4626383.svg)](https://doi.org/10.5281/zenodo.4626383)

[<img src="../assets/logo.svg" height="120" width="240" align="right"/>]()

# Getting started
gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations with GROMACS 
files.

Most of the documentation below is found in the [Amber manual][1], we will point out what is new or different. 
Neither of these should be considered as a “black-box”, and users should be familiar with Amber and MM/PB(GB)SA 
method before at-tempting these sorts of calculations. These scripts automate a series of calculations, and cannot 
trap all the types of errors that might occur. 

!!! important
    We do not intend to replace the original [MMPBSA.py][2]; instead, we have implemented and improved some 
    functionalities, and what is most important, made this valuable tool available for GROMACS users. Most of the 
    documentation below is found in the [Amber manual][1], we will point out what is new or different.


  [1]: https://ambermd.org/doc12/Amber20.pdf#chapter.34
  [2]: https://pubs.acs.org/doi/10.1021/ct300418h
