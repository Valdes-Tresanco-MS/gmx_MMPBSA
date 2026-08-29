---
template: main.html
title: Q&A - Installation
---

# Installation
This page describes commonly reported installation problems and possible solutions.

???+ example "I cannot find the `amber.python` executable"

    #### **Solution:**
    
    1. Make sure that you have installed `Ambertools20` and sourced the `amber.sh(zch)` file
    2. If you installed `AmberTools20` from conda, use the Python executable from that environment

???+ example "I get an error related to MPI when I try to install gmx_MMPBSA"
    If you get an error like this:    

        error: Cannot compile MPI programs. Check your configuration!!!
    
    #### **Solution:**

    Try installing the OpenMPI library:

         sudo apt install openmpi-bin libopenmpi-dev openssh-client

    or reinstall it:

        sudo apt install --reinstall openmpi-bin libopenmpi-dev openssh-client

???+ example "I cannot find the gmx_MMPBSA executable"
    
    #### **Solution:**

    1. Make sure that you have installed gmx_MMPBSA ([See here][1])
    2. Check whether the Miniconda `bin` folder is in `PATH` ([see here][2])
    3. Check whether the gmx_MMPBSA executable has execute permission
    
???+ example "When I run gmx_MMPBSA I get this error `ModuleNotFoundError: No module named 'parmed'`"

    #### **Solution:**

    See [issue 2][3] for the solution.
    
    




  [1]: ../installation.md#installing-gmx_mmpbsa
  [2]: ../installation.md#troubleshooting-after-installation
  [3]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/2
