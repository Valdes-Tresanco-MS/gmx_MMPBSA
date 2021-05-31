---
template: main.html
title: Q&A - Installation
---

# Installation
Here we describe a series of more frequent reported problems related mainly to the installation process and their 
possible solutions.

???+ example "I don't find `amber.python` executable"

    #### **Solution:**
    
    1. Make sure that you have installed `Ambertools20` and sourced the `amber.sh(zch)` file
    2. If you installed `Ambertools20` from conda, use that python executable

???+ example "I get an error related to MPI when I try to install gmx_MMPBSA"
    If you get an error like this:    

        error: Cannot compile MPI programs. Check your configuration!!!
    
    #### **Solution:**

    Please try installing/reinstalling the OpenMPI library like this:

         sudo apt install openmpi-bin libopenmpi-dev openssh-client

    or this way

        sudo apt install --reinstall openmpi-bin libopenmpi-dev openssh-client

???+ example "I don't find the gmx_MMPBSA executable"
    
    #### **Solution:**

    1. Make sure that you have installed gmx_MMPBSA ([See here][1])
    2. Check if the miniconda bin folder are in the PATH ([See here][2])
    3. Check if the gmx_MMPBSA application has permission to run as a program
    
???+ example "When I run gmx_MMPBSA I get this error `ModuleNotFoundError: No module named 'parmed'`"

    #### **Solution:**

    Please see this [issue][3] to see the solution
    
    




  [1]: ../installation.md#installation
  [2]: ../installation.md#after-install
  [3]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/2