---
template: main.html
title: Q&A - OS issues
---

# Operative System compatibility issues
Here we describe a series of more frequent reported problems related mainly to the compatibility issues and 
supported OS and its possible solutions.

??? "gmx_MMPBSA fail in macOS < BigSur"
    
    We have only tested gmx_MMPBSA on macOS BigSur and it works fine. Since a problem related to a version lower 
    than BigSur was reported, we assume that it is not compatible. Please check this [thread][1]

??? "gmx_MMPBSA_ana fail in Debian, Centos or HPC"

    **I am using Debian or Centos as OS**
    : Since Debian and Centos are Linux distributions for servers, their main focus is stability, which is why they 
    keep old versions of most libraries. In addition, its GUI is very basic, as it is generally not very used. 
    `gmx_MMPBSA_ana` uses the latest PyQt5 version which requires updated graphics libraries. That is why we recommend 
    the use of more desktop-focused distributions such as Ubuntu, Linux Mint, Fedora, OpenSuse, Majaro, etc.
    
        !!! warning
            Please note that we cannot support the huge number of Linux distributions out there today.
    
    **I am running gmx_MMPBSA on HPC**
    : By default, HPCs have the GUI disabled, therefore any application that depends on the graphics libraries will 
    end in error. `gmx_MMPBSA` has the `-nogui` option to avoid this type of errors since it prevents `gmx_MMBSA_ana` from 
    being executed at the end of the calculations
    
    **I get an error related to Qt plugins**
    : If you get the following or similar error:

            qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
    
        A possible solution may be to reinstall these libraries
            
            sudo apt install --reinstall libxcb-xinerama0


  [1]: https://groups.google.com/g/gmx_mmpbsa/c/bk-PZl4hZzo