---
template: main.html
title: Q&A - OS issues
---

# Operating-system compatibility issues
This page describes commonly reported operating-system compatibility problems and possible solutions.

???+ example "gmx_MMPBSA fails on macOS versions earlier than Big Sur"

    #### **Solution:**    
    We tested gmx_MMPBSA on macOS Big Sur. Problems have been reported on earlier versions; see this [thread][1]
    for details.

???+ example "gmx_MMPBSA_ana fails in non-native Linux distribution (_i.e._, Windows Subsystem for Linux (WSL), Debian, Centos or HPC)"

    #### **Solution:**
    **I am using Windows Subsystem for Linux (WSL)**
    : A WSL installation without Linux GUI support cannot run `gmx_MMPBSA_ana`. Use `gmx_MMPBSA` with `-nogui` in
    that environment.

    **I am using Debian or Centos as OS**
    : Server-oriented distributions may provide older graphics libraries. `gmx_MMPBSA_ana` requires graphics libraries
    compatible with its PyQt version. A desktop-oriented distribution may therefore be easier to configure.
    
    **I am running gmx_MMPBSA on HPC**
    : HPC compute nodes commonly lack a graphical display, so applications that require graphics libraries may fail.
    Add `-nogui` to the gmx_MMPBSA command to prevent `gmx_MMPBSA_ana` from opening after the calculation.
    
???+ example "Error with `qt.qpa.plugin`"

    **I get an error related to Qt plugins**
    : If you get the following or similar error:

            qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
        
        #### **Solution:**    
        One possible solution is to reinstall the following library:
            
            sudo apt install --reinstall libxcb-xinerama0


  [1]: https://groups.google.com/g/gmx_mmpbsa/c/bk-PZl4hZzo