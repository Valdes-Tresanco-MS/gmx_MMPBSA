---
template: main.html
title: Command-line
---

## `gmx_MMPBSA_ana` command-line
<div class="termy">
    ```console
    $ gmx_MMPBSA_ana -h
    
    usage: run_ana.py [-h] [-v] [-f [FILES [FILES ...]]] [-r]
    
    This program is part of gmx_MMPBSA and will show a workspace to analyze the 
    gmx_MMPBSA results
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    Info file:
      -f [FILES [FILES ...]], --files [FILES [FILES ...]]
                            gmx_MMPBSA info files or container folder or list of 
                            them (default: [.] Current working dir)
      -r, --recursive       Search recursively in this folder at depth = 1 
                            (default: False)
    
    gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
    Based on MMPBSA.py (version 16.0) and AmberTools20
    ```
</div>