---
template: main.html
title: Installation
---
# Installing gmx_MMPBSA

## Requirements

`gmx_MMPBSA` requires **AmberTools20 or 21** to be installed in your machine with **Python3**, and the shell 
environment correctly set up for Amber. **The AmberTools suite is free of charge**. You can check [Amber web page][1] 
for a detailed installation guide. Of note, you can have more than one AmberTools installed in your machine. In case 
AmberTools20/21 is not the default Amber in your computer, just make sure to source AmberTools20/21 before 
installing/updating/running `gmx_MMPBSA`. `gmx_MMPBSA` also requires GROMACS (series 4.x.x or 5.x.x or 20xx.x) to be 
installed in your computer, and the shell environment correctly set up for GROMACS. `gmx_MMPBSA` has been tested 
with GROMACS 4.6.7, 5.1.2, 2018.3, and 2020.4, although it should run smoothly with any GROMACS present in the PATH 
and that is compatible with the files you are using.

### Dependencies
| Dependency      |     gmx_MMPBSA                                                    |   gmx_MMPBSA_ana   |  gmx_MMPBSA_test   |
|:----------------|:-----------------------------------------------------------------:|:------------------:|:------------------:|
| Python3         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium }        | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |
| AmberTools20/21 | :octicons-check-circle-fill-16:{ .req .scale_icon_medium }        | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |
| PyQt5           |                                            | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |
| Matplotlib      |                                            | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |
| Pandas          |                                            | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |
| Seaborn         |                                            | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |
| Scipy        |                                            | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |
| mpi4py          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |                    |
| openmpi-bin     | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |                    |
| libopenmpi-dev  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |                    |
| openssh-client  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |                    |                    |
| Git             | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium } |                    | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |


:octicons-check-circle-fill-16:{ .req } -> Required -- :octicons-check-circle-fill-16:{ .req_optrec } -> 
Optional, but recommended -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

!!! danger "Install Dependencies"
    
    === "MPI"
        Make sure that you install the OpenMPI library
        
            sudo apt install openmpi-bin libopenmpi-dev openssh-client
        
        !!! note ""
            `mpi4py` is installed with `gmx_MMPBSA`
    
    === "PyQt5"
        `gmx_MMPBSA_ana` requires PyQt5 which is not installed automatically. That way, `gmx_MMPBSA_ana` cannot be 
        opened in HPC. 

            amber.python -m pip install PyQt5

        !!! warning
            * Note that if you don't install PyQt5 you won't be able to open `gmx_MMPBSA_ana`
            * **Valid for versions > 1.4.0**

    === "Other Python Packages"
        `matplotlib`, `pandas`, (`mpi4py`), `seaborn` and `scipy` are installed automatically, so the user does not 
        have to worry about installing them

    



## Installation

=== "`stable`"
    You can install `gmx_MMPBSA` from the `stable` version on [PYPI][3]:

        amber.python -m pip install gmx_MMPBSA

    !!! danger
        The latest version of miniconda/anaconda seems to have a problem with `pip`. This is because, when `conda` is 
        updated, the latest version of `setuptools` (v.57.x) is installed, which is incompatible with `pip` (21.2.4), so 
        the latter was removed. If you find that `amber.pip`, `amber.python -m pip` or `python -m pip` (with conda 
        environment active) don't work, you must install it from `conda` as follows:

            conda install pip 
        or 

            amber.conda install pip

[comment]: <> (=== "`development`")

[comment]: <> (    or the `development` version from [GitHub][4]:)
    
[comment]: <> (        amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA)
    
[comment]: <> (    !!! warning)

[comment]: <> (        Make sure that you have `git` installed. If not you can install it as follows:)
    
[comment]: <> (            sudo apt install git)

## Update

If you already have installed a previous `gmx_MMPBSA` version, you can update it as follows:

=== "`stable`"
    `stable` version (recommended):
    
        amber.python -m pip install gmx_MMPBSA --upgrade

[comment]: <> (=== "`development`")

[comment]: <> (    `development` version from GitHub:)
    
[comment]: <> (        amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA --upgrade )
    
[comment]: <> (    !!! warning)

[comment]: <> (        Make sure that you have `git` installed. If not you can install it as follows:)

[comment]: <> (            sudo apt install git)

## Troubleshooting after installation

Once the installation is completed, the following warning may appear:

    WARNING: The scripts gmx_MMPBSA, gmx_MMPBSA_ana and gmx_MMPBSA_test 
    are installed in '/home/user/path_to_amber_install/amber20/miniconda/bin'
    which is not on PATH.

This warning is because `pip` installs the executables (`gmx_MMPBSA`, `gmx_MMPBSA_ana` and `gmx_MMPBSA_test`) in 
`installation_path/amber20/miniconda/bin`.

You have two options to solve this:

* Add this folder (*/amber20/miniconda/bin) to PATH:

        export PATH="/path_to_amber_install/amber20/miniconda/bin:$PATH"
    
    !!! tip
        * This option is more permanent and is recommended if you don't want to activate and deactivate the conda 
        environment
        * Make sure to update **path_to_amber_install** in the PATH variable

* Initializing the environment of conda amber:

        amber.conda init bash

    You can deactivate like this:
    
        conda deactivate

!!! note
    After using one of the above options, you should be able to run `gmx_MMPBSA`, `gmx_MMPBSA_ana` and `gmx_MMPBSA_test` 
    through the terminal

If when running `gmx_MMPBSA`, you get an error like this:

    ModuleNotFoundError: No module named 'parmed'

please see the following [issue][2] to see the solution

  [1]: https://ambermd.org/GetAmber.php#ambertools
  [2]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/2
  [3]: https://pypi.org/project/gmx-MMPBSA
  [4]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA


## Autocompletion script
Since `gmx_MMPBSA` has many flags, we believe that this autocompletion can significantly improve productivity, be
more user-friendly and reduce the number of unforced errors. That is why we created this script, which manages the
autocompletion of the `gmx_MMPBSA`, `gmx_MMPBSA_ana` and `gmx_MMPBSA_test`.


**Execution:**
Enter the following command in the terminal:
    
    source /path/to/ambertools/lib/python3.x/site-packages/GMXMMPBSA/GMXMMPBSA.sh

!!! tip
    If you want it to be activated automatically, add that command to your .bashrc

!!! warning
    * This script requires that `gmx_MMPBSA`, `gmx_MMPBSA_ana` and `gmx_MMPBSA_test` be accessible in PATH
    * If the command-line above end in error, please make sure the file has executed permissions. 
        
        On Ubuntu, Debian, Linux Mint or related:
        
        * GUI:

            * `Right-click on GMXMMPBSA.sh file` >

            * `Properties` > 

            * `Permissions` > 

            * `Mark the checkbox "Allow to execute the file as a program"`
        
        * Terminal:
            
                chmod 755 /path/to/ambertools/lib/python3.x/site-packages/GMXMMPBSA/GMXMMPBSA.sh
    
        
**Once you make the source of GMXMMPBSA.sh you can check its operation as follows:**

_All you have to do is enter the name of the program in the terminal and press the tab key twice:_
    
    gmx_MMPBSA <tab> <tab>

## Testing the operation of gmx_MMPBSA
After preparing everything to run `gmx_MMPBSA`, it only remains to check its correct operation. To know how to do it, 
consult the documentation of [`gmx_MMPBSA_test`](command-line.md#running-gmx_mmpbsa_test)
