---
template: main.html
title: Installation
---
# Installing gmx_MMPBSA

## Requirements
`gmx_MMPBSA` requires **[GROMACS](https://manual.gromacs.org/)** and
**[AmberTools](https://ambermd.org/AmberTools.php)** to be installed on your computer with **Python 3**.
For conda installations, Python `>=3.11,<3.13`, AmberTools `>=24.8,<27`, and GROMACS `>=2022,<2027`
are the recommended dependency boundaries. This keeps the environment compatible with the tested Python 3.12
stack without pinning users to one AmberTools or GROMACS release. `gmx_MMPBSA` supports a broad range of
GROMACS versions and should run with any GROMACS in the `PATH` that is compatible with the files you are using.

gmx_MMPBSA can be installed in two ways:

`Conda environment`
:   **Recommended, especially if you want to keep older versions of gmx_MMPBSA**. A conda environment provides a clean and efficient installation. It also allows you to keep
different versions of gmx_MMPBSA in isolated environments, reducing the possibility of incompatibility with
other packages. Installation time is also less since it does not require the compilation of AmberTools or GROMACS.

`AmberTools compilation`
:   This method assumes that AmberTools is compiled on your computer and that you want to use gmx_MMPBSA without
activating or deactivating a conda environment. You must also compile GROMACS, which increases installation time.
Because the installed packages must remain compatible, dependency errors are more common with this method.

!!! info "Installation"
    === "Conda environment"

        ??? Info "Install Miniconda"
            Install Miniconda on your computer:

            <div class="termy">

            ```bash
            $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

            $ chmod +x Miniconda3-latest-Linux-x86_64.sh

            $ ./Miniconda3-latest-Linux-x86_64.sh

            Successful Miniconda installation
            ```

            </div>

        === "`*.yml file`"    
            Assuming miniconda is installed, we can install gmx_MMPBSA using a yml file. 

            :material-file-download-outline:{:.heart } Download **[env.yml](env.yml)** file
    
            <div class="termy">

            ```console
            // Create a new environment and use the *.yml file to install dependencies
            $ conda env create --file env.yml

            // To use gmx_MMPBSA, just activate the environment
            $ conda activate gmxMMPBSA
            ```
                
            </div>

        === "`pip`"    
            Installing dependencies
    
            <div class="termy">

            ```console
            // Update conda
            $ conda update conda
            
            // Create a new environment and activate it
            $ conda create -n gmxMMPBSA python=3.12 -y -q
            $ conda activate gmxMMPBSA
            
            // Install mpi4py and AmberTools
            $ conda install -c conda-forge "mpi4py>=4.0.1,<5" "ambertools>=24.8,<27" -y -q

            // Install plotting dependencies
            $ conda install -c conda-forge "numpy>=1.26.4,<2" "matplotlib>=3.8,<4" "scipy>=1.14.1,<2" "pandas>=2.2,<3" "seaborn>=0.13,<0.14" -y -q

            // Install PyQt6 required to use the GUI analyzer tool (gmx_MMPBSA_ana). Not needed for HPC
            $ conda install -c conda-forge pyqt6 -y -q

            // (Optional) Install GROMACS
            $ conda install -c conda-forge "gromacs>=2022,<2027" pocl -y -q

            // Install gmx_MMPBSA
            $ python -m pip install gmx_MMPBSA
            ```
                
            </div>
    
    === "AmberTools compilation"
        [Follow the official AmberTools installation instructions for your OS](https://ambermd.org/Installation.php)
        !!! note
            We assume that AmberTools and its shell environment are configured correctly.
    
        **Installation**
        <div class="termy">
        ```console
        // Install gmx_MMPBSA
        $ amber.python -m pip install gmx_MMPBSA                                               
        ```
        </div>
    
        !!! danger
            If you get an error related to installing `mpi4py`, you may want to install this package manually from 
            `conda-forge` as follows:
    
            ```
            amber.conda install -c conda-forge mpi4py=3.1.3
            ```
            
            If you get an error related to `pip`, you may want to install this package manually as follows:
            
            ```
            amber.conda install pip
            ```

### Extra Dependencies
Some features require additional dependencies. Your operating system may also require one or more of the packages
listed below.

`ParmEd`
:  The current version of ParmEd implemented in AmberTools has some limitations that have been resolved in the [GitHub 
repository](https://github.com/ParmEd/ParmEd/tree/16fb2364c284f7c1dd716ee912c5c674b5d31e46) by its author Jason 
Swails and others with our help.

    Some of these limitations are:

    - Error reading topology when it has insertion codes
    - Error processing topologies generated with the Amber ff19SB force field
    - New PBRadii sets for GAFF and CHARMM force fields

!!! danger
    The gmx_MMPBSA installation process has been optimized to be as straightforward as possible. In rare cases, a 
    few extra dependencies may be needed.

`pip`
:   In some cases, the miniconda environment created in the AmberTools compilation does not have the `pip` module, so
any installation that depends on this package will fail. Required only if you did not install gmx_MMPBSA via `conda`

    ```
    amber.conda install pip    
    ```

`Git`
:   Used by **gmx_MMPBSA_test** to download the GitHub repository to get the examples' folder when running in the
default clone mode, or to install the development version. If you already have a local checkout, you can skip `git`
for testing by passing `--examples-dir /path/to/examples` (or setting `GMXMMPBSA_TEST_EXAMPLES_DIR`).

    ```
    conda install -c anaconda git
    ```
    or 
    ```
    sudo apt install git 
    ```

`mpi`
:  In some cases it is necessary to install the MPI dependencies. Required only if you did not install 
gmx_MMPBSA via `conda`
   ```
   sudo apt install openmpi-bin libopenmpi-dev openssh-client
   ```

`libxcb`
: If you get an error related to Qt plugins:
    ```
    sudo apt install --reinstall libxcb-xinerama0
    ```    
---

## Troubleshooting after installation

Once the installation is completed, the following warning may appear:

    WARNING: The scripts gmx_MMPBSA, gmx_MMPBSA_ana and gmx_MMPBSA_test 
    are installed in '/home/user/path_to_amber_install/amber20/miniconda/bin'
    which is not on PATH.

This warning appears because `pip` installs the executables (`gmx_MMPBSA`, `gmx_MMPBSA_ana`, and `gmx_MMPBSA_test`) in
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

see [issue 2][2] for the solution.

  [1]: https://ambermd.org/GetAmber.php#ambertools
  [2]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/2
  [3]: https://pypi.org/project/gmx-MMPBSA
  [4]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA


## Autocompletion script
Because `gmx_MMPBSA` has many options, shell autocompletion can improve productivity and reduce typing errors. The
provided script adds autocompletion for `gmx_MMPBSA`, `gmx_MMPBSA_ana`, and `gmx_MMPBSA_test`.


**Execution:**
Enter the following command in the terminal:
    
    source /path/to/ambertools/lib/python3.x/site-packages/GMXMMPBSA/GMXMMPBSA.sh

!!! tip
    If you want it to be activated automatically, add that command to your .bashrc

!!! warning
    * This script requires that `gmx_MMPBSA`, `gmx_MMPBSA_ana` and `gmx_MMPBSA_test` be accessible in PATH
    * If the command above fails, make sure the file has execute permission.
        
        On Ubuntu, Debian, Linux Mint or related:
        
        * GUI:

            * `Right-click on GMXMMPBSA.sh file` >

            * `Properties` > 

            * `Permissions` > 

            * `Select "Allow executing file as program"`
        
        * Terminal:
            
                chmod 755 /path/to/ambertools/lib/python3.x/site-packages/GMXMMPBSA/GMXMMPBSA.sh
    
        
**After sourcing `GMXMMPBSA.sh`, check that it works as follows:**

_All you have to do is enter the name of the program in the terminal and press the tab key twice:_
    
    gmx_MMPBSA <tab> <tab>

## Testing the operation of gmx_MMPBSA
After installing `gmx_MMPBSA`, verify the installation by following the
[`gmx_MMPBSA_test` instructions](examples/gmx_MMPBSA_test.md#running-gmx_mmpbsa_test).
