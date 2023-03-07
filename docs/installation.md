---
template: main.html
title: Installation
---
# Installing gmx_MMPBSA

## Requirements
`gmx_MMPBSA` requires **[GROMACS](https://manual.gromacs.org/) (series `4.x.x` or `5.x.x` or `20xx.x`)** and 
**[AmberTools20 or 21](https://ambermd.org/AmberTools.php)** to be installed in your machine with **Python3**.
`gmx_MMPBSA` has been tested with GROMACS `4.6.7`, `5.1.2`, `2018.3`, `2020.4` and `2022.4`, although it should run 
smoothly with any GROMACS present in the `PATH` and that is compatible with the files you are using.

## Installing gmx_MMPBSA v1.5.x

!!! danger
    gmx_MMPBSA v1.5.x includes a number of new functionalities and parts of the code have been completely rewritten, 
    hence it is incompatible with previous versions.

Currently, gmx_MMPBSA can be installed using two ways:

`conda environment`
:   The conda environment provides a clean and efficient way of installing gmx_MMPBSA. It also allows to have 
different versions of gmx_MMPBSA in isolated environments, thus reducing the possibility of incompatibility with 
other packages. Installation time is also less since it does not require the compilation of AmberTools or GROMACS. 
(**Recommended, especially if you want to keep older versions of gmx_MMPBSA**)

`AmberTools compilation`
:   In this way, we assume that you have AmberTools compiled on your machine and that you want to do an installation 
without worrying about enabling or disabling conda environments. It also involves user compilation of GROMACS, which 
takes considerable installation time. This way also requires installed packages to be compatible and installation 
errors are more frequent.

!!! info "Installation"
    === "conda environment"

        !!! Info "Important"
            Make sure to have conda installed in your computer. Check the third tab "Miniconda installation" for more 
            info. 

        === "`*.yml file`"    
            Installing gmx_MMPBSA using a yml file. 

            :material-file-download-outline:{:.heart } Download **[env.yml](env.yml)** file
    
            <div class="termy">

            ```console
            // Create a new environment and use the *.yml file to install dependencies
            $ conda env create -n gmxMMPBSA --file env.yml

            // To use gmx_MMPBSA, just activate the environment
            $ conda activate gmxMMPBSA
            ```
                
            </div>

            ??? note "Copy described intructions"     

                ``` bash 
                conda env create -n gmxMMPBSA --file env.yml                                    # (1)
                conda activate gmxMMPBSA                                                        # (2)
               
                ```
            
                1. Create the `gmxMMPBSA` environment and use the *.yml file to install dependencies
                2. Activate `gmxMMPBSA` environment


        === "`pip`"    
            Installing dependencies
    
            <div class="termy">

            ```console
            // Update conda
            $ conda update conda
            
            // Create a new environment and activate it
            $ conda create -n gmxMMPBSA python=3.9 -y -q 
            $ conda activate gmxMMPBSA
            
            // Install mpi4py and AmberTools
            $ conda install -c conda-forge mpi4py=3.1.3 ambertools=21.12 compilers=1.2.0 -y -q
            
            // Install updated version of ParmEd
            $ python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            
            // Install PyQt5 required to use the GUI analyzer tool (gmx_MMPBSA_ana). Not needed for HPC
            $ python -m pip install pyqt5

            // (Optional) Install GROMACS
            $ conda install -c conda-forge gromacs==2022.4 -y -q
            ```
                
            </div>

            ??? note "Copy described intructions"     

                ``` bash 
                conda update conda
                conda create -n gmxMMPBSA python=3.9 -y -q                                      # (1)
                conda activate gmxMMPBSA                                                        # (2)
                conda install -c conda-forge mpi4py=3.1.3 ambertools=21.12 compilers=1.2.0 -y -q      # (3)
                python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4 # (4)
                python -m pip install pyqt5                                                     # (5)
                # Optional
                conda install -c conda-forge gromacs==2022.4 -y -q                                 # (6)
               
                ```
            
                1. Create `gmxMMPBSA` environment
                2. Activate `gmxMMPBSA` environment
                3. Install dependencies
                4. Install ParmEd
                5. Install PyQt5 if you will use gmx_MMPBSA_ana
                6. (Optional) Install GROMACS if GROMACS is not installed in your machine

            === "Rolling/stable release"
                
                **INSTALLATION**
                <div class="termy">
                ```console
                // INSTALLATION
                $ python -m pip install gmx_MMPBSA
                ```
                </div>

                **UPDATE**
                <div class="termy">
                ```console
                // UPDATE
                $ python -m pip install gmx_MMPBSA -U
                ```
                </div>

                !!! info 
                    Install/update gmx_MMPBSA from PyPI. PyPI has the latest version of *gmx_MMPBSA* including stable 
                    and beta versions.
            
            === "development version" 

                **INSTALLATION**
                <div class="termy">
                ```console
                // INSTALLATION
                $ python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
                ```
                </div>

                **UPDATE**
                <div class="termy">
                ```console
                // UPDATE
                $ python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA -U
                ```
                </div>

                !!! warning
                    Install gmx_MMPBSA from the master branch of GitHub repository. This is only recommended 
                    for testing new versions or temporary solutions to reported bugs.
    
[comment]: <> (        === "`conda`")

[comment]: <> (            Stable version only)
            
[comment]: <> (            ``` bash )

[comment]: <> (            conda update conda)

[comment]: <> (            conda create -n gmxMMPBSA python=3.9 -y -q                                      # &#40;1&#41;        )

[comment]: <> (            conda activate gmxMMPBSA                                                        # &#40;2&#41;                        )

[comment]: <> (            conda install -c conda-forge mpi4py=3.1.3 ambertools=21.12 compilers -y -q                  # &#40;3&#41;)

[comment]: <> (            python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236          # &#40;4&#41;)

[comment]: <> (            python -m pip install pyqt5                                                     # &#40;5&#41;)

[comment]: <> (            # Optional)

[comment]: <> (            conda install -c conda-forge gromacs==2022.4 -y -q                                 # &#40;6&#41;)

[comment]: <> (            ```)
            
[comment]: <> (            1. Create `gmxMMPBSA` environment)

[comment]: <> (            2. Activate `gmxMMPBSA` environment)

[comment]: <> (            3. Install dependencies)

[comment]: <> (            4. Install ParmEd)

[comment]: <> (            5. Install PyQt5 if you will use gmx_MMPBSA_ana)

[comment]: <> (            6. &#40;Optional&#41; Install GROMACS if GROMACS is not installed in your machine)
    
[comment]: <> (            **INSTALLATION**)

[comment]: <> (            ```bash        )

[comment]: <> (            conda install -c conda-forge gmx_mmpbsa                                         # &#40;1&#41;    )

[comment]: <> (            ```)
            
[comment]: <> (            1. Install gmx_MMPBSA from conda-forge. This package will install all dependencies automatically)
            
[comment]: <> (            **UPDATE**)

[comment]: <> (            ```bash        )

[comment]: <> (            conda update gmx_mmpbsa    )

[comment]: <> (            ```)
    
        [Miniconda]: https://docs.conda.io/en/latest/miniconda.html
    
    === "AmberTools compilation"
        [Follow the oficial AmberTools installation according to your OS](https://ambermd.org/Installation.php)
        !!! note
            We asume that AmberTools and their shell environment are correctly configured
    
        === "Rolling/stable release"
            **INSTALLATION**
            <div class="termy">
            ```console
            // Install uodated ParmEd
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            // Install gmx_MMPBSA
            $ amber.python -m pip install gmx_MMPBSA                                               
            ```
            </div>

            **UPDATE**
            <div class="termy">
            ```console
            // Update gmx_MMPBSA
            $ amber.python -m pip install gmx_MMPBSA -U
            ```
            </div>
    
            !!! info 
                Install gmx_MMPBSA from PyPI PyPI has the latest version of *gmx_MMPBSA* including stable and beta
                versions.
            
        === "development version" 
            **INSTALLATION**
            <div class="termy">
            ```bash
            // Install updated ParmEd
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            // Install gmx_MMPBSA
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA     
            ```
            </div>

            **UPDATE**
            <div class="termy">
            ```bash
            amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA -U
            ```
            </div>

            !!! warning
                Install/update gmx_MMPBSA from the master branch of GitHub repository. This version is only recommended 
                to test a new version or to try temporary solutions to reported bugs.
    
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
            
    
    === "Miniconda Installation"
    
        Download and install [Miniconda]

        <div class="termy">

        ```bash
        $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        ---> 100%

        $ chmod +x Miniconda3-latest-Linux-x86_64.sh

        $ ./Miniconda3-latest-Linux-x86_64.sh
        ---> 100%

        Successful miniconda intallation
        ```

        </div>

        ??? note "Copy described intructions"     

            ``` bash 
            curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh   # (1)
            chmod +x Miniconda3-latest-Linux-x86_64.sh                                      # (2)
            ./Miniconda3-latest-Linux-x86_64.sh                                             # (3) 
           
            ```
        
            1. Download Miniconda installer
            2. Change permissions for the installer
            3. Execute and install miniconda

### Extra Dependencies
gmx_MMPBSA uses some dependencies for other functions independent of calculations or in some cases they may be 
necessary due to the nature of your OS.

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
:   Used by **gmx_MMPBSA_test** to download the GitHub repository to get the examples' folder or to install the 
development version.

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
: If you get an error related to Qt plugins
    ```
    sudo apt install --reinstall libxcb-xinerama0
    ```    
---

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
consult the documentation of [`gmx_MMPBSA_test`](examples/gmx_MMPBSA_test.md#running-gmx_mmpbsa_test)
