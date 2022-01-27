---
template: main.html
title: 
---

## Running gmx_MMPBSA_ana
In order to analyze multiple systems in the same section and implement the correlation between them, we improved the 
file input to gmx_MMPBSA_ana. Currently, gmx_MMPBSA_ana supports either several info files, a folder that contains the 
info file, or a list of folders containing info files.

Make sure to use the following structure when working with several folders:

=== "Valid folder structure"
    
    Here all folders (5 systems) will be processed by `gmx_MMPBSA_ana`
    ```
    Defined folder
      ├── System-1
      │      └── _GMXMMPBSA_info
      ├── System-2
      │      └── _GMXMMPBSA_info
      ├── System-3
      │      └── _GMXMMPBSA_info
      ├── System-4
      │      └── _GMXMMPBSA_info
      └── System-5
             └── _GMXMMPBSA_info
    ```

=== "Valid folders + files structure"

    Here all systems (4) will be processed by `gmx_MMPBSA_ana`
    ```
    System-1      
      └──_GMXMMPBSA_info
    Folder-1
      ├── System-2.1
      │      └── _GMXMMPBSA_info
      └── System-2.2
             └── _GMXMMPBSA_info
    System-3 (_GMXMMPBSA_info)
    ```
    
    !!! tip ""
        The command-line for this is:

            gmx_MMPBSA_ana -f /path/to/System-1 /path/to/System-3/_GMXMMPBSA_info /path/to/Folder-1 -r
        
        This estructure only work if recursive option was defined. See the examples below
    
    !!! note ""
        Note that:

        * `System-1` is defined as folder that contain a _GMXMMPBSA_info
        * Folder-1 contain two folder (Systems), each containing a _GMXMMPBSA_info
        * System-3 is defined as a _GMXMMPBSA_info file


=== "Wrong folder structure"
    
    Here only  4 (systems) folders will be processed by gmx_MMPBSA_ana. The systems in the `Internal folder` will be 
    ignored 
    ```
    Defined folder
      ├── System-1
      │      └── _GMXMMPBSA_info
      ├── System-2
      │      └── _GMXMMPBSA_info
      ├── Internal folder
      │      ├─X─ System-3.1
      │      │      └── _GMXMMPBSA_info
      │      └─X─ System-3.2
      │             └── _GMXMMPBSA_info
      ├── System-4
      │      └── _GMXMMPBSA_info
      └── System-5
             └── _GMXMMPBSA_info
    ```

!!! examples
        
    === "One file"
        
        Passing a _GMXMMPBSA_info file as input:
        
        * Current directory
        
                gmx_MMPBSA_ana -f _GMXMMPBSA_info
        
        * other location  :material-new-box:{: .medium .heart } Version: 1.4.0

                gmx_MMPBSA_ana -f /path/to/_GMXMMPBSA_info

    === "One folder" 
        :material-new-box:{: .medium .heart } Version: 1.4.0

        Passing a folder as input:  
        
        * Current directory
        
                gmx_MMPBSA_ana -f .
        
        * other location

                gmx_MMPBSA_ana -f /path/to/folder

        !!! warning "Remember"
            This folder must contain a valid `_GMXMMPBSA_info` file
    
    === "Files + Folders"
        :material-new-box:{: .medium .heart } Version: 1.4.0
        
            gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/_GMXMMPBSA_info-1 /path/to/folder-2
        
        !!! warning "Remember"
            * All defined folders must contain a valid `_GMXMMPBSA_info` file
            * All `_GMXMMPBSA_info` files defined must be valid
    
    === "Recursive option"
        :material-new-box:{: .medium .heart } Version: 1.4.0

        Passing a folder as input with `recursive` option:  
        
        * Current directory
        
                gmx_MMPBSA_ana -f . -r
        
        * other location

                gmx_MMPBSA_ana -f /path/to/folder --recursive

        * combine multiple folders

                gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/folder-2 /path/to/folder-3  -r

            !!! warning ""
            Folders can contain one or more systems

        * combine multiple folders and files

                gmx_MMPBSA_ana -f /path/to/folder-1 /path/to/_GMXMMPBSA_info-1 /path/to/folder-3  -r

            !!! warning ""
            * Folders can contain one or more systems
            * Note that if you remove the option -r, each folder must contain a valid _GMXMMPBSA_info file.