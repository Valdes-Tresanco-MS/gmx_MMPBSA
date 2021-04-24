---
template: main.html
title: Q&A - Calculations
---

# `gmx_MMPBSA` Calculations
Here we describe a series of frequent issues related to calculations and their possible solutions.

!!! note 
    Most of the errors noted here are the result of inconsistent input files. Please read the documentation and make 
    sure your files are consistent.



??? example "Inconsistent energy"
    **I get this error: `ValueError: could not convert string to float: '*************'`**
    : This error has two possible causes
        
        1. The structure defined in `-cs`, `-ps`, or `-ls` options is inconsistent, or the trajectory
        has not been fitted (remove PBC) properly. This is the most common error. Many times it is because the 
        system under study is longer than some edges of the box.
        
            **Solution:**
            
            ??? tip "Check for structure consistency"
                
                Visualize the structure contained in the structure input file given in the `-cs`, `-rs`, or `-ls` 
                options and make sure it is consistent (as shown in Fig 1, right panel). On the other hand, if the 
                structure is "broken" (as shown in Fig 1, left panel) this will generate inconsistent results.
                
                **Generate the structure from tpr file**
                    
                    gmx editconf -f md.tpr -o md.pdb
    
            ??? tip "Make sure you have fitted the trajectory"

                Visualize the trajectory given in the `-ct`, `-rt`, or `-lt` options and make sure the PBC has been 
                removed (as shown in Fig 2, right panel). On the other hand, if the trajectory has not been fitted (as 
                shown in Fig 2, left panel) this will generate inconsistent results.
                
                Steps:
    
                1. Generate a group that contains both molecules
                    
                        gmx make_ndx -n index.ndx
                input: 1 | 12
                input: q

                    _Assuming 1 is the receptor and 12 is the ligand. This creates a new group (number 20 in this example)_
                
                2. remove the PBC
                    
                        gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center -n -ur compact
                center: 20 (created group)
                output: 0
                
                3. remove the rotation and translation with respect to the reference structure (optional)
                    
                        gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_fit.xtc -n -fit rot+trans
                fit: 20 (created group)
                output: 0
                    
                4. Visualization
                    
                    Make sure that the trajectory is consistent (as shown in Fig 2, right panel)

                <figure markdown="1">
                [![overview][2]][2]
                  <figcaption markdown="1" style="margin-top:0;">
                **Figure 2.** Vizualization of two different input trajectory files. Left: Trajectory with PBC; 
                Right: Trajectory centered, fitted and with PBC removed.
                  </figcaption>
                </figure>
                
                  [2]:../assets/images/q_a/traj_comp.gif
        
        2. You are trying to calculate the energetic contribution of a too-large group. Technically, the energy 
        value should not exceed 7 digits, so if you get a value higher than this, this error will occur. Although 
        `gmx_MMPBSA` can handle very large systems as described in example [Ribosomal50S_Mycalamide_A][9], it cannot 
        determine certain energetic terms. This is a limitation of Sander when writing the output file.
           
            !!! note
            The error could be solved by recompiling Sander with some modifications in the output function. 
            However, this is not recommended since the error can be large. Another possible solution could be modifying 
            the parameters of the calculation (solvent model, internal dielectric constant).

??? example "NMODE calculation finish in error"
    The only error reported is probably related to RAM saturation. NMODE calculations require a considerable amount 
    of RAM depending on the number of atoms in your system. The amount of total RAM consumed during the calculation 
    will be: `RAM for 1 frame * number of threads`    



  [9]: ../examples/Ribosomal50S_Mycalamide_A/README.md