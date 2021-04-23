---
template: main.html
title: Q&A - Calculations
---

# `gmx_MMPBSA` Calculations
Here we describe a series of more frequent reported problems related mainly to the preparation, execution, 
calculations, results issues its possible solutions.

!!! note 
    Most of the errors noted here are the result of inconsistent input files. Please read the documentation and make 
    sure your files are consistent. Although we have optimized gmx_MMPBSA to handle a wide variety of systems, as well 
    as the efficient and consistent handling of structures, it will not be able to recognize unforced errors that are 
    the responsibility of the user.


??? "Inconsistent energy"
    **I get this error: `ValueError: could not convert string to float: '*************'`**
    : This error has two possible causes
        
        1. The structure defined in `-cs`, `-ps`, or `-ls` options are inconsistent, or the trajectory has not been fitted 
        (remove PBC). This is the most common error. Probably the structure in the tpr file (pdb or gro) is 
        inconsistent. Many times it is because the system under study is longer than some edges of the box (See 
        this thread).
        
            **Solution:**
            
            ??? "Check the structure consistecy"
                
                Visualice the structure contained in the file that you are going to define in the `-cs`, `-ps`, or `-ls` 
                options and make sure it is consistent
                
                **Get the structure from tpr file**
                    
                    gmx editconf -f md.tpr -o md.pdb
    
            ??? "Make sure you have fitted the trajectory"
                
                Steps:
    
                1. Generate a group that contains both molecules
                    
                        gmx make_ndx -n index.ndx
                input: 1 | 12
                input: q
                
                2. remove the PBC
                    
                        gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center -n
                center: 20 (created group)
                output: 0
                
                3. remove the rotation and translation with respect to the reference structure (optional)
                    
                        gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_fit.xtc -n -fit rot + trans
                fit: 20 (created group)
                output: 0
                    
                4. Visualization
                    
                    Make sure that the trajectory is consistent
        
        2. You are trying to calculate the energetic contribution of a too-large group. Technically, the energy 
        value should not exceed 7 digits, so if you get a value higher than this, this error will occur. Although 
        `gmx_MMPBSA` can handle very large systems as described in example [Ribosomal50S_Mycalamide_A][1], it cannot 
        determine certain energetic terms. This limitation is not typical of `gmx_MMPBSA` but of Sander when writing 
        the output file.
           
            !!! note
            The error could be solved by recompiling Sander with some modifications in the output function. 
            However, this is not recommended since the error can be large. In any case, if you need to, we 
            could try to solve it.

    **I get positive energy value**
        

??? "NMODE calculation finish in error"
    The only error reported is probably related to RAM saturation. NMODE calculations require a considerable amount 
    of RAM depending on the number of atoms in your system. The amount of total RAM consumed during the calculation 
    will be: `RAM for 1 frame * number of threads`    



  [1]: ../examples/Ribosomal50S_Mycalamide_A/README.md