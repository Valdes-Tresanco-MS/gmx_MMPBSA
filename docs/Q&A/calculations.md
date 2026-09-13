---
template: main.html
title: Q&A - Calculations
---

# `gmx_MMPBSA` Calculations
This page describes common calculation problems and possible solutions.

!!! note 
    Most of the errors noted here are the result of inconsistent input files. Please read the documentation and make 
    sure your files are consistent.



???+ example "ValueError: could not convert string to float: '*************'"
    : This error has two possible causes:
        
        1. The structure supplied with `-cs`, `-rs`, or `-ls` is inconsistent, or the trajectory has not been fitted
        and processed to remove PBC artifacts correctly. This is the most common cause and often occurs when the system
        is longer than one or more edges of the simulation box.
        
            #### **Possible solutions:**
            
            ??? tip "Check for structure consistency"
                
                Visualize the structure contained in the structure input file given in the `-cs`, `-rs`, or `-ls` 
                options and make sure it is intact and centered (Figure 1, right). A "broken" structure (Figure 1,
                left) can produce inconsistent results.
                
                **Generate the structure from a TPR file**
                    
                    gmx editconf -f md.tpr -o md.pdb

                <figure markdown="1">
                [![overview][1]][1]
                  <figcaption markdown="1" style="margin-top:0;">
                **Figure 1.** Visualization of two input structures. Left: "broken" structure; right: centered structure
                  </figcaption>
                </figure>
                
                  [1]:../assets/images/q_a/inconsistent_str.png
    
            ??? tip "Make sure you have fitted the trajectory"

                Visualize the trajectory supplied with `-ct`, `-rt`, or `-lt` and make sure PBC artifacts have been
                removed (Figure 2, right). An unfitted or broken trajectory (Figure 2, left) can produce inconsistent
                results.
                
                Steps:
    
                1. Generate a group that contains both molecules
                    
                        gmx make_ndx -n index.ndx
                
                        >1 | 12
                
                        >q

                    _Assuming 1 is the receptor and 12 is the ligand. This creates a new group (number 20 in this example)_
                
                2. Remove PBC artifacts
                    
                        gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center -n -ur compact
                        center: 20 (created group)
                        output: 0
                
                3. Remove rotation and translation relative to the reference structure (optional)
                    
                        gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_fit.xtc -n -fit rot+trans
                        fit: 20 (created group)
                        output: 0
                    
                4. Inspect the processed trajectory
                    
                    Make sure that the trajectory is intact and centered (Figure 2, right).

                5. If the process is unsuccessful, consider another option such as `-pbc nojump` (as suggested [here][4]).

                <figure markdown="1">
                [![overview][2]][2]
                  <figcaption markdown="1" style="margin-top:0;">
                **Figure 2.** Visualization of two input trajectories. Left: trajectory with PBC artifacts;
                right: centered and fitted trajectory with PBC artifacts removed.
                  </figcaption>
                </figure>
                
                  [2]:../assets/images/q_a/traj_comp.gif
                  [4]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/353
        
        2. You are trying to calculate the energetic contribution of a very large group. Technically, the energy 
        value should not exceed 7 digits, so if you get a value higher than this, this error will occur. Although 
        `gmx_MMPBSA` can handle very large systems, it cannot 
        determine certain energetic terms. This is a `sander` limitation when writing the output file.
           
            #### **Possible solutions:**
    
            : The error could be solved by recompiling `sander` with some modifications in the output function. 
            However, this is not recommended since the error can be large. Another possible solution could be modifying 
            the parameters of the calculation (solvent model, internal dielectric constant) or just performing the 
            calculation for a part of the system (sub-system).

    
???+ example "I get high values for the solvation energy when using PB model"    

    : When using the PB model, `inp=1` is the default. The total non-polar solvation free energy is modeled as a
    single term linearly proportional to the solvent-accessible surface area. To use the two-term cavity plus
    dispersion model, set `inp=2` explicitly. The dispersion term is computed with a surface-based integration
    method closely related to the PCM solvent for quantum chemical programs.

        #### **Possible solutions:**
    
        :  You may want to try inp=1 and avoid the EDISPER contribution. This way, the total non-polar solvation 
        free energy will be modeled as a single term linearly proportional to the solvent-accessible surface area. Just 
        add `inp=1` in the `&pb` namelist variables in the input file. See example below:

            ```
            &general
            startframe=5, endframe=100, interval=5, verbose=2, 
            /
            &pb
            istrng=0.15, fillratio=4.0, inp=1
            /
            ```

        : A legacy post-processing workaround can remove the stored EDISPER column from a rewritten report, but it is
        not an alternative PB calculation. Work on a copy of the complete result bundle, preserve the original
        `_GMXMMPBSA_info` and output files, change the copied value of `INPUT['pb']['inp']` to 1, and run:

            ```
            gmx_MMPBSA --rewrite-output
            ```

        :  `--rewrite-output` reparses the energies already stored in the copied result; it does not rerun PB or
        recompute the alternate `inp=1` non-polar model. The rewritten report therefore reports the existing ENPOLAR
        term while omitting EDISPER from the displayed model metadata. Do not present it as a recalculated `inp=1`
        result, and do not edit the original result in place.

        
        !!! info
            **The deliberate `inp=1` rerun and this post-processing workaround are different operations and can yield
            different values for the non-polar component of the solvation energy.
            (see [here](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/273#issuecomment-1207144247)).** 
            Use one or another depending on your interest.

    : Check this [publication][11] and see the drawbacks of modeling the total non-polar solvation free energy with 
    two terms, _i.e._, the cavity term and the dispersion term. Sometimes there are imbalances in the 
    cancellation of error between the two components and this can produce unrealistic non-polar energy values.

???+ example "The NMODE calculation ends with an error"
    This error is often caused by insufficient RAM. NMODE calculations can require a considerable amount of memory,
    depending on the number of atoms in the system. Estimate the total memory requirement as
    `RAM for one frame × number of threads`.



  [10]: ../input_file.md#pb-namelist-variables
  [11]: https://pubs.acs.org/doi/full/10.1021/jp073399n
