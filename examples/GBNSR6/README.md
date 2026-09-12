---
template: main.html
title: GBNSR6
---

# Binding free energy calculation with the GBNSR6 model

This example calculates the binding free energy of a protein-protein complex with the single-trajectory protocol
and the GBNSR6 implicit-solvent model. It processes ten frames and uses an ionic strength of 0.15 M.

<div class="example-card-grid" markdown>

-   **Method**

    GBNSR6

-   **System**

    Protein-protein complex

-   **Protocol**

    Single trajectory

-   **Bundled test**

    `gmx_MMPBSA_test -t 24`

</div>

## Before you begin

The manual workflow uses the following files and selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **GROMACS system**

    Structure `com.tpr` (`-cs`) and topology `topol.top` (`-cp`). Keep any `*.itp` files referenced by the topology
    in the same directory.

-   **Trajectory**

    PBC-corrected and fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index.ndx` (`-ci`) and receptor/ligand group names or zero-based group numbers (`-cg`)

</div>

A complex reference structure without hydrogens may also be supplied with `-cr`. It is optional but recommended
when you need specific chain IDs or residue numbering. See the [complete command-line reference][1] for all options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 24
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the GBNSR6 example as a ZIP archive][6].

Extract the archive, change to the `GBNSR6` directory, and choose either the serial or MPI command. You can also
[view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg 3 4 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg 3 4 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the minimal `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input gbnsr6` and then adapted with the example-specific values. The concise block is the
runnable starting point; the generated block exposes additional options and defaults, so the two blocks are not
textually identical.

For this example, `-cp topol.top` supplies the GROMACS topology parameters used by the calculation. The generated
`forcefields` line is retained as an all-options/default field and should not be read as evidence that both input
blocks use identical topology-preparation paths.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for GB calculation using GBNSR6 model
    # This sample input is intended only to demonstrate that gmx_MMPBSA works.
    # Although it follows the recommendations in the Amber manual, some parameters
    # have been adjusted to keep the computational cost reasonable. Modify them as
    # appropriate for your system.

    &general
    sys_name="GBNSR6",
    startframe=1,
    endframe=10,
    /

    &gbnsr6
    istrng=0.15
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gbnsr6"
    Input block generated for the 1.7.0 release with --create_input gbnsr6 and adapted for this example.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "GBNSR6"                       # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 10                                     # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "leaprc.protein.ff14SB"      # Force fields; e.g. "leaprc.protein.ff14SB"
      ions_parameters                = 1                                      # Ion params; e.g. 1
      PBRadii                        = 4                                      # PB radii set; 1-7
      temperature                    = 298.15                                 # Temperature (K); e.g. 298.15
      qh_entropy                     = 0                                      # Legacy QH output reader; new calculations reject 1
      interaction_entropy            = 0                                      # Run IE entropy; 0/1
      ie_segment                     = 25                                     # IE tail diagnostic only (%); not primary IE; e.g. 25
      c2_entropy                     = 0                                      # Run C2 entropy; 0/1
      assign_chainID                 = 0                                      # Assign chain IDs; 0/1
      exp_ki                         = 0.0                                    # Experimental Ki (nM); e.g. 0.0
      full_traj                      = 0                                      # Write full trajectory; 0/1
      gmx_path                       = ""                                     # GROMACS path; e.g. "/usr/bin"
      keep_files                     = 2                                      # Files to keep; 0-2
      netcdf                         = 0                                      # Use NetCDF; 0/1
      solvated_trajectory            = 1                                      # Clean solvated traj.; 0/1
      explicit_waters                = 0                                      # Explicit waters; e.g. 10
      explicit_waters_mask           = ""                                     # Water reference; e.g. ":1-10", "within 4", "dASA"
      explicit_waters_group          = ""                                     # Solvent group; e.g. "TIP3"
      explicit_waters_dasa_cutoff    = 0.5                                    # dASA cutoff; e.g. 0.5
      explicit_waters_as             = "receptor"                             # Water owner; e.g. "receptor"
      explicit_waters_extra_points   = "error"                                # Virtual sites; "error" or "strip"
      verbose                        = 1                                      # Output verbosity; 0-2
    /

    # GBNSR6 namelist variables
    &gbnsr6
      b                              = 0.028                                  # GBNSR6 offset; e.g. 0.028
      alpb                           = 1                                      # Use ALPB; 0/1
      epsin                          = 1.0                                    # Solute dielectric; e.g. 1.0
      epsout                         = 78.5                                   # Solvent dielectric; e.g. 78.5
      istrng                         = 0.15                                    # Ionic strength (M); e.g. 0.150
      rs                             = 0.52                                   # Boundary shift; e.g. 0.52
      dprob                          = 1.4                                    # Probe radius (A); e.g. 1.4
      space                          = 0.5                                    # Grid spacing (A); e.g. 0.5
      arcres                         = 0.2                                    # Arc resolution; e.g. 0.2
      radiopt                        = 0                                      # Radii option; e.g. 0
      chagb                          = 0                                      # Use CHAGB; 0/1
      roh                            = 1                                      # RzOH value; e.g. 1
      tau                            = 1.47                                   # CHAGB tau; e.g. 1.47
      cavity_surften                 = 0.005                                  # Cavity surften; e.g. 0.005
    /

    ```

!!! info "Keep in mind"
    This input provides a practical starting point and can serve as the basis for production calculations. Review the
    available [input-file options][2], their accepted values, and adjust settings that depend on your system or protocol.
    Additional sample inputs are available [here][3].

## How this example works

The single-trajectory approximation generates the receptor and ligand structures and trajectories from the complex.
In this protein-protein system, the second protein is treated as the ligand. The command selects index groups `3`
and `4` as the receptor and ligand, respectively.

The input processes ten frames with the GBNSR6 model and an ionic strength of 0.15 M.

!!! note "About the GBNSR6 model"
    - GBNSR6 computes effective Born radii numerically through R6 integration over the solute molecular surface
      ([reference][222]).
    - Unlike most practical GB models, GBNSR6 is parameter-free in the same sense as the numerical PB framework.
      Consequently, its accuracy relative to the PB standard is largely unaffected by the selected input atomic radii.
    - `gmx_MMPBSA` automatically prepares temporary GBNSR6 topology copies for the calculation while preserving the
      original complex, receptor, and ligand topologies for output parsing. No additional input option is required.
    - Chapter [5 of the Amber 2021 Reference Manual](https://ambermd.org/doc12/Amber21.pdf#chapter.5) provides a more
      detailed description of the model and its parameters.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the plain-text energy summary and statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/GBNSR6&fileName=gmx_MMPBSA-GBNSR6&rootDirectory=GBNSR6
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/GBNSR6
  [222]: https://pubs.acs.org/doi/abs/10.1021/ct200786m
