---
template: main.html
title: nmode Entropy
---

# Normal-mode entropy calculations

This example estimates the entropic contribution to protein-protein binding with normal-mode analysis (NMODE). It
combines a ten-frame MM/GBSA calculation using GB-Neck2 at 303.15 K with a frame-matched ten-frame NMODE
demonstration.

<div class="example-card-grid" markdown>

-   **Entropy method**

    Normal-mode analysis

-   **MM/GBSA model**

    GB-Neck2 (`igb=8`)

-   **NMODE frames**

    Frames 1-10

-   **Bundled test**

    `gmx_MMPBSA_test -t 17`

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
gmx_MMPBSA_test -t 17
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the NMODE example as a ZIP archive][6].

Extract the archive, change to the `nmode` directory, and choose either the serial or MPI command. You can also
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
`gmx_MMPBSA --create_input gb nmode` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same MM/GBSA and NMODE calculation.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for entropy calculations (nmode)
    # This sample input is intended only to demonstrate that gmx_MMPBSA works.
    # Although it follows the recommendations in the Amber manual, some parameters
    # have been adjusted to keep the computational cost reasonable. Modify them as
    # appropriate for your system.

    &general
    sys_name="NMode",
    startframe=1,
    endframe=10,
    PBRadii=4,
    temperature=303.15,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    # nmode processes the same ten frames selected in &general. The intentionally
    # relaxed maxcyc and drms settings reduce the runtime and are not appropriate
    # for a production calculation.
    &nmode
    nmstartframe=1,
    nmendframe=10,
    nminterval=1,
    nmode_igb=1, nmode_istrng=0.150,
    maxcyc=5, drms=2,
    /
    ```

=== "Generated input - all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb nmode"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "NMode"                       # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 10                                     # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1


      PBRadii                        = 4                                      # PB radii set; 1-7
      temperature                    = 303.15                                 # Temperature (K); e.g. 298.15
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

    # (AMBER) Generalized-Born namelist variables
    &gb
      igb                            = 8                                      # GB model, e.g. 2 or 8
      intdiel                        = 1.0                                    # Internal dielectric; e.g. 1.0
      extdiel                        = 78.5                                   # External dielectric; e.g. 78.5
      saltcon                        = 0.150                                  # Salt conc. (M); e.g. 0.150
      surften                        = 0.0072                                 # Surface tension; e.g. 0.0072
      surfoff                        = 0.0                                    # Surface offset; e.g. 0.0
      molsurf                        = 0                                      # Use molsurf; 0/1
      msoffset                       = 0.0                                    # Molsurf offset; e.g. 0.0
      probe                          = 1.4                                    # Probe radius (A); e.g. 1.4
      ifqnt                          = 0                                      # Enable QM/MM; 0/1
      qm_theory                      = "PM6-DH+"                              # QM theory; e.g. "PM6-DH+"
      qm_residues                    = ""                                     # QM residues; e.g. ":1-5"
      com_qmmask                     = ""                                     # Complex QM mask; e.g. ":1-5"
      rec_qmmask                     = ""                                     # Receptor QM mask; e.g. ":1-5"
      lig_qmmask                     = ""                                     # Ligand QM mask; e.g. ":1"
      qmcharge_com                   = 0                                      # Complex QM charge; e.g. 0
      qmcharge_lig                   = 0                                      # Ligand QM charge; e.g. 0
      qmcharge_rec                   = 0                                      # Receptor QM charge; e.g. 0
      qmcut                          = 9999.0                                 # QM cutoff (A); e.g. 9999
      scfconv                        = 1e-08                                  # SCF convergence; e.g. 1.0e-8
      itrmax                         = 1000                                   # Maximum SCF iterations; e.g. 5000
      # ndiis_attempts                 = None                                 # Maximum DIIS attempts per SCF cycle; e.g. 700
      peptide_corr                   = 0                                      # Peptide correction; 0/1
      writepdb                       = 1                                      # Write QM PDB; 0/1
      verbosity                      = 0                                      # QM/MM verbosity; 0-5
      alpb                           = 0                                      # Use ALPB; 0/1
      arad_method                    = 1                                      # ALPB size method; e.g. 1
    /

    # Normal Modes Entropy namelist variables
    &nmode
      nmstartframe                   = 1                                      # First NM frame; e.g. 1
      nmendframe                     = 10                                # Last NM frame; e.g. 100
      nminterval                     = 1                                      # NM frame stride; e.g. 1
      nmode_igb                      = 1                                      # GB model, e.g. 1
      nmode_istrng                   = 0.150                                    # NM ionic strength (M); e.g. 0.0
      dielc                          = 1.0                                    # NM dielectric; e.g. 1.0
      drms                           = 2                                  # Min. gradient cutoff; e.g. 0.001
      maxcyc                         = 5                                  # Max minimization cycles; e.g. 10000
    /

    ```

!!! info "Keep in mind"
    This input provides a practical starting point and can serve as the basis for production calculations. Review the
    available [input-file options][2], their accepted values, and adjust settings that depend on your system or protocol.
    Additional sample inputs are available [here][3].

## How this example works

The single-trajectory approximation generates the receptor and ligand Amber-format topologies and trajectories from
the complex. In this protein-protein system, the second protein is treated as the ligand. The command selects index
groups `3` and `4` as the receptor and ligand, respectively.

The MM/GBSA calculation processes ten frames with GB-Neck2 (`igb=8`), the matching mbondi3 radii (`PBRadii=4`),
a salt concentration of 0.15 M, and a temperature of 303.15 K. These settings now match the IE and C2 examples.
NMODE also processes frames 1 through 10, giving a frame-matched ten-frame comparison.

The NMODE engine has a separate implicit-solvent control and supports only vacuum (`nmode_igb=0`) or its built-in
HCT treatment (`nmode_igb=1`), not GB-Neck2. This example therefore retains `nmode_igb=1` while matching the ionic
strength at `nmode_istrng=0.150` M.

## Interpreting this demonstration

!!! warning
    Normal-mode calculations are substantially more computationally expensive than IE or C2. This example sets
    `maxcyc=5` and `drms=2` solely to keep the bundled test short. These values provide much less minimization than
    the current defaults (`maxcyc=10000`, `drms=0.001`) and must not be treated as production-quality convergence
    settings.

For production work, use convergence criteria appropriate for the system and verify that the selected NMODE frames
minimize successfully.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/GBSA summary and NMODE entropy results.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy and entropy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][4] for usage details.

  [1]: ../../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../input_file.md#the-input-file
  [3]: ../../../input_file.md#sample-input-files
  [4]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Entropy_calculations/nmode&fileName=gmx_MMPBSA-NMODE&rootDirectory=nmode
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Entropy_calculations/nmode
