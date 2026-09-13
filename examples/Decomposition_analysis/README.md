---
template: main.html
title: Decomposition
---

# Per-residue decomposition analysis

This example decomposes the binding energy of a protein-protein complex into per-residue contributions. It uses the
same solvated simulation, frame range, and GB-Neck2 model as the alanine-scanning and stability examples.

<div class="example-card-grid" markdown>

-   **Analysis**

    Per-residue decomposition

-   **System**

    Protein-protein complex

-   **Solvent model**

    GB-Neck2 (`igb=8`)

-   **Bundled test**

    `gmx_MMPBSA_test -t 14`

</div>

## Before you begin

The manual workflow uses the following files and selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **GROMACS system**

    Structure `com.tpr` (`-cs`) and topology `topol.top` (`-cp`). Keep the `toppar` directory containing the
    referenced `*.itp` files beside `topol.top`.

-   **Trajectory**

    PBC-corrected and fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index.ndx` (`-ci`) with the `SOLU_chain1` and `SOLU_chain2` groups (`-cg`)

</div>

A complex reference structure may also be supplied with `-cr` when specific chain IDs or residue numbering must be
preserved in the decomposition selection. See the [complete command-line reference][1] for all options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 14
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the Decomposition example as a ZIP archive][6].

Extract the archive, change to the `Decomposition_analysis` directory, and choose either the serial or MPI command.
You can also
[view the example files on GitHub][10] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg SOLU_chain1 SOLU_chain2 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv \
      -do FINAL_DECOMP_MMPBSA.dat \
      -deo FINAL_DECOMP_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg SOLU_chain1 SOLU_chain2 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv \
      -do FINAL_DECOMP_MMPBSA.dat \
      -deo FINAL_DECOMP_MMPBSA.csv
    ```

## Configure the calculation

The example uses the minimal `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input gb decomp` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same MM/GBSA decomposition calculation.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file with decomposition analysis
    # This input provides a practical starting point for decomposition analysis.
    # Review the selected residues and decomposition output settings for your system.

    &general
    sys_name="Decomposition",
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    # Include at least one residue from both the receptor and ligand in print_res.
    # The "within" selection satisfies this requirement automatically.
    &decomp
    idecomp=2, dec_verbose=1,
    print_res="within 4",
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb decomp"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "Decomposition"                       # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 10                                     # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1


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
      explicit_waters_mask           = "dASA"                                     # Water reference; e.g. ":1-10", "within 4", "dASA"
      explicit_waters_group          = "automatic"                                     # Solvent group; e.g. "TIP3" or "automatic"
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

    # Decomposition namelist variables
    &decomposition
      idecomp                        = 2                                      # Decomp mode; 0-4
      dec_verbose                    = 1                                      # Decomp verbosity; 0-3
      print_res                      = "within 4"                             # Residues to print; e.g. "all", "within 6", "A/2-10"
      csv_format                     = 1                                      # Write CSV output; 0/1
    /

    ```

!!! info "Keep in mind"
    This input provides a practical starting point and can serve as the basis for production calculations. Review the
    available [input-file options][2], their accepted values, and adjust settings that depend on your system or protocol.
    Additional sample inputs are available [here][3].

## How this example works

The single-trajectory calculation selects the 608-atom protein-protein solute from the solvated source trajectory.
It processes frames 1 through 10 with GB-Neck2 (`igb=8`), the matching mbondi3 radii (`PBRadii=4`), and a salt
concentration of 0.15 M.

`idecomp=2` reports per-residue contributions with 1-4 electrostatic terms included in EEL and 1-4 van der Waals
terms included in VDW. `dec_verbose=1` reports the delta total, side-chain, and backbone contributions. The
`print_res="within 4"` selection includes residues from both components that lie within 4 Å of the interface.

!!! warning "Residue selection"
    A manual `print_res` selection must include at least one receptor residue and one ligand residue. Distance-based
    `within` selections enforce this requirement automatically. Large selections increase parsing time and output
    size, particularly for pairwise decomposition (`idecomp=3` or `4`).

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/GBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms (or the filename supplied with `-eo`).
- `FINAL_DECOMP_MMPBSA.dat`: the per-residue decomposition summary requested with `-do`.
- `FINAL_DECOMP_MMPBSA.csv`: the per-residue decomposition data (or the filename supplied with `-deo`).

## Visualize residue contributions

Decomposition values can be mapped onto the generated structures in PyMOL, VMD, or Chimera. In PyMOL,
`set cartoon_side_chain_helper, 1` can hide backbone atoms or display terminal residues incompletely; use
`set cartoon_side_chain_helper, 0` if that occurs. See the [VMD][8] and [Chimera][9] demonstrations for alternatives.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Decomposition_analysis&fileName=gmx_MMPBSA-Decomposition&rootDirectory=Decomposition_analysis
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://www.youtube.com/watch?v=PeboM8KE5SA
  [9]: https://www.youtube.com/watch?v=jKA4fuYuKps
  [10]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Decomposition_analysis
