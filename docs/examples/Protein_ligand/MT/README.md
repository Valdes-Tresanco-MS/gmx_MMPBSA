---
template: main.html
title: Protein-ligand (MT)
---

# Protein-ligand binding free energy: multiple trajectories

This example demonstrates the files and command-line options used by the multiple-trajectory (MT) workflow for a
protein-small-molecule complex. The complex also serves as the source for any receptor or ligand component that is
not supplied separately.

<div class="example-card-grid" markdown>

-   **Protocol**

    Multiple trajectories

-   **System**

    Protein-small molecule

-   **Solvent model**

    GB-Neck2 (`igb=8`)

-   **Bundled test**

    `gmx_MMPBSA_test -t 16`

</div>

## Before you begin

The bundled command supplies a separate structure, trajectory, index, and topology for each state:

<div class="example-card-grid" markdown>

-   **Complex**

    `com.tpr`, `com_traj.xtc`, `index.ndx`, and `topol.top`

-   **Receptor**

    `rec.pdb`, `rec_traj.xtc`, `rec_index.ndx`, and `rec.top`

-   **Ligand**

    `lig.pdb`, `lig_traj.xtc`, `lig_index.ndx`, and `lig.top`

-   **Shared parameters**

    The `toppar` directory referenced by all three topology files

</div>

The complex index contains the `receptor` and `ligand` groups. The receptor and ligand indexes each use their
`System` group. See the [complete command-line reference][1] for all options.

!!! info "Supplying receptor and ligand components"
    The receptor and ligand inputs are independent of one another. If either component is not provided separately,
    `gmx_MMPBSA` extracts its structure and trajectory from the complex and generates its topology from the selected
    complex group. You can therefore supply both components, only the receptor, only the ligand, or neither one.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 16
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the protein-ligand MT example as a ZIP archive][6].

Extract the archive, change to the `MT` directory, and choose either the serial or MPI command. You can also
[view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg receptor ligand \
      -cp topol.top \
      -rs rec.pdb \
      -rt rec_traj.xtc \
      -ri rec_index.ndx \
      -rg System \
      -rp rec.top \
      -ls lig.pdb \
      -lt lig_traj.xtc \
      -li lig_index.ndx \
      -lg System \
      -lp lig.top \
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
      -cg receptor ligand \
      -cp topol.top \
      -rs rec.pdb \
      -rt rec_traj.xtc \
      -ri rec_index.ndx \
      -rg System \
      -rp rec.top \
      -ls lig.pdb \
      -lt lig_traj.xtc \
      -li lig_index.ndx \
      -lg System \
      -lp lig.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the minimal `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input gb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same calculation.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for GB calculation
    # This sample input is intended only to demonstrate that gmx_MMPBSA works.
    # Although it follows the recommendations in the Amber manual, some parameters
    # have been adjusted to keep the computational cost reasonable. Modify them as
    # appropriate for your system.

    &general
    sys_name="Prot-Lig-MT",
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    ```

=== "Generated input - all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "Prot-Lig-MT"                       # System name; e.g. "complex"
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

    ```

!!! info "Keep in mind"
    This input provides a practical starting point and can serve as the basis for production calculations. Review the
    available [input-file options][2], their accepted values, and adjust settings that depend on your system or protocol.
    Additional sample inputs are available [here][3].

## How this example works

The complex inputs (`-cs`, `-ct`, `-ci`, `-cg`, and `-cp`) define the bound system and the selections used for any
component that must be generated from it. The `-rs`/`-rt`/`-ri`/`-rg`/`-rp` options override that source for the
receptor, while `-ls`/`-lt`/`-li`/`-lg`/`-lp` do the same for the ligand. Omitting one component leaves the other
component's independently supplied inputs unchanged.

The calculation processes frames 1 through 10 with GB-Neck2 (`igb=8`), the matching mbondi3 radii (`PBRadii=4`),
and 0.15 M salt.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/GBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][4] for usage details.

  [1]: ../../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../input_file.md#the-input-file
  [3]: ../../../input_file.md#sample-input-files
  [4]: ../../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand/MT&fileName=gmx_MMPBSA-Protein-Ligand-MT&rootDirectory=MT
  [7]: ../../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand/MT
