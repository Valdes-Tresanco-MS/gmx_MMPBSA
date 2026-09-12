---
template: main.html
title: 3D-RISM
---

# Protein-protein binding free energy calculations with MM/3D-RISM

This example calculates the binding free energy of a protein-protein complex with the single-trajectory protocol
and the MM/3D-RISM method. It uses four frames and the Kovalenko-Hirata closure to provide a short, reproducible
demonstration rather than a production calculation.

<div class="example-card-grid" markdown>

-   **Method**

    MM/3D-RISM

-   **System**

    Protein-protein complex

-   **Protocol**

    Single trajectory

-   **Bundled test**

    `gmx_MMPBSA_test -t 18`

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
gmx_MMPBSA_test -t 18
```

See the [`gmx_MMPBSA_test` documentation][6] for download, selection, and cleanup options.

### Run it manually

[Download the 3D-RISM example as a ZIP archive][5].

Extract the archive, change to the `3D-RISM` directory, and choose either the serial or MPI command. You can also
[view the example files on GitHub][7] before downloading them.

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
`gmx_MMPBSA --create_input rism` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same calculation; the generated version also documents every available `&general` and `&rism` variable.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for MM/3D-RISM
    # This sample input is intended only to demonstrate that gmx_MMPBSA works. Although
    # it follows the recommendations in the Amber manual, some parameters have been adjusted
    # to keep the computational cost reasonable. Modify them as appropriate for your system.

    &general
    sys_name="3D-RISM",
    startframe=5,
    endframe=8,
    /
    &rism
    polardecomp=0, tolerance=0.001, rism_verbose=2, closure="kh"
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input rism"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "3D-RISM"                # System name; e.g. "complex"
      startframe                     = 5                                      # First frame; e.g. 1
      endframe                       = 8                                           # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "oldff/leaprc.ff99SB,leaprc.gaff"      # Force fields; e.g. "leaprc.protein.ff14SB"
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

    # 3D-RISM namelist variables
    &rism
      closure                        = "kh"                                   # Closure equation; e.g. "kh"
      gfcorrection                   = 0                                      # GF correction; 0/1
      pcpluscorrection               = 0                                      # PC+ correction; 0/1
      noasympcorr                    = 1                                      # Disable asymptotic corr.; 0/1
      buffer                         = 14.0                                   # Grid buffer (A); e.g. 14
      solvcut                        = -1.0                                   # Solvent cutoff (A); e.g. -1
      grdspc                         = 0.5,0.5,0.5                            # Grid spacing; e.g. 0.5,0.5,0.5
      ng                             = -1,-1,-1                               # Grid points; e.g. -1,-1,-1
      solvbox                        = -1,-1,-1                               # Solvent box; e.g. -1,-1,-1
      tolerance                      = 0.001                                  # Convergence tol.; e.g. 1.0e-5
      ljTolerance                    = -1.0                                   # LJ tolerance; e.g. -1.0
      asympKSpaceTolerance           = -1.0                                   # K-space tolerance; e.g. -1.0
      treeDCF                        = 1                                      # Use DCF treecode; 0/1
      treeTCF                        = 1                                      # Use TCF treecode; 0/1
      treeCoulomb                    = 0                                      # Use Coulomb treecode; 0/1
      treeDCFMAC                     = 0.1                                    # DCF MAC; e.g. 0.1
      treeTCFMAC                     = 0.1                                    # TCF MAC; e.g. 0.1
      treeCoulombMAC                 = 0.1                                    # Coulomb MAC; e.g. 0.1
      treeDCFOrder                   = 2                                      # DCF tree order; e.g. 2
      treeTCFOrder                   = 2                                      # TCF tree order; e.g. 2
      treeCoulombOrder               = 2                                      # Coulomb tree order; e.g. 2
      treeDCFN0                      = 500                                    # DCF leaf size; e.g. 500
      treeTCFN0                      = 500                                    # TCF leaf size; e.g. 500
      treeCoulombN0                  = 500                                    # Coulomb leaf size; e.g. 500
      mdiis_del                      = 0.7                                    # MDIIS step size; e.g. 0.7
      mdiis_nvec                     = 5                                      # MDIIS vectors; e.g. 5
      mdiis_restart                  = 10.0                                   # MDIIS restart; e.g. 10.0
      maxstep                        = 10000                                  # Max iterations; e.g. 10000
      npropagate                     = 5                                      # Propagation history; e.g. 5
      polardecomp                    = 0                                      # Polar decomposition; 0/1
      entropicdecomp                 = 0                                      # Entropic decomposition; 0/1
      rism_verbose                   = 2                                      # RISM verbosity; 0-2
    /

    ```

!!! info "Keep in mind"
    This input provides a practical starting point and can serve as the basis for production calculations. Review the
    available [input-file options][2], their accepted values, and adjust settings that depend on your system or protocol.
    Additional sample inputs are available [here][3].

## How this example works

The single-trajectory approximation generates the receptor and ligand Amber topologies and trajectories from the
complex. In this protein-protein system, the second protein is treated as the ligand. The command selects index
groups `3` and `4` as the receptor and ligand, respectively.

The input processes four frames with the Kovalenko-Hirata closure. Its convergence tolerance is `0.001`, increased
from the default of `0.00001` to keep the runtime practical for a test calculation.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the plain-text energy summary and statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Troubleshooting

!!! warning "AmberTools/Fortran runtime compatibility"
    Some conda AmberTools builds linked with newer Fortran runtime libraries can stop before the 3D-RISM calculation
    starts. This is an AmberTools/runtime compatibility problem, not an input-preparation error in `gmx_MMPBSA`.

??? example "Show the error and a tested workaround"

    The affected runtime can report:

    ```text
    Fortran runtime error: Missing comma between descriptors
    amber_rism_interface.F90
    ```

    One tested workaround uses `gmx_MMPBSA` 1.6.4 with Python 3.9, AmberTools 23, and compatible GCC runtime
    libraries:

    ```bash
    conda create -n gmxMMPBSA_rism -c conda-forge python=3.9 ambertools=23 "libgfortran5<13" "libgcc-ng<13"
    conda activate gmxMMPBSA_rism
    python -m pip install "gmx_MMPBSA==1.6.4"
    ```

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][4] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [4]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [5]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/3D-RISM&fileName=gmx_MMPBSA-3D-RISM&rootDirectory=3D-RISM
  [6]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [7]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/3D-RISM
