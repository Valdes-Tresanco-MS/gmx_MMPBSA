---
template: main.html
title: Nonlinear PB
---

# Binding free energy calculation with the nonlinear PB equation

This example calculates the binding free energy of a protein-protein complex with the single-trajectory protocol
and the nonlinear Poisson-Boltzmann equation (NLPBE). It processes ten frames at an ionic strength of 0.15 M.
Nonlinear PB calculations have been available in `gmx_MMPBSA` since version 1.5.0.

<div class="example-card-grid" markdown>

-   **Method**

    Nonlinear PB (NLPBE)

-   **System**

    Protein-protein complex

-   **Protocol**

    Single trajectory

-   **Bundled test**

    `gmx_MMPBSA_test -t 21`

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
gmx_MMPBSA_test -t 21
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the nonlinear PB example as a ZIP archive][6].

Extract the archive, change to the `NonLinear_PB_solver` directory, and choose either the serial or MPI command.
You can also [view the example files on GitHub][9] before downloading them.

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
`gmx_MMPBSA --create_input pb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same nonlinear PB calculation.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for nonlinear PB calculation
    # This sample input is intended only to demonstrate that gmx_MMPBSA works.
    # Although it follows the recommendations in the Amber manual, some parameters
    # have been adjusted to keep the computational cost reasonable. Modify them as
    # appropriate for your system.

    &general
    sys_name="NonLinear_PB",
    startframe=1,
    endframe=10,
    /
    &pb
    npbopt=1,
    indi=1.0, istrng=0.15,
    radiopt=0,
    eneopt=1, cutnb=8.0,
    /
    # Check these threads
    # http://archive.ambermd.org/201203/0191.html
    # http://archive.ambermd.org/201610/0114.html
    # for more information about NLPB.
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input pb"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "NonLinear_PB"                       # System name; e.g. "complex"
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

    # (AMBER) Poisson-Boltzmann namelist variables
    &pb
      ipb                            = 2                                      # PB model; e.g. 2
      inp                            = 1                                      # Nonpolar method; 1 or 2
      indi                           = 1.0                                    # Internal dielectric; e.g. 1.0
      exdi                           = 78.5                                   # External dielectric; e.g. 78.5
      emem                           = 4.0                                    # Membrane dielectric; e.g. 4.0
      smoothopt                      = 1                                      # Dielectric smoothing; 0-2
      istrng                         = 0.15                                    # Ionic strength (M); e.g. 0.150
      radiopt                        = 0                                      # Use optimized radii; 0/1
      prbrad                         = 1.4                                    # Probe radius (A); e.g. 1.4
      iprob                          = 2.0                                    # Ion probe (A); e.g. 2.0
      sasopt                         = 0                                      # PB surface option; 0/1
      arcres                         = 0.25                                   # Arc resolution (A); e.g. 0.25
      memopt                         = 0                                      # Use membrane PB; 0/1
      mprob                          = 2.7                                    # Membrane probe (A); e.g. 2.7
      mthick                         = "automatic"                            # Membrane thickness (A), or automatic
      mctrdz                         = "automatic"                            # Membrane Z offset (A), or automatic
      membrane_atoms                 = "P"                                    # Atom names for automatic membrane parameters; semicolon-separated
      poretype                       = 1                                      # Pore type; 1 or 2
      npbopt                         = 1                                      # Use nonlinear PB; 0/1
      solvopt                        = 1                                      # PB solver; e.g. 1
      accept                         = 0.001                                  # Convergence; e.g. 0.001
      linit                          = 1000                                   # SCF iterations; e.g. 1000
      fillratio                      = 4.0                                    # Grid fill ratio; e.g. 4
      scale                          = 2.0                                    # Grid scale; e.g. 2
      nbuffer                        = 0.0                                    # Grid buffer; e.g. 0
      nfocus                         = 2                                      # Focus levels; e.g. 2
      fscale                         = 8                                      # Focus scale; e.g. 8
      npbgrid                        = 1                                      # Grid update freq.; e.g. 1
      bcopt                          = 5                                      # Boundary condition; e.g. 5
      eneopt                         = 1                                      # Energy option; e.g. 2
      frcopt                         = 0                                      # Force output; e.g. 0
      scalec                         = 0                                      # Reaction field option; e.g. 0
      cutfd                          = 5.0                                    # FD cutoff (A); e.g. 5
      cutnb                          = 8.0                                    # Nonbonded cutoff (A); e.g. 0
      nsnba                          = 1                                      # Pairlist frequency; e.g. 1
      decompopt                      = 2                                      # Decomp scheme; 1 or 2
      use_rmin                       = 1                                      # Use Rmin radii; 0/1
      sprob                          = 1.4                                    # SASA probe (A); e.g. 1.4
      vprob                          = 1.3                                    # Volume probe (A); e.g. 1.3
      rhow_effect                    = 1.129                                  # Water density; e.g. 1.129
      use_sav                        = 1                                      # Use SAV cavity; 0/1
      cavity_surften                 = 0.005                                  # Cavity surften; e.g. 0.005
      cavity_offset                  = 0.0                                    # Cavity offset; e.g. 0.0
      maxsph                         = 400                                    # Max surface dots; e.g. 400
      maxarcdot                      = 1500                                   # Max arc dots; e.g. 1500
      npbverb                        = 0                                      # PB verbosity; 0/1
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

The input processes ten frames with the nonlinear PB solver (`npbopt=1`), an internal dielectric constant of 1.0,
and an ionic strength of 0.15 M. It uses topology radii (`radiopt=0`) and a nonbonded cutoff of 8.0 Å.

## Interpreting nonlinear PB energies

!!! warning
    With `eneopt=1`, total electrostatic energies and forces are computed with the particle-particle
    particle-mesh (P3M) procedure described by [Lu and Luo][8]. The output therefore reports `EPB` as zero and
    combines the reaction-field and Coulombic energies in `EEL`. The van der Waals energy is evaluated together
    with the particle-particle contribution to the Coulombic energy.

    This setting requires a nonzero `cutnb` (`8.0` Å here) and `bcopt=5`, which is the default boundary
    condition used by this input.

    Because `EPB` and `EEL` are combined in the gas-phase term, ΔGGAS and ΔGSOLV are not separately meaningful
    for this calculation. ΔTOTAL remains the relevant combined result because it includes both contributions.

The comments at the end of `mmpbsa.in` link to additional Amber mailing-list discussions of
[nonlinear PB settings from 2012](http://archive.ambermd.org/201203/0191.html) and
[2016](http://archive.ambermd.org/201610/0114.html).

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the plain-text energy summary and statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][4] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [4]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/NonLinear_PB_solver&fileName=gmx_MMPBSA-Nonlinear-PB&rootDirectory=NonLinear_PB_solver
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://doi.org/10.1063/1.1622376
  [9]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/NonLinear_PB_solver
