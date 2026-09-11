---
template: main.html
title: Protein-ligand with CHARMM
---

# Protein-ligand binding with a CHARMM force field

This example calculates the binding free energy of the `JZ4` ligand to a protein using a system prepared with a
CHARMM force field and the single-trajectory approximation.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **Force field**

    CHARMM

-   **Solvent model**

    Linear PB with CHARMM radii

-   **Bundled test**

    `gmx_MMPBSA_test -t 10`

</div>

!!! info "Representative system"
    This protein-ligand complex demonstrates topology conversion for a GROMACS system prepared with CHARMM. CHARMM
    support is not limited to protein-ligand systems; the same topology-based workflow applies to other
    receptor-ligand compositions supported by `gmx_MMPBSA`. Method-specific restrictions and CHARMM conversion
    limitations still apply.

!!! warning "CHARMM CMAP conversion"
    The current GROMACS-to-AMBER topology conversion omits CHARMM CMAP terms and reports this during setup. This
    example exercises the complete workflow, but quantitative CHARMM applications should assess the effect of the
    missing CMAP contribution before interpreting binding energies.

!!! note "CHARMM PB radii"
    `PBRadii=7` selects the `charmm_radii` set, which is intended only for systems prepared with CHARMM force fields.
    Its protein radii draw on work by [Nina, Belogv, and Roux][9], nucleic-acid radii on [Banavali and Roux][10],
    and additional elements on [Fortuna and Costa][11]. With `radiopt=0`, PBSA reads these radii from the generated
    AMBER topologies.

## Before you begin

The manual workflow uses the following files and selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **GROMACS system**

    Structure `com.tpr` (`-cs`) and topology `topol.top` (`-cp`). Keep the `toppar` directory containing the
    referenced CHARMM `*.itp` files beside `topol.top`.

-   **Trajectory**

    PBC-corrected and fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index.ndx` (`-ci`) with receptor `Protein` and ligand `JZ4` (`-cg`)

</div>

The complete solvated system contains 46,407 atoms. The selected binding system contains the 2,614-atom protein and
22-atom ligand. A complex reference structure without hydrogens may also be supplied with `-cr` when specific chain
IDs or residue numbering are required. See the [complete command-line reference][1] for all options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 10
```

This is a slow test because it performs PB calculations. See the [`gmx_MMPBSA_test` documentation][7] for download,
selection, and cleanup options.

### Run it manually

[Download the protein-ligand CHARMM example as a ZIP archive][6].

Extract the archive, change to the `Protein_ligand_CHARMMff` directory, and choose either the serial or MPI command.
You can also [view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg Protein JZ4 \
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
      -cg Protein JZ4 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input pb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks describe the same
linear PB calculation.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for PB calculation
    # This input provides a practical starting point for protein-ligand calculations
    # with CHARMM topologies. Review the PB settings for your system and protocol.

    &general
    sys_name="Prot-Lig-CHARMM",
    startframe=1,
    endframe=4,
    PBRadii=7,
    /
    &pb
    radiopt=0, istrng=0.150, fillratio=4.0,
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input pb"
    Input block generated for the 1.7.0 release; the generator's development-version header is omitted from this documentation.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "Prot-Lig-CHARMM"                    # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 4                                      # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "leaprc.protein.ff14SB"               # Force fields; e.g. "leaprc.protein.ff14SB"
      ions_parameters                = 1                                      # Ion params; e.g. 1
      PBRadii                        = 7                                      # PB radii set; 1-7
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

    # (AMBER) Poisson-Boltzmann namelist variables
    &pb
      ipb                            = 2                                      # PB model; e.g. 2
      inp                            = 1                                      # Nonpolar method; 1 or 2
      indi                           = 1.0                                    # Internal dielectric; e.g. 1.0
      exdi                           = 78.5                                   # External dielectric; e.g. 78.5
      emem                           = 4.0                                    # Membrane dielectric; e.g. 4.0
      smoothopt                      = 1                                      # Dielectric smoothing; 0-2
      istrng                         = 0.150                                  # Ionic strength (M); e.g. 0.150
      radiopt                        = 0                                      # Use optimized radii; 0/1
      prbrad                         = 1.4                                    # Probe radius (A); e.g. 1.4
      iprob                          = 2.0                                    # Ion probe (A); e.g. 2.0
      sasopt                         = 0                                      # PB surface option; 0/1
      arcres                         = 0.25                                   # Arc resolution (A); e.g. 0.25
      memopt                         = 0                                      # Use membrane PB; 0-3
      mprob                          = 2.7                                    # Membrane probe (A); e.g. 2.7
      mthick                         = "automatic"                            # Membrane thickness (A), or automatic
      mctrdz                         = "automatic"                            # Membrane Z offset (A), or automatic
      membrane_atoms                 = "P"                                    # Atom names for automatic membrane parameters; semicolon-separated
      poretype                       = 1                                      # Pore type; 1 or 2
      npbopt                         = 0                                      # Use nonlinear PB; 0/1
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
      eneopt                         = 2                                      # Energy option; e.g. 2
      frcopt                         = 0                                      # Force output; e.g. 0
      scalec                         = 0                                      # Reaction field option; e.g. 0
      cutfd                          = 5.0                                    # FD cutoff (A); e.g. 5
      cutnb                          = 0.0                                    # Nonbonded cutoff (A); e.g. 0
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
    This input provides a practical starting point for CHARMM protein-ligand PB calculations. Review sampling, PB
    radii, dielectric treatment, grid convergence, and the CMAP limitation for the intended system. Additional
    [input-file options][2] may be needed for a production protocol.

## How this example works

The single-trajectory approximation extracts `Protein` and `JZ4` from the same four complex frames. Solvent and ions
remain present in the source TPR and trajectory but are not retained in the final complex, receptor, or ligand
calculation topologies.

`PBRadii=7` assigns CHARMM-specific radii during topology conversion, while `radiopt=0` instructs PBSA to use those
topology radii. Because `topol.top` is provided, no `forcefields` variable is needed in the concise input; the CHARMM
bonded and nonbonded parameters are read from the topology include tree. CMAP terms are the stated exception.

The calculation processes frames 1 through 4 with the linear PB equation (`npbopt=0`), an ionic strength of 0.15 M,
and a grid fill ratio of 4.0.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/PBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand_CHARMMff&fileName=gmx_MMPBSA-Protein-Ligand-CHARMM&rootDirectory=Protein_ligand_CHARMMff
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand_CHARMMff
  [9]: https://doi.org/10.1021/jp970736r
  [10]: https://doi.org/10.1021/jp025852v
  [11]: https://doi.org/10.1021/acs.jcim.1c00177
