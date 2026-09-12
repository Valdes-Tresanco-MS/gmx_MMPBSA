---
template: main.html
title: Protein-protein with OPLS
---

# Protein-protein binding with OPLS

This example calculates the binding free energy of a protein-protein complex prepared with the OPLS force field. It
uses the single-trajectory approximation and the linear Poisson–Boltzmann model.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **Force field**

    OPLS

-   **Solvent model**

    Linear PB

-   **Frames**

    1–8

</div>

!!! info "Representative system"
    This protein-protein complex demonstrates how a supplied GROMACS OPLS topology is converted and analyzed. OPLS
    support is not limited to protein-protein systems; other receptor-ligand compositions can use the same
    topology-based workflow when compatible with `gmx_MMPBSA`. Method-specific restrictions and implicit-solvent
    validation still apply.

!!! caution "OPLS and implicit-solvent radii"
    The OPLS bonded and nonbonded parameters are read from the supplied GROMACS topology, while the PB calculation
    uses mbondi3 radii assigned during conversion to AMBER topology format. This combination is supported by
    `gmx_MMPBSA`, but it has not been extensively benchmarked as a force-field/implicit-solvent pairing. The example
    provides a reproducible starting point; validate the PB and radius choices for quantitative applications.

## Before you begin

The manual workflow uses the following files and selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **GROMACS system**

    Complex coordinates `com.pdb` (`-cs`) and OPLS topology `topol.top` (`-cp`). Keep the `toppar` directory
    containing the referenced `*.itp` files beside `topol.top`.

-   **Trajectory**

    PBC-corrected and fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index.ndx` (`-ci`) with the `Protein_chain1` and `Protein_chain2` groups (`-cg`)

</div>

A separate reference structure may be supplied with `-cr` when specific chain IDs or residue numbering are needed.
See the [complete command-line reference][1] for all options.

## Run the example

[Download the OPLS protein-protein example as a ZIP archive][6].

Extract the archive, change to the `OPLS/protein_protein` directory, and choose either the serial or MPI command. You
can also [view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.pdb \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg Protein_chain1 Protein_chain2 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.pdb \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg Protein_chain1 Protein_chain2 \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input pb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same linear-PB calculation.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for PB calculation
    # This input provides a practical starting point for OPLS protein-protein calculations.
    # Review the PB settings and radius model for your system.

    &general
    sys_name="OPLS_Support",
    startframe=1,
    endframe=8,
    PBRadii=4,
    /
    &pb
    radiopt=0, istrng=0.150,
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input pb"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "OPLS_Support"                         # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 8                                      # Last frame; e.g. 100
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
      memopt                         = 0                                      # Use membrane PB; 0/1
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
    This input provides a practical starting point and can serve as the basis for related calculations. Review the
    available [input-file options][2], their accepted values, and settings that depend on the system or intended
    analysis. In particular, assess the PB-radius treatment and sampling convergence before quantitative use.

## How this example works

The index groups define a 478-atom receptor (`Protein_chain1`) and a 130-atom protein ligand (`Protein_chain2`). The
ST approximation extracts both components from the same complex coordinates for every selected frame. Although the
source OPLS topology also contains solvent and ions, the selections restrict the calculation to the 608-atom protein
complex.

The bundled trajectory contains 11 frames; this input processes frames 1 through 8. With `npbopt=0`, the calculation
uses the linear PB equation and an ionic strength of 0.15 M.

## PB radii and OPLS parameters

Because `topol.top` is supplied, the concise input does not need a `forcefields` setting: the OPLS bonded, charge, and
Lennard-Jones parameters are read from the GROMACS topology. During conversion, `PBRadii=4` assigns mbondi3 radii to
the generated AMBER topologies. `radiopt=0` then instructs PBSA to use those stored topology radii rather than its
optimized PB radii.

This separation is why the radius choice deserves independent validation. Changing `PBRadii` changes the radii
assigned during topology conversion; changing `radiopt` controls whether PBSA uses those radii.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/PBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the result with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../docs/input_file.md#the-input-file
  [5]: ../../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/OPLS/protein_protein&fileName=gmx_MMPBSA-OPLS-Protein-Protein&rootDirectory=protein_protein
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/OPLS/protein_protein
