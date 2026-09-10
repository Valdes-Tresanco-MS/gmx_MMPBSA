---
template: main.html
title: QM/MMGBSA
---

# QM/MMGBSA binding free energy

This example combines a semiempirical QM/MM energy model with GB solvation in a single-trajectory binding free-energy
calculation. Residues close to the ligand form the QM region; the remainder of the system is treated classically.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **QM method**

    PM6-DH+

-   **QM region**

    Residues within 4 Å

-   **Bundled test**

    `gmx_MMPBSA_test -t 23`

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

    Index `index.ndx` (`-ci`) with the `receptor` and `ligand` groups (`-cg`)

</div>

A complex reference structure without hydrogens may also be supplied with `-cr`. It is optional but recommended
when you need specific chain IDs or residue numbering. See the [complete command-line reference][1] for all options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 23
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the QM/MMGBSA example as a ZIP archive][6].

Extract the archive, change to the `QM_MMGBSA` directory, and choose either the serial or MPI command. You can also
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
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input gb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks describe the
same calculation.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for QM/MMGBSA calculation
    # This input provides a practical starting point for QM/MMGBSA calculations.
    # Review the QM region, charge, Hamiltonian, and convergence settings for your system.

    &general
    sys_name="QM/MMGBSA",
    startframe=1,
    endframe=10,
    PBRadii=2,
    /
    &gb
    igb=1, saltcon=0.150,
    ifqnt=1, qm_theory=PM6-DH+,

    # Residues to be treated with QM can be selected using different approaches. Make sure to include at least
    # one residue from both the receptor and ligand in the qm_residues mask when using 'ifqnt'. This requirement is
    # automatically fulfilled when using the within keyword https://groups.google.com/g/gmx_mmpbsa/c/GNb4q4YGCH8

    # Residue selection by distance (recommended)
    qm_residues="within 4",

    ## Explicit residue selection
    #qm_residues="A/40-41,44,47,78,81-82,85,88,115,118,122,215,218-220,232 B/241"

    # Residue selection with amber masks
    #com_qmmask="(:44,47,85,88,218&!@N,H,CA,HA,C,O) | :241"
    #rec_qmmask="(:44,47,85,88,218&!@N,H,CA,HA,C,O)"
    #lig_qmmask=":1"
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb"
    Input block generated for the 1.7.0 release; the generator's development-version header is omitted from this documentation.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "QM/MMGBSA"                            # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 10                                     # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "oldff/leaprc.ff99SB,leaprc.gaff"      # Force fields; e.g. "leaprc.protein.ff14SB"
      ions_parameters                = 1                                      # Ion params; e.g. 1
      PBRadii                        = 2                                      # PB radii set; 1-7
      temperature                    = 298.15                                 # Temperature (K); e.g. 298.15
      qh_entropy                     = 0                                      # Legacy QH output reader; new calculations reject 1
      interaction_entropy            = 0                                      # Run IE entropy; 0/1
      ie_segment                     = 25                                     # IE segment length (%); e.g. 25
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
      igb                            = 1                                      # GB model, e.g. 2 or 8
      intdiel                        = 1.0                                    # Internal dielectric; e.g. 1.0
      extdiel                        = 78.5                                   # External dielectric; e.g. 78.5
      saltcon                        = 0.150                                  # Salt conc. (M); e.g. 0.150
      surften                        = 0.0072                                 # Surface tension; e.g. 0.0072
      surfoff                        = 0.0                                    # Surface offset; e.g. 0.0
      molsurf                        = 0                                      # Use molsurf; 0/1
      msoffset                       = 0.0                                    # Molsurf offset; e.g. 0.0
      probe                          = 1.4                                    # Probe radius (A); e.g. 1.4
      ifqnt                          = 1                                      # Enable QM/MM; 0/1
      qm_theory                      = "PM6-DH+"                              # QM theory; e.g. "PM6-DH+"
      qm_residues                    = "within 4"                             # QM residues; e.g. ":1-5"
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
    QM/MM results can be sensitive to the QM-region boundary, net charge, Hamiltonian, and SCF convergence. Inspect
    `qmmm_region.pdb`, confirm that the selected residues form a chemically sensible region, and validate the setup
    before applying it to production calculations.

## How this example works

The ST approximation extracts the 893-atom receptor and 30-atom ligand from each selected complex frame. With
`qm_residues="within 4"`, `gmx_MMPBSA` identifies residues from both components that lie within 4 Å of their interface.
For this bundle, the selection contains five receptor residues plus the ligand. The corresponding QM charges are
calculated from the topology and assigned automatically when explicit QM masks and user-defined charges are absent.

The calculation processes frames 1 through 10 using the GB-HCT model (`igb=1`), mbondi radii (`PBRadii=2`), a salt
concentration of 0.15 M, and PM6-DH+. The alternative selections retained in the concise input show how to specify
residues directly or provide separate Amber masks. Every QM region must contain atoms from both receptor and ligand.
If `qm_theory` is omitted, `gmx_MMPBSA` uses PM6-DH+ by default.

## References for PM6-DH+

PM6-DH+ adds dispersion and hydrogen-bond corrections that are important for many biomolecular noncovalent
interactions. Relevant method and application studies include:

1. Řezáč and Hobza, *J. Chem. Theory Comput.* **2009**, 5, 1749–1760. [doi:10.1021/ct9000922](https://doi.org/10.1021/ct9000922)
2. Korth, *J. Chem. Theory Comput.* **2010**, 6, 3808–3816. [doi:10.1021/ct100408b](https://doi.org/10.1021/ct100408b)
3. Grimme and Brandenburg, *Front. Chem.* **2015**, 3, 8. [PMC4881564](https://pmc.ncbi.nlm.nih.gov/articles/PMC4881564/)
4. Thapa *et al.*, *J. Phys. Chem. B* **2018**, 122, 7866–7878. [doi:10.1021/acs.jpcb.8b03655](https://doi.org/10.1021/acs.jpcb.8b03655)
5. *Commun. Biol.* **2025**. [doi:10.1038/s42003-025-09143-z](https://doi.org/10.1038/s42003-025-09143-z)
6. Muddana and Gilson, *J. Chem. Theory Comput.* **2012**, 8, 2868–2880. [doi:10.1021/ct3002738](https://doi.org/10.1021/ct3002738)

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the QM/MMGBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.
- `qmmm_region.pdb`: the selected QM region for visual inspection.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/QM_MMGBSA&fileName=gmx_MMPBSA-QM-MMGBSA&rootDirectory=QM_MMGBSA
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/QM_MMGBSA
