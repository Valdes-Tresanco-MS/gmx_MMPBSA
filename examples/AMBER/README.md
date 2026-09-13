---
template: main.html
title: Native AMBER input files
---

# Binding free energy from native AMBER files

This example calculates the binding free energy of the RAS–RAF protein-protein complex directly from AMBER
topologies and a trajectory. It uses `amber_MMPBSA`, so no GROMACS structure, index, trajectory, or topology is
required.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **Input format**

    Native AMBER files

-   **Solvent model**

    GB-HCT (`igb=1`)

-   **Bundled test**

    `gmx_MMPBSA_test -t 25`

</div>

!!! info "Representative system"
    This tutorial uses a protein-protein complex to demonstrate the native AMBER workflow. Native AMBER support is
    not limited to protein-protein systems: other receptor-ligand compositions can be analyzed when their AMBER
    topologies, trajectories, and component selections are compatible with `amber_MMPBSA`. Individual calculation
    methods may impose additional restrictions.

## Before you begin

The manual workflow uses the following files and residue selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **Complex**

    Topology `ras-raf_complex.prmtop` (`-cp`) and trajectory `prod_complex.mdcrd` (`-ct`)

-   **Component topologies**

    Receptor `ras.prmtop` (`-rp`) and ligand `raf.prmtop` (`-lp`)

-   **Component selections**

    Receptor mask `:1-166` and ligand mask `:167-242` (`-cm`)

</div>

The folder also contains `ras-raf_complex.inpcrd`, a coordinate snapshot matching the complex topology. It is useful
for inspection and provenance but is not required by the command below. See the complete
[`amber_MMPBSA` reference][1] for all supported inputs and options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 25
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the native AMBER example as a ZIP archive][6].

Extract the archive, change to the `AMBER` directory, and choose either the serial or MPI command. You can also
[view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    amber_MMPBSA -O \
      -i mmpbsa.in \
      -cp ras-raf_complex.prmtop \
      -ct prod_complex.mdcrd \
      -rp ras.prmtop \
      -lp raf.prmtop \
      -cm ":1-166" ":167-242" \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 amber_MMPBSA -O \
      -i mmpbsa.in \
      -cp ras-raf_complex.prmtop \
      -ct prod_complex.mdcrd \
      -rp ras.prmtop \
      -lp raf.prmtop \
      -cm ":1-166" ":167-242" \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`amber_MMPBSA --create_input gb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same native AMBER MM/GBSA calculation.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Input file for a short AMBER input files test
    # This input provides a practical starting point for native AMBER calculations.
    # Review the model settings and frame range for your system.

    &general
    sys_name="AMBER",
    startframe=1,
    endframe=5,
    /
    &gb
    igb=1, saltcon=0.100,
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "AMBER"                                # System name; e.g. "complex"
      startframe                     = 1                                      # First frame; e.g. 1
      endframe                       = 5                                      # Last frame; e.g. 100
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
      igb                            = 1                                      # GB model, e.g. 2 or 8
      intdiel                        = 1.0                                    # Internal dielectric; e.g. 1.0
      extdiel                        = 78.5                                   # External dielectric; e.g. 78.5
      saltcon                        = 0.100                                  # Salt conc. (M); e.g. 0.150
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
    This five-frame calculation provides a practical, runnable introduction to native AMBER inputs. For production
    work, verify sampling convergence and adjust the frame range, implicit-solvent model, salt concentration, and
    other system-dependent settings. See the available [input-file options][2] before extending the protocol.

## How this example works

The ST approximation evaluates the complex, receptor, and ligand using coordinates from the same five trajectory
frames. The complex masks assign residues 1–166 to RAS and residues 167–242 to RAF. The supplied `ras.prmtop` and
`raf.prmtop` provide the component parameters; if `-rp` or `-lp` is omitted, `amber_MMPBSA` instead extracts that
component topology from the complex.

The three bundled topologies store `mbondi` radii. This example therefore uses the conventionally matched GB-HCT
model (`igb=1`) with a salt concentration of 0.10 M. Native AMBER topologies already contain atomic parameters,
charges, radii, and screening values. The legacy topology-preparation settings do not rebuild or replace those data
during the normal native AMBER workflow.

## Native topology radii

`amber_MMPBSA` preserves the `RADII`, `SCREEN`, and `RADIUS_SET` data stored in each input topology. If the selected
`igb` does not conventionally match the stored radius set, the program warns but does not alter the topology. Choose
the desired radius set when creating the topology in `tleap`; see the [`amber_MMPBSA` reference][1] for conventional
`igb`/radius pairings and precedence details.

## Multiple-trajectory calculations

For the MT approximation, provide the unbound receptor topology, mask, and trajectory with `-rp`, `-rm`, and `-rt`,
and the corresponding ligand inputs with `-lp`, `-lm`, and `-lt`. The selected complex, receptor, and ligand
trajectories must contain the same number of frames after applying `startframe`, `endframe`, and `interval`.

When multiple files are supplied to one trajectory option, they are concatenated in command-line order and analyzed
as one pooled trajectory. They are not treated as independent replicas.

## Mask selections

Native AMBER masks in this workflow select whole residues. Non-contiguous residue ranges are supported; for example,
the following selection assigns residues 1–120 and 181–260 to the receptor and residues 121–180 to the ligand:

```bash
-cm ":1-120,181-260" ":121-180"
```

Atom-name and boolean mask expressions are not accepted for these component selections.

## Explicit receptor waters

For an explicit-water AMBER calculation, supply a matching solvated `-cp`/`-ct` pair, omit `-rp` and `-lp`, and set
`explicit_waters`, `explicit_waters_mask`, and `solvated_trajectory=1` in the input. The `-cm` masks must select only
the solute components. Selected waters are assigned to the receptor, while the ligand topology remains dry. See the
[explicit-water details in the command reference][9] for supported models and water-selection behavior.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/GBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the result with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/amber_MMPBSA.md
  [2]: ../../docs/input_file.md#the-input-file
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/AMBER&fileName=gmx_MMPBSA-Native-AMBER&rootDirectory=AMBER
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/AMBER
  [9]: ../../docs/amber_MMPBSA.md#input-files
