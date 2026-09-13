---
template: main.html
title: Correlation analysis
---

# Correlation analysis

This example compares calculated binding energies with experimental inhibition constants for a wild-type
protein-protein complex and six alanine mutants. All seven calculations use the same solvated simulation and
GB-Neck2 settings so that the resulting systems can be loaded together in `gmx_MMPBSA_ana`.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **Series**

    Wild type and six mutants

-   **Solvent model**

    GB-Neck2 (`igb=8`)

-   **Analysis**

    Recursive correlation

</div>

## Before you begin

The seven calculation directories share the following files and selections from the `Correlation` directory:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    One `mmpbsa.in` file in each of `WT`, `T13`, `H15`, `Y23`, `Y25`, `S32`, and `W34`

-   **GROMACS system**

    Structure `com.tpr` (`-cs`) and topology `topol.top` (`-cp`). Keep the `toppar` directory containing the
    referenced `*.itp` files beside `topol.top`.

-   **Trajectory**

    PBC-corrected and fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index.ndx` (`-ci`) with the `SOLU_chain1` and `SOLU_chain2` groups (`-cg`)

</div>

The molecular bundle is the same protein-protein system used by the protein-protein binding example. See the
[complete command-line reference][1] for all options.

## Run the calculations

[Download the Correlation example as a ZIP archive][6].

Extract the archive and change to the `Correlation` directory. You can also [view the example files on GitHub][8]
before downloading them. The explicit directory list ensures that only the seven systems are processed; the
`toppar` dependency directory is not treated as a calculation.

=== "Serial"

    ```bash
    for system in WT T13 H15 Y23 Y25 S32 W34; do
      (
        cd "$system"
        gmx_MMPBSA -O \
          -i mmpbsa.in \
          -cs ../com.tpr \
          -ct ../com_traj.xtc \
          -ci ../index.ndx \
          -cg SOLU_chain1 SOLU_chain2 \
          -cp ../topol.top \
          -o FINAL_RESULTS_MMPBSA.dat \
          -eo FINAL_RESULTS_MMPBSA.csv \
          -nogui
      )
    done
    ```

=== "With MPI"

    ```bash
    for system in WT T13 H15 Y23 Y25 S32 W34; do
      (
        cd "$system"
        mpirun -np 2 gmx_MMPBSA -O \
          -i mmpbsa.in \
          -cs ../com.tpr \
          -ct ../com_traj.xtc \
          -ci ../index.ndx \
          -cg SOLU_chain1 SOLU_chain2 \
          -cp ../topol.top \
          -o FINAL_RESULTS_MMPBSA.dat \
          -eo FINAL_RESULTS_MMPBSA.csv \
          -nogui
      )
    done
    ```

## Configure the calculations

Each directory contains a concise input. The wild-type and H15 inputs are shown below to illustrate the two
patterns. The all-options version was generated with `gmx_MMPBSA --create_input gb ala` and then updated with the
same H15-specific values. It therefore describes the same mutant calculation as the concise H15 input.

=== "Wild-type input"

    ```yaml linenums="1" title="WT/mmpbsa.in"
    Sample input file for correlation analysis
    # This input provides a practical starting point for correlation calculations.
    # Review the model settings and experimental affinity for your system.

    &general
    sys_name="WT",
    exp_ki=18,
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    ```

=== "Mutant input"

    ```yaml linenums="1" title="H15/mmpbsa.in"
    Sample input file for correlation analysis
    # This input provides a practical starting point for correlation calculations.
    # Review the mutation, model settings, and experimental affinity for your system.

    &general
    sys_name="H15",
    exp_ki=15,
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    &alanine_scanning
    mutant="ALA", mutant_res="A:15", mutant_only=1,
    /
    ```

=== "Generated input - all options"

    ```yaml linenums="1" title="H15/mmpbsa.in generated with --create_input gb ala"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "H15"                                  # System name; e.g. "complex"
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
      exp_ki                         = 15                                     # Experimental Ki (nM); e.g. 0.0
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

    # Alanine scanning namelist variables
    &alanine_scanning
      mutant_res                     = "A:15"                                 # Residue to mutate; e.g. "A/23"
      mutant                         = "ALA"                                  # Mutation target; "ALA" or "GLY"
      mutant_only                    = 1                                      # Mutant energies only; 0/1
      cas_intdiel                    = 0                                      # Set intdiel by residue; 0/1
      intdiel_nonpolar               = 1                                      # Nonpolar intdiel; e.g. 1
      intdiel_polar                  = 3                                      # Polar intdiel; e.g. 3
      intdiel_positive               = 5                                      # Positive intdiel; e.g. 5
      intdiel_negative               = 5                                      # Negative intdiel; e.g. 5
    /
    ```

!!! info "Keep in mind"
    These inputs provide a practical starting point and can serve as the basis for comparable calculations across a
    related series. Review the available [input-file options][2], their accepted values, and the experimental values
    assigned to each system before interpreting a correlation.

## Systems in the series

The seven inputs differ only in the system name, experimental `Ki`, and-for the mutants-the selected residue:

- `WT`: wild type, `exp_ki=18` nM.
- `T13`: Thr13 to Ala, `exp_ki=20` nM.
- `H15`: His15 to Ala, `exp_ki=15` nM.
- `Y23`: Tyr23 to Ala, `exp_ki=50` nM.
- `Y25`: Tyr25 to Ala, `exp_ki=30` nM.
- `S32`: Ser32 to Ala, `exp_ki=37` nM.
- `W34`: Trp34 to Ala, `exp_ki=30` nM.

The values above are example inhibition constants bundled with this repository. The repository does not record their
source, assay conditions, or temperature, so treat them as illustrative inputs rather than a quantitative validation
set. `exp_ki` is expressed in nM and is used by the analyzer to derive the experimental binding energy for the
correlation.

## How this example works

The ST approximation extracts the two protein components from the shared complex trajectory for every calculation.
All systems process frames 1 through 10 with GB-Neck2 (`igb=8`), mbondi3 radii (`PBRadii=4`), and a salt
concentration of 0.15 M.

The wild-type directory evaluates the original complex. Each mutant input enables alanine scanning and sets
`mutant_only=1`, so its output contains the selected mutant rather than recalculating the shared wild type. Using the
same coordinates, frames, and implicit-solvent parameters across the series isolates the changes required for this
comparative example.

## Expected outputs

Each successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/GBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the correlation

After all seven calculations finish, remain in-or return to-the `Correlation` directory and recursively load the
results:

```bash
gmx_MMPBSA_ana -r
```

The analyzer reads the system name and experimental `Ki` stored with each result. Select the systems in the
Correlation panel to compare calculated and experimental binding energies. See the
[`gmx_MMPBSA_ana` documentation][5] for more details.

The following video demonstrates the correlation workflow in `gmx_MMPBSA_ana`:

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/QyaUTjmfYvc" frameborder="0" allowfullscreen></iframe>
</div>

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Correlation&fileName=gmx_MMPBSA-Correlation&rootDirectory=Correlation
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Correlation
