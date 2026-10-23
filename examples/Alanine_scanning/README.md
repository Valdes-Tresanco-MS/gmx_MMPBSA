---
template: main.html
title: Alanine scanning
---

# Alanine scanning

This example estimates the effect of mutating receptor residue Tyr23 to alanine in a protein-protein complex. It
uses the same solvated simulation, frame range, and GB-Neck2 model as the decomposition and stability examples.
The [multiple-mutant example](#run-a-multiple-alanine-or-glycine-mutant) below extends this calculation to two
receptor residues mutated together.

<div class="example-card-grid" markdown>

-   **Analysis**

    Computational alanine scanning

-   **Mutation**

    Single: Tyr23 to Ala

    Multiple: Thr13 + Tyr23 to Ala or Gly

-   **Solvent model**

    GB-Neck2 (`igb=8`)

-   **Bundled test**

    `gmx_MMPBSA_test -t 12`

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

The folder includes the solvated system snapshot `com.pdb`, matching the corresponding files in the decomposition
and stability examples. The 608-atom `reference.pdb` structure (`-cr`) preserves the chain ID and residue number used
by `mutant_res="R:23"`. See the [complete command-line reference][1] for all options.

## Run the example

### Run the bundled test

The quickest way to reproduce this example is through the test runner:

```bash
gmx_MMPBSA_test -t 12
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the Alanine Scanning example as a ZIP archive][6].

Extract the archive, change to the `Alanine_scanning` directory, and choose either the serial or MPI command. You can also
[view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs com.tpr \
      -ct com_traj.xtc \
      -ci index.ndx \
      -cg SOLU_chain1 SOLU_chain2 \
      -cp topol.top \
      -cr reference.pdb \
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
      -cg SOLU_chain1 SOLU_chain2 \
      -cp topol.top \
      -cr reference.pdb \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the minimal `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input gb ala` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks therefore
describe the same MM/GBSA alanine-scanning calculation.

=== "Minimal input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for Alanine scanning
    # This input provides a practical starting point for alanine scanning.
    # Verify the mutation selection and model settings for your system.

    &general
    sys_name="Alanine_Scanning",
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150,
    /
    &alanine_scanning
    mutant='ALA', mutant_res='R:23', cas_intdiel=1,
    /
    ```

=== "Generated input — all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input gb ala"
    Input block generated for the 1.7.0 release.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "Alanine_Scanning"                       # System name; e.g. "complex"
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

    # Alanine scanning namelist variables
    &alanine_scanning
      mutant_res                     = "R:23"                                     # Residue to mutate; e.g. "A/23"
      mutant                         = "ALA"                                  # Mutation target; "ALA" or "GLY"
      mutant_only                    = 0                                      # Mutant energies only; 0/1
      cas_intdiel                    = 1                                      # Set intdiel by residue; 0/1
      intdiel_nonpolar               = 1                                      # Nonpolar intdiel; e.g. 1
      intdiel_polar                  = 3                                      # Polar intdiel; e.g. 3
      intdiel_positive               = 5                                      # Positive intdiel; e.g. 5
      intdiel_negative               = 5                                      # Negative intdiel; e.g. 5
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

`mutant_res="R:23"` identifies Tyr23 in receptor chain `R`, and `mutant="ALA"` creates the corresponding alanine
variant without requiring a separately prepared mutant topology. With `mutant_only=0`, the calculation evaluates
both the original and mutant systems. Because tyrosine is classified as polar, `cas_intdiel=1` applies
`intdiel_polar=3` to the relevant GB calculations.

## Run a multiple alanine or glycine mutant

To estimate the combined effect of two mutations, select Thr13 and Tyr23 together with `mutant_res="R/13,23"`.
Both residues are present in chain `R` of the supplied `reference.pdb` and belong to the receptor. With
`mutant="ALA"`, this produces the **T13A/Y23A multiple alanine mutant**. With `mutant="GLY"`, it produces the
**T13G/Y23G multiple glycine mutant**.

Make a separate copy of the example folder for each variant so its inputs and results stay together. Save one of
the following inputs as `mmpbsa_multiple.in` in that folder. Both variants use the same trajectory frames, radii,
salt concentration, and internal dielectric.

=== "Multiple alanine mutant"

    ```yaml linenums="1" title="mmpbsa_multiple.in"
    &general
    sys_name="Multiple_Alanine_Mutant",
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150, intdiel=3,
    /
    &alanine_scanning
    mutant="ALA", mutant_res="R/13,23",
    cas_intdiel=0, mutant_only=0,
    /
    ```

=== "Multiple glycine mutant"

    ```yaml linenums="1" title="mmpbsa_multiple.in"
    &general
    sys_name="Multiple_Glycine_Mutant",
    startframe=1,
    endframe=10,
    PBRadii=4,
    /
    &gb
    igb=8, saltcon=0.150, intdiel=3,
    /
    &alanine_scanning
    mutant="GLY", mutant_res="R/13,23",
    cas_intdiel=0, mutant_only=0,
    /
    ```

Multiple mutations require `cas_intdiel=0`; choose the dielectric explicitly in the solvent-model namelist.
Here, `intdiel=3` matches the dielectric applied to the single Tyr23 mutation above. For comparisons with other
single mutants, use the same explicit dielectric and the same frame range in every calculation. All selected
residues must belong to one component and use the same target, either `ALA` or `GLY`.

From the folder containing the chosen input, run:

```bash
gmx_MMPBSA -O \
  -i mmpbsa_multiple.in \
  -cs com.tpr \
  -ct com_traj.xtc \
  -ci index.ndx \
  -cg SOLU_chain1 SOLU_chain2 \
  -cp topol.top \
  -cr reference.pdb \
  -o MULTIPLE_MUTANT_RESULTS.dat \
  -eo MULTIPLE_MUTANT_RESULTS.csv
```

For MPI, prepend `mpirun -np 2` to the same command.

### Interpret the multiple-mutant result

With `mutant_only=0`, the output includes the original system, the multiple mutant, and their difference. The
mutation effect is the change in binding energy for **both substitutions together**, relative to the original
system. This run does not calculate separate T13A and Y23A effects (or T13G and Y23G effects); run each selection
separately if those individual effects are needed. The combined effect need not equal the sum of the individual
effects.

The program generates mutant topologies and trajectory coordinates from the original system. Check that both
selected residues have become ALA or GLY in the generated mutant structures before interpreting the result.

## Verify the mutation selection

Residue selection depends on stable chain IDs and residue numbers. Supplying `reference.pdb` with `-cr` makes the
mapping explicit while `com.pdb` retains the complete solvated snapshot for consistency. Confirm the selected
mutation in `_GMXMMPBSA_COM_FIXED.pdb` and the generated mutant structures before interpreting results.

One run can mutate a single residue or change several residues together to produce a multiple alanine mutant or
a multiple glycine mutant. It does not independently scan each selected residue. `mutant_res` also accepts `CHAIN/RESNUM` notation and insertion codes as described in the
[alanine-scanning input options][9].

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the original, mutant, and alanine-scanning energy summaries.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../docs/input_file.md#the-input-file
  [3]: ../../docs/input_file.md#sample-input-files
  [5]: ../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Alanine_scanning&fileName=gmx_MMPBSA-Alanine-Scanning&rootDirectory=Alanine_scanning
  [7]: ../../docs/examples/gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Alanine_scanning
  [9]: ../../docs/input_file.md#alanine_scanning-namelist-variables
