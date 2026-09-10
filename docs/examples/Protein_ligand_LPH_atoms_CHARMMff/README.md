---
template: main.html
title: Halogenated CHARMM ligand with LPH sites
---

# Halogenated CHARMM ligand with LPH virtual sites

This example shows how a CHARMM protein-ligand system containing lone-pair halogen (`LPH`) virtual sites can be
converted into an input that gmx_MMPBSA can process. The calculation uses the single-trajectory approximation and
the linear PB model.

<div class="example-card-grid" markdown>

-   **Protocol**

    Single trajectory

-   **Force field**

    CHARMM with LPH sites removed

-   **Solvent model**

    Linear PB with CHARMM radii

-   **Bundled test**

    `gmx_MMPBSA_test -t 22`

</div>

!!! info "Representative system"
    The protein-ligand complex is a representative CHARMM system. The CHARMM topology workflow is not limited to
    protein-ligand complexes, although the manual LPH-removal steps on this page apply specifically to unsupported
    ligand virtual sites.

!!! danger "Removing LPH sites changes the electrostatic model"
    gmx_MMPBSA does not support the massless LPH sites in the original topology. When removing an LPH site, transfer
    its charge to the parent halogen instead of simply deleting it. In this example, each removed LPH site carries
    `+0.05 e`; therefore, the charge of each parent bromine is changed from `-0.180 e` to `-0.130 e`. This preserves
    the original ligand charge of `-1.00 e`.

    Charge conservation does not restore the off-center positive sites or their directional halogen-bonding
    electrostatics. For quantitative work, a validated site-free ligand parameterization remains preferable to
    treating LPH removal and charge transfer as a complete reparameterization.

!!! warning "CHARMM CMAP conversion"
    The current GROMACS-to-AMBER topology conversion also omits CHARMM CMAP terms and reports this during setup.
    Quantitative CHARMM applications should assess this additional approximation before interpreting binding energies.

!!! note "Halogen and CHARMM PB radii"
    `PBRadii=7` selects `charmm_radii`, including radii of 1.86 Å for Cl, 1.98 Å for Br, and 2.24 Å for I. These
    radii are intended for CHARMM systems without explicit halogen extra-point charges. With `radiopt=0`, PBSA reads
    the assigned radii from the generated AMBER topologies. See the underlying CHARMM radii sources for
    [proteins][9], [nucleic acids][10], and [additional elements][11].

LPH sites are positively charged virtual particles placed near halogens to represent the anisotropic electrostatic
potential associated with halogen bonding. See the [LPH parameterization study][12] for background.

## Before you begin

The runnable, LPH-free calculation uses the following files and selections:

<div class="example-card-grid" markdown>

-   **Calculation settings**

    `mmpbsa.in` (`-i`)

-   **Prepared system**

    LPH-free structure `str_noLP.pdb` (`-cs`) and modified topology `topol.top` (`-cp`). Keep the `toppar` directory
    containing the referenced CHARMM `*.itp` files beside `topol.top`.

-   **Prepared trajectory**

    LPH-free, fitted trajectory `com_traj.xtc` (`-ct`)

-   **Molecular selections**

    Index `index_mod_gromacs.ndx` (`-ci`) with receptor `Protein` and LPH-free ligand `lig` (`-cg`)

</div>

The prepared complex contains 5,610 atoms: a 5,580-atom protein and a 30-atom ligand. The original ligand and
solvated system contain 32 and 70,483 atoms, respectively. See the [complete command-line reference][1] for all
options.

## Prepare an LPH-free input

The runnable files are already included with the example. The following steps document how they were derived from
`com.tpr`, `traj_fit.xtc`, and the original LPH-containing ligand topology.

### 1. Create LPH-free index groups

Start `make_ndx` with the original TPR:

```bash
gmx make_ndx -f com.tpr -o index_mod_gromacs.ndx
```

In the interactive prompt, split the original 32-atom ligand group (`13`), select its two LPH sites, exclude them,
and combine the resulting 30-atom ligand with the protein:

```text
splitat 13
47|48
13&!49
name 50 lig
1|50
del 17-49
q
```

After the intermediate groups are deleted and the remaining groups are renumbered, `lig` is group 17 and
`Protein_lig` is group 18. The latter contains 5,610 atoms.

### 2. Strip the LPH sites from the coordinates

Use the `Protein_lig` group to create a matching structure and trajectory:

```bash
echo 18 | gmx trjconv \
  -s com.tpr \
  -f traj_fit.xtc \
  -dump 0 \
  -n index_mod_gromacs.ndx \
  -o str_noLP.pdb

echo 18 | gmx trjconv \
  -s com.tpr \
  -f traj_fit.xtc \
  -n index_mod_gromacs.ndx \
  -o com_traj.xtc
```

Inspect `str_noLP.pdb` and confirm that it contains 5,610 atoms and no `LP1` or `LP2` records.

### 3. Prepare a matching topology

The example retains the source topology as `toppar/HETA_original_with_LPH_info.itp` and uses the modified
`toppar/HETA.itp` in `topol.top`. Relative to the source file, prepare the modified file manually as follows:

- Remove atoms 31 and 32 (`LP1` and `LP2`) from `[ atoms ]`.
- Remove every `[ pairs ]` entry involving atoms 31 or 32.
- Delete the `[ virtual_sites3 ]` definitions for the two sites.
- Delete the corresponding `[ exclusions ]` records.

Then transfer each removed site charge to its parent halogen in `[ atoms ]`:

```text
BR1: -0.180 + 0.050 = -0.130
BR2: -0.180 + 0.050 = -0.130
```

Finally, sum the ligand charges and confirm that they remain equal to the original molecular charge (`-1.00 e` in
this example). Apply the same accounting to the actual LPH charges and parent halogens in your topology; do not
assume that every LPH model uses `+0.05 e`.

Coordinate removal and topology removal must be performed together so that the atom order and count remain
consistent. A reusable scientific model also requires validation of the resulting site-free charge distribution.

## Run the example

### Run the bundled test

The quickest way to reproduce the prepared example is through the test runner:

```bash
gmx_MMPBSA_test -t 22
```

See the [`gmx_MMPBSA_test` documentation][7] for download, selection, and cleanup options.

### Run it manually

[Download the LPH CHARMM example as a ZIP archive][6].

Extract the archive, change to the `Protein_ligand_LPH_atoms_CHARMMff` directory, and choose either the serial or
MPI command. You can also [view the example files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs str_noLP.pdb \
      -ct com_traj.xtc \
      -ci index_mod_gromacs.ndx \
      -cg Protein lig \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs str_noLP.pdb \
      -ct com_traj.xtc \
      -ci index_mod_gromacs.ndx \
      -cg Protein lig \
      -cp topol.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input pb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks describe the same
linear PB calculation using the already prepared LPH-free files.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for PB calculation with a halogenated CHARMM ligand
    # This input demonstrates the PB setup for the bundled topology from which the
    # LPH virtual sites have been removed. Review the approximation before reuse.

    &general
    sys_name="Prot-Lig-LPH-CHARMM",
    startframe=5,
    endframe=9,
    solvated_trajectory=0,
    PBRadii=7,
    /
    &pb
    radiopt=0, istrng=0.150, fillratio=1.25, inp=1,
    cavity_surften=0.005, cavity_offset=0.0,
    /
    ```

=== "Generated input - all options"

    ```yaml linenums="1" title="mmpbsa.in generated with --create_input pb"
    Input block generated for the 1.7.0 release; the generator's development-version header is omitted from this documentation.
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
      sys_name                       = "Prot-Lig-LPH-CHARMM"                # System name; e.g. "complex"
      startframe                     = 5                                      # First frame; e.g. 1
      endframe                       = 9                                      # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "leaprc.protein.ff14SB"               # Force fields; e.g. "leaprc.protein.ff14SB"
      ions_parameters                = 1                                      # Ion params; e.g. 1
      PBRadii                        = 7                                      # PB radii set; 1-7
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
      solvated_trajectory            = 0                                      # Clean solvated traj.; 0/1
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
      fillratio                      = 1.25                                   # Grid fill ratio; e.g. 4
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
      sprob                          = 0.557                                  # SASA probe (A); e.g. 0.557
      vprob                          = 1.3                                    # Volume probe (A); e.g. 1.3
      rhow_effect                    = 1.129                                  # Water density; e.g. 1.129
      use_sav                        = 1                                      # Use SAV cavity; 0/1
      cavity_surften                 = 0.005                                  # Cavity surften; e.g. 0.0378
      cavity_offset                  = 0.0                                    # Cavity offset; e.g. -0.5692
      maxsph                         = 400                                    # Max surface dots; e.g. 400
      maxarcdot                      = 1500                                   # Max arc dots; e.g. 1500
      npbverb                        = 0                                      # PB verbosity; 0/1
    /
    ```

!!! warning "Interpretation"
    The PB settings and CHARMM halogen radii are internally documented, but they do not restore the directional LPH
    electrostatics. The manual charge transfer preserves the ligand's `-1.00 e` total charge, but a validated
    site-free ligand model is still preferable before drawing quantitative conclusions.

## How this example works

The prepared trajectory already contains only `Protein` and the 30-atom `lig`, so
`solvated_trajectory=0` prevents an unnecessary solvent-stripping step. The calculation processes frames 5 through
9 with the linear PB equation, an ionic strength of 0.15 M, and the CHARMM-specific PB radii.

Because `topol.top` is provided, no `forcefields` variable is needed in the concise input; the CHARMM parameters are
read from the topology include tree. The topology has already been modified to match the LPH-free structure and
trajectory, with its total charge preserved and the directional-electrostatics and CMAP limitations stated above.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/PBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand_LPH_atoms_CHARMMff&fileName=gmx_MMPBSA-Protein-Ligand-LPH-CHARMM&rootDirectory=Protein_ligand_LPH_atoms_CHARMMff
  [7]: ../gmx_MMPBSA_test.md#gmx_mmpbsa_test-command-line
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Protein_ligand_LPH_atoms_CHARMMff
  [9]: https://doi.org/10.1021/jp970736r
  [10]: https://doi.org/10.1021/jp025852v
  [11]: https://doi.org/10.1021/acs.jcim.1c00177
  [12]: https://doi.org/10.1016/j.bmc.2016.06.034
