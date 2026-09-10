---
template: main.html
title: Protein-protein from PSF/DCD files
---

# Protein-protein binding from PSF/DCD files

This example converts a solvated CHARMM protein-protein simulation stored as PSF, CRD, and DCD files into the
GROMACS-compatible structure, trajectory, topology, and index required by `gmx_MMPBSA`. It then calculates the
binding free energy with the single-trajectory approximation and a linear PB solvent model.

<div class="example-card-grid" markdown>

-   **Source format**

    CHARMM PSF/CRD and DCD

-   **Conversion tools**

    AmberTools `cpptraj`, ParmEd, and GROMACS

-   **Protocol**

    Protein-protein, single trajectory

-   **Solvent model**

    Linear PB with CHARMM radii

</div>

!!! info "Representative system"
    This tutorial uses a protein-protein complex to demonstrate PSF/DCD conversion. The preparation workflow is not
    limited to protein-protein systems: other molecular compositions can be processed when conversion produces a
    compatible topology, structure, trajectory, index, and receptor/ligand selections. Individual calculation
    methods may impose additional restrictions.

!!! important "PSF and DCD are preparation inputs"
    `gmx_MMPBSA` does not consume the PSF and DCD files directly. The workflow first converts them to
    `gromacs.pdb`, `traj.xtc`, `gromacs.top`, and `index.ndx`. These converted files must describe the same atoms in
    the same order.

!!! warning "CHARMM CMAP conversion"
    The source PSF contains 280 CMAP cross-terms. The current GROMACS-to-AMBER topology conversion omits CHARMM CMAP
    terms and reports this during setup. This example exercises the complete conversion workflow, but quantitative
    CHARMM applications should assess the missing CMAP contribution before interpreting binding energies.

!!! note "CHARMM PB radii"
    `PBRadii=7` selects the `charmm_radii` set, which is intended only for systems prepared with CHARMM force fields.
    Its protein radii draw on work by [Nina, Belogv, and Roux][9], nucleic-acid radii on [Banavali and Roux][10],
    and additional elements on [Fortuna and Costa][11]. With `radiopt=0`, PBSA reads these radii from the generated
    AMBER topologies.

## Before you begin

The example contains the following source files:

<div class="example-card-grid" markdown>

-   **Topology and coordinates**

    `step3_input.psf` and its matching `step3_input.crd`

-   **Trajectory**

    `traj.dcd`

-   **CHARMM parameters**

    The parameter files under `toppar/`

-   **Conversion script**

    `script.py`, which creates `gromacs.pdb` and `gromacs.top`

</div>

Install `gmx_MMPBSA` in a dedicated environment containing AmberTools, ParmEd, and GROMACS before running the
conversion. See the [installation instructions][12].

The solvated source system contains 59,505 atoms, and `traj.dcd` contains 17 coordinate sets. `PROA` contains atoms
1–3,220 and acts as the receptor; `PROB` contains atoms 3,221–4,124 and acts as the protein ligand. The remaining
atoms are solvent and ions.

## Prepare the gmx_MMPBSA files

### 1. Convert the trajectory

Use `cpptraj` to remove water and ions from the PSF/DCD system and write the complete dry trajectory:

```bash
cpptraj -p step3_input.psf <<'EOF'
trajin traj.dcd
strip :POT,CLA,TIP3,LIT,SOD,RUB,CES,BAR
trajout traj.xtc
run
exit
EOF
```

This creates `traj.xtc` with 17 frames and the 4,124 protein atoms in the original PSF order.

### 2. Convert the structure and topology

Run the included ParmEd script:

```bash
python script.py
```

The script performs five operations:

```python
import parmed as pmd

psf = pmd.load_file('step3_input.psf')
psf.coordinates = pmd.load_file('step3_input.crd').coordinates
psf.strip(':POT, CLA, TIP3, LIT, SOD, RUB, CES, BAR')

chain_map = {'PROA': 'A', 'PROB': 'B'}
chain_residue_numbers = {chain: 0 for chain in chain_map.values()}
for residue in psf.residues:
    residue.chain = chain_map[residue.segid]
    chain_residue_numbers[residue.chain] += 1
    residue.number = chain_residue_numbers[residue.chain]

for number, atom in enumerate(psf.atoms, start=1):
    atom.number = number

pmd.formats.PDBFile.write(psf, 'gromacs.pdb', renumber=False)

params = pmd.charmm.CharmmParameterSet(
    'toppar/par_all36_carb.prm',
    'toppar/par_all36_cgenff.prm',
    'toppar/par_all36_lipid.prm',
    'toppar/par_all36m_prot.prm',
    'toppar/par_all36_na.prm',
    'toppar/par_interface.prm',
    'toppar/toppar_water_ions.str',
)
psf.load_parameters(params)
psf.save('gromacs.top', overwrite=True)
```

The solvent and ion mask must match the removal performed with `cpptraj`. The conversion maps the PSF segments
`PROA` and `PROB` to the valid one-character PDB chain IDs `A` and `B`. It also renumbers residues sequentially within
each derived chain because the source `PROB` segment begins with residue numbers -3 through 0. The source files are
not modified. Likewise, the parameter list must include every CHARMM parameter file required by the PSF.

!!! note "Active ATOMS sections"
    Some CHARMM-GUI NAMD parameter files comment out their `ATOMS`/`MASS` records with `!`. ParmEd requires those
    records when loading parameters. Use files with active `ATOMS` sections or uncomment the required records before
    running `script.py`. The parameter files bundled with this example are already prepared accordingly.

### 3. Create named molecular selections

Create an index containing the two proteins using their known atom ranges:

```bash
gmx select \
  -s gromacs.pdb \
  -on index.ndx \
  -select '"PROA" atomnr 1 to 3220; "PROB" atomnr 3221 to 4124'
```

The resulting `PROA` and `PROB` groups are used directly with `-cg`; their numerical group positions do not need to
be tracked. See the [GROMACS selection syntax][13] for additional ways to define static index groups.

### 4. Check the converted files

Confirm the structure and trajectory atom counts before starting the calculation:

```bash
gmx check -f gromacs.pdb
gmx check -f traj.xtc
```

Both files must report 4,124 atoms. If they differ, revisit the solvent/ion stripping masks before continuing.

## Run the example

[Download the PSF/DCD protein-protein example as a ZIP archive][6].

Extract the archive, complete the conversion steps above, and choose either the serial or MPI command. You can also
[view the source files on GitHub][8] before downloading them.

=== "Serial"

    ```bash
    gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs gromacs.pdb \
      -ct traj.xtc \
      -ci index.ndx \
      -cg PROA PROB \
      -cp gromacs.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

=== "With MPI"

    ```bash
    mpirun -np 2 gmx_MMPBSA -O \
      -i mmpbsa.in \
      -cs gromacs.pdb \
      -ct traj.xtc \
      -ci index.ndx \
      -cg PROA PROB \
      -cp gromacs.top \
      -o FINAL_RESULTS_MMPBSA.dat \
      -eo FINAL_RESULTS_MMPBSA.csv
    ```

See the [complete command-line reference][1] for all options.

## Configure the calculation

The example uses the concise `mmpbsa.in` shown first below. The all-options version was generated with
`gmx_MMPBSA --create_input pb` and then adapted with the same example-specific values. The concise block is the runnable starting point; the generated block includes additional options and defaults, so the two blocks are not textually identical. Both blocks describe the same
linear PB calculation.

=== "Concise input"

    ```yaml linenums="1" title="mmpbsa.in"
    Sample input file for PB calculation from converted PSF/DCD files
    # This input provides a practical starting point for a CHARMM protein-protein
    # system. Review the PB settings and the CMAP limitation before reuse.

    &general
    sys_name="PSF-DCD-Prot-Prot",
    startframe=5,
    endframe=15,
    PBRadii=7,
    assign_chainID=1,
    solvated_trajectory=0,
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
      sys_name                       = "PSF-DCD-Prot-Prot"                  # System name; e.g. "complex"
      startframe                     = 5                                      # First frame; e.g. 1
      endframe                       = 15                                     # Last frame; e.g. 100
      interval                       = 1                                      # Frame interval; e.g. 1
      forcefields                    = "leaprc.protein.ff14SB"               # Force fields; e.g. "leaprc.protein.ff14SB"
      ions_parameters                = 1                                      # Ion params; e.g. 1
      PBRadii                        = 7                                      # PB radii set; 1-7
      temperature                    = 298.15                                 # Temperature (K); e.g. 298.15
      qh_entropy                     = 0                                      # Legacy QH output reader; new calculations reject 1
      interaction_entropy            = 0                                      # Run IE entropy; 0/1
      ie_segment                     = 25                                     # IE segment length (%); e.g. 25
      c2_entropy                     = 0                                      # Run C2 entropy; 0/1
      assign_chainID                 = 1                                      # Assign chain IDs; 0/1
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
      sprob                          = 0.557                                  # SASA probe (A); e.g. 0.557
      vprob                          = 1.3                                    # Volume probe (A); e.g. 1.3
      rhow_effect                    = 1.129                                  # Water density; e.g. 1.129
      use_sav                        = 1                                      # Use SAV cavity; 0/1
      cavity_surften                 = 0.0378                                 # Cavity surften; e.g. 0.0378
      cavity_offset                  = -0.5692                                # Cavity offset; e.g. -0.5692
      maxsph                         = 400                                    # Max surface dots; e.g. 400
      maxarcdot                      = 1500                                   # Max arc dots; e.g. 1500
      npbverb                        = 0                                      # PB verbosity; 0/1
    /
    ```

!!! info "Keep in mind"
    This input provides a practical starting point for the converted CHARMM protein-protein system. Review sampling,
    PB radii, dielectric treatment, grid convergence, and the CMAP limitation for the intended application. Additional
    [input-file options][2] may be needed for a production protocol.

## How this example works

The source PSF/CRD/DCD files are used only during preparation. `cpptraj` and ParmEd independently remove the same
solvent and ion residues so that `gromacs.pdb`, `traj.xtc`, and `gromacs.top` retain an identical 4,124-atom ordering.
The derived PDB uses chain `A` for PSF segment `PROA` and chain `B` for segment `PROB`; residue numbering restarts at
1 in each chain. The index then assigns the first 3,220 atoms to receptor `PROA` and the remaining 904 atoms to ligand
`PROB`.

The single-trajectory approximation extracts both proteins from every selected complex frame. Because the converted
trajectory is already dry, `solvated_trajectory=0` avoids a redundant stripping step. The explicit chain identifiers
in `gromacs.pdb` are therefore retained throughout the calculation.

The calculation processes 11 of the 17 available frames (frames 5 through 15) with the linear PB equation, an ionic
strength of 0.15 M, and CHARMM-specific topology radii. Since `gromacs.top` is supplied, no `forcefields` variable is
needed in the concise input; the bonded and nonbonded parameters come from the converted CHARMM topology. CMAP terms
are the stated exception.

## Expected outputs

A successful calculation produces:

- `FINAL_RESULTS_MMPBSA.dat`: the MM/PBSA summary and binding-energy statistics.
- `FINAL_RESULTS_MMPBSA.csv`: the per-frame energy terms requested with `-eo`.

## Analyze the results

Open the results with `gmx_MMPBSA_ana` for interactive inspection and plotting. See the
[`gmx_MMPBSA_ana` documentation][5] for usage details.

  [1]: ../../../docs/gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../../docs/input_file.md#the-input-file
  [5]: ../../../docs/analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/psf_dcd/protein_protein&fileName=gmx_MMPBSA-PSF-DCD-Protein-Protein&rootDirectory=protein_protein
  [8]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/psf_dcd/protein_protein
  [9]: https://doi.org/10.1021/jp970736r
  [10]: https://doi.org/10.1021/jp025852v
  [11]: https://doi.org/10.1021/acs.jcim.1c00177
  [12]: ../../../docs/installation.md
  [13]: https://manual.gromacs.org/current/onlinehelp/selections.html
