---
template: main.html
title: Application note draft - Advances from v1.4.3 to v1.7.0
---

# gmx_MMPBSA advances from v1.4.3 to v1.7.0: structured report for a JCIM application note

!!! note "Draft status — updated through v1.7.0 and subsequent development"
    This structured report uses `gmx_MMPBSA` v1.4.3 as the published-era
    baseline and the published v1.7.0 tag (`8d7aad33`) as the release endpoint. It distinguishes the
    cumulative v1.5.x-v1.6.x development from additions and behavior changes
    introduced in v1.7.0. A separate section records subsequent development
    verified through commit `ac637eba` on 2026-10-03; those changes are not
    attributed to the published tag. Users reproducing or comparing calculations should
    also consult the [changelog](changelog.md), the
    [v1.6.5-to-v1.7.0 migration guide](compatibility.md#migrating-from-165-to-170),
    and the current [input](input_file.md) and [output](output.md) references.

## Abstract / overview

Since the original publication of `gmx_MMPBSA` and the v1.4.3 release, the
software has evolved from a GROMACS-oriented interface to AmberTools end-state
free-energy calculations into a broader calculation and analysis platform for
MM/PB(GB)SA workflows. The v1.5.x series redesigned the calculation, output,
input-processing, analyzer, and Python API layers. The v1.6.x series added
GBNSR6, expanded system compatibility, and strengthened analysis and runtime
support. Version 1.7.0 consolidates those developments into a more explicit and
reproducible workflow: GROMACS calculations now require the original complex
topology, native AMBER files have a dedicated `amber_MMPBSA` entry point,
continuum-radius assignment is recorded, and output, logging, error reporting,
testing, and dependency boundaries are more rigorously defined.

The scientific scope has expanded through nonlinear PB, C2 entropy, ALPB,
PC+-corrected 3D-RISM, GBNSR6, improved QM/MMGBSA, composite alanine/glycine
mutations, decomposition analysis, membrane workflows, explicit receptor
waters, and experimental normal-mode support for CHARMM topologies. Version
1.7.0 also changes how several results must be interpreted: it corrects the
Interaction Entropy (IE) estimator, replaces C2 uncertainty resampling with
block-based diagnostics, rejects new quasi-harmonic (QH) calculations, and
changes the omitted implicit-solvent defaults. These changes improve the
scientific contract of new calculations but mean that v1.7.0 results should not
be assumed numerically identical to v1.6.5 results.

## Comparison framework

### Published-era baseline: v1.4.3

Version 1.4.3 already provided GROMACS users with AmberTools-based MM/PBSA and
MM/GBSA calculations, QH, normal-mode and IE entropy calculations, decomposition, alanine/glycine scanning, a graphical analyzer, examples, and an early Python API. The v1.4.x
series also included multi-system analysis, correlation analysis, improved
decomposition selection, PyMOL visualization, and an established documentation
set. Basic membrane PB controls, output rewriting, and cleanup were also
already available. Version 1.4.3 is the publication-era baseline used here;
the first software release was v1.0.0.

### Intermediate development: v1.5.x-v1.6.5

The v1.5.x cycle was a breaking redesign. It reorganized input variables and
calculation processing, added C2 entropy and new implicit-solvent controls,
introduced a new analyzer/API data model, and improved support for additional
force fields and preparation workflows. Versions 1.5.1-1.5.7 then expanded
QM/MM, OPLS, PSF/DCD, 3D-RISM, analyzer performance, correlation, logging, and
output correctness. Version 1.6.0 added GBNSR6 with enthalpy and decomposition,
while the remaining v1.6.x releases focused on compatibility, correctness, and
experimental CHARMM normal-mode support.

### Current endpoint: v1.7.0

Version 1.7.0 is both a feature release and a reproducibility boundary. Its new
capabilities include native-AMBER execution, explicit receptor waters,
composite mutations, automatic membrane parameter detection, GBNSR6 topology
handling, a modernized API, and stronger
validation infrastructure. At the same time, it changes topology requirements,
implicit-solvent defaults, the IE estimator, uncertainty reporting, output
paths, and supported QH behavior. Consequently, the most useful comparison is
not simply “more methods than v1.4.3,” but “broader methods plus a clearer
scientific and provenance contract in v1.7.0.”

## Scientific advances

### Entropy calculations and uncertainty

Version 1.5.0 introduced C2 entropy and expanded the output and diagnostic
handling for IE and C2. Subsequent releases corrected rewrite-output behavior,
alanine-scanning entropy differences, duplicated analyzer items, and
`ie_segment` handling.

Version 1.7.0 makes a more fundamental correction. IE now uses one shared
running ensemble mean for each trajectory prefix and a numerically stable
log-sum-exp evaluation. The full selected ensemble is the primary IE estimate;
`ie_segment` is retained only as a tail-convergence diagnostic. IE and C2 also
report deterministic, nonoverlapping block diagnostics. Frame-based SD and SEM
remain available for compatibility, while Block SD and Block SEM describe
between-block variation and should be reported with the block size and number
of blocks. Short trajectories necessarily provide weak block evidence.

C2 retains the central estimator `sigma(E)^2 / (2RT)` used in v1.6.5, with a
more precise molar gas constant. Its principal statistical change replaces
random frame resampling and trimmed bootstrap summaries with deterministic
block estimates, Block SD, Block SEM, and block percentiles. The C2 warning
threshold for internal-energy SD is also corrected to approximately
6.0 kcal/mol (25 kJ/mol). These are distinct from the correction to the IE
running estimator. Nonoverlapping blocks do not by themselves establish
statistical independence; block-size sensitivity and adequate sampling remain
necessary, and block percentiles are not a confidence interval for the
full-ensemble estimate.

The entropy portfolio in v1.7.0 includes normal-mode, IE, and C2 approaches.
Experimental normal-mode support for CHARMM topologies was introduced in
v1.6.5, and v1.7.0 leaves unconverged normal-mode frames as missing values
rather than replacing them with the mean of converged frames. New QH
calculations are unsupported in v1.7.0; historical QH results remain readable
during this compatibility release, with QH support scheduled for removal
afterward.

### PB, GB, ALPB, and GBNSR6

Compared with v1.4.3, the current release exposes a substantially broader
implicit-solvent portfolio. Version 1.5.0 enabled nonlinear PB through `sander`
and exposed additional PBSA controls. Version 1.5.2 added the Analytical
Linearized Poisson-Boltzmann (ALPB) approximation, and v1.6.0 introduced
GBNSR6 enthalpy, per-residue decomposition, and pairwise decomposition.

Version 1.7.0 changes the defaults used when an input omits key
implicit-solvent settings: `igb` changes from 5 to 8, `PBRadii` from 3 to 4
(`mbondi3`), and PB `exdi` from 80.0 to 78.5. Explicit values remain unchanged.
These defaults align the GB model and radius pairing and update the PB external
dielectric, but they can change calculated energies. Reproduction of a v1.6.5
calculation therefore requires the historical values to be specified explicitly.

The release also adds automatic GBNSR6 topology compaction and post-processing
for legacy AmberTools atom-type limits and corrects parser/frame-term merging.
The legacy `sander_apbs` execution route is no longer available for new
calculations; v1.7.0 uses the built-in PBSA solver while retaining the old field
only for reading archived results.

### 3D-RISM

The v1.5.x series exposed the documented 3D-RISM variables, bundled
precalculated `.xvv` files, moved execution from `rism3d.snglpnt` to `sander`,
and added PC+ correction support. These changes integrate standard, Gaussian
fluctuation, and PC+-corrected 3D-RISM results into the common output and API
model.

The application note should distinguish package functionality from external
runtime compatibility. Some AmberTools builds can fail before a 3D-RISM
calculation because of their Fortran runtime linkage. This is an external
AmberTools/runtime compatibility problem rather than a change in the
`gmx_MMPBSA` input model, and the example documentation provides a diagnostic
path for it.

### QM/MMGBSA, explicit receptor waters, and mutations

The v1.5.x cycle expanded QM/MMGBSA inputs, automatic charge calculation, and
residue-selection behavior. Version 1.7.0 changes the default QM theory to
`PM6-DH+` and treats an unconverged QM/MM SCF calculation as a hard error rather
than allowing an ambiguous result to continue.

Version 1.7.0 also adds single-trajectory explicit receptor-water workflows for
GB, GBNSR6, PB, 3D-RISM, normal mode, and QM/MMGBSA. Support is method- and
trajectory-mode-specific: it should not be described as unrestricted retention
of arbitrary solvent in every calculation type.

Alanine/glycine scanning has accumulated fixes for CHARMM systems, terminal
residues, THR-to-ALA mutation, mutant-normal output, and entropy differences.
In v1.7.0, `mutant_res` can select several residues for one composite mutation.
All selected residues must belong to the same component and share one mutation
target. The reported result is the combined effect of that composite mutant,
not an independent scan of each selected residue.

### Decomposition and membrane systems

Decomposition analysis was restructured across the output, analyzer, and API
layers after v1.4.3. The current implementation supports improved per-residue
and pairwise handling, residue selection, term-level tables and plots, and
GBNSR6 decomposition. Version 1.7.0 additionally prevents inactive `&decomp`
template defaults from leaking into ordinary `sander` or GBNSR6 input files.

Basic PB membrane controls were already exposed in v1.4.3. Later development
expanded the documented membrane-protein workflows, including CHARMM systems
and uniform or heterogeneous slab-like dielectric profiles. Version 1.7.0 adds
automatic membrane thickness (`mthick`) and center (`mctrdz`) calculation from
the selected frames of the original, unstripped complex trajectory. Atom names
are selected with `membrane_atoms`; either parameter can also be set numerically.
The `--create_input pb_mem` template provides membrane-oriented settings.

Automatic detection retains `GMXMMPBSA_membrane_parameters.csv` and
`GMXMMPBSA_membrane_parameters.png` regardless of `keep_files`, recording
resolved parameters and frame-wise diagnostics. The trajectory must already
be continuous across periodic boundaries and oriented with the membrane normal
along z; detection does not unwrap or reorient it. These controls improve
preparation and provenance while retaining the method-specific limitations of
an end-state implicit-solvent treatment.

## Topology, system, and file support

### GROMACS topology as required provenance

Post-v1.4.3 releases expanded support for Amber-, CHARMM-, and OPLS-derived
GROMACS systems; CHARMM-GUI preparation; additional histidine variants;
lone-pair detection; terminal atoms; and topology/coordinate consistency
checks. Version 1.7.0 makes the GROMACS topology contract explicit: every
`gmx_MMPBSA` calculation requires the complex topology through `-cp`, including
its referenced `*.itp` files. For unbound multiple-trajectory inputs, matching
receptor or ligand topologies are also required.

This replaces the historical structure-only route that rebuilt Amber
topologies with `tleap`. A ligand `mol2` file is no longer a substitute for the
original GROMACS topology tree. The change makes force-field parameters and
topology provenance part of the required calculation input instead of relying
on reconstruction from extracted coordinates.

### Native AMBER and PSF/DCD preparation routes

The PSF/DCD tutorial documents how a protein-protein system prepared for NAMD,
OpenMM, GENESIS, or a related workflow can be converted into the inputs used by
`gmx_MMPBSA`; PSF and DCD are preparation sources rather than direct
`gmx_MMPBSA` calculation arguments.

Version 1.7.0 adds the independent `amber_MMPBSA` command for native AMBER
topology, trajectory, and mask workflows. It bypasses GROMACS topology
conversion and shares the calculation and analysis stack where supported. It
is not a full feature-parity replacement for `gmx_MMPBSA`: single- versus
multiple-trajectory support, explicit waters, QM/MM, entropy, and ligand
restrictions remain defined by the native-AMBER support matrix.

### Radius provenance

Radius handling differs deliberately between the two topology routes. Native
AMBER workflows preserve the `RADII`, `SCREEN`, and `RADIUS_SET` data in the
input topology. GROMACS conversion applies the selected `PBRadii` through
ParmEd. Version 1.7.0 warns when the effective radius set is not conventionally
paired with the selected GB model and writes `GMXMMPBSA_radii.json` for every
calculation. The file records the requested and effective radius set,
assignment route, model and force-field context, and checksums of the final
radius and screening arrays. Optional `radii_audit=1` adds per-atom audit data.

## Analyzer, API, and user experience

### Analyzer redesign

The v1.5.x `gmx_MMPBSA_ana` redesign remains one of the largest visible changes
from v1.4.3. It uses the modern Python API backend and is not generally
compatible with v1.4.3 result files. The redesigned analyzer supports:

- multi-system loading from individual info files or directory searches;
- selection of calculation types, subsystems, components, mutants, and
  decomposition data;
- line and bar plots, heatmaps, PyMOL residue-energy visualization, tables, and
  output-file viewers;
- configurable chart themes, labels, dimensions, resolution, and formats;
- frame-range, interval, and time-conversion controls; and
- multi-system correlation and regression workflows using ΔG or ΔΔG.

The v1.5.5 documentation reports large loading-time improvements after earlier
v1.5.2 performance problems. Those values should be cited as historical
project benchmarks, not as independently reproduced v1.7.0 benchmarks.
Version 1.7.0 adds theme and plotting robustness updates rather than claiming a
new analyzer performance benchmark.

### Python API and result outputs

The current API is also the data layer used by `gmx_MMPBSA_ana`. It loads
`_GMXMMPBSA_info` or compact `.mmxsa` results into pandas-based structures for
scripted analysis:

```python
from GMXMMPBSA import API

api = API.load("COMPACT_MMXSA_RESULTS.mmxsa")
```

The API exposes metadata, input namelists, file information, energy and entropy
data, and decomposition data. Version 1.7.0 modernizes the loader and includes a
runnable API example, strengthening the API as a supported route for
reproducible downstream analysis rather than only an analyzer implementation
detail.

Per-frame CSV output is now derived automatically from the summary filename,
with explicit output options taking precedence. Filename collisions are
rejected before output is opened. This makes summary and per-frame data easier
to discover while requiring older automation to avoid assuming that both can
share a path.

## Reproducibility, diagnostics, and validation

Version 1.7.0 improves MPI-safe logging, gives rank 0 ownership of the main log,
adds record-based warning and error totals, and provides selectable progress
styles. Warning severity now distinguishes expected automatic actions from
scientific approximations, fallbacks, incomplete convergence, and user-value
mismatches. Logging text is not a stable machine-readable interface; scripts
should use result files and exit status.

Failed calculations create collision-safe diagnostic bundles by default. These
can include logs, setup files, generated intermediates, and a limited number of
trajectory frames, so users should inspect them before sharing. Bundle creation
can be disabled when scientific inputs must not be archived.

The `gmx_MMPBSA_test` infrastructure now uses a manifest-driven example suite,
with expanded unit tests, README parity checks, and local/Colab notebook
validation. Version 1.7.0 supports Python 3.11-3.12 and documents tested
dependency ranges for AmberTools, GROMACS, NumPy, pandas, Matplotlib, SciPy,
Seaborn, `mpi4py`, ParmEd, and Rich. These ranges define the validated release
environment without claiming that every other external-program version is
categorically incompatible.

## Subsequent development after the published v1.7.0 tag

The following updates are present through commit `ac637eba`, beyond the
published tag `8d7aad33`. They extend or repair the v1.7.0 workflows rather than
representing additional methods in the tagged release:

- **Explicit-water setup:** the default `explicit_waters_mask` is now `dASA`,
  and `explicit_waters_group="automatic"` searches common solvent and water-model
  group names. Custom group names remain selectable. The dASA preparation path
  also extracts a first trajectory frame matching the selected complex topology
  before passing it to cpptraj, avoiding mismatches when other atoms, such as
  ions, have been excluded from that topology.
- **Topology preparation and logging:** elapsed progress is reported for
  structure loading, validation, residue mapping, topology conversion, and dASA
  preparation. Interactive terminals show a spinner and elapsed time; completion
  or interruption durations are logged. Routine console messages and radius
  warnings are further deduplicated across MPI ranks and topology records.
- **Structure validation:** duplicate residue IDs are counted once per chain
  rather than repeatedly scanning the residue list, removing the quadratic
  duplicate-detection step. This is an implementation improvement; no measured
  end-to-end speedup is claimed here.
- **Examples and maintenance:** the topology-backed `Comp_receptor` example now
  includes `topol.top` and its referenced `*.itp` files and is restored to tester
  coverage, resolving the missing-topology limitation of the published release.
  Example guidance, analyzer documentation, and Colab runtime checks have also
  been refreshed, and cleanup now matches multi-digit trajectory filenames.

Calculations using these updates should record their exact commit or installed
build in addition to the displayed version. Published-release validation does
not automatically validate every subsequent implementation change.

## Interpreting v1.6.5-to-v1.7.0 comparisons

A v1.7.0 calculation is not numerically equivalent to a v1.6.5 calculation
merely because it uses the same visible input file. A controlled comparison
must match the topology route, trajectory frames, GB model, radius assignment,
PB dielectric, and entropy convention. It must also account for three accepted
sources of change:

- the corrected full-ensemble IE estimator, revised C2 uncertainty calculation
  and gas constant, and new block diagnostics;
- corrected GBNSR6 parsing and frame/term merging; and
- omission of CHARMM CMAP component terms that cancel in the
  single-trajectory binding difference but not necessarily in component totals
  or multiple-trajectory calculations.

Historical v1.6.5 environments, inputs, outputs, and `_GMXMMPBSA_info` files
should therefore be preserved. Production migration should begin with a short
calculation in a copied working directory and inspection of the result schema,
warnings, logs, automatic CSV paths, and radius-provenance file.

## Current limitations and outlook

The v1.4.3-to-v1.7.0 trajectory shows a shift toward a broader and more
auditable free-energy analysis platform. Important boundaries remain:

- native-AMBER and GROMACS-derived routes do not have identical topology or
  radius semantics, and `amber_MMPBSA` does not yet have complete feature
  parity with `gmx_MMPBSA`;
- experimental CHARMM normal-mode and multiple-trajectory entropy workflows
  require cautious interpretation;
- new QH calculations are unsupported, with historical reading retained only
  for the v1.7.0 compatibility window;
- some 3D-RISM failures depend on the AmberTools/Fortran runtime rather than the
  Python package; and
- the validated dependency and example matrices are evidence for the tested
  workflows, not a guarantee for every force field, topology generator, or
  external-program combination.

## Table 1. Version timeline from v1.4.3 to v1.7.0

| Version | Date in changelog | Development theme | Main advances |
| --- | --- | --- | --- |
| v1.4.3 | 2021-05-26 | Published-era baseline | Established MM/PB(GB)SA workflows, entropy, scanning, decomposition, analyzer, correlation, and early API support. |
| v1.5.0 | 2022-02-22 | Breaking redesign | Reworked inputs and processing; C2 entropy; nonlinear PB; expanded PBSA/3D-RISM controls; analyzer redesign; input generation; stronger structure checks. |
| v1.5.1-v1.5.2 | 2022-03-10 to 2022-03-23 | System and solvent expansion | OPLS and PSF/DCD preparation workflows, QM/MM additions, bundled `.xvv`, ALPB, PC+ correction, and `sander`-based 3D-RISM. |
| v1.5.5-v1.5.7 | 2022-06-10 to 2022-09-10 | Analyzer/API and correctness | Modern API foundations, analyzer concurrency and correlation, compact results, output/entropy/decomposition fixes, progress and MPI logging. |
| v1.6.0 | 2023-02-19 | GBNSR6 and compatibility | GBNSR6 enthalpy and decomposition, named/numbered index groups, GROMACS 2023 compatibility, and entropy fixes. |
| v1.6.1-v1.6.5 | 2023-04-04 to 2026-05-22 | Maintenance and CHARMM support | PB/decomposition and analyzer fixes, additional residue handling, automatic CMAP omission during conversion, bounded dependencies, and experimental CHARMM normal mode. |
| v1.7.0 | 2026-09-11 | Reproducibility and workflow expansion | Required GROMACS topology, native `amber_MMPBSA`, explicit receptor waters, composite mutations, corrected IE and revised C2 uncertainty with block diagnostics, automatic membrane parameters, new solvent defaults, radius provenance, GBNSR6/QM/MM improvements, modernized API/logging/testing, and QH compatibility-only status. |

The tables below compare the published release endpoints; the updates in
[Subsequent development](#subsequent-development-after-the-published-v170-tag)
are recorded separately.

## Table 2. Capability comparison

| Capability | v1.4.3 baseline | Development through v1.6.5 | v1.7.0 endpoint |
| --- | --- | --- | --- |
| GB/PB | Established AmberTools GB and linear-PB workflows, including basic membrane controls. | Nonlinear PB, ALPB, expanded controls, radii sets, membrane workflow improvements, and GBNSR6. | New `igb=8`, `PBRadii=4`, and `exdi=78.5` defaults; radius provenance; GBNSR6 compaction; automatic membrane parameters and diagnostics; legacy APBS route rejected. |
| 3D-RISM | Available within the solvent-model portfolio. | Expanded variables, bundled `.xvv`, PC+, and `sander` execution. | Integrated provenance/output contract and explicit-water ST support; external AmberTools runtime caveat remains. |
| Entropy | QH, normal mode, and IE already implemented. | C2 added; IE/C2 output and analyzer handling improved; experimental CHARMM normal mode. | Corrected full-ensemble IE, stable evaluation, revised C2 uncertainty, block diagnostics, missing unconverged NMODE frames, and no new QH calculations. |
| QM/MM and scanning | QM/MMGBSA and single-residue alanine/glycine scanning. | Automatic charges, selection improvements, CHARMM/terminal fixes, and term-level mutation differences. | `PM6-DH+` default, hard SCF failure, explicit receptor waters, and multi-residue composite mutations. |
| Topology routes | GROMACS-centered conversion with structure-assisted reconstruction. | Broader Amber/CHARMM/OPLS and preparation-workflow support. | Original GROMACS complex topology required; separate native-AMBER command added. |
| Analyzer/API | Functional GUI and early dict-like API. | Redesigned high-capacity analyzer, pandas-backed API, compact results, correlation, tables, and documented performance gains. | Modernized canonical loader and runnable API example; theme/plot robustness updates. |
| Reproducibility | Basic logs, outputs, examples, and version reporting. | Better command reconstruction, checks, tester concurrency, and dependency guidance. | Radius manifest, automatic CSV naming, collision checks, MPI-safe logs, diagnostic bundles, manifest-driven tests, and bounded validated environment. |

## Table 3. Public interface changes to emphasize

| Interface area | Important changes from v1.4.3 to v1.7.0 |
| --- | --- |
| Command line | Added `--create_input` (including `pb_mem`), progress-style and error-bundle controls, improved `gmx_MMPBSA_test`, and the new `amber_MMPBSA` entry point. Output rewriting and `--clean` already existed in v1.4.3; their behavior and compatibility handling subsequently evolved. GROMACS calculations require `-cp`. |
| Input namelists | Reworked variables across GB, GBNSR6, PB, RISM, decomposition, normal mode, alanine scanning, and QM/MM; changed solvent defaults; automatic membrane parameters; composite mutations; QH and APBS restrictions. |
| Output and provenance | Compact `.mmxsa` results, automatic per-frame CSV paths, block statistics, `GMXMMPBSA_radii.json`, optional per-atom radius audit, retained membrane CSV/PNG diagnostics, and diagnostic error bundles. |
| Analyzer | New backend and result model with multi-system plots, tables, PyMOL, correlation, and frame/time controls; historical files should be read with their matching environment when reproducibility matters. |
| Python API | Canonical `GMXMMPBSA.API.load()` interface with pandas-based energy, entropy, decomposition, metadata, input, and file access. |
| Scientific comparison | v1.7.0 changes defaults, the IE estimator, and C2 uncertainty handling; cross-version comparisons require explicit model, topology, frame, radius, and uncertainty matching. |

## Suggested figures

### Figure 1. Development and workflow map

Use four linked blocks:

1. **Preparation and provenance**: GROMACS topology/trajectory inputs or native
   AMBER inputs; Amber/CHARMM/OPLS workflows; required topology and radius
   provenance.
2. **Calculation**: GB, PB, ALPB, 3D-RISM/PC+, GBNSR6, QM/MMGBSA,
   explicit-water and automatic membrane workflows, mutations, decomposition, and entropy.
3. **Results and diagnostics**: summary and per-frame outputs, block statistics,
   compact results, logs, radius manifest, membrane diagnostics, and error bundles.
4. **Analysis**: `gmx_MMPBSA_ana`, PyMOL, correlation, tables, and
   `API.load()`-based scripted workflows.

The existing `docs/assets/images/workflow.svg` can provide the visual
foundation, extended to show the separate GROMACS and native-AMBER topology
routes introduced by the v1.7.0 workflow.

### Figure 2. Scientific and reproducibility changes by release family

Use a three-stage visual comparison: v1.4.3 baseline, v1.5.x-v1.6.5 expansion,
and the v1.7.0 endpoint. Separate scientific-method additions from changes that
affect reproducibility or numerical interpretation. This avoids presenting a
corrected estimator or changed default as if it were merely another optional
feature.

## Source audit checklist

- `docs/changelog.md`: release chronology and v1.7.0 behavior changes.
- `docs/compatibility.md`: migration requirements and accepted numerical
  differences from v1.6.5.
- `docs/analyzer.md` and `docs/api.md`: analyzer/API architecture and historical
  performance evidence.
- `docs/amber_MMPBSA.md`: native-AMBER command and support boundaries.
- `docs/input_file.md` and `docs/output.md`: current defaults, method controls,
  statistics, outputs, and provenance files.
- `docs/examples/README.md`: representative systems, preparation routes, and
  validated example coverage.
- [GitHub v1.7.0 release metadata](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/releases/tag/1.7.0):
  published tag and release limitations; compare subsequent updates with the
  explicitly recorded development commit.

## Citation targets

At minimum, the manuscript should cite:

- Original `gmx_MMPBSA` publication: Valdes-Tresanco et al., *Journal of
  Chemical Theory and Computation* 2021, 17, 6281-6291,
  DOI: 10.1021/acs.jctc.1c00645.
- MMPBSA.py: Miller et al., *Journal of Chemical Theory and Computation* 2012,
  8, 3314-3321, DOI: 10.1021/ct300418h.
- Primary method papers for IE, C2, ALPB, GBNSR6, PBSA, 3D-RISM/PC+,
  QM/MMGBSA, computational alanine scanning, and membrane PB as required by the
  final manuscript.

## Concise conclusion

The advances from v1.4.3 to v1.7.0 justify an updated application note because
they change both the scientific scope and the reproducibility model of
`gmx_MMPBSA`. The package now covers more solvent models, entropy estimators,
topology routes, force-field workflows, mutations, and biomolecular systems.
Equally important, v1.7.0 requires topology provenance, records continuum-radius
assignment, corrects the IE estimator, revises C2 uncertainty, exposes block
diagnostics, strengthens
failure reporting, and defines a validated software environment. The result is
a more capable and auditable MM/PB(GB)SA platform for GROMACS, native-AMBER, and
adjacent molecular-simulation workflows.
