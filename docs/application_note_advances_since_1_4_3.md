---
template: main.html
title: Application note draft - Advances since v1.4.3
---

# gmx_MMPBSA advances since v1.4.3: structured report for a JCIM application note

!!! note "Draft status"
    This document is a structured report draft for an application note. It uses
    `gmx_MMPBSA` v1.4.3 as the baseline and summarizes advances through v1.6.5.
    Claims are grounded in the local project changelog, compatibility notes,
    analyzer/API documentation, examples, and GitHub release metadata.

## Abstract / Overview

Since the original publication of `gmx_MMPBSA` and the v1.4.3 release, the
software has evolved from a GROMACS-oriented wrapper around AmberTools
end-state free-energy calculations into a broader analysis platform for
MM/PB(GB)SA workflows. The post-v1.4.3 development cycle introduced a breaking
v1.5.x redesign of the calculation, output, analyzer, and input-processing
layers, followed by v1.6.x additions focused on expanded solvent models,
compatibility with newer molecular simulation software, richer system support,
and more reproducible analysis interfaces.

The main scientific advances include broader control of PBSA and 3D-RISM
parameters, nonlinear PB calculations, C2 entropy, the ALPB approximation,
PC+ corrections in 3D-RISM, GBNSR6 support, improved QM/MMGBSA workflows,
enhanced alanine/glycine scanning, improved decomposition analysis, membrane
protein workflows, and experimental normal-mode entropy support for CHARMM
topologies. In parallel, `gmx_MMPBSA_ana` was redesigned into a higher-capacity
graphical analysis environment with multi-system analysis, correlation tools,
table and figure export, PyMOL visualization, frame/time controls, and major
performance improvements. The Python API was also modernized into a pandas-based
data layer shared with the analyzer, improving programmatic reuse of results.

## Baseline at v1.4.3

Version 1.4.3 represented the published-era baseline used here. At that point,
`gmx_MMPBSA` already provided GROMACS users with access to AmberTools-based
MM/PBSA and MM/GBSA calculations, entropy approximations, decomposition,
alanine/glycine scanning, a graphical analyzer, examples, and early API support.
The v1.4.x series also introduced multi-system analyzer capabilities,
correlation analysis, improved decomposition selection, PyMOL visualization,
and a more mature documentation set. The post-publication changes described
below therefore represent extensions, redesigns, and stabilization beyond an
already functional v1.4.3 application.

## Expanded scientific capabilities

### Entropy calculations

The post-v1.4.3 releases substantially expanded and stabilized entropy support.
Version 1.5.0 introduced the C2 entropy method, added warnings for problematic
interaction-entropy and C2 conditions, and improved output of interaction
entropy data. Later releases corrected C2 and interaction-entropy calculations
in rewrite-output workflows, fixed Delta Delta entropy values for alanine
scanning, removed duplicated entropy items in analyzer output, and improved
the handling of `ie_segment` changes inside `gmx_MMPBSA_ana`.

The current entropy portfolio includes quasi-harmonic, normal-mode,
interaction-entropy, and C2 approximations. Normal-mode calculations remain
available through AmberTools, and v1.6.5 adds experimental normal-mode support
for CHARMM topologies. The report should clearly distinguish this CHARMM nmode
work as experimental, while presenting C2 and interaction entropy as established
post-v1.4.3 additions that are now integrated into output, analyzer, and API
workflows.

### PB, GB, ALPB, and GBNSR6

Version 1.5.0 enabled the nonlinear PB solver in `sander` and exposed the full
set of relevant PBSA options, giving users more control over Poisson-Boltzmann
calculations than v1.4.3. The same release added modified PB radii sets for
GAFF and CHARMM force fields, improving compatibility with mixed and
nonstandard molecular systems.

Version 1.5.2 added the Analytical Linearized Poisson-Boltzmann approximation
for GB calculations. Version 1.6.0 added GBNSR6 support, including enthalpy,
per-residue decomposition, and pairwise decomposition workflows. In practical
terms, this expanded `gmx_MMPBSA` beyond the original GB/PB/RISM emphasis into
a more diverse implicit-solvent platform where users can compare traditional GB
models, PB variants, ALPB-corrected GB calculations, and GBNSR6 within a common
workflow.

### 3D-RISM updates

The v1.5.0 release exposed all documented 3D-RISM variables, and v1.5.1 added
precalculated `.xvv` files. Version 1.5.2 changed 3D-RISM execution to use
`sander` instead of `rism3d.snglpnt`, and added PC+ correction support.
Together these changes improved both control and integration of RISM-based
solvation calculations inside the broader `gmx_MMPBSA` execution model.

The report should avoid claiming that 3D-RISM is universally stable across all
AmberTools/Python combinations. The current documentation notes compatibility
warnings and a known workaround involving older `gmx_MMPBSA`/Python versions for
some 3D-RISM environments. Those caveats belong in the reproducibility and
outlook sections.

### QM/MMGBSA and residue selection

The v1.5.x cycle improved QM/MMGBSA by adding new input variables, automatic
charge calculation for complex/receptor/ligand QM regions, and better residue
selection logic. The `print_res` function was also improved for decomposition
and QM workflows, reducing friction when defining residues in biologically
meaningful selections rather than low-level topology indexes.

These changes should be framed as usability and reproducibility improvements:
they do not merely add flags, but reduce user intervention in charge assignment,
selection syntax, and the consistency of residue mapping across GROMACS and
AmberTools representations.

### Alanine/glycine scanning and decomposition

Post-v1.4.3 development fixed and improved alanine scanning in several areas:
CHARMM force-field cases, terminal residues, THR-to-ALA mutation, mutant-normal
output, Delta Delta entropy reporting, and compact binary loading for alanine
scanning results. Version 1.5.0 also began reporting energy differences for
every term in alanine scanning, improving interpretability beyond total binding
energy changes.

Decomposition analysis was restructured across output, analyzer, and API layers.
The releases fixed inconsistent decomposition output, improved per-residue and
pairwise handling, added or repaired decomposition tables and plots, improved
residue selection, and later added GBNSR6 decomposition support. The application
note should present this as a key example of the project maturing from a
calculation launcher into an analysis platform that preserves detailed
component-level data.

### Membrane systems

The documentation now includes explicit examples for membrane protein
calculations, including CHARMM force-field membrane workflows. The input
documentation describes PB-based implicit membrane calculations using uniform
or heterogeneous slab-like dielectric profiles. In the application note, this
should be presented as an important expansion of biological scope, especially
because membrane protein binding calculations remain challenging for routine
end-state free-energy analysis.

## Broader system and file support

### Force fields and preparation workflows

Compared with v1.4.3, current `gmx_MMPBSA` supports a broader set of preparation
workflows and force-field conventions. The post-v1.4.3 releases added support
for CHARMM-GUI files generated for Amber force fields, improved CHARMM
topology handling, introduced OPLS force-field support, and added modified PB
radii sets for GAFF and CHARMM systems. The current documentation presents
Amber, OPLS, and CHARMM as supported force-field families for users preparing
systems in GROMACS-compatible workflows.

Residue and atom handling also improved. Notable changes include support for
additional histidine variants, stricter detection of lone-pair atoms, improved
terminal hydrogen handling, checks that receptor and ligand masks do not
overlap, and structure/topology consistency checks. These features should be
discussed as practical reliability improvements for real-world biomolecular
systems, where residue naming, terminal patches, alternate protonation states,
and topology/coordinate mismatches are common sources of failure.

### psf/dcd and AMBER-native inputs

The documentation now includes tutorials for `psf`/`dcd` files used by NAMD,
OpenMM, GENESIS, and related workflows. These examples extend the user base
beyond standard GROMACS file sets while preserving the same end-state
free-energy analysis logic.

Recent development also added `amber_MMPBSA`, an independent command-line module
within the package for AMBER-native topology, coordinate, trajectory, and mask
inputs. This module removes the GROMACS topology-conversion step for users who
already have AMBER files while keeping the same calculation and analysis
features where supported. This is a notable architectural broadening: the
package is no longer only a GROMACS bridge, but also provides a shared analysis
stack around both GROMACS-derived and AMBER-native workflows.

## Analyzer and user experience

### Redesign of gmx_MMPBSA_ana

The v1.5.x analyzer redesign is one of the most visible changes since v1.4.3.
The current `gmx_MMPBSA_ana` keeps the general goal of the original analyzer but
uses a different backend API and is not compatible with old result files from
previous versions. The redesign added:

- multi-system loading from info files, folders, and recursive folder searches;
- selection of systems, calculation types, subsystems, components, mutants, and
  decomposition data;
- line plots, bar plots, heatmaps, PyMOL residue-energy visualization, summary
  tables, and output-file viewers;
- chart property controls for fonts, palettes, themes, figure sizes, dpi,
  formats, labels, rotations, and plot-specific settings;
- frame-range, interval, and time-conversion controls;
- table views that can be copied into spreadsheet tools;
- per-system settings and better start-dialog controls.

This should be framed as a shift from viewing individual outputs to managing
large, multi-calculation result workspaces.

### Correlation workflows

Correlation analysis was present in v1.4.x, but post-v1.4.3 development made it
more usable and better integrated. The analyzer gained normal correlations using
Delta G, mutant correlations using Delta Delta G, energy-value tables,
interactive system selection, an experimental Ki editor, regression plots, and
chart options for regression/scatter/distribution views. The input file and
system-selection workflows now expose fields such as `sys_name`, `exp_ki`, and
`temperature`, making multi-system comparison more reproducible.

### Performance and large result sets

The v1.5.5 analyzer documentation reports major performance improvements after
the v1.5.2 performance problems. The documented benchmark shows that systems
that failed or required hundreds to thousands of seconds in v1.5.2-era behavior
could be loaded in seconds after the redesign. The documented examples include
up to 56x improvement for one energy system with 200,000 frames and up to 265x
for four comparable energy systems. For energy plus pairwise decomposition, the
reported improvements were 21x and 69x for one and four systems, respectively.

The report should quote these numbers as documented analyzer benchmarks, not as
new independent performance measurements. The underlying implementation changes
include multiprocessing for reading systems, multithreading for data processing,
optimized output parsing, pandas-backed storage, removal of redundant data and
pop-ups, nonblocking GUI processing, waiting indicators, and an option to store
temporary data on disk instead of memory.

## Developer architecture and reproducibility

### Input and calculation layer

Version 1.5.0 was intentionally incompatible with previous versions because of
major changes in calculation and processing modules. This release removed
deprecated variables from v1.4.x, reorganized input variables, exposed more
calculation controls, changed `EnergyVector` into an `ndarray` subclass, and
made output/statistical handling more consistent.

The input interface also became more self-documenting. The `--create_input`
command can generate default input files for calculation classes such as GB,
PB, RISM, alanine scanning, decomposition, nmode, GBNSR6, and all combined
templates. The command-line interface also supports `--rewrite-output` for
reparsing previous calculations and rewriting outputs, and `--clean` for
removing temporary files.

### Data model, compact results, and Python API

The current API is the same data layer used by `gmx_MMPBSA_ana`. It loads
`_GMXMMPBSA_info` files or compact `.mmxsa` result files and returns pandas-based
data structures for custom analysis. The modern loader is:

```python
from GMXMMPBSA import API

api = API.load("COMPACT_MMXSA_RESULTS.mmxsa")
```

The API exposes metadata through `get_info()`, input namelists through
`get_input()`, file metadata through `get_files()`, energy data through
`get_energy()`, entropy data through `get_entropy()`, and decomposition data
through `get_decomp_energy()`. The older `load_gmxmmpbsa_info()` function remains
available but now returns a loaded `MMPBSA_API` object and emits a deprecation
warning.

This API modernization matters for a JCIM application note because it opens
reproducible post-processing pipelines beyond the GUI: users can extract
per-frame energy terms, summaries, entropy values, correlation-ready data, and
decomposition tables directly into Python workflows.

### Testing, logging, and installation

The `gmx_MMPBSA_test` tool has been improved repeatedly since v1.4.3, including
parallel/concurrent example execution, lazy executable discovery, reuse of test
folders, better frame counting, and clearer failure status reporting. The
examples folder was moved to the repository root, and the documentation now
organizes examples by system class, force-field family, file type, and analysis
type.

Logging and debugging also improved. The releases added reproducible command
lines to log files, version information, debug-level input-file logging, better
external-program discovery, better error parsing, and explicit checks for
structure/topology mismatches. These changes are less visible than new solvent
models but are important for reproducibility, support, and long-running
high-throughput workflows.

The development environment also changed. Version 1.6.5 updated dependency
guidance for Python, AmberTools, ParmEd, and GROMACS, added a conda/pip
installation script, and restricts the package to Python 3.11 in `setup.py`.
The application note should present this as active maintenance around a complex
scientific software stack, while also acknowledging that version coupling across
AmberTools, GROMACS, Python, and Fortran runtimes remains a practical concern.

## Outlook and current limitations

The post-v1.4.3 trajectory of `gmx_MMPBSA` shows a clear movement toward a
broader, more integrated free-energy analysis platform. Near-term areas to
highlight cautiously include:

- GROMACS 2026 support, which has a documented workaround but is not yet fully
  automated in the current workflow;
- experimental CHARMM nmode support, which should be presented as emerging
  rather than mature;
- continued stabilization of 3D-RISM compatibility across AmberTools and Python
  versions;
- continued consolidation of the Python API as the common backend for scripted
  and graphical analysis;
- documentation and examples as key infrastructure for scientific
  reproducibility.

## Table 1. Version timeline since v1.4.3

| Version | Date in changelog | Development theme | Main advances |
| --- | --- | --- | --- |
| v1.4.3 | 2021-05-26 | Baseline | Published-era baseline; correlation p-values; QM/MMGBSA and LPH atom tutorials; improved force-field parsing and logging. |
| v1.5.0 | 2022-02-22 | Breaking redesign | New input format and variables; C2 entropy; nonlinear PB; expanded PBSA and 3D-RISM variables; CHARMM-GUI Amber workflows; experimental H5 output; analyzer redesign; `--create_input`; improved structure checks. |
| v1.5.0.1-v1.5.0.3 | 2022-02-22 to 2022-02-26 | Early stabilization | Group-name selection, analyzer launch fixes, nmode/statistical output fixes, and updated tutorials. |
| v1.5.1 | 2022-03-10 | Expanded system support | OPLS support, psf/dcd tutorials, PyQt6 support, QM/MM variable additions, improved `print_res`, and bundled `.xvv` files. |
| v1.5.2 | 2022-03-23 | Solvent-model expansion | ALPB approximation, PC+ 3D-RISM correction, `sander`-based 3D-RISM, and analyzer tooltip/chart fixes. |
| v1.5.5 | 2022-06-10 | Analyzer/API performance | API methods for enthalpy, entropy, binding, decomposition, and analyzer data; multiprocessing/threading; correlation workflow; large performance gains; improved IE/C2 handling. |
| v1.5.6 | 2022-07-06 | Output correctness | Fixed decomposition output, IE/C2 rewrite-output failures, Delta Delta entropy values, compact binary loading, mutant-normal analyzer data, and PyMOL energy consistency. |
| v1.5.7 | 2022-09-10 | Runtime usability | Progress bars, API/logging verbosity, MPI logging, compact result fixes, QH/decomposition fixes, and automatic `mpi4py`/fake-MPI selection. |
| v1.6.0 | 2023-02-19 | New GB model and compatibility | GBNSR6 implementation with enthalpy and decomposition, group names/numbers in index selection, GROMACS 2023 compatibility, entropy fixes, and default `inp = 1`. |
| v1.6.1-v1.6.4 | 2023-04-04 to 2024-04-11 | Maintenance and compatibility | PB decomposition fixes, error parsing, minor fixes, and compatibility updates. |
| v1.6.5 | 2026-05-22 | Maintenance plus experimental CHARMM nmode | Experimental nmode for CHARMM topologies, additional histidine variants, stricter LP atom detection, entropy/analyzer fixes, terminal hydrogen handling, Python 3.11 requirement, updated dependencies and installation script. |

!!! note
    GitHub releases list v1.6.5 as released on 2026-05-23, while the local
    changelog lists 2026-05-22. Use one date consistently in the manuscript or
    state the discrepancy in internal notes only.

## Table 2. Capability matrix

| Capability | v1.4.3 baseline | Current state through v1.6.5 |
| --- | --- | --- |
| GB calculations | Established MM/GBSA workflows through AmberTools. | Added ALPB, improved dielectric/radii controls, GBNSR6, QM/MMGBSA improvements, and better output/API integration. |
| PB calculations | Linear PB workflows available. | Nonlinear PB solver enabled; expanded PBSA controls; membrane PB workflows documented; PB decomposition fixes. |
| 3D-RISM | Available as part of solvent-model portfolio. | Full variable exposure, bundled `.xvv` files, PC+ correction, and `sander`-based execution; compatibility caveats remain. |
| GBNSR6 | Not available. | Implemented in v1.6.0 with enthalpy and per-residue/pairwise decomposition. |
| Entropy | QH, nmode, and interaction entropy available or emerging. | C2 entropy added; IE/C2 output and analyzer integration improved; experimental CHARMM nmode added. |
| Alanine/glycine scanning | Available, including v1.4.x improvements. | Better CHARMM and terminal-residue handling, term-level energy differences, mutant-normal data, and entropy correction fixes. |
| Decomposition | Available with analyzer support. | Restructured output/API/analyzer handling; improved residue selection; GBNSR6 decomposition; fixes for inconsistent outputs. |
| Correlation | Analyzer correlation and p-values available. | More integrated multi-system correlation, Ki editing, regression plots, and correlation-ready API summaries. |
| Analyzer | Functional GUI with v1.4.x multi-system improvements. | Redesigned high-capacity GUI with chart customization, tables, PyMOL, frame/time controls, compact results, and documented performance gains. |
| Python API | Earlier dict-like API documented. | Modern `API.load()` interface backed by pandas structures and shared with `gmx_MMPBSA_ana`. |
| File/force-field support | GROMACS-centered workflows with Amber/CHARMM improvements. | Broader Amber, CHARMM, OPLS, psf/dcd, CHARMM-GUI, and AMBER-native `amber_MMPBSA` workflows. |

## Table 3. Public interface changes to emphasize

| Interface area | Important post-v1.4.3 changes |
| --- | --- |
| Command line | `--create_input`, `--rewrite-output`, `--clean`, improved `gmx_MMPBSA_test`, `amber_MMPBSA`, group names or numbers for group selection, and expanded input-generation modes. |
| Input namelists | New/reworked variables in `&general`, `&gb`, `&gbnsr6`, `&pb`, `&rism`, `&decomp`, `&nmode`, and QM/MMGBSA controls; removed deprecated v1.4.x variables. |
| Output files | Continued use of `FINAL_RESULTS_MMPBSA.dat` and `FINAL_DECOMP_MMPBSA.dat`; CSV export for energy and decomposition terms; compact `.mmxsa` result files for portable analysis. |
| Analyzer | v1.5.x result files are not compatible with v1.4.3 analyzer files; new analyzer backend uses the modern API and supports multi-system, multi-plot, table, PyMOL, and correlation workflows. |
| Python API | New canonical `GMXMMPBSA.API.load()` loader; pandas-based data; metadata, input, file, energy, entropy, decomposition, and correlation accessors. |
| Reproducibility | Improved logging, version reporting, command reconstruction, structure/topology consistency checks, test-runner controls, environment documentation, and dependency pinning. |

## Suggested figures

### Figure 1. Post-publication development map

Use a workflow diagram with four grouped blocks:

1. **Preparation**: GROMACS, psf/dcd, AMBER-native inputs; Amber/CHARMM/OPLS;
   topology conversion or direct AMBER topology use; structure consistency
   checks.
2. **Calculation**: GB, PB, nonlinear PB, ALPB, 3D-RISM, PC+, GBNSR6,
   QM/MMGBSA, membrane PB, alanine/glycine scanning, decomposition, entropy.
3. **Analysis**: redesigned `gmx_MMPBSA_ana`, plots, tables, PyMOL,
   correlation, frame/time controls, compact result loading.
4. **Programmatic access**: `API.load()`, pandas DataFrames, metadata,
   summaries, decomposition extraction, reproducible downstream scripts.

The existing workflow figure in `docs/assets/images/workflow.svg` can be used
as the visual foundation, with a new overlay or expanded version for the
post-publication architecture.

### Figure 2. Analyzer evolution

Use the existing analyzer overview image
`docs/assets/images/gmx_mmpbsa_ana_overview.png` and pair it with a small
benchmark panel based on the documented v1.5.5 analyzer performance table:

| Scenario | Previous behavior | Redesigned analyzer | Documented improvement |
| --- | --- | --- | --- |
| 1 energy system, 200,000 frames | v1.5.2 failed; v1.5.2+20 took 960 s | 17 s | 56x |
| 4 energy systems, 200,000 frames | v1.5.2 failed; v1.5.2+20 took 4515 s | 17 s | 265x |
| 1 energy + pairwise decomposition system, 11 frames | v1.5.2 failed; v1.5.2+20 took 192 s | 9 s | 21x |
| 4 energy + pairwise decomposition systems, 11 frames | v1.5.2 failed; v1.5.2+20 took 960 s | 14 s | 69x |

## Source audit checklist

- `docs/changelog.md`: release chronology, additions, fixes, and behavior
  changes.
- `docs/compatibility.md`: v1.4.3 to v1.5.0 breaking changes and input
  variable migration.
- `docs/analyzer.md`: analyzer redesign, performance benchmarks, plots, tables,
  PyMOL, and user workflow.
- `docs/api.md`: current Python API and deprecation of the historical v1.4.x
  dict-like interface.
- `docs/examples/README.md`: supported systems, force fields, file types, and
  calculation/tutorial coverage.
- `docs/amber_MMPBSA.md`: AMBER-native command-line workflow.
- `docs/input_file.md`: detailed input namelists and model-specific settings.
- GitHub releases: external verification of current latest release metadata.

## Citation targets

At minimum, the manuscript should cite:

- Original `gmx_MMPBSA` publication: Valdes-Tresanco et al.,
  *Journal of Chemical Theory and Computation* 2021, 17, 6281-6291,
  DOI: 10.1021/acs.jctc.1c00645.
- MMPBSA.py: Miller et al., *Journal of Chemical Theory and Computation* 2012,
  8, 3314-3321, DOI: 10.1021/ct300418h.
- Method papers for interaction entropy, C2 entropy, ALPB, GBNSR6, PBSA,
  3D-RISM/PC+, QM/MMGBSA, computational alanine scanning, and membrane PB as
  required by the final prose.

## Concise conclusion

The advances after v1.4.3 justify an application note because they affect both
the scientific scope and the engineering model of `gmx_MMPBSA`. Scientifically,
the package now covers more solvent models, entropy corrections, force-field
families, file formats, and complex biomolecular systems. Architecturally, the
software now exposes a shared analyzer/API data layer, compact result files,
improved testing, richer logging, and stronger compatibility guidance. The
result is a more complete and reusable MM/PB(GB)SA platform for GROMACS,
AMBER-native, and adjacent molecular simulation workflows.
