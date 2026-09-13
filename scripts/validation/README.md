# Calculation validation scripts

These scripts run calculations in physical copies of `examples/` under an
external timestamped results directory. They do not write generated outputs to
the repository.

Use the explicit environment interpreters. The full current calculation pass
is:

```bash
python scripts/validation/run_calculation_matrix.py \
  --python /path/to/gmxMMPBSA_dev/bin/python \
  --label current-all
```

The default selector `101` means all 23 currently runnable manifest examples.
The manifest defines 23 runnable examples, including `Comp_receptor` (case 9),
whose required `topol.top` topology is now shipped. Selectors can also
be individual IDs such as `4 12 23` or aliases such as `gbnsr6`.

The legacy comparison runner uses a conservative set of examples shared by the
1.6.5-era and current calculation interfaces. It archives tag `1.6.5` outside
the checkout unless `--baseline-source` is supplied:

```bash
python scripts/validation/compare_legacy_matrix.py \
  --current-python /path/to/gmxMMPBSA_dev/bin/python \
  --baseline-python /path/to/gmxMMPBSA_165/bin/python
```

For parity runs, copied inputs explicitly set the legacy `PBRadii=3`,
`igb=5`, and PB `exdi=80` values. This is not applied to canonical examples.
CSV numeric values are compared with `atol=1e-3` and `rtol=1e-6` by default.
Known compatibility boundaries are recorded in
`legacy_comparison_policy.json`.  CMAP-sensitive STP cases compare the
`Delta Energy Terms` section while retaining raw component differences in the
report.  Known baseline gaps and intentional behavior changes remain non-zero
in strict mode; use `--allow-known-exceptions` only when accepting the
documented exceptions for an exploratory extended run.

For the extended cases exercised during the first calculation pass:

```bash
python scripts/validation/compare_legacy_matrix.py \
  6 8 10 17 18 25 26 \
  --current-python /path/to/gmxMMPBSA_dev/bin/python \
  --baseline-python /path/to/gmxMMPBSA_165/bin/python \
  --allow-known-exceptions
```

## Calculation combination matrix

`calculation_combination_matrix.json` defines the next current-version
calculation gate. It uses pairwise and one-factor-at-a-time rows so a failure
can be attributed to a solver, setting, trajectory mode, analysis mode, or
advanced workflow. It deliberately does not form a Cartesian product of every
input variable.

The matrix currently contains 44 planned rows in three executable tiers:
`core`, `analysis`, and `advanced`. Use the planning command to inspect rows
without copying examples or launching calculations:

```bash
python scripts/validation/plan_combination_matrix.py core
python scripts/validation/plan_combination_matrix.py full --json
```

The planner is non-executing by design. The executor creates one physical
`/tmp` copy per row, applies only that row's overrides, records the final input
and command, and writes logs/results outside Git. For the complete current
pass:

```bash
python scripts/validation/run_combination_matrix.py full \
  --python /path/to/gmxMMPBSA_dev/bin/python \
  --label combination-full-current --execute
```

Rows marked `legacy_candidate` can be rerun with `--legacy-settings` against
both environments. Compare their isolated `run.json` files with:

```bash
python scripts/validation/compare_combination_matrix.py \
  /tmp/current/run.json /tmp/baseline/run.json \
  --output /tmp/current/combination-comparison.json
```

The comparison uses `atol=1e-3` and `rtol=1e-6`; a current pass with an old
baseline failure is reported as `BASELINE-UNSUPPORTED`, not numerical parity.
Negative/error,
GUI/analyzer/API, MPI/concurrency, cleanup/restart, and packaging cases remain
separate round-two suites.

The first calculation round intentionally excluded analyzer/GUI, API,
diagnostic bundles, negative cases, cleanup/restart behavior, and concurrency
stress. The API and analyzer sections below are the first staged parts of the
second round; the remaining surfaces still follow separately.

## Python API validation

The API runner consumes existing outputs without launching calculations. Give it
one or more external result directories; it discovers both `_GMXMMPBSA_info`
and `COMPACT_MMXSA_RESULTS.mmxsa` files and writes its report outside Git:

```bash
python scripts/validation/run_api_validation.py \
  /tmp/gmx_MMPBSA-validation/current/examples \
  /tmp/gmx_MMPBSA-validation/baseline/examples \
  --output /tmp/gmx_MMPBSA-validation/api-validation.json
```

It checks the public loader and accessors, legacy loader compatibility,
energy/entropy/binding/decomposition routing, summary statistics, and
in-memory analyzer preparation. It does not launch calculations or require a
parquet engine. Disk-backed analyzer mode requires `pyarrow` or `fastparquet`
and should be validated separately in an environment that provides one.

To compare paired API data from the current legacy-settings run and the
1.6.5 baseline, use:

```bash
python scripts/validation/compare_api_results.py \
  --current-root /tmp/gmx_MMPBSA-validation/current/examples \
  --baseline-root /tmp/gmx_MMPBSA-validation/baseline/examples \
  --output /tmp/gmx_MMPBSA-validation/api-comparison.json
```

The comparison includes raw numeric cells and derived summaries. Its explicit
policy classifies documented IE/C2, GBNSR6, and component-level CHARMM-CMAP
changes as `EXPECTED-DIFFERENCE`; any other mismatch remains a failure.

## Analyzer intake validation

The analyzer intake runner validates the real file-discovery code and Qt
initialization dialog in offscreen mode. It does not enter the GUI event loop:

```bash
QT_QPA_PLATFORM=offscreen python scripts/validation/run_analyzer_validation.py \
  /tmp/gmx_MMPBSA-validation/current/examples/Protein_ligand/ST \
  --output /tmp/gmx_MMPBSA-validation/analyzer-validation.json
```

Use `--recursive` for the analyzer's documented one-directory-level recursive
search, or pass individual result files when validating a deeper tree.

## MPI, concurrency, and cleanup validation

The MPI prerequisite report checks the selected environment's launcher,
`mpi4py` world size, and GROMACS build mode. It does not claim that the
application itself passed MPI:

```bash
python scripts/validation/run_mpi_validation.py \
  --python /path/to/gmxMMPBSA_dev/bin/python \
  --output /tmp/gmx_MMPBSA-validation/mpi.json
```

For application-level MPI evidence, run a calculation row with multiple
ranks. The following serial/MPI pair was used for the current gate:

```bash
python scripts/validation/run_combination_matrix.py GB-01 \
  --python /path/to/gmxMMPBSA_dev/bin/python --ranks 2 --execute
```

Independent calculation rows can be launched concurrently; every row receives
its own physical example copy and log:

```bash
python scripts/validation/run_concurrency_validation.py GB-01 GB-02 \
  --python /path/to/gmxMMPBSA_dev/bin/python --execute
```

Cleanup validation creates only synthetic files under the external results
root. It checks that minimal cleanup preserves metadata, compact results, and
membrane diagnostics, while full cleanup removes generated results and keeps
the diagnostics:

```bash
python scripts/validation/run_cleanup_validation.py \
  --output /tmp/gmx_MMPBSA-validation/cleanup.json
```

## Negative/error validation

The negative runner executes an invalid GB setting in two isolated copies. It
expects a diagnostic error bundle by default and no bundle when
`--no-error-bundle` is supplied:

```bash
python scripts/validation/run_negative_validation.py \
  --python /path/to/gmxMMPBSA_dev/bin/python \
  --output /tmp/gmx_MMPBSA-validation/negative.json --execute
```

## Disk-backed analyzer validation

This check accepts either a successful parquet-backed run or the intentional
dependency guard when neither parquet engine is installed:

```bash
python scripts/validation/run_disk_analyzer_validation.py \
  /tmp/gmx_MMPBSA-validation/current/COMPACT_MMXSA_RESULTS.mmxsa \
  --output /tmp/gmx_MMPBSA-validation/disk-analyzer.json
```

## Packaging/install validation

The packaging runner builds an sdist and wheel from a physical temporary copy,
installs the wheel without dependencies into another external directory, and
checks imports plus bundled data files:

```bash
python scripts/validation/run_packaging_validation.py \
  --source /path/to/gmx_MMPBSA \
  --python /path/to/gmxMMPBSA_dev/bin/python \
  --output /tmp/gmx_MMPBSA-validation/packaging.json
```

## Coverage vs Amber kcal/mol goldens

Unit tests cover IE/C2 math and energy **parsing/aggregation** from fixture
mdouts (`tests/test_energy_parse_golden.py`). They do **not** replace an
end-to-end AmberTools regression that regenerates GB/PB totals from topologies
and trajectories.

What exists today:

| Layer | What it checks |
| --- | --- |
| `tests/test_calculation.py` | IE/C2 equations and diagnostics |
| `tests/test_energy_parse_golden.py` | GB mdout parse → ΔTOTAL/ΔGGAS/ΔGSOLV means (fixture) |
| `scripts/validation/run_calculation_matrix.py` | Full example runs under current env |
| `scripts/validation/compare_legacy_matrix.py` | Current vs archived 1.6.5 (legacy defaults) |
| `scripts/validation/run_combination_matrix.py` | Isolated current calculation-combination rows |
| `scripts/validation/compare_combination_matrix.py` | CSV comparison for matched combination rows |

To add a true Amber kcal/mol golden later: archive a short AmberTools `sander`
(or `MMPBSA.py`) reference CSV for one example, pin AmberTools version in the
fixture README, and compare `FINAL_RESULTS_MMPBSA.csv` component means within a
stated absolute tolerance. Keep those binaries out of the default unit-test
path so CI stays dependency-light.
