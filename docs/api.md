---
template: main.html
title: Python API
---

# Python API

The `GMXMMPBSA.API` module provides programmatic access to parsed `gmx_MMPBSA`
results. The current API is the same data layer used by `gmx_MMPBSA_ana`: it
loads `_GMXMMPBSA_info` files or compact `COMPACT_MMXSA_RESULTS.mmxsa` files and
returns pandas-based data structures for custom analysis.

!!! note
    Older documentation described a dict-like `mmpbsa_data` object for
    `gmx_MMPBSA` 1.4.x. That interface is no longer the canonical API. The
    historical `load_gmxmmpbsa_info()` function is still available, but it now
    returns a loaded `MMPBSA_API` object and emits a deprecation warning.

!!! warning
    When loading an `_GMXMMPBSA_info` file, the topology, output, and related
    files referenced by that info file must still be available. The compact
    `.mmxsa` file is usually the most portable input for API use.

## Loading Results

Use `load()` for new code:

```python
from GMXMMPBSA import API

api = API.load("COMPACT_MMXSA_RESULTS.mmxsa")
```

The same loader accepts an info file:

```python
from GMXMMPBSA import API

api = API.load("_GMXMMPBSA_info")
```

If you need frame indexes expressed as simulation time, pass the same arguments
accepted by `MMPBSA_API.setting_time()`:

```python
from GMXMMPBSA import API

api = API.load(
    "COMPACT_MMXSA_RESULTS.mmxsa",
    time={"timestep": 10, "timeunit": "ps"},
)
```

The lower-level class remains public:

```python
from GMXMMPBSA.API import MMPBSA_API

api = MMPBSA_API()
api.setting_time(timestep=10, timeunit="ps")
api.load_file("COMPACT_MMXSA_RESULTS.mmxsa")
```

For a complete runnable example that loads a real compact result and writes a
CSV file, see the [Python API extraction example](examples/API/README.md).

## Metadata

The loaded object exposes calculation metadata:

```python
info = api.get_info()
input_options = api.get_input()
files = api.get_files()

print(info["numframes"])
print(input_options["general"]["temperature"])
print(files.output_file)
```

`get_info()` returns a dictionary with values such as frame counts, output text,
mutation information, and the complex PDB text. `get_input()` returns the parsed
input namelists. `get_files()` returns the file namespace stored in the result.

## Energy Data

`get_energy()` returns a dictionary with four main entries:

* `map`: available result keys after filtering.
* `data`: pandas DataFrames containing per-frame values plus `Average`, `SD`,
  and `SEM` rows.
* `summary`: pandas DataFrames containing only `Average`, `SD`, and `SEM`.
* `correlation`: compact energy summary used by correlation analysis.

Example:

```python
energy = api.get_energy()

gb_delta = energy["data"]["normal"]["gb"]["delta"]
total = gb_delta["TOTAL"]

print(total.head())
print(energy["summary"]["normal"]["gb"]["delta"]["TOTAL"])
```

You can restrict the returned data:

```python
energy = api.get_energy(
    etype=("normal",),
    model=("gb",),
    mol=("complex", "delta"),
    term=("VDWAALS", "EEL", "TOTAL"),
    startframe=1,
    endframe=100,
    interval=2,
)
```

Common `etype` values are `normal`, `mutant`, and `mutant-normal`. Common model
values include `gb`, `pb`, `rism std`, `rism gf`, `rism pcplus`, and `gbnsr6`,
depending on the calculation.

## Key Reference

The keys available through the API depend on the calculation that produced the
result file. For example, `pb` keys are only present when PB was run, `mutant`
keys are only present for alanine scanning, and decomposition keys are only
present when decomposition was enabled.

Use the returned `map` object to inspect what is actually available:

```python
api = API.load("examples/API/COMPACT_MMXSA_RESULTS.mmxsa")

energy = api.get_energy()
entropy = api.get_entropy()

print(energy["map"])
print(entropy["map"])

# For a decomposition result:
# decomp = api.get_decomp_energy()
# print(decomp["map"])
```

### Result Types

| Key | Meaning |
| --- | --- |
| `normal` | Results for the original system. |
| `mutant` | Results for the mutated system when alanine or glycine scanning was run. |
| `mutant-normal` | Difference between mutant and normal results. |
| `decomp_normal` | Decomposition data for the original system. Used as an `etype` argument to `get_decomp_energy()`. |
| `decomp_mutant` | Decomposition data for the mutated system. Used as an `etype` argument to `get_decomp_energy()`. |

### Molecular Components

| Key | Meaning |
| --- | --- |
| `complex` | Complex energy or entropy. Present for stability and binding calculations. |
| `receptor` | Receptor contribution. Present for binding calculations. |
| `ligand` | Ligand contribution. Present for binding calculations. |
| `delta` | Binding difference, usually `complex - receptor - ligand`. |

### Energy Models

| Key | Meaning |
| --- | --- |
| `gb` | Generalized Born energy. |
| `gbnsr6` | GBNSR6 energy. |
| `pb` | Poisson-Boltzmann energy. |
| `rism std` | Standard 3D-RISM energy. |
| `rism gf` | Gaussian fluctuation 3D-RISM energy. |
| `rism pcplus` | PC+ corrected 3D-RISM energy. |

### Energy Terms

| Key | Meaning |
| --- | --- |
| `BOND` | Bond energy. |
| `ANGLE` | Angle energy. |
| `DIHED` | Dihedral energy. |
| `VDWAALS` | van der Waals energy. |
| `EEL` | Electrostatic energy. |
| `1-4 VDW` | 1-4 van der Waals energy. |
| `1-4 EEL` | 1-4 electrostatic energy. |
| `UB` | Urey-Bradley term for CHARMM-style topologies. |
| `IMP` | Improper dihedral term for CHARMM-style topologies. |
| `CMAP` | Correction map term for CHARMM-style topologies. |
| `ESCF` | Semiempirical correction term for supported QM/MM calculations. |
| `EGB` | GB polar solvation energy. |
| `ESURF` | GB non-polar solvation energy. |
| `EPB` | PB polar solvation energy. |
| `ENPOLAR` | PB non-polar cavity/dispersion-style contribution, depending on PB settings. |
| `EDISPER` | PB dispersion contribution. |
| `POLAR SOLV` | 3D-RISM polar solvation energy. |
| `APOLAR SOLV` | 3D-RISM apolar solvation energy. |
| `ERISM` | Combined 3D-RISM solvation energy when polar/apolar decomposition is not available. |
| `GGAS` | Composite gas-phase energy. |
| `GSOLV` | Composite solvation energy. |
| `TOTAL` | Total energy. |

### Entropy Models And Terms

| Key | Meaning |
| --- | --- |
| `nmode` | Normal mode entropy. |
| `qh` | Quasi-harmonic entropy. |
| `ie` | Interaction entropy. |
| `c2` | C2 entropy. |
| `TRANSLATIONAL` | Translational entropy term for NMODE/QH outputs. |
| `ROTATIONAL` | Rotational entropy term for NMODE/QH outputs. |
| `VIBRATIONAL` | Vibrational entropy term for NMODE/QH outputs. |
| `TOTAL` | Total NMODE/QH entropy. |
| `AccIntEnergy` | Accumulated interaction energy series used in IE output. |
| `ie` | Interaction entropy values in IE DataFrames. |
| `sigma` | IE/C2 sigma value. |
| `c2` | C2 entropy value in C2 DataFrames. |

### Decomposition Keys

| Key | Meaning |
| --- | --- |
| `TDC` | Total decomposition contribution. |
| `SDC` | Sidechain decomposition contribution. |
| `BDC` | Backbone decomposition contribution. |
| Residue key | Per-residue identifier from the parsed decomposition output. |
| Pair residue key | Partner residue identifier for pairwise decomposition. |
| `int` | Internal contribution. |
| `vdw` | van der Waals contribution. |
| `eel` | Electrostatic contribution. |
| `pol` | Polar solvation contribution. |
| `sas` | Non-polar solvation contribution. |
| `tot` | Total decomposition contribution for the residue or pair. |

### Return Object Keys

| Key | Returned by | Meaning |
| --- | --- | --- |
| `map` | `get_energy()`, `get_entropy()`, `get_decomp_energy()`, `get_binding()` | Nested dictionary/list of available keys after filtering. |
| `data` | `get_energy()`, `get_entropy()`, `get_decomp_energy()`, `get_binding()` | Main pandas result objects. |
| `summary` | `get_energy()`, `get_entropy()` | `Average`, `SD`, and `SEM` summaries. |
| `correlation` | `get_energy()`, `get_binding()`, optional `get_ana_data()` | Compact values prepared for correlation analysis. |
| `keys` | `get_ana_data()` | Analyzer chart-data keys. |

## Entropy And Binding

`get_entropy()` combines available entropy methods:

```python
entropy = api.get_entropy(model=("ie", "c2", "nmode", "qh"))
print(entropy["map"])
```

Specific helpers are also available:

```python
ie = api.get_ie_entropy(ie_segment=25)
c2 = api.get_c2_entropy()
nmode = api.get_nmode_entropy()
qh = api.get_qh_entropy()
```

To compute binding tables that combine enthalpy and entropy summaries:

```python
energy = api.get_energy()
entropy = api.get_entropy()

binding = api.get_binding(
    energy_summary=energy["summary"],
    entropy_summary=entropy["summary"],
)

print(binding["data"]["normal"]["gb"]["ie"])
```

## Decomposition Data

`get_decomp_energy()` returns decomposition data when decomposition was enabled
in the original calculation:

```python
decomp = api.get_decomp_energy(
    model=("gb",),
    mol=("complex", "delta"),
    contribution=("TDC", "SDC"),
    res_threshold=0.5,
)

print(decomp["map"])
print(decomp["data"]["normal"]["gb"].head())
```

For pairwise decomposition, use `res1` and `res2` to restrict residue pairs.
For per-residue decomposition, use `res1` and `term`.

## Analyzer Data

`get_ana_data()` prepares the same hierarchical data consumed by
`gmx_MMPBSA_ana`. Public API use defaults to in-memory pandas objects:

```python
ana_data = api.get_ana_data(
    energy_options={"etype": ("normal",)},
    entropy_options={"etype": ("normal",)},
    decomp_options={"res_threshold": 0.5},
    performance_options={"energy_memory": True, "decomp_memory": True},
    correlation=True,
)
```

Disk-backed analyzer data uses parquet files. If you request disk-backed data,
install either `pyarrow` or `fastparquet`:

```python
api = API.load("COMPACT_MMXSA_RESULTS.mmxsa", data_on_disk=True)
ana_data = api.get_ana_data(
    energy_options={},
    performance_options={"energy_memory": False, "decomp_memory": True},
)
```

## Compatibility Loader

`load_gmxmmpbsa_info()` is kept for older imports:

```python
from GMXMMPBSA import API

api = API.load_gmxmmpbsa_info("_GMXMMPBSA_info")
energy = api.get_energy()
```

This compatibility function intentionally returns the modern `MMPBSA_API`
object. Scripts written for the old dict-like API should be migrated to
`get_energy()["data"]`, `get_entropy()["data"]`, and `get_decomp_energy()["data"]`.

## Validation Status

This API modernization is intentionally staged on a separate branch. The current
minimum validation target is:

```bash
python3 -m py_compile GMXMMPBSA/API.py
```

Import and fixture-loading smoke tests should be run in the supported Python
3.11 environment with the project dependencies installed. Full API regression
coverage is deferred until the public interface is reviewed.
