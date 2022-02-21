---
template: main.html
title: Python API
---

# Python API

The aim of the `gmx_MMPBSA` API is to provide you with direct access to the raw data produced during a `gmx_MMPBSA`
calculation. By default, gmx_MMPBSA calculates an average, standard deviation, and standard error of the mean for 
all the generated data sets, but it does not support custom analyses. The API reads an _GMXMMPBSA_info file, from which 
it will determine what kind of calculation was performed, then automatically parse the output files and load the 
data into arrays.

!!! warning
    The topology files you used in the `gmx_MMPBSA` calculation must also be available in the location specified in the 
    _`GMXMMPBSA_info` file.

## Using the API

We have derived a new API to reorganize the data so that it is arranged more hierarchically. This makes easier to
transform the data into graphs in the `gmx_MMPBSA_ana`. **The original and the current API only differ in the name of
the callable function, the disposition of the data in Per-wise decomposition analysis and in the new 'delta' key.

The function `load_gmxmmpbsa_info` takes the name of a `gmx_MMPBSA` info file (typically `_GMXMMPBSA_info`)
and returns a populated `mmpbsa_data` instance with all the parsed data. An example code snippet that creates
a `mmpbsa_data` instance from the information in _`GMXMMPBSA_info` is shown below.

!!! important
    Unlike MMPBSA.py, `load_gmxmmpbsa_info` does not need to be located in the folder that contains the 
    `_GMXMMPBSA_info` file.

    _New in v1.4.0_

```python
from GMXMMPBSA import API as gmxMMPBSAapi

data = gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA_info")
```


## Properties of mmpbsa_data

The `mmpbsa_data` class is a nested dictionary structure (`mmpbsa_data` is actually derived from `dict`). The various
attributes of `mmpbsa_data` are described below followed by the defined operators.

### Attributes

If the numpy package is installed and available, all data arrays will be `numpy.ndarray` instances. Otherwise, all data
arrays will be `array.array` instances with the ’d’ data type specifier (for a double precision float). The data is
organized in an `mmpbsa_data` instance in the following manner:

```python
mmpbsa_data_instance['calc_key']['system_component']['energy_term']
```

In this example, `calc_key` is a `dict` key that is paired to another `dict` (`mmpbsa_data_instance` is the
first-level `dict`, in this case) (Table 2). The keys of these second-level dict instances
(`system_component`) pair to another `dict` (Table 3).

Table 2. List and description of `calc_key` dict keys that may be present in instances of the `mmpbsa_data`
class.

| Dictionary Key (calc_key) | Calculation Type                     |
|:-------------------------:|:-------------------------------------|
|           `gb`            | Generalized Born Results             |
|           `pb`            | Poisson-Boltzmann Results            |
|         `rism gf`         | Gaussian Fluctuation 3D-RISM Results |
|        `rism std`         | Standard 3D-RISM Results             |
|           `ie`            | Interaction Entropy Results          |
|          `nmode`          | Normal Mode Analysis Results         |
|           `qh`            | Quasi-harmonic Approximation Results |

Table 3. List and description of system_component keys that may be present in instances of the mmpbsa_data class.

| Dictionary Key (system_component) | Description                                      |
|:---------------------------------:|:-------------------------------------------------|
|             `complex`             | Data sets for the complex. (Stability & Binding) |
|            `receptor`             | Data sets for the receptor. (Binding only)       |
|             `ligand`              | Data sets for the ligand. (Binding only)         |
|              `delta`              | Data sets for the delta. (Binding only)          |

The keys of these inner-most (third-level) `dict` instances are paired with the data arrays for that energy term (Table
4). The various dictionary keys are listed below for each level. If alanine scanning was performed, the
`mmpbsa_data_instance` also has a `"mutant"` attribute that contains the same dictionary structure as
`mmpbsa_data` does for the normal system. If not, the mutant attribute is None. The only difference is that the data is
accessed as follows:

```python
mmpbsa_data_instance.mutant['calc_key']['system_component']['energy_term']
```

!!! warning
    All keys are case-sensitive, and if a space appears in the key, it must be present in your program. Also, if
    polar/non-polar decomposition is not performed for `3D-RISM`, then the `POLAR SOLV` and `APOLAR SOLV` keys are 
    replaced with the single key `ERISM`

Table 4. List and description of energy_term keys that may be present in instances of the mmpbsa_data class. The
allowed values of energy_term depend on the value of calc_key above in Table 2. The `energy_term` keys are listed 
for each `calc_key` enumerated above, accompanied by a description. The RISM keys are the same for both `rism gf` 
and `rism std` although the value of `POLAR SOLV` and `APOLAR SOLV` will differ depending on the method chosen. 
Those keys marked with * are specific to the CHARMM force field used through chamber. Those arrays are all 0 for 
normal Amber topology files.

| Description                 |   `gb`    |   `pb`    |    `RISM`     |
|:----------------------------|:---------:|:---------:|:-------------:|
| Bond energy                 |  `BOND`   |  `BOND`   |    `BOND`     |
| Angle energy                |  `ANGLE`  |  `ANGLE`  |    `ANGLE`    |
| Dihedral Energy             |  `DIHED`  |  `DIHED`  |    `DIHED`    |
| Urey-Bradley*               |   `UB`    |   `UB`    |       —       |
| Improper Dihedrals*         |   `IMP`   |   `IMP`   |       —       |
| Correction Map*             |  `CMAP`   |  `CMAP`   |       —       |
| 1-4 van der Waals energy    | `1-4 VDW` | `1-4 VDW` |   `1-4 VDW`   |
| 1-4 Electrostatic energy    | `1-4 EEL` | `1-4 EEL` |   `1-4 EEL`   |
| van der Waals energy        | `VDWAALS` | `VDWAALS` |   `VDWAALS`   |
| Electrostatic energy        |   `EEL`   |   `EEL`   |     `EEL`     |
| Polar solvation energy      |   `EGB`   |   `EPB`   | `POLAR SOLV`  |
| Non-polar solvation energy  |  `ESURF`  | `ENPOLAR` | `APOLAR SOLV` |
| Total solvation free energy | `G solv`  | `G solv`  |   `G solv`    |
| Total gas phase free energy |  `G gas`  |  `G gas`  |    `G gas`    |
| Total energy                |  `TOTAL`  |  `TOTAL`  |    `TOTAL`    |

Table 5. Same as Table 4 for the entropy (nmode and qh) data.

| Description           |     `nmode`     |      `qh`       |
|:----------------------|:---------------:|:---------------:|
| Translational entropy | `Translational` | `Translational` |
| Rotational entropy    |  `Rotational`   |  `Rotational`   |
| Vibrational entropy   |  `Vibrational`  |  `Vibrational`  |
| Total entropy         |     `Total`     |     `Total`     |

Table 6. Same as Table 5 for the Interaction Entropy data.

| Description                                  |   `IE`   |
|:---------------------------------------------|:--------:|
| Data per-frame                               |  `data`  |
| Mean of the selected interval                |  `value` |
| Star and End frames of the selected interval | `frames` |

### Defined operators

In-place addition: It extends all the arrays that are common to both `mmpbsa_data` instances. This is useful if, for
instance, you run two `gmx_MMPBSA` calculations, and you use -prefix <new_prefix> for the second simulation. Assuming
that <new_prefix> is `_GMXMMPBSA2_` for the second `gmx_MMPBSA` calculation, the following pseudo-code will generate a
mmpbsa_data instance with all the data in concatenated arrays. The pseudo-code assumes GMXMMPBSA.API was imported as
demonstrated below.

```python
data = gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA_info")
data += gmxMMPBSAapi.load_gmxmmpbsa_info("_GMXMMPBSA2_info")
```


## Example API Usage

In many cases, the autocorrelation function of the energy can aid in the analysis of MM/PBSA data, since it provides a
way of determining the statistical independence of your data points. For example, 1000 correlated snapshots provide less
information, and therefore less statistical certainty, than 1000 uncorrelated snapshots. The standard error of the mean
calculation performed by `gmx_MMPBSA` assumes a completely uncorrelated set of snapshots, which means that it is a lower
bound of the true standard error of the mean, and a plot of the autocorrelation function may help determine the actual
value.

The example program below will calculate the autocorrelation function of the total energy (complex only for both the
normal and alanine mutant systems) from a GB calculation and plot the resulting code using matplotlib.

```python
import os
import sys

# append AMBERHOME/bin to sys.path
sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))
# Now import the MMPBSA API
from GMXMMPBSA import API as gmxMMPBSAapi
import matplotlib.pyplot as plt
import numpy as np

data = gmxMMPBSAapi.load_gmxmmpbsa_info('_GMXMMPBSA_info')
total = data['gb']['complex']['TOTAL'].copy()
data = gmxMMPBSAapi.load_gmxmmpbsa_info('_GMXMMPBSA_info')
total_mut = data.mutant['gb']['complex']['TOTAL'].copy()
# Create a second copy of the data set. The np.correlate function does not
# normalize the correlation function, so we modify total and total2 to get
# that effect
total -= total.mean()
total /= total.std()
total2 = total.copy() / len(total)
acor = np.correlate(total, total2, 'full')
total_mut -= total_mut.mean()
total_mut /= total_mut.std()
total2_mut = total_mut.copy() / len(total_mut)
acor_mut = np.correlate(total_mut, total2_mut, 'full')
# Now generate the 'lag' axis
xdata = np.arange(0, len(total))
# The acor data set is symmetric about the origin, so only accept the
# positive lag times. Graph the result
plt.plot(xdata, acor[len(acor) // 2:], xdata, acor_mut[len(acor) // 2:])
plt.show()
```

## Decomposition Data

When performing decomposition analysis, the various decomp data is stored in a separate tree of dicts referenced with
the `decomp` key. The key sequence is similar to the sequence for the `normal` data described above, where `decomp` is
followed by the solvent model (`GB` or `PB`), followed by the species (`complex`, `receptor`, or `ligand`) 
(additionally, we include `delta` key), followed by the decomposition components (total, backbone, or sidechain), 
followed by the residue number (and residue pair for pairwise decomposition), finally followed by the contribution 
(`internal`, `van der Waals`, `electrostatics`, etc.) The available keys are shown in Figure 1 (and each key is 
described afterwards).

**Decomp Key Descriptions:**

* `gb` All Generalized Born results
* `pb` All Poisson-Boltzmann results
* `complex` All results from the complex trajectory
* `receptor` All results from the receptor trajectory
* `ligand` All results from the ligand trajectory
* `delta` All results from delta total decomposition [ complex TDC - (receptor TDC + ligand TDC) ]
* `TDC` All results from the total decomposition
* `SDC` All results from the sidechain decomposition
* `BDC` All results from the backbone decomposition
* `#` All data from residue number `#` in per-residue and per-wise decomposition (same residue numbering scheme as in
  each respective topology file)
    * `##` All interaction energies between residues `##` and their respective pair `#` in per-wise decomposition (same
      residue numbering scheme as in each respective topology file)
* `int` Internal energy contributions (see the idecomp variable description above)
* `vdw` van der Waals energy contributions
* `eel` Electrostatic energy contributions
* `pol` Polar solvation free energy contributions
* `sas` Non-polar solvation free energy contributions
* `tot` Total free energy contributions (sum of previous 5).


<figure markdown="1">
[![decomp][1]][1]
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 1**. Tree of `dict` keys following the `decomp` key in a `mmpbsa_data` instance.
  </figcaption>
</figure>

[1]: assets/images/decomp_dict_keys.png

