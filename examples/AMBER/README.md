---
template: main.html
title: AMBER input files
---

# AMBER input files binding free energy calculations

!!! info
    This example can be found in the [examples/AMBER][6] directory in the repository folder. If you didn't
    use gmx_MMPBSA_test before, use [downgit](https://downgit.github.io/#/home) to download the specific folder from
    gmx_MMPBSA GitHub repository.

## Requirements

In this case, `gmx_MMPBSA_amber` requires:

| Input File required            | Required | Type | Description |
|:-------------------------------|:--------:|:----:|:------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `in` | Input file containing all calculation specifications |
| Complex topology file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `prmtop` | AMBER topology file for the complex |
| Complex structure file         | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `inpcrd` | AMBER coordinate file for the complex |
| Complex trajectory file        | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `mdcrd` | AMBER trajectory file for the complex |
| Receptor and ligand masks      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `masks` | AMBER masks identifying the receptor and ligand in the complex |
| Receptor topology file         | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium } | `prmtop` | AMBER topology file for the receptor |
| Ligand topology file           | :octicons-check-circle-fill-16:{ .req_opt .scale_icon_medium } | `prmtop` | AMBER topology file for the ligand |

:octicons-check-circle-fill-16:{ .req } -> Must be defined -- :octicons-check-circle-fill-16:{ .req_opt } -> Optional

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line

Once you are in the folder containing all files, the command-line will be as follows:

=== "Serial"

        gmx_MMPBSA_amber -O -i mmpbsa.in -cp ras-raf_complex.prmtop -cs ras-raf_complex.inpcrd -ct prod_complex.mdcrd -rp ras.prmtop -lp raf.prmtop -cm ":1-166" ":167-242" -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv

=== "gmx_MMPBSA_test"

        gmx_MMPBSA_test -t 25

where the `mmpbsa.in` input file is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for AMBER input files"
Input file for a short AMBER input files test
&general
   sys_name="AMBER",
   startframe=1,
   endframe=5,
   interval=1,
   verbose=1,
/
&gb
   igb=2, saltcon=0.100,
/
```

!!! info "Keep in mind"
    This example is meant only to show that `gmx_MMPBSA_amber` works with native AMBER input files. It uses a short
    five-frame trajectory and preserves the radii already stored in the AMBER topology files.

## Considerations

In this case, a single trajectory (ST) approximation is followed. The receptor and ligand are selected from the
complex with AMBER masks (`:1-166` and `:167-242`). Receptor and ligand topology files are provided directly, so
`gmx_MMPBSA_amber` can build the working topologies without using GROMACS files.

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/AMBER
