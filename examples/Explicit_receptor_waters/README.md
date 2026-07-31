---
template: main.html
title: Explicit receptor waters
---

# ST MM/PB(GB)SA with explicit receptor waters

!!! info
    This example can be found in the [examples/Explicit_receptor_waters][6] directory in the repository folder. It
    reuses the topology, structure, trajectory, and index files from the [Protein-protein][8] example, so no additional
    coordinate fixture is required.

## Requirements

In this case, `gmx_MMPBSA` requires:

| Input File required            | Required |           Type             | Description |
|:-------------------------------|:--------:|:--------------------------:|:-------------------------------------------------------------------------------------------------------------|
| Input parameters file          | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `in`          | Input file containing all the specifications regarding the type of calculation that is going to be performed |
| The MD Structure+mass(db) file | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |    `tpr` `pdb`    | Structure file containing the system coordinates |
| An index file                  | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |          `ndx`    | File containing the receptor and ligand in separated groups |
| Receptor and ligand group      | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |        `integers`       | Group numbers in the index files |
| A trajectory file              | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } | `xtc` `pdb` `trr` | Final GROMACS MD trajectory |
| A topology file                | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |           `top`         | GROMACS topology file. The `*.itp` files defined in the topology must be in the same folder |
| PyMOL executable               | :octicons-check-circle-fill-16:{ .req .scale_icon_medium } |       executable        | Required only when `explicit_waters_mask="pymol"` |

:octicons-check-circle-fill-16:{ .req } -> Must be defined

_See a detailed list of all the flags in gmx_MMPBSA command line [here][1]_

## Command-line

Run from the `examples/Protein_protein` folder so the topology can resolve its local `toppar/` includes:

=== "Serial"

        cd examples/Protein_protein
        gmx_MMPBSA -O -i ../Explicit_receptor_waters/mmpbsa_explicit_waters_pymol.in -cs com.tpr -ct com_traj.xtc -ci index.ndx -cg 3 4 -cp topol.top -o FINAL_RESULTS_EXPLICIT_WATERS_PYMOL.dat -eo FINAL_RESULTS_EXPLICIT_WATERS_PYMOL.csv

where the `mmpbsa_explicit_waters_pymol.in` input file is a text file containing the following lines:

``` yaml linenums="1" title="Sample input file for ST GB calculation with explicit receptor waters"
Sample input file for ST GB calculation with explicit receptor waters
# This input keeps 10 waters closest to the PyMOL/dASA interface.
# PyMOL is required only when explicit_waters_mask="pymol".

&general
sys_name="Prot-Prot-ExpWat",
startframe=1,
endframe=10,
forcefields="leaprc.protein.ff14SB",
explicit_waters=10,
explicit_waters_mask="pymol",
explicit_waters_pymol_cutoff=0.5,
/
&gb
igb=2, saltcon=0.150,
/
```

!!! info "Keep in mind"
    See a detailed list of all the options in `gmx_MMPBSA` input file [here][2] as well as several [examples][3].
    This example is meant to show the explicit-water workflow. It is recommended to review the input variables and
    adapt the number of waters, interface definition, GB model, and frame selection to your system.

## Considerations

This mode keeps a fixed number of explicit water molecules in the working complex topology and assigns those waters to
the receptor. It is currently supported for single-trajectory GB or PB calculations only. RISM, GBNSR6, entropy, and
multi-trajectory inputs are not supported with `explicit_waters > 0`.

The PyMOL interface mode identifies interface residues from the first extracted complex PDB using a dASA cutoff. Then
`cpptraj closest` selects the closest waters to that static interface mask in each trajectory frame. This means the
interface residue mask is static, while the water identities can change frame by frame.

If PyMOL is not available, define `explicit_waters_mask` as an Amber residue mask or as a `within <distance>` selection
instead. The selected waters are assigned to the receptor internally, and the ligand topology remains dry.

By default, the explicit-water setup looks for common solvent index groups such as `SOLV`, `SOL`, `Water`, `WAT`,
`TP3`, and `OPC`. If the solvent group has a custom name, set `explicit_waters_group` in `&general`.

Extra-point water models such as OPC or TIP4P can fail in `sander` because of their virtual-site atoms. By default,
`gmx_MMPBSA` stops when these atoms are found. Set `explicit_waters_extra_points="strip"` only if you intentionally
want to remove the virtual sites and use the result as an approximate relative comparison.

Useful generated files for checking the setup are:

* `_GMXMMPBSA_explicit_waters_interface.dat`: PyMOL/dASA interface residues used to build the water reference mask
* `_GMXMMPBSA_explicit_waters_closest_0.dat`: water molecules selected by `cpptraj closest` in each frame
* `COM.prmtop`: complex topology with the requested number of water residues
* `REC.prmtop`: receptor topology with those water residues assigned to the receptor
* `LIG.prmtop`: dry ligand topology

A plain text output file with all the statistics and a CSV-format output file containing all energy terms for every
frame in every calculation will be saved.

!!! note
    Once the calculation is done, the results can be analyzed in `gmx_MMPBSA_ana` if `-nogui` was not used in the
    command-line. Please, check the [gmx_MMPBSA_ana][5] section for more information.

  [1]: ../../gmx_MMPBSA_command-line.md#gmx_mmpbsa-command-line
  [2]: ../../input_file.md#the-input-file
  [3]: ../../input_file.md#sample-input-files
  [5]: ../../analyzer.md#gmx_mmpbsa_ana-the-analyzer-tool
  [6]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Explicit_receptor_waters
  [8]: ../Protein_protein/README.md
