"""Automatic membrane center and thickness estimation.

The estimator follows the Amber MMPBSA.py membrane implementation: selected
head-group atom z coordinates are pooled over the selected complex trajectory
frames, the membrane center is their mean, and thickness is the distance
between the means of the coordinates above and below that center.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from GMXMMPBSA.exceptions import InputError, MMPBSA_Error


AUTOMATIC = 'automatic'
DEFAULT_ATOMS = ('P',)


def normalize_parameter(value, name):
    """Return a membrane parameter as a float or the ``automatic`` marker."""
    if isinstance(value, str):
        value = value.strip()
        if value.lower() == AUTOMATIC:
            return AUTOMATIC
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(
            f'{name} must be a number or "{AUTOMATIC}"; got {value!r}'
        ) from exc


def parse_atom_names(value):
    """Parse the simple semicolon-separated atom-name selection."""
    if value is None or not str(value).strip():
        raise InputError('membrane_atoms must contain at least one atom name')
    names = tuple(dict.fromkeys(
        item.strip() for item in str(value).replace(',', ';').split(';') if item.strip()
    ))
    if not names:
        raise InputError('membrane_atoms must contain at least one atom name')
    if any(any(char.isspace() for char in name) for name in names):
        raise InputError('membrane_atoms accepts atom names separated by semicolons, not masks')
    if any(any(char in name for char in ':@*?[](),/\\') for name in names):
        raise InputError('membrane_atoms accepts atom names separated by semicolons, not masks')
    return names


def needs_automatic_parameters(input_data):
    """Whether membrane coordinates are needed for the current PB settings."""
    pb = input_data.get('pb')
    if pb is None or pb.get('memopt', 0) <= 0:
        return False
    return pb['memopt'] > 0 and (pb['mctrdz'] == AUTOMATIC or pb['mthick'] == AUTOMATIC)


def _read_pdb_z(path):
    coordinates = []
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        try:
            coordinates.append(float(line[46:54]))
        except (IndexError, ValueError):
            try:
                coordinates.append(float(line[30:54].split()[2]))
            except (IndexError, ValueError) as exc:
                try:
                    coordinates.append(float(line.split()[7]))
                except (IndexError, ValueError):
                    raise MMPBSA_Error(
                        f'Could not read a z coordinate from membrane atom file {path}'
                    ) from exc
    return coordinates


def read_extracted_coordinates(directory, atom_names):
    """Read ``maskpdb`` output into an ordered frame-to-z-coordinate mapping."""
    frames = defaultdict(list)
    for atom_name in atom_names:
        base = Path(directory) / f'{atom_name}.pdb'
        files = sorted(base.parent.glob(f'{base.name}*'))
        if not files:
            raise MMPBSA_Error(
                f'No coordinates were found for membrane atom name {atom_name!r}. '
                'Check membrane_atoms against the atom names in -ct.'
            )
        found = 0
        for path in files:
            suffix = path.name[len(base.name):]
            if suffix and suffix[1:].isdigit():
                frame = int(suffix[1:])
            elif not suffix:
                frame = 1
            else:
                continue
            coordinates = _read_pdb_z(path)
            if coordinates:
                frames[frame].extend(coordinates)
                found += len(coordinates)
        if found == 0:
            raise MMPBSA_Error(
                f'Membrane atom name {atom_name!r} was not present in the selected -ct frames.'
            )
    if not frames:
        raise MMPBSA_Error('No membrane atom coordinates were extracted from -ct')
    return dict(sorted(frames.items()))


def calculate_parameters(frame_coordinates, center_setting, thickness_setting):
    """Calculate resolved membrane parameters and per-frame diagnostics."""
    frame_arrays = [np.asarray(values, dtype=float) for values in frame_coordinates.values()]
    z_coordinates = np.concatenate(frame_arrays)

    if center_setting == AUTOMATIC:
        center = float(np.mean(z_coordinates))
    else:
        center = float(center_setting)

    if thickness_setting == AUTOMATIC:
        upper = z_coordinates[z_coordinates > center]
        lower = z_coordinates[z_coordinates <= center]
        if not len(upper) or not len(lower):
            raise MMPBSA_Error(
                'Automatic membrane thickness requires atoms on both sides of the membrane center.'
            )
        thickness = float(np.mean(upper) - np.mean(lower))
    else:
        thickness = float(thickness_setting)

    if not math.isfinite(center) or not math.isfinite(thickness) or thickness <= 0:
        raise MMPBSA_Error(
            f'Invalid membrane parameters calculated from -ct: center={center}, thickness={thickness}'
        )

    diagnostics = []
    for frame, values in zip(frame_coordinates, frame_arrays):
        upper = values[values > center]
        lower = values[values <= center]
        frame_thickness = float(np.mean(upper) - np.mean(lower)) if len(upper) and len(lower) else math.nan
        diagnostics.append({
            'frame': frame,
            'n_atoms': len(values),
            'n_above': len(upper),
            'n_below': len(lower),
            'mean_z_A': float(np.mean(values)),
            'thickness_A': frame_thickness,
        })
    return center, thickness, diagnostics


def diagnostic_paths(prefix=''):
    """Return retained output paths, outside any temporary-file prefix."""
    return (Path('GMXMMPBSA_membrane_parameters.csv'), Path('GMXMMPBSA_membrane_parameters.png'))


def write_diagnostics(frame_diagnostics, csv_path, png_path, atom_names,
                      center_setting, thickness_setting, center, thickness):
    """Write an auditable CSV and a compact diagnostic plot."""
    with Path(csv_path).open('w', newline='') as handle:
        handle.write(f'# selected_atoms={";".join(atom_names)}\n')
        handle.write(f'# center_setting={center_setting}\n')
        handle.write(f'# thickness_setting={thickness_setting}\n')
        handle.write(f'# resolved_center_A={center:.6f}\n')
        handle.write(f'# resolved_thickness_A={thickness:.6f}\n')
        writer = csv.DictWriter(handle, fieldnames=(
            'frame', 'n_atoms', 'n_above', 'n_below', 'mean_z_A', 'thickness_A'
        ))
        writer.writeheader()
        writer.writerows(frame_diagnostics)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    frames = [row['frame'] for row in frame_diagnostics]
    means = [row['mean_z_A'] for row in frame_diagnostics]
    frame_thickness = [row['thickness_A'] for row in frame_diagnostics]
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(frames, means, marker='o', markersize=2, linewidth=1)
    axes[0].axhline(center, color='tab:red', linestyle='--', label=f'resolved center = {center:.3f} Å')
    axes[0].set_ylabel('Frame mean z (Å)')
    axes[0].legend(loc='best')
    axes[1].plot(frames, frame_thickness, marker='o', markersize=2, linewidth=1)
    axes[1].axhline(thickness, color='tab:red', linestyle='--', label=f'resolved thickness = {thickness:.3f} Å')
    axes[1].set_xlabel('Selected -ct frame')
    axes[1].set_ylabel('Frame thickness (Å)')
    axes[1].legend(loc='best')
    fig.suptitle(f'Membrane parameters from {"; ".join(atom_names)} atoms')
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
