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


def _set_context_ylim(axis, values, reference):
    """Set an honest, readable y-scale for a stability trace.

    Keep a minimum one-percent-of-reference context around small variations so
    absolute changes are not visually overstated, while retaining extra room
    around the observed range.
    """
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if not len(finite_values):
        return
    observed_min = float(np.min(finite_values))
    observed_max = float(np.max(finite_values))
    observed_span = observed_max - observed_min
    minimum_span = max(abs(float(reference)) * 0.01, 0.1)
    visible_span = max(observed_span * 1.2, minimum_span)
    observed_midpoint = (observed_min + observed_max) / 2.0
    axis.set_ylim(
        observed_midpoint - visible_span / 2.0,
        observed_midpoint + visible_span / 2.0,
    )


def write_diagnostics(frame_diagnostics, csv_path, png_path, atom_names,
                      center_setting, thickness_setting, center, thickness,
                      frame_coordinates=None):
    """Write an auditable CSV and a leaflet-resolved diagnostic plot.

    ``frame_coordinates`` is optional for compatibility with callers that
    only have the summary rows. When supplied, the plot shows the selected
    membrane atom coordinates for each frame instead of reducing each frame
    to a single mean value.
    """
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
    import seaborn as sns
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FormatStrFormatter, MaxNLocator

    sns.set_theme(style='white')

    frames = [row['frame'] for row in frame_diagnostics]
    means = [row['mean_z_A'] for row in frame_diagnostics]
    frame_thickness = [row['thickness_A'] for row in frame_diagnostics]
    if frame_coordinates is None:
        frame_coordinates = {
            row['frame']: [row['mean_z_A']] for row in frame_diagnostics
        }

    fig, axes = plt.subplots(
        3, 1, figsize=(9, 9.6), sharex=True,
        gridspec_kw={'height_ratios': (2.4, 1, 1)},
    )
    coordinate_axis, center_axis, thickness_axis = axes

    # Plot a deterministic, lightly jittered view of the actual selected
    # atoms. Limit only the rendered view; the CSV retains every diagnostic
    # frame. Consolidating points into two scatter collections also keeps
    # rendering fast for long trajectories.
    lower_color = '#4C78A8'
    upper_color = '#C27C5D'
    center_color = '#8C6BB1'
    thickness_color = '#4FAF8F'
    boundary_color = thickness_color
    slab_color = '#A8DCC6'
    max_plot_points = 20_000
    frame_values = []
    total_points = 0
    coordinate_min = math.inf
    coordinate_max = -math.inf
    for frame in frames:
        values = np.asarray(frame_coordinates.get(frame, ()), dtype=float)
        values = values[np.isfinite(values)]
        frame_values.append(values)
        if len(values):
            total_points += len(values)
            coordinate_min = min(coordinate_min, float(np.min(values)))
            coordinate_max = max(coordinate_max, float(np.max(values)))

    point_stride = max(1, math.ceil(total_points / max_plot_points))
    lower_x, lower_y, upper_x, upper_y = [], [], [], []
    plotted_coordinates = []
    for x_position, values in enumerate(frame_values, start=1):
        if not len(values):
            continue
        if point_stride > 1:
            selected = np.linspace(
                0, len(values) - 1,
                max(1, math.ceil(len(values) / point_stride)), dtype=int,
            )
            values = np.sort(values)[selected]
        plotted_coordinates.append(values)
        jitter = ((np.arange(len(values)) * 0.61803398875) % 1.0 - 0.5) * 0.58
        lower = values <= center
        upper = ~lower
        lower_x.extend((x_position + jitter[lower]).tolist())
        lower_y.extend(values[lower].tolist())
        upper_x.extend((x_position + jitter[upper]).tolist())
        upper_y.extend(values[upper].tolist())

    coordinate_axis.scatter(
        lower_x, lower_y, color=lower_color, s=8, alpha=0.55,
        linewidths=0, rasterized=True,
    )
    coordinate_axis.scatter(
        upper_x, upper_y, color=upper_color, s=8, alpha=0.55,
        linewidths=0, rasterized=True,
    )

    slab_low = center - thickness / 2.0
    slab_high = center + thickness / 2.0
    coordinate_axis.axhspan(
        slab_low, slab_high, color=slab_color, alpha=0.16,
        label=f'resolved slab = {thickness:.1f} Å',
    )
    coordinate_axis.axhline(
        center, color=center_color, linestyle='--', linewidth=1.4,
        label=f'center = {center:.3f} Å',
    )
    coordinate_axis.axhline(slab_low, color=boundary_color, linestyle=':', linewidth=1)
    coordinate_axis.axhline(slab_high, color=boundary_color, linestyle=':', linewidth=1)
    dimension_x = 1 + 0.92 * max(len(frames) - 1, 1)
    coordinate_axis.annotate(
        '', xy=(dimension_x, slab_high), xytext=(dimension_x, slab_low),
        arrowprops={
            'arrowstyle': '<->', 'color': boundary_color,
            'linewidth': 1.4, 'shrinkA': 0, 'shrinkB': 0,
        },
    )
    coordinate_axis.text(
        dimension_x - 0.015 * max(len(frames), 1),
        center + 0.10 * thickness,
        f'thickness = {thickness:.1f} Å', color=thickness_color,
        ha='right', va='bottom', fontsize=9,
        bbox={'facecolor': 'white', 'alpha': 0.75, 'edgecolor': 'none', 'pad': 2},
    )
    coordinate_axis.text(
        0.06, center + 0.04 * thickness,
        f'center = {center:.3f} Å', color=center_color,
        ha='left', va='bottom', fontsize=10,
        transform=coordinate_axis.get_yaxis_transform(),
        bbox={'facecolor': 'white', 'alpha': 0.75, 'edgecolor': 'none', 'pad': 2},
    )
    coordinate_axis.set_ylabel(f'{"; ".join(atom_names)} atoms z coordinate (Å)')
    coordinate_axis.set_title(f'Selected membrane atoms {"; ".join(atom_names)} resolve into two leaflets')
    coordinate_axis.set_xlim(0.5, len(frames) + 0.5)
    coordinate_axis.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    if plotted_coordinates:
        coordinate_min = min(coordinate_min, slab_low)
        coordinate_max = max(coordinate_max, slab_high)
        coordinate_padding = max((coordinate_max - coordinate_min) * 0.06, 0.5)
        coordinate_axis.set_ylim(
            coordinate_min - coordinate_padding,
            coordinate_max + coordinate_padding,
        )
    coordinate_axis.legend(handles=(
        Line2D([], [], marker='o', linestyle='None', color=lower_color,
               markersize=5, label='lower leaflet'),
        Line2D([], [], marker='o', linestyle='None', color=upper_color,
               markersize=5, label='upper leaflet'),
        Line2D([], [], color=center_color, linestyle='--',
               label=f'center = {center:.3f} Å'),
        Line2D([], [], color=boundary_color, linestyle=':',
               label='slab boundaries'),
        Line2D([], [], color=thickness_color, linestyle='-', linewidth=1.4,
               label=f'thickness = {thickness:.1f} Å'),
        Patch(facecolor=slab_color, alpha=0.16,
              label=f'resolved slab = {thickness:.1f} Å'),
    ), loc='lower center', bbox_to_anchor=(0.5, 1.1), borderaxespad=0,
        ncol=3, frameon=True)
    # coordinate_axis.grid(axis='y', alpha=0.25)

    center_axis.plot(frames, means, color=center_color, marker='o', markersize=3, linewidth=1.0)
    center_axis.axhline(
        center, color=center_color, linestyle='--', linewidth=1.2,
        label=f'resolved center = {center:.3f} Å',
    )
    center_axis.set_ylabel('Frame mean z (Å)')
    center_axis.set_title('Frame-to-frame center stability', loc='left', fontsize=10)
    center_axis.ticklabel_format(axis='y', style='plain', useOffset=False)
    center_axis.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    _set_context_ylim(center_axis, means, center)
    center_axis.text(
        0.99, 0.06, f'observed range = {max(means) - min(means):.3f} Å',
        transform=center_axis.transAxes, ha='right', va='bottom',
        fontsize=8, color='dimgray',
    )
    center_axis.legend(loc='best')
    # center_axis.grid(axis='y', alpha=0.25)

    thickness_axis.plot(
        frames, frame_thickness, color=thickness_color, marker='o', markersize=3,
        linewidth=1.0,
    )
    thickness_axis.axhline(
        thickness, color=thickness_color, linestyle='--', linewidth=1.2,
        label=f'resolved thickness = {thickness:.1f} Å',
    )
    thickness_axis.set_xlabel('Selected -ct frame')
    thickness_axis.set_ylabel('Thickness (Å)')
    thickness_axis.set_title('Frame-to-frame thickness stability', loc='left', fontsize=10)
    thickness_axis.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    _set_context_ylim(thickness_axis, frame_thickness, thickness)
    finite_thickness = [value for value in frame_thickness if math.isfinite(value)]
    if finite_thickness:
        thickness_axis.text(
            0.99, 0.06,
            f'observed range = {max(finite_thickness) - min(finite_thickness):.1f} Å',
            transform=thickness_axis.transAxes, ha='right', va='bottom',
            fontsize=8, color='dimgray',
        )
    thickness_axis.legend(loc='best')
    # thickness_axis.grid(axis='y', alpha=0.25)

    fig.align_ylabels(axes)
    # fig.suptitle(
    #     f'Membrane parameters from {"; ".join(atom_names)} atoms',
    #     fontsize=14, x=0.55, y=0.97,
    # )
    fig.subplots_adjust(
        left=0.12, right=0.98, bottom=0.09, top=0.84, hspace=0.3,
    )
    fig.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
