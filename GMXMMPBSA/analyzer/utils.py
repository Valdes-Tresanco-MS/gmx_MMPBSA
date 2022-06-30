# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################
try:
    from PyQt6.QtCore import *
except:
    from PyQt5.QtCore import *

import logging
import math
from typing import Union

from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import pandas as pd
import numpy as np
from queue import Queue
import multiprocessing
from pathlib import Path

R = 0.001987

ncpu = multiprocessing.cpu_count()

style = Path(__file__).joinpath('style')


def bar_label(ax, container, labels=None, *, fmt="%g", label_type="edge", padding=0, **kwargs):
    """
    Fixed bar_label for inverted axes. matplotlib >= 3.5.2 ???

    PR: https://github.com/matplotlib/matplotlib/pull/22447

    Label a bar plot.

    Adds labels to bars in the given `.BarContainer`.
    You may need to adjust the axis limits to fit the labels.

    Parameters
    ----------
    container : `.BarContainer`
        Container with all the bars and optionally errorbars, likely
        returned from `.bar` or `.barh`.

    labels : array-like, optional
        A list of label texts, that should be displayed. If not given, the
        label texts will be the data values formatted with *fmt*.

    fmt : str, default: '%g'
        A format string for the label.

    label_type : {'edge', 'center'}, default: 'edge'
        The label type. Possible values:

        - 'edge': label placed at the end-point of the bar segment, and the
          value displayed will be the position of that end-point.
        - 'center': label placed in the center of the bar segment, and the
          value displayed will be the length of that segment.
          (useful for stacked bars, i.e.,
          :doc:`/gallery/lines_bars_and_markers/bar_label_demo`)

    padding : float, default: 0
        Distance of label from the end of the bar, in points.

    **kwargs
        Any remaining keyword arguments are passed through to
        `.Axes.annotate`. The alignment parameters (
        *horizontalalignment* / *ha*, *verticalalignment* / *va*) are
        not supported because the labels are automatically aligned to
        the bars.

    Returns
    -------
    list of `.Text`
        A list of `.Text` instances for the labels.
    """
    from matplotlib import _api
    import itertools

    for key in ['horizontalalignment', 'ha', 'verticalalignment', 'va']:
        if key in kwargs:
            raise ValueError(
                f"Passing {key!r} to bar_label() is not supported.")

    # want to know whether to put label on positive or negative direction
    # cannot use np.sign here because it will return 0 if x == 0
    a, b = ax.yaxis.get_view_interval()
    y_inverted = a > b
    c, d = ax.xaxis.get_view_interval()
    x_inverted = c > d

    def sign(x):
        return 1 if x >= 0 else -1

    _api.check_in_list(['edge', 'center'], label_type=label_type)

    bars = container.patches
    errorbar = container.errorbar
    datavalues = container.datavalues
    orientation = container.orientation

    if errorbar:
        # check "ErrorbarContainer" for the definition of these elements
        lines = errorbar.lines  # attribute of "ErrorbarContainer" (tuple)
        barlinecols = lines[2]  # 0: data_line, 1: caplines, 2: barlinecols
        barlinecol = barlinecols[0]  # the "LineCollection" of error bars
        errs = barlinecol.get_segments()
    else:
        errs = []

    if labels is None:
        labels = []

    annotations = []

    for bar, err, dat, lbl in itertools.zip_longest(
            bars, errs, datavalues, labels
    ):
        (x0, y0), (x1, y1) = bar.get_bbox().get_points()
        xc, yc = (x0 + x1) / 2, (y0 + y1) / 2

        if orientation == "vertical":
            extrema = max(y0, y1) if dat >= 0 else min(y0, y1)
            length = abs(y0 - y1)
        elif orientation == "horizontal":
            extrema = max(x0, x1) if dat >= 0 else min(x0, x1)
            length = abs(x0 - x1)

        if err is None:
            endpt = extrema
        elif orientation == "vertical":
            endpt = err[:, 1].max() if dat >= 0 else err[:, 1].min()
        elif orientation == "horizontal":
            endpt = err[:, 0].max() if dat >= 0 else err[:, 0].min()

        if label_type == "center":
            value = sign(dat) * length
        elif label_type == "edge":
            value = extrema

        if label_type == "center":
            xy = xc, yc
        elif label_type == "edge" and orientation == "vertical":
            xy = xc, endpt
        elif label_type == "edge" and orientation == "horizontal":
            xy = endpt, yc

        if orientation == "vertical":
            y_direction = -1 if y_inverted else 1
            xytext = 0, y_direction * sign(dat) * padding
        else:
            x_direction = -1 if x_inverted else 1
            xytext = x_direction * sign(dat) * padding, 0

        if label_type == "center":
            ha, va = "center", "center"
        elif label_type == "edge":
            if orientation == "vertical":
                ha = 'center'
                if y_inverted:
                    va = 'top' if dat > 0 else 'bottom'  # also handles NaN
                else:
                    va = 'top' if dat < 0 else 'bottom'  # also handles NaN
            elif orientation == "horizontal":
                if x_inverted:
                    ha = 'right' if dat > 0 else 'left'  # also handles NaN
                else:
                    ha = 'right' if dat < 0 else 'left'  # also handles NaN
                va = 'center'

        if np.isnan(dat):
            lbl = ''

        annotation = ax.annotate(fmt % value if lbl is None else lbl,
                                 xy, xytext, textcoords="offset points",
                                 ha=ha, va=va, **kwargs)
        annotations.append(annotation)

    return annotations


def com2str(com_pdb_str):
    import parmed
    import tempfile
    fp = tempfile.TemporaryFile(mode='w+t')
    fp.writelines(com_pdb_str)
    fp.seek(0)
    return parmed.read_PDB(fp)


def multiindex2dict(p: Union[pd.MultiIndex, dict]) -> dict:
    """
    Converts a pandas Multiindex to a nested dict
    :parm p: As this is a recursive function, initially p is a pd.MultiIndex, but after the first iteration it takes
    the internal_dict value, so it becomes to a dictionary
    """
    internal_dict = {}
    end = False
    for x in p:
        # Since multi-indexes have a descending hierarchical structure, it is convenient to start from the last
        # element of each tuple. That is, we start by generating the lower level to the upper one. See the example
        if isinstance(p, pd.MultiIndex):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used. Only for 2 levels
            # pd.MultiIndex
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = [x[-1]]
            else:
                internal_dict[t].append(x[-1])
        elif isinstance(x, tuple):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = {x[-1]: p[x]}
            else:
                internal_dict[t][x[-1]] = p[x]
    if end:
        return internal_dict
    return multiindex2dict(internal_dict)


def calculatestar(args):
    return run_process(*args)


def run_process(func, args):
    results = func(args[1])
    return (args, results)


class worker(QThread):
    job_finished = pyqtSignal()

    def __init__(self):
        super(worker, self).__init__()

    def define_dat(self, function, queue: Queue, result_q: Queue, jobs: int = 1):
        self.fn = function
        self.queue = queue
        self.result_queue = result_q
        self.jobs = jobs

    def run(self):
        size = self.queue.qsize()
        TASKS = [[self.fn, self.queue.get_nowait()] for _ in range(size)]
        self.jobs = min(self.jobs, len(TASKS))
        with multiprocessing.Pool(self.jobs) as pool:
            imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)
            for result in imap_unordered_it:
                self.job_finished.emit()
                self.result_queue.put(result)


def energy2pdb_pml(residue_list, options, pml_path: Path, pdb_path: Path):
    pymol_colors = []
    with open(pml_path, 'w') as bf:
        bf.write(f'load {pdb_path}\n')
        # backwards compatibility
        bf.write("hide lines\n")
        bf.write(f"set cartoon_oval_length, {options['cartoon_oval_length']}\n")
        bf.write(f"set cartoon_rect_length, {options['cartoon_rect_length']}\n")
        bf.write(f"set cartoon_rect_width, {options['cartoon_rect_width']}\n")
        bf.write(f"set cartoon_side_chain_helper, {options['cartoon_side_chain_helper']}\n")
        bf.write(f"set light_count, {options['light_count']}\n")
        for c, color in enumerate(options['colors']):
            bf.write(f'set_color gmxc{c} = {color}\n')
            pymol_colors.append(f'gmxc{c}')
        bf.write(f"set bg_rgb, {options['bg_rgb']}\n")

        minimum = 99999
        maximum = -99999
        select = {}
        for res, energy in residue_list.items():
            icode = ''
            if res.count(':') == 3:
                if res.split(':')[0] in ['R', 'L']:
                    id, chain, name, number = res.split(':')
                else:
                    chain, name, number, icode = res.split(':')
            elif res.count(':') == 4:
                id, chain, name, number, icode = res.split(':')
            else:
                chain, name, number = res.split(':')

            if chain not in select:
                select[chain] = []
            select[chain].append(f'{number}{icode}')
            if energy < minimum:
                minimum = math.floor(energy)
            if energy > maximum:
                maximum = math.ceil(energy)

        chain_sele = []
        for chain in select:
            text = f"(chain {chain} and resi {'+'.join(select[chain])})"
            chain_sele.append(text)
        select_text = ' or '.join(chain_sele)

        if abs(minimum) > abs(maximum):
            maximum = abs(minimum)
        else:
            minimum = -abs(maximum)

        bf.write("show cartoon, (all)\n")
        bf.write(f'select sele_residues, {select_text}\n')
        if options['representation'] not in ['lines+mesh', 'sticks+mesh', 'lines+dots', 'sticks+dots']:
            bf.write(f"show {options['representation']}, sele_residues\n")
        else:
            bf.write(f"show {options['representation'].split('+')[0]}, sele_residues\n")
            bf.write(f"show {options['representation'].split('+')[1]}, sele_residues\n")
        bf.write('remove (h. and (e. c extend 1))\n')
        bf.write(f'spectrum b, {" ".join(pymol_colors)}, minimum={minimum}, maximum={maximum}\n')
        bf.write(f'ramp_new energybar, none, {np.linspace(minimum, maximum, len(pymol_colors)).tolist()},'
                 f' {pymol_colors}\n')
        bf.write('center sele_residues\n')
        bf.write('orient sele_residues\n')
        bf.write('zoom sele_residues, 15\n')


def ki2energy(ki, temp):
    # deltaG (inhibition) = R * T * ln ( Ki )
    return R * temp * math.log(ki * 1e-9) if ki else np.nan


def make_corr_DF(corr_data: dict) -> pd.DataFrame:
    data = []
    for x, value in corr_data.items():
        if x == 'mutant':
            continue
        for m in value['ΔG']:
            curr = [x] + [np.nanmean(corr_data[x]['ΔG'][m][d]) if type(corr_data[x]['ΔG'][m][d]) == np.ndarray
                          else np.nan for d in corr_data[x]['ΔG'][m]] + [m, corr_data[x]['Exp.Energy']]
            data.append(curr)
    return pd.DataFrame(
        data=data,
        columns=[
            'System',
            'ΔH',
            'ΔH+IE',
            'ΔH+NMODE',
            'ΔH+QH',
            'MODEL',
            'Exp.Energy',
        ],
    )


def get_files(parser_args):
    info_files = []
    if not parser_args.files:
        parser_args.files = [Path('.')]
    for cf in parser_args.files:
        if not cf.exists():
            logging.warning(f'{cf} not found')
            continue
        if cf.is_dir():
            recursive_files = []
            if parser_args.recursive:
                if files := list(cf.absolute().glob('*/*_info')):
                    recursive_files.extend(files)
                else:
                    recursive_files.extend(cf.absolute().glob('*/COMPACT_MMXSA_RESULTS.mmxsa'))
            elif files := list(cf.absolute().glob('*_info')):
                recursive_files.extend(files)
            else:
                recursive_files.extend(cf.absolute().glob('COMPACT_MMXSA_RESULTS.mmxsa'))
            for rf in recursive_files:
                if rf in info_files:
                    logging.warning(f'{rf} is duplicated and will be ignored')
                    continue
                info_files.append(rf)
        else:
            if cf in info_files:
                logging.warning(f'{cf} is duplicated and will be ignored')
                continue
            info_files.append(cf)
    if not len(info_files):
        GMXMMPBSA_ERROR('No info files found!')
    return info_files
