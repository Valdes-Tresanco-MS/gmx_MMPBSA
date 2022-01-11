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

import math
from typing import Union

from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, GMXMMPBSA_WARNING
import pandas as pd
import numpy as np
from queue import Queue
from PyQt5.QtCore import *
import multiprocessing
from pathlib import Path

R = 0.001987

ncpu = multiprocessing.cpu_count()

style = Path(__file__).joinpath('style')


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


def energy2pdb_pml(residue_list, colors, pml_path: Path, pdb_path: Path):
    pymol_colors = []
    with open(pml_path, 'w') as bf:
        bf.write(f'load {pdb_path}\n')
        bf.write('set cartoon_oval_length, 1.0\n')
        bf.write('set cartoon_rect_length, 1.2\n')
        bf.write('set cartoon_rect_width, 0.3\n')
        bf.write('set cartoon_side_chain_helper, 1\n')
        bf.write('set light_count, 0\n')
        for c, color in enumerate(colors):
            bf.write(f'set_color gmxc{c} = {color}\n')
            pymol_colors.append(f'gmxc{c}')
        bf.write('set bg_rgb, gray50\n')

        minimum = 999
        maximum = -999
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

        bf.write(f'show sticks, {select_text}\n')
        bf.write('remove (h. and (e. c extend 1))\n')
        bf.write(f'spectrum b, {" ".join(pymol_colors)}, minimum={minimum}, maximum={maximum}\n')
        bf.write(f'ramp_new colorbar, none, {np.linspace(minimum, maximum, len(colors)).tolist()}, {pymol_colors}\n')
        bf.write(f'center {select_text}\n')
        bf.write(f'orient {select_text}\n')
        bf.write(f'zoom {select_text}, 15\n')


def ki2energy(ki, temp):
    # deltaG (inhibition) = R * T * ln ( Ki )
    if not ki:
        return np.nan
    return R * temp * math.log(ki * 1e-9)


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
            GMXMMPBSA_WARNING(f'{cf} not found')
            continue
        if cf.is_dir():
            recursive_files = []
            if parser_args.recursive:
                recursive_files.extend(cf.absolute().glob('*/*_info'))
            else:
                recursive_files.extend(cf.absolute().glob('*_info'))
            for rf in recursive_files:
                if rf in info_files:
                    GMXMMPBSA_WARNING(f'{rf} is duplicated and will be ignored')
                    continue
                info_files.append(rf)
        else:
            if cf in info_files:
                GMXMMPBSA_WARNING(f'{cf} is duplicated and will be ignored')
                continue
            info_files.append(cf)
    if not len(info_files):
        GMXMMPBSA_ERROR('No info files found!')
    return info_files
