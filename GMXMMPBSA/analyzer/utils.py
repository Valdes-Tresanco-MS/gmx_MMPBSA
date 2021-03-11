# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
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
from pathlib import Path
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import logging
import pandas as pd
import numpy as np
from queue import Queue
from PyQt5.QtCore import *


R = 0.001987


class worker(QThread):
    job_finished = pyqtSignal()
    def __init__(self):
        super(worker, self).__init__()
        # self.parent = parent
    def define_dat(self, function, queue: Queue, result_q: Queue):
        self.fn = function
        self.queue = queue
        self.result_queue = result_q
        self.running = True

    def run(self):
        while self.running:
            if self.queue.qsize() == 0:
                return
            items = self.queue.get()
            try:
                results = self.fn(items)
                if len(results) == 1:
                    self.result_queue.put(results)
                else:
                    self.result_queue.put([x for x in results])
            finally:
                self.queue.task_done()
                self.job_finished.emit()


def energy2pdb_pml(residue_list, pml_path: Path, pdb_path: Path):
    with open(pml_path, 'w') as bf:
        bf.write(f'load {pdb_path}\n')
        bf.write(f'set cartoon_oval_length, 1.0\n')
        bf.write(f'set cartoon_rect_length, 1.2\n')
        bf.write(f'set cartoon_rect_width, 0.3\n')
        bf.write(f'set cartoon_side_chain_helper, 1\n')
        bf.write(f'set light_count, 0\n')
        bf.write('set_color gmxc1 = (0.0, 0.0, 0.3)\n')
        bf.write('set_color gmxc2 = (0.0, 0.0, 1.0)\n')
        bf.write('set_color gmxc3 = (1.0, 1.0, 1.0)\n')
        bf.write('set_color gmxc4 = (1.0, 0.0, 0.0)\n')
        bf.write('set_color gmxc5 = (0.5, 0.0, 0.0)\n')
        bf.write('set bg_rgb, gray50\n')

        minimum = 999
        maximum = -999
        selections = []
        for res, energy in residue_list.items():
            if res.count(':') == 3:
                chain, name, number, icode = res.split(':')
            else:
                chain, name, number = res.split(':')
            sele = f'(chain {chain} and resi {number})'
            selections.append(sele)
            bf.write(f'show sticks, {sele}\n')
            if energy < minimum:
                minimum = math.floor(energy)
            if energy > maximum:
                maximum = math.ceil(energy)
        if abs(minimum) > abs(maximum):
            maximum = abs(minimum)
        else:
            minimum = -abs(maximum)
        bf.write(f'remove (h. and (e. c extend 1))\n')
        bf.write(f'spectrum b, gmxc1 gmxc2 gmxc3 gmxc4 gmxc5, minimum={minimum}, maximum={maximum}\n')
        bf.write(f'ramp_new colorbar, none, [{minimum}, 0, {maximum}], [gmxc1, gmxc2, gmxc3, gmxc4, gmxc5]\n')
        bf.write(f'center {" or ".join(selections)}\n')
        bf.write(f'orient {" or ".join(selections)}\n')
        bf.write(f'zoom {" or ".join(selections)}, 15\n')


def ki2energy(ki, temp):
    # deltaG (inhibition) = R * T * ln ( Ki )
    if not ki:
        return np.nan
    return R * temp * math.log(ki * 1e-9)


def make_corr_DF(corr_data: dict) -> pd.DataFrame:
    data = []
    for x in corr_data:
        if x == 'mutant':
            continue
        for m in corr_data[x]['ΔG']:
            curr = [x] + [corr_data[x]['ΔG'][m][d] for d in corr_data[x]['ΔG'][m]] + [m, corr_data[x]['Exp.Energy']]
            data.append(curr)
    df = pd.DataFrame(data=data, columns=['System', 'ΔH', 'ΔH+IE', 'ΔH+NMODE', 'ΔH+QH', 'MODEL', 'Exp.Energy'])
    return df

def get_files(parser_args):
    info_files = []
    if not parser_args.files:
        parser_args.files = [Path('.')]
    for cf in parser_args.files:
        if not cf.exists():
            logging.warning(f'{cf} not found')
            continue
        if cf.is_dir():
            if parser_args.recursive:
                info_files.extend(cf.glob('*/*_info'))
            else:
                info_files.extend(cf.glob('*_info'))
        else:
            info_files.append(cf)
    if not len(info_files):
        GMXMMPBSA_ERROR('No info files found!')
    return info_files
