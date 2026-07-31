# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdés-Tresanco and Mario E. Valdés-Tresanco   #
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

import os
import logging
import subprocess
import shutil
import multiprocessing
import shlex
import sys
import re
from pathlib import Path
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import time


def calculatestar(args):
    return run_process(*args)


def _parse_test_entry(entry):
    path, name = entry[0], entry[1]
    options = entry[2] if len(entry) > 2 else {}
    return path, name, options


def run_process(work_dir, display_name, sys_name, args, log_file):
    time.sleep(0.1)
    logging.info(f"{display_name:60}{'RUNNING':>10}")
    os.chdir(work_dir)
    with open(log_file, 'a') as system_log:
        g_p = subprocess.Popen(args, stdout=system_log, stderr=system_log)
        if g_p.wait():
            return sys_name, True
    return sys_name, False


def _has_known_rism_fortran_runtime_error(log_file: Path):
    try:
        log_text = log_file.read_text(errors='replace')
    except OSError:
        return False
    return (
        'Fortran runtime error: Missing comma between descriptors' in log_text
        and 'amber_rism_interface.F90' in log_text
    )


def _get_frames(input_file: Path):
    values = {'startframe': 1, 'endframe': 21, 'interval': 1}
    assignment = re.compile(r'\b(startframe|endframe|interval)\s*=\s*(\d+)\b')
    with input_file.open() as ifile:
        for line in ifile:
            line = line.split('#', 1)[0]
            for key, value in assignment.findall(line):
                values[key] = int(value)
    if values['interval'] < 1:
        GMXMMPBSA_ERROR(f'Invalid interval in {input_file}. It must be greater than 0')
    if values['endframe'] < values['startframe']:
        GMXMMPBSA_ERROR(f'Invalid frame range in {input_file}. endframe must be greater than or equal to startframe')
    return ((values['endframe'] - values['startframe']) // values['interval']) + 1


def _find_executable(executable: str):
    exe_path = shutil.which(executable)
    if exe_path is None:
        GMXMMPBSA_ERROR(f'Please make sure {executable} is in the PATH and try again...')
    return exe_path


def run_test(parser):
    git_path = _find_executable('git')

    if not parser.folder.exists():
        GMXMMPBSA_ERROR(f'{parser.folder} does not exist or is inaccessible. Please define a new folder and try again...')

    gmx_mmpbsa_test_folder = parser.folder.joinpath('gmx_MMPBSA_test').absolute()

    if not gmx_mmpbsa_test_folder.exists() and parser.reuse:
        GMXMMPBSA_ERROR(f'The examples directory {gmx_mmpbsa_test_folder.exists()} does not exist. To use the -r '
                        f'option you must first have cloned the repository')
    if not gmx_mmpbsa_test_folder.exists():
        clonning = True
    elif gmx_mmpbsa_test_folder.exists() and not parser.reuse:
        shutil.rmtree(gmx_mmpbsa_test_folder)
        clonning = True
    else:
        clonning = False

    if clonning:
        logging.info(f'Cloning gmx_MMPBSA repository in {gmx_mmpbsa_test_folder}')
        git_p = subprocess.Popen(['git', 'clone', '--depth', '1', 'https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
                                  gmx_mmpbsa_test_folder.as_posix()])
        if git_p.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('git failed when trying to clone the gmx_MMPBSA repository')

        logging.info('Cloning gmx_MMPBSA repository...Done.')

    examples = gmx_mmpbsa_test_folder.joinpath('examples')
    test_sys = {
        3: [examples.joinpath('Protein_ligand', 'ST'), 'Protein-Ligand (Single trajectory approximation)'],
        4: [examples.joinpath('Protein_protein'), 'Protein-Protein'],
        5: [examples.joinpath('Protein_DNA'), 'Protein-DNA'],
        6: [examples.joinpath('Protein_membrane'), 'Protein-Membrane'],
        7: [examples.joinpath('Protein_glycan'), 'Protein-Glycan'],
        8: [examples.joinpath('Metalloprotein_ligand'), 'Metalloprotein-ligand'],
        9: [examples.joinpath('Comp_receptor'), 'Comp_receptor'],
        10: [examples.joinpath('Protein_ligand_CHARMMff'), 'Protein-Ligand (CHARMM force field)'],
        11: [examples.joinpath('Protein_membrane_CHARMMff'), 'Protein-ligand complex in membrane with CHARMMff'],
        12: [examples.joinpath('Alanine_scanning'), 'Alanine Scanning'],
        13: [examples.joinpath('Stability'), 'Stability calculation'],
        14: [examples.joinpath('Decomposition_analysis'), 'Decomposition Analysis'],
        15: [examples.joinpath('Entropy_calculations', 'Interaction_Entropy'), 'Interaction Entropy approximation'],
        16: [examples.joinpath('Protein_ligand', 'MT'), 'Protein-Ligand (Multiple trajectory approximation)'],
        17: [examples.joinpath('Entropy_calculations', 'nmode'), 'Entropy calculation using Normal Mode '
                                                                 'approximation'],
        18: [examples.joinpath('3D-RISM'), 'Calculations using 3D-RISM approximation'],
        19: [examples.joinpath('Entropy_calculations', 'C2_Entropy'), 'C2 Entropy approximation'],
        20: [examples.joinpath('Linear_PB_solver'), 'LPB Calculation'],
        21: [examples.joinpath('NonLinear_PB_solver'), 'NLPB Calculation'],
        22: [examples.joinpath('Protein_ligand_LPH_atoms_CHARMMff'), 'Protein-Ligand_LPH (CHARMM force field)'],
        23: [examples.joinpath('QM_MMGBSA'), 'QM/MMGBSA Calculation'],
        24: [examples.joinpath('GBNSR6'), 'GBNSR6 Calculation'],
        25: [examples.joinpath('AMBER'), 'AMBER input files'],
        26: [examples.joinpath('Explicit_receptor_waters'), 'ST MM/PB(GB)SA with explicit receptor waters']
    }

    if parser.test == [0]:
        key_list = list(range(3, 27))
    elif parser.test == [1]:
        key_list = [3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15]
    elif parser.test == [2]:
        key_list = [3, 4, 5, 7, 9, 12, 13, 14, 15]
    elif parser.test == [101]:
        key_list = list(range(3, 27))
    else:
        key_list = parser.test
    if not key_list:
        GMXMMPBSA_ERROR('No test was selected. Please define at least one test number')

    test_meta = {}
    req_cpus = {}
    for x in key_list:
        path, name, options = _parse_test_entry(test_sys[x])
        input_file = options.get('input_file', path.joinpath('mmpbsa.in'))
        work_dir = options.get('cwd', path)
        req_cpus[x] = _get_frames(input_file)
        test_meta[x] = {
            'path': path,
            'name': name,
            'work_dir': work_dir,
            'log_file': work_dir.joinpath(f'{x}.log'),
        }

    if parser.num_processors > multiprocessing.cpu_count():
        logging.warning(f'The number cpus defined {parser.num_processors} is greater than the system cpu'
                        f' {multiprocessing.cpu_count()}. All the cpus will be used...')
        parser.num_processors = multiprocessing.cpu_count()

    if parser.num_concurrent < 1:
        GMXMMPBSA_ERROR('The number of concurrent examples must be greater than 0')
    if parser.num_concurrent > len(key_list):
        parser.num_concurrent = len(key_list)
    if parser.num_processors * parser.num_concurrent > multiprocessing.cpu_count():
        logging.warning(f'The requested test concurrency can use up to '
                        f'{parser.num_processors * parser.num_concurrent} MPI ranks across '
                        f'{parser.num_concurrent} examples, which is greater than the system cpu '
                        f'{multiprocessing.cpu_count()}. Consider reducing -n or -j...')

    TASKS = []
    for x in key_list:
        path, name, options = _parse_test_entry(test_sys[x])
        meta = test_meta[x]
        command_line = options.get('command')
        if not command_line:
            with open(path.joinpath('README.md')) as readme:
                for line in readme:
                    if 'gmx_MMPBSA -O -i mmpbsa.in' in line or 'amber_MMPBSA -O -i mmpbsa.in' in line:
                        command_line = line.strip()
                        break
        if not command_line:
            GMXMMPBSA_ERROR(f'No runnable command found for test {x} ({name})')

        executable = _find_executable('amber_MMPBSA' if 'amber_MMPBSA' in command_line else 'gmx_MMPBSA')
        command = (['mpirun', '-np',
                    f'{req_cpus[x] if req_cpus[x] <= parser.num_processors else parser.num_processors}']
                   + [executable] + shlex.split(command_line)[1:] + ['-nogui'])
        TASKS.append((meta['work_dir'], name, x, command, meta['log_file']))

    result_list = []
    logging.info(f"{'Example':^60}{'STATE':>10}")
    print(80 * '-')
    any_failed = False
    c = 1
    with multiprocessing.Pool(parser.num_concurrent) as pool:
        imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)
        for x in imap_unordered_it:
            sys_name, failed = x
            if failed:
                any_failed = True
                log_file = test_meta[sys_name]['log_file']
                logging.error(f"{test_meta[sys_name]['name']:55}[{c:2}/{len(key_list):2}]{'ERROR':>8}\n"
                              f"           Please, check the test log\n"
                              f"           ({log_file})")
                if sys_name == 18 and _has_known_rism_fortran_runtime_error(log_file):
                    logging.warning('The 3D-RISM log matches a known AmberTools/Fortran runtime issue. '
                                    'Conda AmberTools builds linked with newer libgfortran can abort before '
                                    'the RISM calculation starts. A known working workaround is gmx_MMPBSA 1.6.4 '
                                    'with Python 3.9/3.10, AmberTools 23, and libgfortran5/libgcc-ng 12.x, or a '
                                    'patched AmberTools build.')
            else:
                logging.info(f"{test_meta[sys_name]['name']:55}[{c:2}/{len(key_list):2}]{'DONE':>8}")
                result_list.append(test_meta[sys_name]['work_dir'])

            c += 1

    if any_failed:
        sys.exit(1)

    if not parser.nogui:
        gmx_mmpbsa_ana_path = _find_executable('gmx_MMPBSA_ana')
        print(80 * '-')
        logging.info('Opening gmx_MMPBSA_ana...')
        g_p = subprocess.Popen([gmx_mmpbsa_ana_path, '-f'] + result_list, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        if g_p.wait():
            error = g_p.stderr.read().decode("utf-8")
            sys.stderr.write(error)
