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

import os
import logging
import subprocess
import shutil
import multiprocessing
from pathlib import Path
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import time


def calculatestar(args):
    return run_process(*args)


def run_process(system, sys_name, args):
    time.sleep(0.1)
    logging.info(f"{system[1]:60}{'RUNNING':>10}")
    os.chdir(system[0])
    system_log = open(f'{sys_name}.log', 'a')
    g_p = subprocess.Popen(args, stdout=system_log, stderr=system_log)
    if g_p.wait():
        return sys_name, True
    return sys_name, False


def _get_frames(input_file: Path):
    startframe = 1
    endframe = 21
    with input_file.open() as ifile:
        for line in ifile:
            line = line.strip('\n').strip(',')
            if line.startswith('startframe'):
                startframe = int(line.split('=')[1])
            elif line.startswith('endframe'):
                endframe = int(line.split('=')[1])
    return 1 + endframe - startframe


def run_test(parser):
    # get the git path
    git_path = [os.path.join(path, 'git') for path in os.environ["PATH"].split(os.pathsep)
                if os.path.exists(os.path.join(path, 'git')) and
                os.access(os.path.join(path, 'git'), os.X_OK)][0]
    gmx_mmpbsa_path = [os.path.join(path, 'gmx_MMPBSA') for path in os.environ["PATH"].split(os.pathsep)
                       if os.path.exists(os.path.join(path, 'gmx_MMPBSA')) and
                       os.access(os.path.join(path, 'gmx_MMPBSA'), os.X_OK)][0]
    gmx_mmpbsa_ana_path = [os.path.join(path, 'gmx_MMPBSA_ana') for path in os.environ["PATH"].split(os.pathsep)
                           if os.path.exists(os.path.join(path, 'gmx_MMPBSA_ana')) and
                           os.access(os.path.join(path, 'gmx_MMPBSA_ana'), os.X_OK)][0]
    if not git_path:
        GMXMMPBSA_ERROR('Git not found. Please install git ( sudo apt install git ) or make sure git is in the PATH '
                        'and try again...')
    if not gmx_mmpbsa_path or not gmx_mmpbsa_ana_path:
        GMXMMPBSA_ERROR('Please make sure gmx_MMPBSA and gmx_MMPBSA_ana are in the PATH and try again...')

    if not parser.folder.exists():
        GMXMMPBSA_ERROR(f'{parser.folder} not exists or is inaccessible. Please define a new folder and try again...')

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
        git_p = subprocess.Popen(['git', 'clone', 'https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
                                  gmx_mmpbsa_test_folder.as_posix()])
        if git_p.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('git failed when try to clone the gmx_MMPBSA repository')

        logging.info('Cloning gmx_MMPBSA repository...Done.')

    examples = gmx_mmpbsa_test_folder.joinpath('examples')
    test_sys = {
        3: [examples.joinpath('Protein_ligand', 'ST'), 'Protein-Ligand (Single trajectory approximation)'],
        4: [examples.joinpath('Protein_protein'), 'Protein-Protein'],
        5: [examples.joinpath('Protein_DNA'), 'Protein-DNA'],
        6: [examples.joinpath('Protein_membrane'), 'Protein-Membrane'],
        7: [examples.joinpath('Protein_glycan'), 'Protein-Glycan'],
        8: [examples.joinpath('Metalloprotein_peptide'), 'Metalloprotein-Peptide'],
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
        23: [examples.joinpath('QM_MMGBSA'), 'QM/MMGBSA Calculation']
    }

    if parser.test == [0]:
        key_list = range(3, 19)
    elif parser.test == [1]:
        key_list = [3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15]
    elif parser.test == [2]:
        key_list = [3, 4, 5, 7, 9, 12, 13, 14, 15]
    elif parser.test == [101]:
        key_list = [3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    else:
        key_list = parser.test

    req_cpus = {x: _get_frames(test_sys[x][0].joinpath('mmpbsa.in')) for x in key_list}

    if parser.num_processors > multiprocessing.cpu_count():
        logging.warning(f'The number cpus defined {parser.num_processors} is greater than the system cpu'
                        f' {multiprocessing.cpu_count()}. All the cpus will be used...')
        parser.num_processors = multiprocessing.cpu_count()

    TASKS = []
    for x in key_list:
        with open(test_sys[x][0].joinpath('README.md')) as readme:
            for line in readme:
                if 'gmx_MMPBSA -O -i mmpbsa.in' in line:
                    command = (['mpirun', '-np',
                                f'{req_cpus[x] if req_cpus[x] <= parser.num_processors else parser.num_processors}']
                               + [gmx_mmpbsa_path, 'MPI'] + line.strip('\n').split()[1:] + ['-nogui'])
                    TASKS.append((test_sys[x], x, command))

    result_list = []
    logging.info(f"{'Example':^60}{'STATE':>10}")
    print(80 * '-')
    exitcode = 0
    c = 1
    with multiprocessing.Pool(1) as pool:
        imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)
        for x in imap_unordered_it:
            sys_name, exitcode = x
            if exitcode:
                logging.error(f"{test_sys[sys_name][1]:55}[{c:2}/{len(key_list):2}]{'ERROR':>8}\n"
                              f"           Please, check the test log\n"
                              f"           ({test_sys[sys_name][0].joinpath(f'{sys_name}')}.log)")
            else:
                logging.info(f"{test_sys[sys_name][1]:55}[{c:2}/{len(key_list):2}]{'DONE':>8}")
                result_list.append(test_sys[sys_name][0])

            c += 1

    if not parser.nogui and not exitcode:
        print(80 * '-')
        logging.info('Opening gmx_MMPBSA_ana...')
        g_p = subprocess.Popen([gmx_mmpbsa_ana_path, '-f'] + result_list, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        if g_p.wait():
            error = g_p.stderr.read().decode("utf-8")
            import sys
            sys.stderr.write(error)
