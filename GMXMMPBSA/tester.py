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
from queue import Queue
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR


def calculatestar(args):
    return run_process(*args)

def run_process(system, sys_name, args):
    logging.info(f"{system[1]:60}{'RUNNING':>10}")
    os.chdir(system[0])
    system_log = open(sys_name + '.log', 'a')
    g_p = subprocess.Popen(args, stdout=system_log, stderr=system_log)
    if g_p.wait():
        return sys_name, False
    return sys_name, True


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

    gmx_mmpbsa_test_folder = parser.folder.joinpath('gmx_MMPBSA_test')
    if gmx_mmpbsa_test_folder.exists():
        shutil.rmtree(gmx_mmpbsa_test_folder)
    logging.info(f'Cloning gmx_MMPBSA repository in {gmx_mmpbsa_test_folder}')
    git_p = subprocess.Popen(['git', 'clone', 'https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
                              gmx_mmpbsa_test_folder.as_posix()])
    if git_p.wait():  # if it quits with return code != 0
        GMXMMPBSA_ERROR('git failed when try to clone the gmx_MMPBSA repository')

    logging.info(f'Cloning gmx_MMPBSA repository...Done.')

    examples = gmx_mmpbsa_test_folder.joinpath('docs', 'examples')
    test = {
        'prot_lig_st': [examples.joinpath('Protein_ligand', 'ST'), 'Protein-Ligand (Single trajectory approximation)'],
        'prot_prot': [examples.joinpath('Protein_protein'), 'Protein-Protein'],
        'prot_dna': [examples.joinpath('Protein_DNA'), 'Protein-DNA'],
        'memb_prot': [examples.joinpath('Protein_membrane'), 'Protein-Membrane'],
        'prot_glycan': [examples.joinpath('Protein_glycan'), 'Protein-Glycan'],
        'metalloprot_pep': [examples.joinpath('Metalloprotein_peptide'), 'Metalloprotein-Peptide'],
        'prot_dna_rna_ions_lig': [examples.joinpath('Protein_DNA_RNA_Ion_ligand'), 'Protein-DNA-RNA-IONs-Ligand'],
        'prot_lig_charmm': [examples.joinpath('Protein_ligand_CHARMMff'), 'Protein-Ligand (CHARMM force field)'],
        'ala_scan': [examples.joinpath('Alanine_scanning'), 'Alanine Scanning'],
        'stability': [examples.joinpath('Stability'), 'Stability calculation'],
        'decomp': [examples.joinpath('Decomposition_analysis'), 'Decomposition Analysis'],
        'prot_lig_mt': [examples.joinpath('Protein_ligand', 'MT'), 'Protein-Ligand (Multiple trajectory '
                                                                   'approximation)'],
        'ie': [examples.joinpath('Entropy_calculations', 'Interaction_Entropy'), 'Interaction Entropy approximation'],
        'nmode': [examples.joinpath('Entropy_calculations', 'nmode'), 'Entropy calculation using Normal Mode '
                                                                      'approximation '],
        '3drism': [examples.joinpath('3D-RISM'), 'Calculations using 3D-RISM approximation']
    }
    all = ['prot_lig_st', 'prot_prot', 'prot_dna', 'memb_prot', 'prot_glycan', 'metalloprot_pep',
           'prot_dna_rna_ions_lig', 'prot_lig_charmm', 'ala_scan', 'stability', 'decomp', 'prot_lig_mt', 'ie',
           'nmode', '3drism']
    minimal = ['prot_lig_st', 'prot_prot', 'prot_dna', 'memb_prot', 'prot_glycan', 'metalloprot_pep',
               'prot_dna_rna_ions_lig', 'prot_lig_charmm', 'ala_scan', 'stability', 'decomp', 'ie']

    if parser.test == 'all':
        key_list = all.copy()
    elif parser.test == 'minimal':
        key_list = minimal.copy()
    else:
        key_list = [parser.test]

    # Create tasks
    TASKS = []
    for x in key_list:
        with open(test[x][0].joinpath('README.md')) as readme:
            for line in readme:
                if 'gmx_MMPBSA -O -i mmpbsa.in' in line:
                    command = line.strip('\n').split() + ['-nogui']
                    TASKS.append((test[x], x, command))

    if parser.num_processors > multiprocessing.cpu_count():
        logging.warning(f'Using all processors')
        parser.num_processors = multiprocessing.cpu_count()

    if len(key_list) < parser.num_processors:
        jobs = len(key_list)
    else:
        jobs = parser.num_processors

    result_list = []

    logging.info(f"{'Example':^60}{'STATE':>10}")
    print(80*'-')
    c = 1
    with multiprocessing.Pool(jobs) as pool:
        imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)
        for x in imap_unordered_it:
            sys_name, result = x
            if result:
                logging.info(f"{test[sys_name][1]:55}[{c:2}/{len(key_list):2}]{'DONE':>8}")
                result_list.append(test[sys_name][0])
            else:
                logging.error(f"{test[sys_name][1]:55}[{c:2}/{len(key_list):2}]{'ERROR':>8}\n"
                              f"           Please, check the test log\n"
                              f"           ({test[sys_name][0].joinpath(f'{sys_name}')}.log)")
            c += 1
    if not parser.nogui:
        print(80 * '-')
        logging.info('Opening gmx_MMPBSA_ana...')
        g_p = subprocess.Popen(['gmx_MMPBSA_ana', '-f'] + result_list, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        if g_p.wait():
            error = g_p.stderr.read().decode("utf-8")
            import sys
            sys.stderr.write(error)
