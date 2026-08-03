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
import sys
import re
from pathlib import Path
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
from GMXMMPBSA.test_manifest import load_manifest, snapshot_outputs, verify_outputs
import time


def calculatestar(args):
    return run_process(*args)


def run_process(work_dir, display_name, sys_name, args, log_file, expected_outputs, skip_output_check):
    time.sleep(0.1)
    logging.info(f"{display_name:60}{'RUNNING':>10}")
    os.chdir(work_dir)
    output_snapshot = snapshot_outputs(work_dir, expected_outputs)
    with open(log_file, 'a') as system_log:
        g_p = subprocess.Popen(args, stdout=system_log, stderr=system_log)
        if g_p.wait():
            return sys_name, True, []

    if skip_output_check:
        return sys_name, False, []

    missing = verify_outputs(work_dir, expected_outputs, output_snapshot)
    return sys_name, bool(missing), missing


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


def _resolve_examples_dir(parser):
    env_dir = os.getenv('GMXMMPBSA_TEST_EXAMPLES_DIR')
    if parser.examples_dir:
        examples = Path(parser.examples_dir).expanduser().resolve()
        if not examples.is_dir():
            GMXMMPBSA_ERROR(f'{examples} does not exist or is inaccessible.')
        return examples, False
    if env_dir:
        examples = Path(env_dir).expanduser().resolve()
        if not examples.is_dir():
            GMXMMPBSA_ERROR(f'{examples} from GMXMMPBSA_TEST_EXAMPLES_DIR does not exist or is inaccessible.')
        return examples, False

    if not parser.folder.exists():
        GMXMMPBSA_ERROR(f'{parser.folder} does not exist or is inaccessible. Please define a new folder and try again...')

    gmx_mmpbsa_test_folder = parser.folder.joinpath('gmx_MMPBSA_test').absolute()

    if not gmx_mmpbsa_test_folder.exists() and parser.reuse:
        GMXMMPBSA_ERROR(f'The examples directory {gmx_mmpbsa_test_folder} does not exist. To use the -r '
                        f'option you must first have cloned the repository')
    if not gmx_mmpbsa_test_folder.exists():
        clonning = True
    elif gmx_mmpbsa_test_folder.exists() and not parser.reuse:
        shutil.rmtree(gmx_mmpbsa_test_folder)
        clonning = True
    else:
        clonning = False

    if clonning:
        _find_executable('git')
        logging.info(f'Cloning gmx_MMPBSA repository in {gmx_mmpbsa_test_folder}')
        git_p = subprocess.Popen(['git', 'clone', '--depth', '1', 'https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
                                  gmx_mmpbsa_test_folder.as_posix()])
        if git_p.wait():
            GMXMMPBSA_ERROR('git failed when trying to clone the gmx_MMPBSA repository')
        logging.info('Cloning gmx_MMPBSA repository...Done.')

    examples = gmx_mmpbsa_test_folder.joinpath('examples')
    return examples, True


def _preflight_executables(test_entries, clone_mode, open_analyzer):
    required = set()
    for test in test_entries:
        required.update(test.requires)
    if clone_mode:
        required.add('git')
    if open_analyzer:
        required.add('gmx_MMPBSA_ana')

    missing = [exe for exe in sorted(required) if shutil.which(exe) is None]
    if missing:
        GMXMMPBSA_ERROR('Missing executables for selected tests: ' + ', '.join(missing))


def run_test(parser):
    manifest = load_manifest()
    selectors = [str(value) for value in parser.test]
    try:
        key_list = manifest.resolve_test_ids(selectors)
    except ValueError as exc:
        GMXMMPBSA_ERROR(str(exc))

    examples, clone_mode = _resolve_examples_dir(parser)
    test_entries = [manifest.get_test(test_id) for test_id in key_list]
    _preflight_executables(test_entries, clone_mode, not parser.nogui)

    for test in test_entries:
        if 'rism_fortran' in test.known_issues:
            logging.warning('Test 18 uses AmberTools 3D-RISM and may fail on some AmberTools/Fortran runtime builds.')

    test_meta = {}
    req_cpus = {}
    for test_id in key_list:
        test = manifest.get_test(test_id)
        work_dir = examples.joinpath(test.workdir)
        input_file = work_dir.joinpath(test.input)
        req_cpus[test_id] = _get_frames(input_file)
        test_meta[test_id] = {
            'name': test.name,
            'work_dir': work_dir,
            'log_file': work_dir.joinpath(f'{test_id}.log'),
            'expected_outputs': test.expected_outputs,
            'known_issues': test.known_issues,
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
    for test_id in key_list:
        test = manifest.get_test(test_id)
        meta = test_meta[test_id]
        executable = _find_executable(test.executable)
        np = req_cpus[test_id] if req_cpus[test_id] <= parser.num_processors else parser.num_processors
        command = (['mpirun', '-np', f'{np}', executable]
                   + test.command_args
                   + ['-nogui'])
        TASKS.append((
            meta['work_dir'],
            meta['name'],
            test_id,
            command,
            meta['log_file'],
            meta['expected_outputs'],
            parser.skip_output_check,
        ))

    result_list = []
    logging.info(f"{'Example':^60}{'STATE':>10}")
    print(80 * '-')
    any_failed = False
    c = 1
    with multiprocessing.Pool(parser.num_concurrent) as pool:
        imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)
        for test_id, failed, missing_outputs in imap_unordered_it:
            if failed:
                any_failed = True
                log_file = test_meta[test_id]['log_file']
                logging.error(f"{test_meta[test_id]['name']:55}[{c:2}/{len(key_list):2}]{'ERROR':>8}\n"
                              f"           Please, check the test log\n"
                              f"           ({log_file})")
                if missing_outputs:
                    logging.error('           Missing expected output files: ' + ', '.join(missing_outputs))
                if test_id == 18 and _has_known_rism_fortran_runtime_error(log_file):
                    logging.warning('The 3D-RISM log matches a known AmberTools/Fortran runtime issue. '
                                    'Conda AmberTools builds linked with newer libgfortran can abort before '
                                    'the RISM calculation starts. A known working workaround is gmx_MMPBSA 1.6.4 '
                                    'with Python 3.9/3.10, AmberTools 23, and libgfortran5/libgcc-ng 12.x, or a '
                                    'patched AmberTools build.')
            else:
                logging.info(f"{test_meta[test_id]['name']:55}[{c:2}/{len(key_list):2}]{'DONE':>8}")
                result_list.append(test_meta[test_id]['work_dir'])

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
