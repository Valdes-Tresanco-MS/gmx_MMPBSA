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

import sys
from os.path import split
import logging
from pathlib import Path

try:
    from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, InputError, CommandlineError
    from GMXMMPBSA.infofile import InfoFile
    from GMXMMPBSA import main
    from GMXMMPBSA.tester import run_test
    from GMXMMPBSA.commandlineparser import anaparser, testparser
    from GMXMMPBSA.utils import create_input_args
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))


def gmxmmpbsa():
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(levelname)-7s] %(message)s")
    stream_handler.setFormatter(formatter)
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("gmx_MMPBSA.log", 'w'),
            stream_handler])
    # Just for compatibility as mpi4py works as serial when run without mpirun
    # (since v1.4.2)
    if len(sys.argv) > 1 and sys.argv[1] in ['MPI', 'mpi']:
        from mpi4py import MPI
    else:
        # If we're not running "gmx_MMPBSA MPI", bring MPI into the top-level namespace
        # (which will overwrite the MPI from mpi4py, which we *want* to do in serial)
        from GMXMMPBSA.fake_mpi import MPI
    # Set up error/signal handlers
    main.setup_run()

    # Instantiate the main MMPBSA_App
    app = main.MMPBSA_App(MPI)

    # Read the command-line arguments
    try:
        app.get_cl_args(sys.argv[1:])
    except CommandlineError as e:
        sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
        sys.exit(1)

    if app.FILES.createinput is not None:
        args_list = create_input_args(app.FILES.createinput)
        app.input_file.print_contents('mmpbsa.in', args_list)
        logging.info(f'Input file creation successful. Path: {Path("mmpbsa.in").absolute()}')
        sys.exit(0)

    # Perform our MMPBSA --clean now
    if app.FILES.clean:
        logging.info('Cleaning temporary files and quitting.\n')
        app.remove(-1)
        sys.exit(0)

    # See if we wanted to print out our input file options
    if app.FILES.infilehelp:
        app.input_file.print_contents(sys.stdout)
        sys.exit(0)

    # If we're not rewriting output do whole shebang, otherwise load info and parms
    # Throw up a barrier before and after running the actual calcs
    if not app.FILES.rewrite_output:
        try:
            app.read_input_file()
        except InputError as e:
            sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
            sys.stderr.write('  Enter `%s --help` for help\n' %
                             (split(sys.argv[0])[1]))
            sys.exit(1)
        app.process_input()
        app.check_for_bad_input()
        app.make_prmtops()
        app.loadcheck_prmtops()
        app.file_setup()
        app.run_mmpbsa()
    # If we are rewriting output, load the info and check prmtops
    else:
        info = InfoFile(app)
        info.read_info()
        app.loadcheck_prmtops()

    # Now we parse the output, print, and finish
    app.parse_output_files()
    app.write_final_outputs()
    app.finalize()


def gmxmmpbsa_ana():
    try:
        from PyQt6.QtWidgets import QApplication
        pyqt = True
    except:
        try:
            from PyQt5.QtWidgets import QApplication
            pyqt = True
        except:
            pyqt = False
    finally:
        if not pyqt:
            GMXMMPBSA_ERROR('Could not import PyQt5/PyQt6. gmx_MMPBSA_ana will be disabled until PyQt5/PyQt6 is '
                            'installed')

    from GMXMMPBSA.analyzer.gui import GMX_MMPBSA_ANA
    from GMXMMPBSA.analyzer.utils import get_files

    app = QApplication(sys.argv)
    app.setApplicationName('gmx_MMPBSA Analyzer (gmx_MMPBSA_ana)')
    try:
        parser = anaparser.parse_args(sys.argv[1:])
    except CommandlineError as e:
        GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
        sys.exit(1)
    ifiles = get_files(parser)
    w = GMX_MMPBSA_ANA(ifiles)
    w.show()
    sys.exit(app.exec())


def gmxmmpbsa_test():
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("gmx_MMPBSA_test.log", 'w'),
            stream_handler])
    try:
        parser = testparser.parse_args(sys.argv[1:])
    except CommandlineError as e:
        GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
        sys.exit(1)
    run_test(parser)

if __name__ == '__main__':


    logging.info('Finished')

    # gmxmmpbsa()
    gmxmmpbsa_ana()