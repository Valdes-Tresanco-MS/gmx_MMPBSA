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
import sys
from os.path import split
import logging
from pathlib import Path

try:
    from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, InputError, CommandlineError
    from GMXMMPBSA.infofile import InfoFile
    from GMXMMPBSA import main
    from GMXMMPBSA.tester import run_test
    from GMXMMPBSA.commandlineparser import anaparser, testparser, amber_parser
    from GMXMMPBSA.error_bundle import create_error_bundle
    from GMXMMPBSA.utils import create_input_args
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))


class WarningSpacingFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._previous_was_warning = False

    def format(self, record):
        message = super().format(record)
        if record.levelno == logging.WARNING:
            prefix = '' if self._previous_was_warning else '\n'
            self._previous_was_warning = True
            return f'{prefix}{message}\n'
        self._previous_was_warning = False
        return message


def _setup_logging(log_file):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
    file_handler = logging.FileHandler(log_file, 'w')
    file_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
    logging.basicConfig(level=logging.DEBUG, handlers=[file_handler, stream_handler], force=True)


def _gmxmmpbsa_base(parser, engine='gmx'):
    _setup_logging("gmx_MMPBSA.log")
    # Just for compatibility as mpi4py works as serial when run without mpirun
    # (since v1.4.2)
    from mpi4py import MPI
    if MPI.COMM_WORLD.size == 1:
        from GMXMMPBSA.fake_mpi import MPI
    # Set up error/signal handlers
    main.setup_run()

    # Instantiate the main MMPBSA_App
    app = main.MMPBSA_App(MPI)
    app.clparser = parser
    app.engine = engine

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

    try:
        # If we're not rewriting output do whole shebang, otherwise load info and parms
        # Throw up a barrier before and after running the actual calcs
        if not app.FILES.rewrite_output:
            try:
                app.read_input_file()
            except InputError as e:
                _maybe_create_error_bundle(app, e)
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
            info = InfoFile(app, True)
            info.read_info()
            app.loadcheck_prmtops()

        # Now we parse the output, print, and finish
        app.parse_output_files()
        app.write_final_outputs()
        app.finalize()
    except SystemExit:
        raise
    except Exception as e:
        _maybe_create_error_bundle(app, e)
        raise


def _maybe_create_error_bundle(app, exc):
    if not getattr(app, 'master', True):
        return
    if getattr(app, '_error_bundle_created', False):
        return
    files = getattr(app, 'FILES', None)
    if files is not None and getattr(files, 'no_error_bundle', False):
        return
    app._error_bundle_created = True
    try:
        bundle = create_error_bundle(app, exc)
    except Exception as bundle_exc:
        sys.stderr.write('\nCould not create gmx_MMPBSA error bundle: %s\n' % bundle_exc)
        return
    sys.stderr.write(
        '\nA gmx_MMPBSA error bundle was created for debugging:\n'
        '  %s\n\n'
        'Please attach this zip file when reporting the issue. It contains logs,\n'
        'input/setup files, generated intermediates, and up to 5 trajectory frames.\n\n' % bundle
    )


def gmxmmpbsa():
    from GMXMMPBSA.commandlineparser import parser
    _gmxmmpbsa_base(parser)

def gmxmmpbsa_amber():
    _gmxmmpbsa_base(amber_parser, 'amber')

def gmxmmpbsa_ana():
    try:
        from PyQt6.QtWidgets import QApplication
        pyqt = True
    except ImportError:
        try:
            from PyQt5.QtWidgets import QApplication
            pyqt = True
        except ImportError:
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
    _setup_logging("gmx_MMPBSA_test.log")
    try:
        from GMXMMPBSA.test_manifest import build_help_text
        from GMXMMPBSA.commandlineparser import test_action
        test_action.help = build_help_text()
        parser = testparser.parse_args(sys.argv[1:])
        if parser.examples_source == 'local' and not parser.examples_dir and not os.getenv('GMXMMPBSA_TEST_EXAMPLES_DIR'):
            testparser.error('--examples-source local requires --examples-dir or GMXMMPBSA_TEST_EXAMPLES_DIR')
    except CommandlineError as e:
        GMXMMPBSA_ERROR('%s: %s' % (type(e).__name__, e))
        sys.exit(1)
    run_test(parser)

if __name__ == '__main__':


    logging.info('Finished')

    # gmxmmpbsa()
    gmxmmpbsa_ana()
