# !/home/mario/programs/amber18/bin/amber.python

"""
MMPBSA.py is a script for performing (M)olecular (M)echanics
(P)oisson (B)oltzmann (S)urface (A)rea to find free energies of
binding. Refer to the AMBER manual and/or relevant literature for a
more thorough overview of the method. This implementation uses AMBER
executables to find energies using either Poisson Boltzmann or
Generalized Born implicit solvent models of a complex of a receptor
with a bound ligand. This script was written by Dwight McGee,
Billy Miller III, and Jason Swails in Adrian Roitberg's research
group at the Quantum Theory Project at the University of Florida.

Last updated: 12/03/2012

                           GPL LICENSE INFO

Copyright (C) 2009 - 2011  Dwight McGee, Billy Miller III, and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

import sys
from os.path import split
try:
    from GMXMMPBSA.exceptions import MMPBSA_Error, InputError, CommandlineError
    from GMXMMPBSA.infofile import InfoFile
    from GMXMMPBSA import main
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))

# Set up the MPI version or not based on whether our executable name ends in MPI
if sys.argv[0].endswith('MPI'):
    try:
        from mpi4py import MPI
    except ImportError:
        raise MMPBSA_Error('Could not import mpi4py package! Use serial version '
                           'or install mpi4py.')
else:
    # If we're not running MMPBSA.py.MPI, bring MPI into the top-level namespace
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
    sys.stderr.write('%s: %s' % (type(e).__name__, e )+ '\n')
    sys.exit(1)

# Perform our MMPBSA --clean now
if app.FILES.clean:
    sys.stdout.write('Cleaning temporary files and quitting.\n')
    app.remove(0)
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
        sys.stderr.write('%s: %s' % (type(e).__name__, e )+ '\n')
        sys.stderr.write('  Enter `%s --help` for help\n' %
                         (split(sys.argv[0])[1]))
        sys.exit(1)
    app.process_input()
    app.check_for_bad_input()
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