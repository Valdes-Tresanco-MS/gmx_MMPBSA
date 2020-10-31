"""
This module contains classes that back up the state of the MM/PBSA calculation
so that future post-processing of the output files can be done without
re-supplying all of the information again.
"""

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

import re
import warnings

class InfoFile(object):
    """
    This is the class responsible for backing up all of the calculation
    information and details so the output files can be re-written or the data
    can be readied for post-processing with external scripts.

    This object is passed a copy of the main app MMPBSA_App class
    """

    # Make a list about which INPUT variables are editable
    EDITABLE_INFO_VARS = ['debug_printlevel', 'verbose',
                          'csv_format', 'dec_verbose']

    def __init__(self, app):
        """ Instantiate me """
        self.app = app

    def write_info(self, name=None):
        """ Writes the info file to the info name """
        if name is None:
            name = self.app.pre + 'info'
        outfile = open(name, 'w')
        # The data we have to write: INPUT, FILES, and the following attributes:
        # numframes, numframes_nmode, mpi_size (also recognize size),
        # input_file_text, mut_str

        # Start with INPUT (and the editable vars). Allow this to recognize INFO
        # files from the last version of gmx_MMPBSA
        outfile.write('# You can alter the variables below to change what info '
                      'is printed out\n')
        for var in InfoFile.EDITABLE_INFO_VARS:
            outfile.write("INPUT['%s'] = %s\n" % (var,
                                                  self.write_var(self.app.INPUT[var])))
        outfile.write('#\n# MODIFY NOTHING BELOW HERE, OR GET WHAT YOU DESERVE\n')
        for var in list(self.app.INPUT.keys()):
            if var in InfoFile.EDITABLE_INFO_VARS: continue
            outfile.write("INPUT['%s'] = %s\n" % (var,
                                                  self.write_var(self.app.INPUT[var])))

        # Now print out the FILES
        for var in dir(self.app.FILES):
            # Skip over __method__ functions and output files
            if var.startswith('_') or var in ('rewrite_output', 'output_file',
                                              'decompout', 'energyout', 'dec_energies', 'overwrite'):
                continue
            outfile.write("FILES.%s = %s\n" % (var,
                                               self.write_var(getattr(self.app.FILES, var))))

        # Now print out the attributes we need to save
        outfile.write('size = %d\n' % self.app.mpi_size)
        outfile.write('numframes = %d\n' % self.app.numframes)
        outfile.write('numframes_nmode = %d\n' % self.app.numframes_nmode)
        outfile.write("mut_str = '%s'\n" % self.app.mut_str)
        outfile.write('using_chamber = %s\n' % self.app.using_chamber)
        outfile.write(self.app.input_file_text)

    def read_info(self, name=None):
        """ Reads a _MMPBSA_info file and populates the app """
        from GMXMMPBSA.commandlineparser import OptionList
        # Give self.app a FILES attribute if it does not have one yet before
        # trying to modify it
        if not hasattr(self.app, 'FILES'):
            self.app.FILES = OptionList()
        if name is None:
            name = self.app.pre + 'info'
        #TODO Add checking for differences in FILES
        inputre = re.compile(r'''INPUT\['(\S+)'\] = ((?:'.*')|(?:".*")|\S+)''')
        filesre = re.compile(r'''FILES\.(\S+) = (.*)''')
        otherre = re.compile(r'(\S+) = (\S+)')
        infile = open(name, 'r')
        input_text = ''
        for line in infile:
            if line.startswith('#'):
                continue
            # Now load each info file line
            if line.startswith('|'):
                input_text += line
                continue
            rematch = inputre.match(line)
            if rematch:
                var, val = rematch.groups()
                val = _determine_type(val)
                self.app.INPUT[var] = val
                continue
            rematch = filesre.match(line)
            if rematch:
                var, val = rematch.groups()
                val = _determine_type(val)
                setattr(self.app.FILES, var, val)
                continue
            rematch = otherre.match(line)
            if rematch:
                var, val = rematch.groups()
                val = _determine_type(val)
                if var == 'size':
                    self.app.mpi_size = val
                else:
                    setattr(self.app, var, val)
                continue
        # Determine stability here:
        self.app.stability = self.app.FILES.stability
        # Set app.pre as prefix
        self.app.pre = self.app.FILES.prefix
        if self.app.FILES.receptor_trajs or self.app.FILES.ligand_trajs:
            self.app.traj_protocol = 'MTP'
        else:
            self.app.traj_protocol = 'STP'
        # Load the input file text
        self.app.input_file_text = input_text

    def write_var(self, var):
        """
        Wrapper to return a string in which str vars are enclosed in quotes and
        numeric types (int and float) are not
        """
        if isinstance(var, str):
            return "'%s'" % var
        else:
            return "%s" % var


def _determine_type(thing):
    """ Determines what type this thing is """
    # If it is a string with quotes, strip off the quotes
    if thing.startswith('"') and thing.endswith('"'):
        return thing[1:len(thing)-1]
    if thing.startswith("'") and thing.endswith("'"):
        return thing[1:len(thing)-1]

    # Check for bool or None
    if thing == 'True':
        return True
    elif thing == 'False':
        return False
    elif thing == 'None':
        return None

    # Check for int, then check for float
    try:
        return int(thing)
    except ValueError:
        pass

    try:
        return float(thing)
    except ValueError:
        pass

    # Check for list
    if thing.startswith('[') and thing.endswith(']'):
        return eval(thing)

    # No idea what else it could be! Return string, but warn
    warnings.warn('Encountered unknown type in info file.\n' +
                  str(thing))
    return thing
