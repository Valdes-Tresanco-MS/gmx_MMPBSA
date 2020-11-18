"""
This module is responsible for finding the necessary programs for all of the
desired calculations in the gmx_MMPBSA calculation.

Methods:
         find_progs(INPUT): Determines which programs are needed and sees if
                            they can be found. Returns dictionary of programs
         which(program): Internal; searches AMBERHOME/bin and optionally
                         PATH to see if a program can be found

Classes: ExternProg: External program class
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

from GMXMMPBSA.exceptions import MMPBSA_Error

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def find_progs(INPUT):
    """ Find the necessary programs based in the user INPUT """
    # List all of the used programs with the conditions that they are needed
    used_progs = { 'cpptraj' : True, 'gmx': True,
                   # look for any available gromacs executable
                   'gmx_mpi': True, 'gmx_d': True, 'gmx_mpi_d': True,
                   # look for gromacs 4.x
                   'make_ndx': True, 'trjconv': True, 'tleap': True, 'parmchk2': True,
                   'mmpbsa_py_energy' : ((INPUT['pbrun'] or INPUT['gbrun'])
                                         and not (INPUT['use_sander'] or
                                                  INPUT['decomprun'])),
                   'sander' : (INPUT['decomprun'] or INPUT['use_sander'] or
                               INPUT['ifqnt'] == 1),
                   'sander.APBS' : INPUT['sander_apbs'] == 1,
                   'mmpbsa_py_nabnmode' : INPUT['nmoderun'],
                   'rism3d.snglpnt' : INPUT['rismrun']
                   }

    # The returned dictionary:
    my_progs = {}

    search_path = True
    force_path = INPUT['gmx_path']
    print(force_path)
    gromacs5x = True
    gromacs5x_d = True
    gromacs5x_mpi = True
    gromacs5x_mpi_d = True
    gromacs4x = True

    for prog in list(used_progs.keys()):
        my_progs[prog] = ExternProg(prog, used_progs[prog], search_path, force_path)
        if used_progs[prog]:
            if not my_progs[prog].full_path:
                if prog == 'gmx':
                    gromacs5x = False
                    continue
                if prog == 'gmx_d':
                    gromacs5x_d = False
                    continue
                if prog == 'gmx_mpi':
                    gromacs5x_mpi = False
                    continue
                if prog == 'gmx_mpi_d':
                    gromacs5x_mpi_d = False
                    continue
                if prog in ['make_ndx', 'trjconv']:
                    gromacs4x = False
                    continue
                raise MMPBSA_Error('Could not find necessary program [%s]' % prog)
            print('%s found! Using %s' % (prog, str(my_progs[prog])))
    if not (gromacs5x or gromacs5x_d or gromacs5x_mpi or gromacs5x_mpi_d) and not gromacs4x:
        raise MMPBSA_Error('Could not find necessary GROMACS program')

    return my_progs

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def which(program, search_path=False, force_path=None):
    """ Searches for a program in $AMBERHOME first, then PATH if we allow
        PATH searching.
    """
    import os

    def is_exe(filename):
        return os.path.exists(filename) and os.access(filename, os.X_OK)

    def get_amberhome():
        ambhome = os.getenv('AMBERHOME')
        if ambhome == None:
            raise MMPBSA_Error('AMBERHOME is not set!')
        return ambhome

    # Check to see that a path was provided in the program name
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program): return program

    # If not (and it generally isn't), then look in AMBERHOME
    amberhome = get_amberhome()

    if is_exe(os.path.join(amberhome, 'bin', program)):
        return os.path.join(amberhome, 'bin', program)

    if force_path:
        exe_file = os.path.join(force_path, program)
        if is_exe(exe_file):
            return exe_file  # if it's executable, return the file

    # If we can search the path, look through it
    if search_path:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file # if it's executable, return the file

    return None  # if program can still not be found... return None

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class ExternProg(object):
    """ An external program """
    def __init__(self, prog_name, needed=False, search_path=False, force_path=None):
        # We need to initialize the program name that we're looking for, whether
        # or not it is necessary by default, and the conditions upon which we
        # must look for it. The full_path is the attribute that stores this
        # program's full PATH name. We also need to know if we are going to search
        # through the whole PATH or just AMBERHOME.
        self.prog_name = prog_name
        self.full_path = None
        self.search_path = search_path
        self.force_path = force_path
        if needed: self.find()

    def find(self):
        """ Find the program """
        self.full_path = which(self.prog_name, self.search_path, self.force_path)

    def __str__(self):
        if self.full_path:
            return self.full_path
        else:
            return self.prog_name

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
