"""
 This package contains all of the functions for gmx_MMPBSA that it
 needs to run smoothly.
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

__all__ = ['alamdcrd', 'amber_outputs', 'calculation', 'commandlineparser',
           'createinput', 'exceptions', 'findprogs', 'fake_mpi', 'infofile',
           'input_parser', 'make_trajs', 'main', 'output_file', 'parm_setup',
           'timer', 'utils', 'API', 'make_top']

__author__ = "Mario S. Valdes Tresanco and Mario E. Valdes Tresanco"
__version__ = "1.1.1"
__license__ = "GPLv3"
__mmpbsa_author__ = "Jason Swails, Dwight McGee, and Bill Miller III"
__mmpbsa_version__ = "16.0"
