"""
This is a module that contains the exceptions thrown by gmx_MMPBSA
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

class MMPBSA_Error(Exception):
    """ Base gmx_MMPBSA error class """

    def __init__(self, msg='gmx_MMPBSA error'):
        self.msg = msg

    def __str__(self):
        return self.msg


class MMPBSA_Warning(Warning):
    """ Base gmx_MMPBSA warning class """

    def __init__(self, msg='gmx_MMPBSA warning'):
        self.msg = msg

    def __str__(self):
        return self.msg


class TimerException(MMPBSA_Error):
    """ Error in timer module """
    pass


class CommandlineError(MMPBSA_Error):
    """ Error parsing the command-line """
    pass


class CalcError(MMPBSA_Error):
    """ Error when running calculations """
    pass


class SetupError(MMPBSA_Error):
    """ Error in standard setup; i.e. environment variables """
    pass


class PrmtopError(MMPBSA_Error):
    """ Error in one of the prmtops """
    pass


class SelectionError(MMPBSA_Error):
    """ Error in which a residue selection is illegal """
    pass


class TrajError(MMPBSA_Error):
    """ Error in trajectory processing """
    pass


class InputError(MMPBSA_Error):
    """ Error in the Input File """
    pass


class IllegalOption(MMPBSA_Error):
    """ Error captured when looking for incompatibilities """
    pass


class InterruptError(MMPBSA_Error):
    """ When the process is interrupted """
    pass


class OutputError(MMPBSA_Error):
    """ Error in parsing the output """
    pass


class InconsistentError1(MMPBSA_Error):
    """ Error when internal potential terms are inconsistent. Specifically
        BOND, ANGLE, and DIHED terms
    """
    pass


class InconsistentError2(MMPBSA_Error):
    """ Error when internal potential terms are inconsistent. Specifically
        1-4 VDW, and 1-4 EEL terms
    """
    pass


class ConvergenceError(MMPBSA_Error):
    """ Error when nmode frame doesn't minimize within tolerance """
    pass


class MutateError(MMPBSA_Error):
    """ Error mutating a trajectory in alanine scanning """
    pass


class MutantResError(MutateError):
    """ Error if we're mutating to a bad residue """
    pass


class CreateInputError(MMPBSA_Error):
    """ Error creating MDIN file """
    pass


class LengthError(MMPBSA_Error):
    """
    Error that occurs when trying to subtract different length EnergyVectors
    """
    pass


class DecompError(MMPBSA_Error):
    """ Error parsing decomp results """
    pass


class InternalError(MMPBSA_Error):
    """ Error from buggy coding """
    pass


class NoFileExists(MMPBSA_Error):
    """ Error if we don't have a file we need """
    pass


class InputWarning(MMPBSA_Warning):
    """ If we have a non-fatal warning """
    pass

class StabilityWarning(MMPBSA_Warning):
    """
    When define stability calculation and protein or ligand
    """
    pass