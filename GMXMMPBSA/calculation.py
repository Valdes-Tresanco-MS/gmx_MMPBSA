"""
This module contains calculation classes that call the necessary programs
for running MM/PBSA calculations.

Methods:
   run_calculations(FILES, INPUT, rank) : Determines which calculations need to
        be run, then sets up the calculations and runs them

Classes:
   Calculation: Base calculation class
   EnergyCalculation: Typical GB/PB binding FE calculations. Handles all sander program calls
   RISMCalculation: RISM binding FE calculation
   NmodeCalc: normal mode entropy calculation
   QuasiHarmCalc: Quasi-harmonic entropy calculation
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
import logging

from GMXMMPBSA.exceptions import CalcError
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import os
import sys
import numpy as np
import math


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class CalculationList(list):
    """ This contains the list of all calculations that need to be run """

    def __init__(self, timer, *args):
        self.timer = timer
        self.timer_keys = []
        self.labels = []
        list.__init__(self)

    def append(self, calc, label='', timer_key=None):
        """ Add a new Calculation instance to the list """
        if not isinstance(calc, Calculation):
            raise TypeError('CalculationList can only take Calculation instances!')

        self.timer_keys.append(timer_key)
        list.append(self, calc)
        self.labels.append(label)

    def extend(self, calcs, labels, timer_keys):
        """ Add a list/iterable of Calculation instances to the list """
        for i, calc in enumerate(calcs):
            CalculationList.append(self, calc, labels[i], timer_keys[i])

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        """ Runs every calculation in the list """
        own_handle = False
        try:
            f = open(stdout, 'w')
            own_handle = True
        except TypeError:
            f = stdout
        try:
            for i, calc in enumerate(self):
                # Start timer, run calculation, then stop the timer
                if self.timer_keys[i] is not None:
                    self.timer.start_timer(self.timer_keys[i])
                if self.labels[i] and rank == 0:
                    logging.info(self.labels[i])
                calc.setup()
                calc.run(rank, stdout=stdout, stderr=stderr)
                if self.timer_keys[i] is not None:
                    self.timer.stop_timer(self.timer_keys[i])
        finally:
            if own_handle: f.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Calculation(object):
    """ Base calculation class. All other calculation classes should be inherited
        from this class.
    """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prog, prmtop, incrd, input_file, output, inptraj, xvv=None):
        self.prmtop = str(prmtop)
        self.incrd = incrd
        self.input_file = input_file
        self.inptraj = inptraj
        self.output = output
        self.program = prog
        self.xvv = xvv

        self.calc_setup = False  # This means that the setup has run successfully

        self.command_args = [self.program]

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        """ Runs the program. All command-line arguments must be set before
            calling this method. Command-line arguments should be set in setup()
        """
        from subprocess import Popen

        # If this has not been set up yet
        # then raise a stink
        if not self.calc_setup:
            raise CalcError('Cannot run a calculation without calling its' +
                            ' its setup() function!')

            # Here, make sure that we could pass a file *OR* a string as stdout/stderr.
        # If they are strings, then open files up with that name, and make sure to
        # close them afterwards. The setup() method should make sure that they are
        # either a file or a string!
        own_handleo = own_handlee = False
        try:
            process_stdout = open(stdout, 'w')
            own_handleo = True
        except TypeError:
            process_stdout = stdout
        try:
            process_stderr = open(stderr, 'w')
            own_handlee = True
        except TypeError:
            process_stderr = stderr

        # The setup() method sets the command-line arguments and makes sure that
        # all of the CL arguments are set. Now all we have to do is start the
        # process and monitor it for success.

        # Popen can only take strings as command-line arguments, so convert
        # everything to a string here. And if it appears to need the rank
        # substituted into the file name, substitute that in here
        try:
            for i in range(len(self.command_args)):
                self.command_args[i] = str(self.command_args[i])
                if '%d' in self.command_args[i]:
                    self.command_args[i] %= rank

            process = Popen(self.command_args, stdin=None, stdout=process_stdout, stderr=process_stderr)

            calc_failed = bool(process.wait())

            if calc_failed:
                raise CalcError(f'{self.program} failed with prmtop {self.prmtop}!')
        finally:
            if own_handleo: process_stdout.close()
            if own_handlee: process_stdout.close()

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def setup(self):
        """ Sets up the Calculation. Finds the program and adds that to the
            first element of the array. Inherited classes should call this
            method first, but then do anything else that is necessary for that
            calculation.
        """
        self.calc_setup = True


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class EnergyCalculation(Calculation):
    """ Uses mmpbsa_py_energy to evaluate energies """

    def __init__(self, prog, prmtop, incrd, input_file, output, restrt=None, inptraj=None, xvv=None):
        Calculation.__init__(self, prog, prmtop, incrd, input_file, output, inptraj, xvv)
        self.restrt = restrt

    def setup(self):
        """
        Sets up the command-line arguments. Sander requires a unique restrt file
        for the MPI version (since one is *always* written and you don't want 2
        threads fighting to write the same dumb file)
        """
        self.command_args.append('-O')  # overwrite flag
        self.command_args.extend(('-i', self.input_file))  # input file flag
        self.command_args.extend(('-p', self.prmtop))  # prmtop flag
        self.command_args.extend(('-c', self.incrd))  # input coordinate flag
        self.command_args.extend(('-o', self.output))  # output file flag
        if self.inptraj:
            self.command_args.extend(('-y', self.inptraj))  # input trajectory flag
        if self.restrt is not None:
            self.command_args.extend(('-r', self.restrt))  # restart file flag
        if self.xvv is not None:
            self.command_args.extend(('-xvv', self.xvv))  # xvv file flag

        # Now test to make sure that the input file exists, since that's the only
        # one that may be absent (due to the use of -use-mdins)
        if not os.path.exists(self.input_file):
            raise IOError("Input file (%s) doesn't exist" % self.input_file)

        self.calc_setup = True


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class RISMCalculation(Calculation):
    """ This class handles RISM calculations """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prog, prmtop, incrd, inptraj, xvvfile, output, INPUT):
        """ Sets up a RISM calculation. It's not as similar to the base class as
            other calculation classes are, but it still inherits useful methods
        """
        # rism3d.snglpnt dumps its output to stdout
        Calculation.__init__(self, prog, prmtop, incrd, None, output, inptraj)

        # Set up instance variables
        self.xvvfile = xvvfile
        self.closure = ','.join(map(str, INPUT['closure']))
        self.polardecomp = INPUT['polardecomp']
        self.ng = ','.join(map(str, INPUT['ng']))
        self.solvbox = ','.join(map(str, INPUT['solvbox']))
        self.buffer = INPUT['buffer']
        self.grdspc = ','.join(map(str, INPUT['grdspc']))
        self.solvcut = INPUT['solvcut']
        self.tolerance = ','.join(map(str, INPUT['tolerance']))
        self.verbose = INPUT['rism_verbose']
        self.solvbox = ','.join(map(str, INPUT['solvbox']))
        self.gf = INPUT['rismrun_gf']

        self.noasympcorr = INPUT['noasympcorr']
        self.mdiis_del = INPUT['mdiis_del']
        self.mdiis_restart = INPUT['mdiis_restart']
        self.mdiis_nvec = INPUT['mdiis_nvec']
        self.maxstep = INPUT['maxstep']
        self.npropagate = INPUT['npropagate']
        # self.centering = INPUT['centering']
        # self.entropicDecomp = INPUT['entropicDecomp']
        # self.pc_plus = INPUT['rismrun_pc+']
        # self.uccoeff = ','.join(map(str, INPUT['uccoeff']))
        self.treeDCF = INPUT['treeDCF']
        self.treeTCF = INPUT['treeTCF']
        self.treeCoulomb = INPUT['treeCoulomb']
        self.treeDCFOrder = INPUT['treeDCFOrder']
        self.treeTCFOrder = INPUT['treeTCFOrder']
        self.treeCoulombOrder = INPUT['treeCoulombOrder']
        self.treeDCFN0 = INPUT['treeDCFN0']
        self.treeTCFN0 = INPUT['treeTCFN0']
        self.treeCoulombN0 = INPUT['treeCoulombN0']
        self.treeDCFMAC = INPUT['treeDCFMAC']
        self.treeTCFMAC = INPUT['treeTCFMAC']
        self.treeCoulombMAC = INPUT['treeCoulombMAC']
        self.asympKSpaceTolerance = INPUT['asympKSpaceTolerance']
        self.ljTolerance = INPUT['ljTolerance']

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def setup(self):
        """ Sets up the RISM calculation. All it has to do is fill in the
            necessary command-line arguments
        """
        # Set up some defaults
        ngflag = self.ng != "-1,-1,-1"
        solvboxflag = self.solvbox != "-1,-1,-1"
        polardecompflag = bool(self.polardecomp)
        gfflag = bool(self.gf)
        # pc_plusflag = bool(self.pc_plus)

        Calculation.setup(self)
        self.command_args.extend(('--xvv', self.xvvfile,
                                  '--closure', self.closure,
                                  '--buffer', self.buffer,
                                  '--grdspc', self.grdspc,
                                  '--solvcut', self.solvcut,
                                  '--tolerance', self.tolerance,
                                  '--verbose', self.verbose,
                                  '--prmtop', self.prmtop,
                                  '--pdb', self.incrd,
                                  '--traj', self.inptraj))
        if ngflag:
            self.command_args.extend(('--ng', self.ng))
        if solvboxflag:
            self.command_args.extend(('--solvbox', self.solvbox))
        if polardecompflag:
            self.command_args.extend(['--polarDecomp'])
        if gfflag:
            self.command_args.extend(['--gf'])
        # if pc_plusflag:
        #     self.command_args.extend(['--pc+'])
        if not os.path.exists(self.xvvfile):
            raise IOError('XVVFILE (%s) does not exist!' % self.xvvfile)

        # additional variables
        var_names = [self.mdiis_del, self.mdiis_restart, self.mdiis_nvec, self.maxstep, self.npropagate,
                     self.treeDCF, self.treeTCF, self.treeCoulomb,
                     self.treeDCFOrder, self.treeTCFOrder, self.treeCoulombOrder, self.treeDCFN0,
                     self.treeTCFN0, self.treeCoulombN0, self.treeDCFMAC, self.treeTCFMAC, self.treeCoulombMAC,
                     self.asympKSpaceTolerance, self.ljTolerance]

        var_input_names = ['mdiis_del', 'mdiis_restart', 'mdiis_nvec', 'maxstep', 'npropagate',
                           'treeDCF', 'treeTCF', 'treeCoulomb',
                           'treeDCFOrder', 'treeTCFOrder', 'treeCoulombOrder', 'treeDCFN0', 'treeTCFN0',
                           'treeCoulombN0', 'treeDCFMAC', 'treeTCFMAC', 'treeCoulombMAC',
                           'asympKSpaceTolerance', 'ljTolerance']

        for i in zip(var_names, var_input_names):
            if i[1] not in ['treeDCF', 'treeTCF', 'treeCoulomb']:
                self.command_args.extend((f'--{i[1]}', str(i[0])))
            elif i[0] != 0:
                self.command_args.extend([f'--{i[1]}'])

        self.calc_setup = True

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def run(self, rank, *args, **kwargs):
        Calculation.run(self, rank, stdout=self.output % rank)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class NmodeCalc(Calculation):
    """ Calculates entropy contribution by normal mode approximation """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prog, prmtop, incrd, inptraj, output, INPUT):
        """ Initializes the nmode calculation. Need to set the options string """
        from math import sqrt
        Calculation.__init__(self, prog, prmtop, incrd, None, output, inptraj)

        kappa = sqrt(0.10806 * INPUT['nmode_istrng'])
        if INPUT['nmode_igb']:
            option_string = ('ntpr=10000, diel=C, kappa=%f, cut=1000, gb=1, ' +
                             'dielc=%f, temp0=%f') % (kappa, INPUT['dielc'], INPUT['temperature'])
        else:
            option_string = ('ntpr=10000, diel=R, kappa=%f, cut=1000, gb=0, ' +
                             'dielc=%f, temp0=%f') % (kappa, INPUT['dielc'], INPUT['temperature'])

        self.option_string = option_string
        self.drms = INPUT['drms']
        self.maxcyc = INPUT['maxcyc']

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def setup(self):
        """ Sets up the simulation """

        self.command_args.extend((self.incrd, self.prmtop, self.maxcyc, self.drms,
                                  self.option_string, self.inptraj))
        self.calc_setup = True

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def run(self, rank, *args, **kwargs):
        Calculation.run(self, rank, stdout=self.output % rank)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class QuasiHarmCalc(Calculation):
    """ Quasi-harmonic entropy calculation class """

    def __init__(self, prog, prmtop, inptraj, input_file, output,
                 receptor_mask, ligand_mask, fnpre):
        """ Initializes the Quasi-harmonic calculation class """
        Calculation.__init__(self, prog, prmtop, None, input_file, output, inptraj)
        self.stability = not bool(receptor_mask) and not bool(ligand_mask)
        self.receptor_mask, self.ligand_mask = receptor_mask, ligand_mask
        self.calc_setup = False
        self.fnpre = fnpre  # file name prefix

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def setup(self):
        """ Sets up a Quasi-harmonic calculation """
        from subprocess import Popen, PIPE

        # Determine the prefix from our input file... hack way to do this
        if self.input_file.startswith(self.fnpre + 'mutant_'):
            prefix = self.fnpre + 'mutant_'
        else:
            prefix = self.fnpre

        # Make sure masks are a list, and that there are enough masks

        # First thing we need is the average PDB as a reference
        ptraj_str = 'trajin %s\naverage %savgcomplex.pdb pdb chainid " "\ngo' % (self.inptraj,
                                                                                 prefix)

        outfile = open(self.fnpre + 'create_average.out', 'w')

        process = Popen([self.program, self.prmtop], stdin=PIPE, stdout=outfile)
        out, err = process.communicate(ptraj_str.encode())

        if process.wait():
            raise CalcError('Failed creating average PDB')

        outfile.close()

        # Now that we have the PDB file

        self.command_args.extend((self.prmtop, self.input_file))

        self.calc_setup = True

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def run(self, rank, *args, **kwargs):
        Calculation.run(self, rank, stdout=self.output)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PBEnergyCalculation(EnergyCalculation):
    """
    Specially handle the PB calculations to extract warnings and errors PBSA
    prints to stdout and redirect them to the user
    """

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        """
        Runs the program. All command-line arguments must be set before calling
        this method. Command-line arguments should be set in setup()
        stdout is ignored here because we need to parse it for errors
        """
        import re
        from subprocess import Popen, PIPE

        # If this has not been set up yet
        # then raise a stink
        if not self.calc_setup:
            raise CalcError('Cannot run a calculation without calling its' +
                            ' its setup() function!')

        errorre = re.compile('(pb (?:bomb)|(?:warning))', re.I)
        # Here, make sure that we could pass a file *OR* a string as stderr.
        own_handle = False
        try:
            process_stderr = open(stderr, 'w')
            own_handle = True
        except TypeError:
            process_stderr = stderr

        # The setup() method sets the command-line arguments and makes sure that
        # all of the CL arguments are set. Now all we have to do is start the
        # process and monitor it for success.

        # Popen can only take strings as command-line arguments, so convert
        # everything to a string here. If rank needs to be substituted in, do that
        # here
        try:
            for i in range(len(self.command_args)):
                self.command_args[i] = str(self.command_args[i])
                if '%d' in self.command_args[i]:
                    self.command_args[i] = self.command_args[i] % rank

            process = Popen(self.command_args, stdin=None, stdout=PIPE,
                            stderr=process_stderr)

            out, err = process.communicate(b'')
            calc_failed = bool(process.wait())
            out = out.decode('utf-8')
            if calc_failed:
                error_list = [s.strip() for s in out.split('\n')
                              if errorre.match(s.strip())]

                GMXMMPBSA_ERROR('%s failed with prmtop %s!\n\t' % (self.program, self.prmtop) +
                                '\n\t'.join(error_list) + '\n' +
                                'If you are using sander and PB calculation, check the *.mdout files to get the sander '
                                'error\n',
                                CalcError)
        finally:
            if own_handle: process_stderr.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SurfCalc(Calculation):
    """
    Base class for a surface area calculation using cpptraj. Children must
    implement _get_instring(self, rank) which returns the string containing the
    necessary cpptraj input
    """

    def __init__(self, prog, prmtop, inptraj, output, probe=1.4, offset=0.0):
        self.prmtop = str(prmtop)
        self.inptraj = inptraj
        self.output = output
        self.program = prog
        self.probe = probe
        self.offset = offset

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        """ Runs the program. All command-line arguments must be set before
            calling this method. Command-line arguments should be set in setup()
        """
        from subprocess import Popen, PIPE

        # If this has not been set up yet
        # then raise a stink
        if not self.calc_setup:
            raise CalcError('Cannot run a calculation without calling its' +
                            ' its setup() function!')

            # Make sure the inptraj and output are rank-substituted
        instring = self._get_instring(rank)

        process = Popen([self.program, self.prmtop], stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)

        out, err = process.communicate(instring.encode())

        calc_failed = bool(process.wait())

        if calc_failed:
            raise CalcError('%s failed with prmtop %s!' % (self.program, self.prmtop))


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class LcpoCalc(SurfCalc):
    """
    Uses LCPO to calculate the surface area
    (Linear Combination of Pairwise Overlaps)
    """

    def _get_instring(self, rank):
        """ Returns the cpptraj input string """
        inptraj = self.inptraj % rank
        output = self.output % rank
        return "trajin %s\nsolvent none\nsurf :* out %s\n" % (inptraj, output)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MolsurfCalc(SurfCalc):
    """ Uses molsurf to calculate the surface area """

    def __init__(self, prog, prmtop, inptraj, output, probe=1.4, offset=0.0):
        SurfCalc.__init__(self, prog, prmtop, inptraj, output)
        self.probe = probe
        self.offset = offset

    def _get_instring(self, rank):
        inptraj = self.inptraj % rank
        output = self.output % rank
        return "trajin %s\nmolsurf :* out %s probe %s offset %s" % (inptraj,
                                                                    output, self.probe, self.offset)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class CopyCalc(Calculation):
    """
    This is for mutant files that are unchanged from 'normal' files (i.e., when
    the mutation is in the receptor, the ligand outputs are copied)
    """

    def __init__(self, orig_name, final_name):
        self.orig_name = orig_name
        self.final_name = final_name

    def run(self, rank, stdout=None, stderr=None):
        from shutil import copy
        # Do rank-substitution if necessary
        if '%d' in self.orig_name:
            orig_name = self.orig_name % rank
        else:
            orig_name = self.orig_name

        if '%d' in self.final_name:
            final_name = self.final_name % rank
        else:
            final_name = self.final_name

        copy(orig_name, final_name)


class MergeOut(Calculation):
    def __init__(self, filename, mdouts):
        self.filename = filename
        self.mdouts = mdouts

    def run(self, rank, stdout=None, stderr=None):
        # Do rank-substitution if necessary
        filename = self.filename % rank if '%d' in self.filename else self.filename
        with open(filename, 'w') as ofile:
            for file in sorted(self.mdouts, key=lambda x: f"{int(x.split('.')[1]):08d}"):
                with open(file % rank) as rf:
                    for line in rf.readlines():
                        ofile.write(line)
                ofile.write('\n')

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PrintCalc(Calculation):
    """
    This is just a way to insert a printed message to the screen during the
    calculation list execution
    """

    def __init__(self, message):
        self.message = message

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        if rank == 0:
            logging.info(self.message)
            # stdout.write(self.message + '\n')


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class InteractionEntropyCalc:
    """
    Class for Interaction Entropy calculation
    :return {IE_key: data}
    """

    def __init__(self, ggas, INPUT, iesegment=None):
        """

        Args:
            ggas: Model GGAS energy
            INPUT: INPUT dict
            iesegment: If not defined, iesegment = INPUT['ie_segment']
        """
        self.ggas = ggas
        self.INPUT = INPUT
        self.isegment = iesegment or INPUT['ie_segment']
        self.data = []

        self._calculate()

    def _calculate(self):
        # boltzmann constant in kcal/(mol⋅K)
        k = 0.001985875
        temperature = self.INPUT['temperature']

        energy_int = np.array([], dtype=np.float)
        a_energy_int = np.array([], dtype=np.float)
        exp_energy_int = np.array([], dtype=np.float)
        self.data = np.array([], dtype=np.float)
        for eint in self.ggas:
            energy_int = np.append(energy_int, eint)
            aeint = energy_int.mean()
            a_energy_int = np.append(a_energy_int, aeint)
            deint = eint - aeint
            if deint > 425:
                logging.warning('The internal energy of your system has very large energy fluctuation so it is not '
                                'possible to continue with the calculations. Please, make sure your system is '
                                'consistent')
                logging.info('The Interaction Entropy will be skipped...')
                self.INPUT['general']['interaction_entropy'] = 0
                break
            eceint = math.exp(deint / (k * temperature))
            exp_energy_int = np.append(exp_energy_int, eceint)
            aeceint = exp_energy_int.mean()
            cts = k * temperature * math.log(aeceint)
            self.data = np.append(self.data, cts)
        numframes = len(self.data)
        self.ie_std = energy_int.std()
        self.ieframes = math.ceil(numframes * (self.isegment / 100))
        self.iedata = self.data[-self.ieframes:]
        self.frames = list(
            range(
                self.INPUT['startframe'],
                self.INPUT['startframe']
                + numframes * self.INPUT['interval'],
                self.INPUT['interval'],
            )
        )

    def save_output(self, filename):
        with open(filename, 'w') as out:
            out.write(f'Calculation for last {self.ieframes} frames:\n')
            out.write(f'Interaction Entropy (-TΔS): {self.iedata.mean():9.4f} +/- {self.iedata.std():7.4f}\n\n')
            out.write('Interaction Entropy per-frame:\n')

            out.write('Frame # | IE value\n')
            for f, d in zip(self.frames, self.data):
                out.write('{:d}  {:.2f}\n'.format(f, d))

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class C2EntropyCalc:
    """
    Class for Interaction Entropy calculation
    :return {IE_key: data}
    """

    def __init__(self, ggas, INPUT):
        self.ggas = ggas
        self.INPUT = INPUT

        self._calculate()

    def _calculate(self):
        # gas constant in kcal/(mol⋅K)
        R = 0.001987
        temperature = self.INPUT['temperature']
        self.ie_std = self.ggas.std()
        self.c2data = (self.ie_std ** 2) / (2 * temperature * R)

        array_of_c2 = np.zeros(2000)
        for i in range(2000):
            idxs = np.random.randint(0, len(self.ggas), len(self.ggas))
            ie_std = self.ggas[idxs].std()
            c2data = (ie_std ** 2) / (2 * temperature * R)
            array_of_c2[i] = c2data

        self.c2_std = np.sort(array_of_c2)[100:1900].std()
        self.c2_ci = np.percentile(np.sort(array_of_c2)[100:1900], [2.5, 97.5])

    def save_output(self, filename):
        with open(filename, 'w') as out:
            out.write(f'C2 Entropy (-TΔS): {self.c2data:.4f}\n\n')


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
