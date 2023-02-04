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
import threading
from pathlib import Path
from tqdm import tqdm
from time import sleep

from GMXMMPBSA.exceptions import CalcError
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
import os
import sys
import numpy as np
import math


TQDM_BAR_FORMAT = '            {l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def pb(output_basename, nframes=1, mpi_size=1, nmode=False):
    pbar = tqdm(total=nframes, ascii=True, bar_format=TQDM_BAR_FORMAT)
    accum_frames = 0
    ctime = 0

    while accum_frames < nframes:
        sleep(ctime)
        frames = 0
        for i in range(mpi_size):
            if 'gbnsr6' in output_basename:
                _output_folder, _output_filename = output_basename.split('/')
                output_folder = Path(_output_folder % i)
                output_filename = Path(_output_filename)
                frames = len(list(output_folder.glob(f"{output_filename.stem}*")))
            else:
                obasename = Path(output_basename % i)
                if not obasename.exists():
                    continue
                with obasename.open() as of:
                    for line in of:
                        if not nmode and line.startswith('                    FINAL RESULTS'):
                            frames += 1
                        elif nmode and line.startswith('Total:'):
                            frames += 1

        if frames - accum_frames:
            pbar.update(frames - accum_frames)
            accum_frames = frames
        else:
            ctime += 1

    pbar.clear()
    pbar.close()


class CalculationList(list):
    """ This contains the list of all calculations that need to be run """

    def __init__(self, timer, *args):
        self.timer = timer
        self.timer_keys = []
        self.labels = []
        self.output_files = []
        self.nframes, self.nmframes, self.mpi_size = args
        list.__init__(self)

    def append(self, calc, label='', timer_key=None, output_basename=None):
        """ Add a new Calculation instance to the list """
        if not isinstance(calc, (Calculation, MultiCalculation)):
            raise TypeError('CalculationList can only take Calculation instances!')

        self.timer_keys.append(timer_key)
        list.append(self, calc)
        self.labels.append(label)
        self.output_files.append(output_basename)

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
                pb_thread = None
                # Start timer, run calculation, then stop the timer
                if self.timer_keys[i] is not None:
                    self.timer.start_timer(self.timer_keys[i])
                if self.labels[i] and rank == 0:
                    logging.info(self.labels[i])
                    if isinstance(calc, (EnergyCalculation, ListEnergyCalculation, NmodeCalc)):
                        if isinstance(calc, (EnergyCalculation, ListEnergyCalculation)):
                            nframes = self.nframes
                            nmode = False
                        else:
                            nframes = self.nmframes
                            nmode = True
                        pb_thread = threading.Thread(target=pb, args=(self.output_files[i], nframes, self.mpi_size,
                                                                      nmode), daemon=True)
                        pb_thread.start()

                calc.setup()
                calc.run(rank, stdout=stdout, stderr=stderr)
                if self.timer_keys[i] is not None:
                    self.timer.stop_timer(self.timer_keys[i])
                    if pb_thread:
                        pb_thread.join()
        finally:
            if own_handle: f.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class MultiCalculation(object):
    def __init__(self):
        self.list_calc = []

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        """ Runs the program. All command-line arguments must be set before
                    calling this method. Command-line arguments should be set in setup()
                """
        from subprocess import Popen

        # If this has not been set up yet
        # then raise a stink
        if not self.calc_setup:
            raise CalcError('Cannot run a calculation without calling its its setup() function!')

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
            for command_args in self.list_calc:
                for i in range(len(command_args)):
                    command_args[i] = str(command_args[i])
                    if '%d' in command_args[i]:
                        command_args[i] %= rank

                process = Popen(command_args, stdin=None, stdout=process_stdout, stderr=process_stderr)

                calc_failed = bool(process.wait())
                # print(command_args, calc_failed)
                if calc_failed:
                    raise CalcError(f'{command_args[0]} failed with prmtop {command_args[1]}!')
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


class Calculation(object):
    """ Base calculation class. All other calculation classes should be inherited
        from this class.
    """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prog, prmtop, incrd, inptraj, input_file, output, xvv=None):
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

    def __init__(self, prog, prmtop, incrd, inptraj, input_file, output, restrt, xvv=None):
        Calculation.__init__(self, prog, prmtop, incrd, inptraj,
                             input_file, output, xvv)
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
        if self.inptraj is not None:
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

class ListEnergyCalculation(MultiCalculation):
    def __init__(self, prog, prmtop, input_file, incrds, outputs, xvv=None):
        super().__init__()
        self.program = prog
        self.prmtop = prmtop
        self.incrds = incrds
        self.input_file = input_file
        self.outputs = outputs
        self.xvv = xvv

    def setup(self):
        """
        Sets up the command-line arguments. Sander requires a unique restrt file
        for the MPI version (since one is *always* written and you don't want 2
        threads fighting to write the same dumb file)
        """
        # print(self.incrds, self.outputs)
        for c, o in zip(self.incrds, self.outputs):
            command_args = [self.program,
                            '-i', self.input_file,  # input file flag
                            '-p', self.prmtop  # prmtop flag
                            ]
            command_args.extend(('-c', c))  # input coordinate flag
            command_args.extend(('-o', o))  # output file flag
            self.list_calc.append(command_args)

            # Now test to make sure that the input file exists, since that's the only
            # one that may be absent (due to the use of -use-mdins)
            # if not os.path.exists(self.input_file):
            #     raise IOError("Input file (%s) doesn't exist" % self.input_file)

        self.calc_setup = True
class RISMCalculation(Calculation):
    """ This class handles RISM calculations """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, prog, prmtop, incrd, inptraj, xvvfile, output, INPUT):
        """ Sets up a RISM calculation. It's not as similar to the base class as
            other calculation classes are, but it still inherits useful methods
        """
        # rism3d.snglpnt dumps its output to stdout
        Calculation.__init__(self, prog, prmtop, incrd, inptraj, None, output)

        # Set up instance variables
        self.xvvfile = xvvfile
        self.closure = ','.join(map(str, INPUT['rism']['closure']))
        self.polardecomp = INPUT['rism']['polardecomp']
        self.ng = ','.join(map(str, INPUT['rism']['ng']))
        self.solvbox = ','.join(map(str, INPUT['rism']['solvbox']))
        self.buffer = INPUT['rism']['buffer']
        self.grdspc = ','.join(map(str, INPUT['rism']['grdspc']))
        self.solvcut = INPUT['rism']['solvcut']
        self.tolerance = ','.join(map(str, INPUT['rism']['tolerance']))
        self.verbose = INPUT['rism']['rism_verbose']
        self.solvbox = ','.join(map(str, INPUT['rism']['solvbox']))
        self.gf = INPUT['rism']['rismrun_gf']

        self.noasympcorr = INPUT['rism']['noasympcorr']
        self.mdiis_del = INPUT['rism']['mdiis_del']
        self.mdiis_restart = INPUT['rism']['mdiis_restart']
        self.mdiis_nvec = INPUT['rism']['mdiis_nvec']
        self.maxstep = INPUT['rism']['maxstep']
        self.npropagate = INPUT['rism']['npropagate']
        # self.centering = INPUT['rism']['centering']
        # self.entropicDecomp = INPUT['rism']['entropicDecomp']
        # self.pc_plus = INPUT['rism']['rismrun_pc+']
        # self.uccoeff = ','.join(map(str, INPUT['rism']['uccoeff']))
        self.treeDCF = INPUT['rism']['treeDCF']
        self.treeTCF = INPUT['rism']['treeTCF']
        self.treeCoulomb = INPUT['rism']['treeCoulomb']
        self.treeDCFOrder = INPUT['rism']['treeDCFOrder']
        self.treeTCFOrder = INPUT['rism']['treeTCFOrder']
        self.treeCoulombOrder = INPUT['rism']['treeCoulombOrder']
        self.treeDCFN0 = INPUT['rism']['treeDCFN0']
        self.treeTCFN0 = INPUT['rism']['treeTCFN0']
        self.treeCoulombN0 = INPUT['rism']['treeCoulombN0']
        self.treeDCFMAC = INPUT['rism']['treeDCFMAC']
        self.treeTCFMAC = INPUT['rism']['treeTCFMAC']
        self.treeCoulombMAC = INPUT['rism']['treeCoulombMAC']
        self.asympKSpaceTolerance = INPUT['rism']['asympKSpaceTolerance']
        self.ljTolerance = INPUT['rism']['ljTolerance']

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
        Calculation.__init__(self, prog, prmtop, incrd, inptraj, None, output)

        kappa = sqrt(0.10806 * INPUT['nmode']['nmode_istrng'])
        if INPUT['nmode']['nmode_igb']:
            option_string = ('ntpr=10000, diel=C, kappa=%f, cut=1000, gb=1, ' +
                             'dielc=%f, temp0=%f') % (kappa, INPUT['nmode']['dielc'], INPUT['general']['temperature'])
        else:
            option_string = ('ntpr=10000, diel=R, kappa=%f, cut=1000, gb=0, ' +
                             'dielc=%f, temp0=%f') % (kappa, INPUT['nmode']['dielc'], INPUT['general']['temperature'])

        self.option_string = option_string
        self.drms = INPUT['nmode']['drms']
        self.maxcyc = INPUT['nmode']['maxcyc']

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
        Calculation.__init__(self, prog, prmtop, None, inptraj,
                             input_file, output)
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
            raise CalcError(f'{self.program} failed with prmtop {self.prmtop}!')


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
    def __init__(self, topology, output_filename, mm_filename, mdout_filenames, idecomp, dec_verbose):
        self.topology = topology
        self.output_filename = output_filename
        self.mm_filename = mm_filename
        self.mdouts = mdout_filenames
        self.idecomp = idecomp
        self.dec_verbose = dec_verbose


    def run(self, rank, stdout=None, stderr=None):
        # Do rank-substitution if necessary
        out_filename = self.output_filename % rank if '%d' in self.output_filename else self.output_filename
        mm_filename = self.mm_filename % rank if '%d' in self.mm_filename else self.mm_filename
        MergeGBNSR6Output(self.topology, out_filename, mm_filename,
                          self.mdouts, self.idecomp, self.dec_verbose)

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
        self.isegment = iesegment or INPUT['general']['ie_segment']
        self.data = []

        self._calculate()

    def _calculate(self):
        # boltzmann constant in kcal/(mol⋅K)
        k = 0.001985875
        temperature = self.INPUT['general']['temperature']

        exp_energy_int = np.array([], dtype=float)
        self.data = np.zeros(self.ggas.size, dtype=float)

        for i in tqdm(range(self.ggas.size), bar_format=TQDM_BAR_FORMAT, ascii=True):
            aeint = self.ggas[:i+1].mean()
            deint = self.ggas[i] - aeint
            try:
                eceint = math.exp(deint / (k * temperature))
            except CalcError:
                logging.warning('The internal energy of your system has very large energy fluctuation so it is not '
                                'possible to continue with the calculations. Please, make sure your system is '
                                'consistent')
                logging.info('The Interaction Entropy will be skipped...')
                self.INPUT['general']['interaction_entropy'] = 0
                break
            exp_energy_int = np.append(exp_energy_int, eceint)
            aeceint = exp_energy_int.mean()
            cts = k * temperature * math.log(aeceint)
            self.data[i] = cts

        numframes = len(self.data)
        self.ie_std = float(self.ggas.std())
        self.ieframes = math.ceil(numframes * (self.isegment / 100))
        self.iedata = self.data[-self.ieframes:]

    def save_output(self, filename):
        frames = list(
            range(
                self.INPUT['general']['startframe'],
                self.INPUT['general']['startframe'] + len(self.data) * self.INPUT['general']['interval'],
                self.INPUT['general']['interval'],
            )
        )
        with open(filename, 'w') as out:
            out.write('| Interaction Entropy results\n')
            out.write(f'IE-frames: last {self.ieframes}\n')
            out.write(f'Internal Energy SD (sigma): {self.ie_std:9.2f}\n')
            out.write(f'| Interaction Entropy (-TΔS): {self.iedata.mean():9.2f} +/- {self.iedata.std():7.2f}\n\n')
            out.write('| Interaction Entropy per-frame:\n')

            out.write('Frame # | IE value\n')
            for f, d in zip(frames, self.data):
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
        temperature = self.INPUT['general']['temperature']
        self.ie_std = float(self.ggas.std())
        self.c2data = (self.ie_std ** 2) / (2 * temperature * R)

        size = self.ggas.size
        array_of_c2 = np.zeros(2000)
        for i in range(2000):
            idxs = np.random.randint(0, size, size)
            ie_std = self.ggas[idxs].std()
            c2data = (ie_std ** 2) / (2 * temperature * R)
            array_of_c2[i] = c2data

        self.c2_std = float(np.sort(array_of_c2)[100:1900].std())
        self.c2_ci = np.percentile(np.sort(array_of_c2)[100:1900], [2.5, 97.5])

    def save_output(self, filename):
        with open(filename, 'w') as out:
            out.write('| C2 Entropy results\n')
            out.write(f'C2 Entropy (-TΔS): {self.c2data:.4f}\n')
            out.write(f'C2 Entropy SD: {self.c2_std:.4f}\n')
            out.write(f'Internal Energy SD (sigma): {self.ie_std:9.2f}\n')
            out.write(f'C2 Entropy CI: {self.c2_ci[0]:.4f} {self.c2_ci[1]:.4f}\n')



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def get_gbnsr6_out(dgij, topology, idecomp=0, dec_verbose=0, res2print=None):
    import parmed
    t = parmed.load_file(topology)
    res_list = {residue.idx + 1: [atm.idx + 1 for atm in residue.atoms] for residue in t.residues}
    if idecomp in [1, 2]:
        pw = {x: {y: {} for y in res_list} for x in res_list if x in res2print}
    else:
        pw = {x: {y: {} for y in res_list if y in res2print} for x in res_list if x in res2print}

    for line in dgij:
        if line.startswith('DGij'):
            kw, at1, at2, energy = line.strip('\n').split()
            res_idx = t.atoms[int(at1) - 1].residue.idx + 1
            res2_idx = t.atoms[int(at2) - 1].residue.idx + 1
            if res_idx not in res2print:
                continue
            if idecomp in [1, 2]:
                pw[res_idx][res2_idx].setdefault((at1, at2), float(energy))
                if res_idx != res2_idx and res2_idx in res2print:
                    pw[res2_idx][res_idx].setdefault((at2, at1), float(energy))
            else:
                if res2_idx not in res2print:
                    continue
                pw[res_idx][res2_idx].setdefault((at1, at2), float(energy))
                if res_idx != res2_idx:
                    pw[res2_idx][res_idx].setdefault((at2, at1), float(energy))
    return _get_decomp(pw, idecomp, dec_verbose, t)


def _get_decomp(pw, idecomp, dec_verbose, t):
    bb = ['CA', 'C', 'O', 'N', 'H', 'OXT', 'H1', 'H2', 'H3']
    decomp = {'TDC': []}
    if dec_verbose in [1, 3]:
        decomp |= {'BDC': [], 'SDC': []}
    for r1, v1 in pw.items():
        if idecomp in [1, 2]:
            TDC = sum(sum(float(x) for x in v2.values()) for r2, v2 in v1.items())
            decomp['TDC'].append(['TDC', r1, TDC])
            if dec_verbose in [1, 3]:
                BDC = sum(sum(v3 for (tr1, tr2), v3 in v2.items() if t.atoms[int(tr1) - 1].name in bb)
                          for r2, v2 in v1.items())
                SDC = TDC - BDC
                decomp['BDC'].append(['BDC', r1, BDC])
                decomp['SDC'].append(['SDC', r1, SDC])
        else:
            for r2, v2 in v1.items():
                TDC = sum(float(x) for x in v2.values())
                decomp['TDC'].append(['TDC', r1, r2, TDC])
                if dec_verbose in [1, 3]:
                    BDC = sum(v3 for (tr1, tr2), v3 in v2.items() if t.atoms[int(tr1) - 1].name in bb)
                    SDC = TDC - BDC
                    decomp['BDC'].append(['BDC', r1, r2, BDC])
                    decomp['SDC'].append(['SDC', r1, r2, SDC])
    return decomp


class MergeGBNSR6Output():
    def __init__(self, topology, output_filename, mm_filename, mdout_filenames, idecomp, dec_verbose):
        self.topology = topology
        self.output_filename = output_filename
        self.mm_filename = mm_filename
        self.mdout_filenames = mdout_filenames
        self.idecomp = idecomp

        self.dec_verbose = dec_verbose
        self.header = '''
          -------------------------------------------------------
          SANDER + GBNSR6
          -------------------------------------------------------\n\n
          '''
        self.resource = ('--------------------------------------------------------------------------------\n'
                         '   1.  ' 'RESOURCE   USE:\n'
                         '--------------------------------------------------------------------------------\n')
        self.control_data = ('--------------------------------------------------------------------------------\n'
                             '   2.  CONTROL  DATA  FOR  THE  RUN\n'
                             '--------------------------------------------------------------------------------\n')
        self.atomic_coor = ('--------------------------------------------------------------------------------\n'
                            '   3.  ATOMIC COORDINATES AND VELOCITIES\n'
                            '--------------------------------------------------------------------------------\n')
        self.results = ('--------------------------------------------------------------------------------\n'
                        '   4.  RESULTS\n'
                        '--------------------------------------------------------------------------------\n')

        self.decomp_labels = {'TDC': 'TOTAL ENERGIES', 'SDC': 'SIDECHAIN ENERGIES', 'BDC': 'BACKBONE ENERGIES'}
        self.decomp_headers = {'pr': '                    PRINT DECOMP - {}\n\n'
                                     '    resid |internal |vdw      |eel      |pol      |sas\n'
                                     '============================================================\n',
                               'pw': '                    PRINT PAIR DECOMP - {}\n\n'
                                     '    resid1 ->resid2 |internal    |vdw         |eel         |pol         |sas\n'
                                     '=============================================================================\n'}

        self.write_output()

    def read_mm_output(self):

        file_assignments = []
        inputfile = []
        resource_section = []
        control_data = []
        atomic_coor_vel = []
        results_section = []
        temp_res2print = None

        with open(self.mm_filename) as mmfile:
            current_section = None
            while line := mmfile.readline():
                if 'File Assignments:' in line:
                    current_section = file_assignments
                    line = mmfile.readline()
                elif line.startswith(' Here is the input file:'):
                    current_section = inputfile
                    line = mmfile.readline()
                if line.startswith('----------------------------------------------------------------------------'):
                    line = mmfile.readline()
                    if line.startswith('   1.  RESOURCE   USE:'):
                        current_section = resource_section
                        mmfile.readline()
                        line = mmfile.readline()
                    elif line.startswith('   2.  CONTROL  DATA  FOR  THE  RUN'):
                        current_section = control_data
                        mmfile.readline()
                        line = mmfile.readline()
                    elif 'ATOMIC COORDINATES AND VELOCITIES' in line:
                        current_section = atomic_coor_vel
                        mmfile.readline()
                        line = mmfile.readline()
                    elif '.  RESULTS' in line:
                        current_section = results_section
                        mmfile.readline()
                        line = mmfile.readline()
                    else:
                        continue
                if current_section is not None:
                    current_section.append(line)
                if line[:4] == 'RES ':
                    temp_res2print = line.split()[1:]

        res2print = []
        if temp_res2print:
            for i in range(0, len(temp_res2print), 2):
                res2print.extend(range(int(temp_res2print[i]), int(temp_res2print[i + 1]) + 1))

        return {'file_assignments': file_assignments, 'inputfile': inputfile, 'resource_section': resource_section,
                'control_data': control_data, 'atomic_coor_vel': atomic_coor_vel,
                'results_section': self._get_energy_decomp(results_section), 'res2print': res2print}

    @staticmethod
    def _get_energy_decomp(results_section):
        energy = {}
        decomp = {}

        c = 0
        while True:
            line = results_section[c]
            if line.startswith('minimizing coord set #'):
                f = int(line.split()[-1])
                energy[f] = {}
                decomp[f] = {}
            if line.startswith(' BOND'):
                words = line.split()
                energy[f][words[0].strip()] = float(words[2])
                energy[f][words[3].strip()] = float(words[5])
                energy[f][words[6].strip()] = float(words[8])
                c += 1
                line = results_section[c]
                words = line.split()
                energy[f][words[0].strip()] = float(words[2])
                energy[f][words[3].strip()] = float(words[5])
                energy[f][words[6].strip()] = float(words[8])
                c += 1
                line = results_section[c]
                words = line.split()
                t = ' '.join([words[0].strip(), words[1].strip()])
                energy[f][t] = float(words[3])
                t = ' '.join([words[4].strip(), words[5].strip()])
                energy[f][t] = float(words[7])
                energy[f][words[8].strip()] = float(words[10])
                c += 1
                # line = results_section[c]
                # words = line.split()
                # energy[f][words[0].strip()] = float(words[2])
            if line[:3] in ['TDC', 'SDC', 'BDC']:
                data = [x.strip().strip('->') for x in line.split()]
                if len(data) == 8:
                    _t, _r1, _r2, _i, _v, _e, _p, _s = data
                    data = [_t, int(_r1), int(_r2), float(_i), float(_v), float(_e), float(_p), float(_s)]
                else:
                    _t, _r1, _i, _v, _e, _p, _s = data
                    data = [_t, int(_r1), float(_i), float(_v), float(_e), float(_p), float(_s)]
                if not decomp[f].get(line[:3]):
                    decomp[f][line[:3]] = [data]
                else:
                    decomp[f][line[:3]].append(data)
            c +=1
            if c == len(results_section):
                break
        return {'energy': energy, 'decomp':decomp}

    def read_gbnsr6_output(self, res2print):
        file_assignments = []
        inputfile = []

        energy = {}
        decomp = {}
        for i, filename in enumerate(self.mdout_filenames, start=1):
            results_section = []
            decomp_section = []
            with open(filename) as mmfile:
                current_section = None
                while line := mmfile.readline():
                    if i == 1 and 'File Assignments:' in line:
                        current_section = file_assignments
                        line = mmfile.readline()
                    elif i == 1 and line.startswith(' Here is the input file:'):
                        current_section = inputfile
                        line = mmfile.readline()
                    if line.startswith('----------------------------------------------------------------------------'):
                        line = mmfile.readline()
                        if '.  RESULTS' in line:
                            current_section = results_section
                            mmfile.readline()
                            line = mmfile.readline()
                        else:
                            current_section = None
                    if current_section is not None:
                        if line.startswith('DGij'):
                            decomp_section.append(line)
                        else:
                            current_section.append(line)
            energy[i] = self._get_energy_gbnsr6(results_section)
            if self.idecomp:
                decomp[i] = get_gbnsr6_out(decomp_section, self.topology, self.idecomp, self.dec_verbose, res2print)

        results = {'energy': energy, 'decomp':decomp}
        return {'file_assignments': file_assignments, 'inputfile': inputfile, 'results_section': results}

    @staticmethod
    def _get_energy_gbnsr6(results_section):
        energy = {}

        store = False
        c = 0
        while True:
            line = results_section[c]
            if "FINAL RESULTS" in line:
                store = True
            if store and line.startswith(' EELEC'):
                words = line.split()
                energy[words[3].strip()] = float(words[5])
                c += 1
                line = results_section[c]
                words = line.split()
                energy[words[0].strip()] = float(words[2])
            c += 1
            if c == len(results_section):
                break
        return energy

    def write_output(self):
        mm = self.read_mm_output()
        gbnsr6 = self.read_gbnsr6_output(mm['res2print'])

        with open(self.output_filename, 'w') as output_file:
            output_file.write(self.header)
            output_file.write('File Assignments:\n')
            output_file.write(' MM:\n')
            for l in mm['file_assignments']:
                output_file.write(l)
            output_file.write(' GBNSR6:\n')
            for l in gbnsr6['file_assignments']:
                output_file.write(l)

            output_file.write(' Here is the input file:\n')
            output_file.write(' MM:\n')
            for l in mm['inputfile']:
                output_file.write(l)
            output_file.write(' GBNSR6:\n')
            for l in gbnsr6['inputfile']:
                output_file.write(l)

            output_file.write(self.resource)
            for l in mm['resource_section']:
                output_file.write(l)

            output_file.write(self.control_data)
            for l in mm['control_data']:
                output_file.write(l)
            output_file.write(self.atomic_coor)
            for l in mm['atomic_coor_vel']:
                output_file.write(l)
            output_file.write(self.results)

            mmenergy, mmdecomp = mm['results_section'].values()
            gbenergy, gbdecomp = gbnsr6['results_section'].values()

            k2print = [['BOND', 'ANGLE', 'DIHED'], ['VDWAALS', 'EEL', 'EGB'], ['1-4 VDW', '1-4 EEL', 'RESTRAINT'],
                       ['ESURF']]
            for i, _ in enumerate(range(len(mmenergy)), start=1):
                output_file.write(f'minimizing coord set #       {i}\n\n')
                mmenergy[i].pop('EGB')
                for e, ev in gbenergy[i].items():
                    mmenergy[i][e] = ev
                for kl in k2print:
                    if len(kl) == 3:
                        f = []
                        for klk in kl:
                            f.extend((klk, mmenergy[i][klk]))
                        output_file.write(' {:8s}={:>14.4f}  {:8s}={:>14.4f}  {:11s}={:>14.4f}\n'.format(*f))
                    else:
                        f = [kl[0], mmenergy[i][kl[0]]]
                        output_file.write(' {:8s}={:>14.4f}\n\n'.format(*f))

                if self.idecomp:
                    if self.idecomp in [1, 2]:
                        for term in mmdecomp[i]:
                            output_file.write(self.decomp_headers['pr'].format(self.decomp_labels[term]))
                            for c, l in enumerate(mmdecomp[i][term]):
                                l[-2] = gbdecomp[i][term][c][-1]
                                output_file.write('{}{:>7d}{:>10.3f}{:>10.3f}{:>10.3f}{:>10.3f}{:>10.3f}\n'.format(*l))
                            output_file.write('\n')
                    else:
                        for term in mmdecomp[i]:
                            output_file.write(self.decomp_headers['pw'].format(self.decomp_labels[term]))
                            for c, l in enumerate(mmdecomp[i][term]):
                                l[-2] = gbdecomp[i][term][c][-1]
                                output_file.write('{}{:>8d}->{:>7d}{:>13.4f}{:>13.4f}{:>13.4f}{:>13.4f}{:>13.4f}\n'.format(*l))
                            output_file.write('\n')