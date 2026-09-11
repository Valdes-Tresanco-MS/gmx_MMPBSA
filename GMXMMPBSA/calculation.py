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
import logging
import threading
from pathlib import Path
from tqdm import tqdm
import json
import re

from GMXMMPBSA.exceptions import CalcError
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
from GMXMMPBSA.qmmm_diagnostics import QMMMDiagnostic, parse_qmmm_diagnostics
from GMXMMPBSA.utils import mdout2json
import os
import sys
import numpy as np
import math

from GMXMMPBSA.progress import monitor_progress, TQDM_BAR_FORMAT


# Backward-compatible name used by callers and third-party integrations.
pb = monitor_progress


class CalculationList(list):
    """ This contains the list of all calculations that need to be run """

    def __init__(self, timer, nframes, nmframes, mpi_size, progress_style='auto'):
        self.timer = timer
        self.timer_keys = []
        self.labels = []
        self.output_files = []
        self.nframes, self.nmframes, self.mpi_size = nframes, nmframes, mpi_size
        self.progress_style = progress_style
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
                        label = self.labels[i].strip().removeprefix('calculating ')
                        # Uppercase only the first character; ``capitalize()``
                        # lowercases the rest, turning ``GB`` into ``Gb``.
                        label = label.removesuffix(' contribution...')
                        label = label[:1].upper() + label[1:]
                        pb_thread = threading.Thread(
                            target=pb,
                            args=(self.output_files[i], nframes, self.mpi_size, nmode),
                            kwargs={'style': self.progress_style, 'label': label},
                            daemon=True,
                        )
                        pb_thread.start()

                if isinstance(calc, (EnergyCalculation, ListEnergyCalculation)):
                    # The master progress monitor reports streaming QM/MM
                    # diagnostics while mdout files grow. Suppress the old
                    # per-rank end-of-run warning path when that monitor is
                    # active; it remains the fallback when progress is off.
                    calc._qmmm_diagnostics_streamed = self.progress_style != 'none'
                calc.setup()
                calc.run(rank, stdout=stdout, stderr=stderr)
                if self.timer_keys[i] is not None:
                    self.timer.stop_timer(self.timer_keys[i])
                    if pb_thread:
                        pb_thread.join()
        finally:
            if own_handle: f.close()


class MultiCalculation(object):
    def __init__(self):
        self.list_calc = []
        self.postprocess_prmtop = None
        self.keep_mdouts = False

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
                if calc_failed:
                    raise CalcError(f'{command_args[0]} failed with prmtop {command_args[1]}!')
                # Each file of gbnsr6 with decomp is huge, so we need to reduce it. Here we transform the file to json
                # to make it small
                if 'gbnsr6' in command_args[0]:
                    postprocess_args = command_args
                    if self.postprocess_prmtop is not None:
                        postprocess_args = command_args.copy()
                        postprocess_args[postprocess_args.index('-p') + 1] = str(self.postprocess_prmtop)
                    mdout2json(postprocess_args, keep_mdout=self.keep_mdouts)
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


class EnergyCalculation(Calculation):
    """ Uses mmpbsa_py_energy to evaluate energies """

    def __init__(self, prog, prmtop, incrd, inptraj, input_file, output, restrt, xvv=None):
        Calculation.__init__(self, prog, prmtop, incrd, inptraj,
                             input_file, output, xvv)
        self.restrt = restrt

    def run(self, rank, stdout=sys.stdout, stderr=sys.stderr):
        try:
            Calculation.run(self, rank, stdout, stderr)
        except CalcError as error:
            diagnostics = self._read_qmmm_diagnostics()
            if diagnostics:
                raise CalcError(self._format_qmmm_error(diagnostics, error)) from error
            raise
        self._check_qmmm_convergence(
            emit_warnings=rank == 0 and not getattr(self, '_qmmm_diagnostics_streamed', False),
        )

    def _check_qmmm_convergence(self, emit_warnings=True):
        """Reject fatal QM/MM diagnostics and report nonfatal warnings."""
        diagnostics = self._read_qmmm_diagnostics()
        fatal = [diagnostic for diagnostic in diagnostics if diagnostic.severity == 'error']
        if emit_warnings:
            for diagnostic in diagnostics:
                if diagnostic.severity == 'warning':
                    logging.warning('QM/MM diagnostic: %s %s', diagnostic.message, diagnostic.remediation)
        if fatal:
            raise CalcError(self._format_qmmm_error(diagnostics))

    def _read_qmmm_diagnostics(self):
        """Read and classify the current QM/MM SANDER output, if applicable."""
        try:
            input_text = Path(self.input_file).read_text(errors='replace').lower()
        except (OSError, TypeError):
            return []
        if '&qmmm' not in input_text:
            return []

        try:
            output_index = self.command_args.index('-o') + 1
            output_file = Path(self.command_args[output_index])
            output_text = output_file.read_text(errors='replace')
        except (OSError, TypeError, ValueError, IndexError):
            return []

        qm_theory = None
        theory_match = re.search(r"\bqm_theory\s*=\s*['\"]?([^'\"\s,]+)", input_text, re.IGNORECASE)
        if theory_match:
            qm_theory = theory_match.group(1)
        return parse_qmmm_diagnostics(output_text, qm_theory=qm_theory)

    def _format_qmmm_error(self, diagnostics, original_error=None):
        """Format fatal QM/MM diagnostics without hiding useful SANDER context."""
        try:
            output_index = self.command_args.index('-o') + 1
            output_file = Path(self.command_args[output_index])
        except (ValueError, IndexError):
            output_file = Path(self.output)

        lines = [f'QM/MM calculation failed; inspect SANDER output: {output_file}']
        for diagnostic in diagnostics:
            prefix = diagnostic.severity.upper()
            lines.append(f'{prefix} [{diagnostic.code}]: {diagnostic.message}')
            lines.append(f'Remediation: {diagnostic.remediation}')
        if original_error is not None:
            lines.append(f'Original SANDER error: {original_error}')
        return ' '.join(lines)

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

        # Verify that the calculation input file exists.
        if not os.path.exists(self.input_file):
            raise IOError("Input file (%s) doesn't exist" % self.input_file)

        self.calc_setup = True


class ListEnergyCalculation(MultiCalculation):
    def __init__(self, prog, prmtop, input_file, incrds, outputs, xvv=None, postprocess_prmtop=None,
                 keep_mdouts=False):
        super().__init__()
        self.program = prog
        self.prmtop = prmtop
        self.incrds = incrds
        self.input_file = input_file
        self.outputs = outputs
        self.xvv = xvv
        self.postprocess_prmtop = postprocess_prmtop
        self.keep_mdouts = keep_mdouts

    def setup(self):
        """
        Sets up the command-line arguments. Sander requires a unique restrt file
        for the MPI version (since one is *always* written and you don't want 2
        threads fighting to write the same dumb file)
        """
        for c, o in zip(self.incrds, self.outputs):
            command_args = [self.program,
                            '-i', self.input_file,  # input file flag
                            '-p', self.prmtop  # prmtop flag
                            ]
            command_args.extend(('-c', c))  # input coordinate flag
            command_args.extend(('-o', o))  # output file flag
            self.list_calc.append(command_args)

            # Input-file validation is handled when the calculation is prepared.
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


class MolsurfCalc(SurfCalc):
    """ Uses molsurf to calculate the surface area """

    def __init__(self, prog, prmtop, inptraj, output, probe=1.4, offset=0.0):
        SurfCalc.__init__(self, prog, prmtop, inptraj, output)
        self.probe = probe
        self.offset = offset

    def _get_instring(self, rank):
        inptraj = self.inptraj % rank
        output = self.output % rank
        return "trajin %s\nmolsurf :* out %s probe %s offset %s\n" % (inptraj,
                                                                       output, self.probe, self.offset)


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
        MergeGBNSR6Output(self.topology, out_filename, mm_filename, self.mdouts, self.idecomp, self.dec_verbose)


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


def _entropy_block_analysis(energies, estimator):
    """Evaluate an entropy estimator over deterministic nonoverlapping blocks."""
    energies = np.asarray(energies, dtype=float)
    numframes = energies.size
    if not numframes:
        return []

    if numframes == 1:
        block_sizes = [1]
    else:
        block_sizes = sorted({max(2, numframes // divisor) for divisor in (16, 8, 4, 2, 1)})

    analysis = []
    for block_size in block_sizes:
        nblocks = numframes // block_size
        if not nblocks:
            continue
        estimates = np.asarray([
            estimator(energies[start:start + block_size])
            for start in range(0, nblocks * block_size, block_size)
        ], dtype=float)
        block_std = float(estimates.std(ddof=1)) if nblocks > 1 else float('nan')
        block_sem = block_std / math.sqrt(nblocks) if nblocks > 1 else float('nan')
        percentile = np.percentile(estimates, [2.5, 97.5])
        analysis.append({
            'block_size': block_size,
            'nblocks': nblocks,
            'used_frames': nblocks * block_size,
            'mean': float(estimates.mean()),
            'std': block_std,
            'sem': block_sem,
            'p025': float(percentile[0]),
            'p975': float(percentile[1]),
        })
    return analysis


def _select_block_diagnostic(block_analysis):
    """Select the largest block size retaining at least eight block estimates."""
    for minimum_blocks in (8, 2, 1):
        eligible = [row for row in block_analysis if row['nblocks'] >= minimum_blocks]
        if eligible:
            return max(eligible, key=lambda row: row['block_size'])
    return None


class InteractionEntropyCalc:
    """
    Class for Interaction Entropy calculation
    :return {IE_key: data}
    """

    # Molar Boltzmann constant (gas constant) in kcal/(mol*K).
    GAS_CONSTANT = 0.00198720425864083

    def __init__(self, ggas, INPUT, method, iesegment=None):
        """

        Args:
            ggas: Model GGAS energy
            INPUT: INPUT dict
            iesegment: If not defined, iesegment = INPUT['general']['ie_segment'].
                Use ``None`` for the input default; ``0`` is a valid diagnostic
                segment (empty IE tail) and must not fall back via truthiness.
        """
        self.ggas = ggas
        self.INPUT = INPUT
        self.method = method
        self.isegment = (
            INPUT['general']['ie_segment'] if iesegment is None else iesegment
        )
        self.data = []

        self._calculate()

    def _calculate(self):
        temperature = self.INPUT['general']['temperature']
        kT = self.GAS_CONSTANT * temperature

        self.data = np.zeros(self.ggas.size, dtype=float)
        running_energy_sum = 0.0
        running_logsumexp = -np.inf

        for i, energy in enumerate(tqdm(self.ggas, bar_format=TQDM_BAR_FORMAT, ascii=True)):
            nframes = i + 1
            running_energy_sum += energy
            mean_energy = running_energy_sum / nframes
            running_logsumexp = np.logaddexp(running_logsumexp, energy / kT)

            # Eq. 6-8 of Duan et al., JACS 2016, 138, 5722-5728:
            # kT * log(mean(exp((E - <E>) / kT))). Rewriting it as
            # kT * (logsumexp(E / kT) - log(N)) - <E> evaluates the same
            # expression for every trajectory prefix without exponential
            # overflow and consistently recenters every included frame.
            self.data[i] = kT * (running_logsumexp - math.log(nframes)) - mean_energy

        # Jensen's inequality makes IE non-negative. Remove only numerical
        # roundoff below zero; genuine negative values indicate an algorithmic
        # error and must not be generated by the expression above.
        self.data = np.maximum(self.data, 0.0)

        numframes = len(self.data)
        self.ie_std = float(self.ggas.std())
        self.ieframes = math.ceil(numframes * (self.isegment / 100))
        # ``data[-0:]`` is the full array in Python; keep an empty tail for 0%.
        self.iedata = self.data[-self.ieframes:] if self.ieframes else np.asarray([], dtype=float)
        self.ie_value = float(self.data[-1]) if numframes else float('nan')
        self.tail_mean = float(self.iedata.mean()) if self.ieframes else float('nan')
        self.tail_std = float(self.iedata.std()) if self.ieframes else float('nan')
        self.block_analysis = _entropy_block_analysis(
            self.ggas, lambda block: self._estimate(block, kT)
        )
        block_diagnostic = _select_block_diagnostic(self.block_analysis)
        self.block_size = block_diagnostic['block_size'] if block_diagnostic else 0
        self.block_nblocks = block_diagnostic['nblocks'] if block_diagnostic else 0
        self.block_std = block_diagnostic['std'] if block_diagnostic else float('nan')
        self.block_sem = block_diagnostic['sem'] if block_diagnostic else float('nan')

    @staticmethod
    def _estimate(energies, kT):
        """Return the full-ensemble IE estimate for one energy block."""
        energies = np.asarray(energies, dtype=float)
        scaled = energies / kT
        maximum = scaled.max()
        estimate = kT * (maximum + math.log(np.exp(scaled - maximum).mean())) - energies.mean()
        return max(0.0, float(estimate))

    def save_output(self, filename):
        frames = list(
            range(
                self.INPUT['general']['startframe'],
                self.INPUT['general']['startframe'] + len(self.data) * self.INPUT['general']['interval'],
                self.INPUT['general']['interval'],
            )
        )
        with open(filename, 'w') as out:
            out.write(f'| Interaction Entropy results for {self.method} calculations\n')
            out.write(f'IE-frames: last {self.ieframes}\n')
            out.write(f'Internal Energy SD (sigma): {self.ie_std:9.2f}\n')
            out.write(f'Full-ensemble Interaction Entropy (-TΔS): {self.ie_value:9.4f}\n')
            out.write(f'Tail convergence mean (last {self.ieframes} prefixes): {self.tail_mean:9.4f}\n')
            out.write(f'Tail convergence SD (last {self.ieframes} prefixes): {self.tail_std:9.4f}\n')
            out.write(
                f'Block diagnostic: size {self.block_size} blocks {self.block_nblocks} '
                f'SD {self.block_std:.4f} SEM {self.block_sem:.4f}\n\n'
            )
            out.write('| Nonoverlapping block analysis: method block_size blocks used_frames mean SD SEM P2.5 P97.5\n')
            for row in self.block_analysis:
                out.write(
                    f"| BLOCK IE {row['block_size']} {row['nblocks']} {row['used_frames']} "
                    f"{row['mean']:.6f} {row['std']:.6f} {row['sem']:.6f} "
                    f"{row['p025']:.6f} {row['p975']:.6f}\n"
                )
            out.write('\n')
            out.write('| Interaction Entropy per-frame:\n')

            out.write('Frame # | IE value\n')
            for f, d in zip(frames, self.data):
                out.write('{:d}  {:.2f}\n'.format(f, d))


class C2EntropyCalc:
    """
    Class for Interaction Entropy calculation
    :return {IE_key: data}
    """

    def __init__(self, ggas, INPUT, method):
        self.ggas = ggas
        self.INPUT = INPUT
        self.method = method

        self._calculate()

    def _calculate(self):
        R = InteractionEntropyCalc.GAS_CONSTANT
        temperature = self.INPUT['general']['temperature']
        self.ie_std = float(self.ggas.std())
        self.c2data = (self.ie_std ** 2) / (2 * temperature * R)
        self.block_analysis = _entropy_block_analysis(
            self.ggas,
            lambda block: float(block.std() ** 2 / (2 * temperature * R)),
        )
        block_diagnostic = _select_block_diagnostic(self.block_analysis)
        self.block_size = block_diagnostic['block_size'] if block_diagnostic else 0
        self.block_nblocks = block_diagnostic['nblocks'] if block_diagnostic else 0
        self.c2_std = block_diagnostic['std'] if block_diagnostic else float('nan')
        self.c2_sem = block_diagnostic['sem'] if block_diagnostic else float('nan')
        self.c2_ci = np.asarray([
            block_diagnostic['p025'], block_diagnostic['p975']
        ]) if block_diagnostic else np.asarray([float('nan'), float('nan')])

    def save_output(self, filename):
        with open(filename, 'w') as out:
            out.write(f'| C2 Entropy results for {self.method} calculations\n')
            out.write(f'C2 Entropy (-TΔS): {self.c2data:.4f}\n')
            out.write(f'C2 Block SD: {self.c2_std:.4f}\n')
            out.write(f'C2 Block SEM: {self.c2_sem:.4f}\n')
            out.write(f'Internal Energy SD (sigma): {self.ie_std:9.2f}\n')
            out.write(f'C2 Block P2.5-P97.5: {self.c2_ci[0]:.4f} {self.c2_ci[1]:.4f}\n')
            out.write(f'Block diagnostic: size {self.block_size} blocks {self.block_nblocks}\n')
            out.write('| Nonoverlapping block analysis: method block_size blocks used_frames mean SD SEM P2.5 P97.5\n')
            for row in self.block_analysis:
                out.write(
                    f"| BLOCK C2 {row['block_size']} {row['nblocks']} {row['used_frames']} "
                    f"{row['mean']:.6f} {row['std']:.6f} {row['sem']:.6f} "
                    f"{row['p025']:.6f} {row['p975']:.6f}\n"
                )


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
    """Merge sander MM (+ optional decomp) with GBNSR6 energies into one mdout-like file.

    Totals replace EEL / 1-4 EEL / EGB from GBNSR6. Decomposition keeps sander
    internal/vdw/eel/sas and overwrites only the polar column from GBNSR6 DGij.
    """
    def __init__(self, topology, output_filename, mm_filename, mdout_filenames, idecomp, dec_verbose):
        self.topology = topology
        self.output_filename = output_filename
        self.mm_filename = mm_filename
        self.mdout_filenames = mdout_filenames
        self.idecomp = idecomp

        self.dec_verbose = dec_verbose
        # Totals: EEL/1-4 EEL/EGB from GBNSR6. Decomp: eel from sander, pol from GBNSR6 DGij.
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
        temp_res2print = []

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
                    temp_res2print.append(line.split()[1:])

        res2print = []
        if temp_res2print:
            for l in temp_res2print:
                for i in range(0, len(l), 2):
                    res2print.extend(range(int(l[i]), int(l[i + 1]) + 1))

        return {'file_assignments': file_assignments, 'inputfile': inputfile, 'resource_section': resource_section,
                'control_data': control_data, 'atomic_coor_vel': atomic_coor_vel,
                'results_section': self._get_energy_decomp(results_section), 'res2print': res2print}

    @staticmethod
    def _get_energy_decomp(results_section):
        import re

        energy = {}
        decomp = {}
        energy_term = re.compile(
            r'([A-Z0-9-]+(?:\s+[A-Z0-9-]+)?)\s*=\s*'
            r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|\*+)'
        )

        # A single MM output can contain several local coordinate sets when
        # frames are distributed over MPI ranks.  The old stateful parser
        # advanced through the energy rows and could leave only the last
        # frame in ``energy``.  Parse each coordinate-set block independently
        # so every frame survives the GBNSR6 merge.
        current_frame = None
        for line in results_section:
            if line.startswith('minimizing coord set #'):
                current_frame = int(line.split()[-1])
                energy[current_frame] = {}
                decomp[current_frame] = {}
                continue

            if current_frame is None:
                continue

            for term, value in energy_term.findall(line):
                # Amber uses asterisks for overflowed values.  They are not
                # usable as energies, but must not prevent later frames from
                # being parsed.
                if '*' not in value:
                    energy[current_frame][term.strip()] = float(value)

            if line[:3] in ['TDC', 'SDC', 'BDC']:
                data = [x.strip().replace('->', '') for x in line.split()]
                if len(data) == 8:
                    _t, _r1, _r2, _i, _v, _e, _p, _s = data
                    data = [_t, int(_r1), int(_r2), float(_i), float(_v), float(_e), float(_p), float(_s)]
                else:
                    _t, _r1, _i, _v, _e, _p, _s = data
                    data = [_t, int(_r1), float(_i), float(_v), float(_e), float(_p), float(_s)]
                if not decomp[current_frame].get(line[:3]):
                    decomp[current_frame][line[:3]] = [data]
                else:
                    decomp[current_frame][line[:3]].append(data)
        return {'energy': energy, 'decomp':decomp}

    def read_gbnsr6_output(self, res2print):

        energy = {}
        decomp = {}
        for i, filename in enumerate(self.mdout_filenames, start=1):
            jsonfilename = Path(filename).with_suffix('.json')
            with open(jsonfilename, 'r') as openfile:
                data = json.load(openfile)
                if i == 1:
                    file_assignments = data['file_assignments']
                    inputfile = data['inputfile']
                energy[i] = data['results_section']['energy']
            if self.idecomp:
                decomp[i] = self.get_decomp(data['results_section']['decomp'], res2print)

        results = {'energy': energy, 'decomp':decomp}
        return {'file_assignments': file_assignments, 'inputfile': inputfile, 'results_section': results}

    def get_decomp(self, decomp, res2print):
        d = {}
        if self.idecomp in [1, 2]:
            for res1, v1 in decomp.items():
                res1 = int(res1)
                if res1 not in res2print:
                    continue
                d[res1] = {'TDC': 0.0}
                if self.dec_verbose in [1, 3]:
                    d[res1]['BDC'] = 0.0
                    d[res1]['SDC'] = 0.0
                for res2, de in v1.items():
                    res2 = int(res2)
                    for t, e in de.items():
                        if t in ['BDC', 'SDC'] and self.dec_verbose not in [1, 3]:
                            continue
                        d[res1][t] += e if res1 == res2 else e/2
        else:
            for res1, v1 in decomp.items():
                res1 = int(res1)
                if res1 not in res2print:
                    continue
                d[res1] = {}
                for res2, de in v1.items():
                    res2 = int(res2)
                    if res2 not in res2print:
                        continue
                    # if res1 == res2:
                    #     continue
                    d[res1][res2] = {'TDC': 0.0}
                    if self.dec_verbose in [1, 3]:
                        d[res1][res2]['BDC'] = 0.0
                        d[res1][res2]['SDC'] = 0.0
                    for t, e in de.items():
                        if t in ['BDC', 'SDC'] and self.dec_verbose not in [1, 3]:
                            continue
                        d[res1][res2][t] = e if res1 == res2 else e/2
        return d

    @staticmethod
    def _get_energy_gbnsr6(results_section):
        energy = {}

        store = False
        c = 0
        while True:
            line = results_section[c]
            if "FINAL RESULTS" in line:
                store = True
            if store and line.startswith(' 1-4 NB'):
                words = line.split()
                energy['1-4 EEL'] = float(words[7])
                c += 1
                line = results_section[c]
                words = line.split()
                energy['EEL'] = float(words[2])
                energy['EGB'] = float(words[5])
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

            frame_ids = sorted(mmenergy)
            if not frame_ids:
                raise CalcError(f'No MM energy frames were found in {self.mm_filename}')
            missing_gb_frames = [frame for frame in frame_ids if frame not in gbenergy]
            if missing_gb_frames:
                raise CalcError(
                    f'Missing GBNSR6 energy frames {missing_gb_frames} for {self.output_filename}'
                )

            k2print = [['BOND', 'ANGLE', 'DIHED']]
            if {'UB', 'IMP', 'CMAP'}.issubset(mmenergy[frame_ids[0]]):
                k2print.append(['UB', 'IMP', 'CMAP'])
            k2print.extend([['VDWAALS', 'EEL', 'EGB'], ['1-4 VDW', '1-4 EEL', 'RESTRAINT'], ['ESURF']])
            for frame in frame_ids:
                output_file.write(f'minimizing coord set #       {frame}\n\n')
                frame_energy = mmenergy[frame].copy()
                frame_energy.pop('EGB', None)
                frame_energy.pop('EEL', None)
                frame_energy.pop('1-4 EEL', None)
                frame_energy.update(gbenergy[frame])
                for kl in k2print:
                    if len(kl) == 3:
                        f = []
                        for klk in kl:
                            f.extend((klk, frame_energy[klk]))
                        output_file.write(' {:8s}={:>14.4f}  {:8s}={:>14.4f}  {:11s}={:>14.4f}\n'.format(*f))
                    else:
                        f = [kl[0], frame_energy[kl[0]]]
                        output_file.write(' {:8s}={:>14.4f}\n\n'.format(*f))

                if self.idecomp:
                    for term in mmdecomp.get(frame, {}):
                        if self.idecomp in [1, 2]:
                            output_file.write(self.decomp_headers['pr'].format(self.decomp_labels[term]))
                            for c, l in enumerate(mmdecomp[frame][term]):
                                r1 = mmdecomp[frame][term][c][1]
                                l[-2] = gbdecomp[frame][r1][term]
                                output_file.write('{}{:>7d}{:>10.3f}{:>10.3f}{:>10.3f}{:>10.3f}{:>10.3f}\n'.format(*l))
                        else:
                            output_file.write(self.decomp_headers['pw'].format(self.decomp_labels[term]))
                            for c, l in enumerate(mmdecomp[frame][term]):
                                r1, r2 = mmdecomp[frame][term][c][1:3]
                                l[-2] = gbdecomp[frame][r1][r2][term]
                                output_file.write('{}{:>8d}->{:>7d}{:>13.4f}{:>13.4f}{:>13.4f}{:>13.4f}{:>13.4f}\n'.format(*l))
                        output_file.write('\n')
