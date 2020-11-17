"""
 This is a module that contains the class of the main gmx_MMPBSA
 Application.
"""

# Import system modules

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

import os
import signal
import sys
import warnings
import numpy as np
from math import exp, log

# Import gmx_MMPBSA modules
from GMXMMPBSA import utils
from GMXMMPBSA.amber_outputs import (QHout, NMODEout, QMMMout, GBout, PBout,
                                     PolarRISM_std_Out, RISM_std_Out,
                                     PolarRISM_gf_Out, RISM_gf_Out,
                                     SingleTrajBinding, MultiTrajBinding)
from GMXMMPBSA.calculation import (CalculationList, EnergyCalculation,
                                   PBEnergyCalculation, RISMCalculation,
                                   NmodeCalc, QuasiHarmCalc, CopyCalc,
                                   PrintCalc, LcpoCalc, MolsurfCalc)
from GMXMMPBSA.commandlineparser import parser
from GMXMMPBSA.createinput import create_inputs
from GMXMMPBSA.exceptions import (MMPBSA_Error, InternalError, InputError,
                                  InputWarning)
from GMXMMPBSA.fake_mpi import MPI as FakeMPI
from GMXMMPBSA.findprogs import find_progs
from GMXMMPBSA.infofile import InfoFile
from GMXMMPBSA.input_parser import input_file as _input_file
from GMXMMPBSA.make_trajs import make_trajectories, make_mutant_trajectories
from GMXMMPBSA.output_file import (write_stability_output,
                                   write_binding_output,
                                   write_decomp_stability_output,
                                   write_decomp_binding_output)
from GMXMMPBSA.parm_setup import MMPBSA_System
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.timer import Timer
from GMXMMPBSA.gui import run as GUI_run

# Global variables for the excepthook replacement at the bottom. Override these
# in the MMPBSA_App constructor and input file reading
_unbuf_stdout = utils.Unbuffered(sys.stdout)  # unbuffered stdout
_unbuf_stderr = utils.Unbuffered(sys.stderr)  # unbuffered stderr
_stdout = sys.stdout
_stderr = sys.stderr
_debug_printlevel = 2
_mpi_size = 1
_rank = 0
_MPI = FakeMPI()


# Main class

class MMPBSA_App(object):
    """ Main MM/PBSA application for driving the entire calculation """
    # The command line parser and input file objects are class attributes here
    clparser = parser
    input_file = _input_file
    debug_printlevel = 2

    def __init__(self, MPI, stdout=None, stderr=None, size=None):
        """
        Sets up the main gmx_MMPBSA driver class. All we set up here is the output
        and error streams (unbuffered by default) and the prefix for the
        intermediate files. Also set up empty INPUT dict
        """
        global _rank, _stdout, _stderr, _mpi_size, _MPI
        _MPI = self.MPI = MPI
        self.pre = '_MMPBSA_'
        self.INPUT = {}
        if stdout is None:
            _stdout = self.stdout = _unbuf_stdout
        else:
            _stdout = self.stdout = stdout

        if stderr is None:
            _stderr = self.stderr = _unbuf_stderr
        else:
            _stderr = self.stderr = stderr

        # MPI-related variables. Squash output for non-master threads
        _rank = self.mpi_rank = self.MPI.COMM_WORLD.Get_rank()
        self.master = self.mpi_rank == 0
        _mpi_size = self.mpi_size = self.MPI.COMM_WORLD.Get_size()
        if not self.master:
            self.stdout = open(os.devnull, 'w')

        # Set up timers
        timers = [Timer() for i in range(self.mpi_size)]
        self.timer = timers[self.mpi_rank]

        # Support possible threading for those that don't use MPI. However, if
        # mpi_size is > 1, just use the MPI mechanism instead
        if size is not None and self.mpi_size == 1:
            self.mpi_size = size

    def file_setup(self):
        """ Sets up the trajectories and input files """
        # If we are rewriting the output file only, bail out here
        if self.FILES.rewrite_output:
            return
        # This work belongs to the 'setup' timer
        self.timer.start_timer('setup')
        if not hasattr(self, 'normal_system'):
            raise InternalError('MMPBSA_App not set up and parms not checked!')
        # Set up some local refs for convenience
        FILES, INPUT, master = self.FILES, self.INPUT, self.master

        # # Now we're getting ready, remove existing intermediate files
        # if master and FILES.use_mdins:
        #     self.remove(-1)
        # elif master and not FILES.rewrite_output:
        #     self.remove(0)

        # Create input files based on INPUT dict
        if master and not FILES.use_mdins:
            create_inputs(INPUT, self.normal_system, self.pre)
        self.timer.stop_timer('setup')
        # Bail out if we only wanted to generate mdin files
        if FILES.make_mdins:
            self.stdout.write('Created mdin files. Quitting.\n')
            sys.exit(0)

        # Now create our trajectory files

        self.timer.add_timer('cpptraj', 'Creating trajectories with cpptraj:')
        self.timer.start_timer('cpptraj')

        if master:
            self.stdout.write('Preparing trajectories for simulation...\n')
            self.numframes, self.numframes_nmode = make_trajectories(INPUT, FILES, self.mpi_size,
                                                                     self.external_progs['cpptraj'].full_path, self.pre)

        self.MPI.COMM_WORLD.Barrier()

        self.timer.stop_timer('cpptraj')

        self.timer.add_timer('muttraj', 'Mutating trajectories:')
        self.timer.start_timer('muttraj')

        if INPUT['alarun']:
            self.stdout.write('Mutating trajectories...\n')
        self.mut_str, mutant_residue = make_mutant_trajectories(INPUT, FILES,
                                                                self.mpi_rank, self.external_progs['cpptraj'].full_path,
                                                                self.normal_system, self.mutant_system, self.pre)

        self.MPI.COMM_WORLD.Barrier()

        if master:
            self.stdout.write(('%d frames were processed by cpptraj for use in '
                               'calculation.\n') % self.numframes)
            if INPUT['nmoderun']:
                self.stdout.write(('%d frames were processed by cpptraj for '
                                   'nmode calculations.\n') % self.numframes_nmode)

        self.timer.stop_timer('muttraj')

        # Add all of the calculation timers
        self.timer.add_timer('calc', 'Total calculation time:')
        if INPUT['gbrun']:
            self.timer.add_timer('gb', 'Total GB calculation time:')
        if INPUT['pbrun']:
            self.timer.add_timer('pb', 'Total PB calculation time:')
        if INPUT['rismrun']:
            self.timer.add_timer('rism', 'Total 3D-RISM calculation time:')
        if INPUT['nmoderun']:
            self.timer.add_timer('nmode', 'Total normal mode calculation time:')
        if INPUT['entropy'] == 1:
            self.timer.add_timer('qh', 'Total quasi-harmonic calculation time:')

        self.sync_mpi()

    def run_mmpbsa(self, rank=None):
        """
        Runs the MM/PBSA analysis. This assumes FILES and INPUT are already set.
        """

        if not hasattr(self, 'external_progs'):
            raise InternalError('external_progs not declared in run_mmpbsa!')

        FILES, INPUT = self.FILES, self.INPUT
        if rank is None:
            rank = self.mpi_rank
        master = rank == 0

        # Load the list of calculations we need to do, then run them.

        if master:
            self.timer.start_timer('calc')

        self.load_calc_list()

        self.stdout.write('\n')

        self.calc_list.run(rank, self.stdout)

        self.sync_mpi()

        if master: self.timer.stop_timer('calc')

        # Write out the info file now
        if master:
            info = InfoFile(self)
            info.write_info(self.pre + 'info')

    def load_calc_list(self):
        """
        Sets up all of the calculations to be run. When adding a new
        calculation type, add a class to calculation.py, import it at the top of
        the file here, then append it to the calc list appropriately
        """
        self.calc_list = CalculationList(self.timer)

        if not self.INPUT['mutant_only']:
            self.calc_list.append(
                PrintCalc('Running calculations on normal system...'),
                timer_key=None)
            self._load_calc_list(self.pre, False, self.normal_system)
        if self.INPUT['alarun']:
            self.calc_list.append(
                PrintCalc('\nRunning calculations on mutant system...'),
                timer_key=None)
            self._load_calc_list(self.pre + 'mutant_', True, self.mutant_system)

    def _load_calc_list(self, prefix, mutant, parm_system):
        """
        Internal routine to handle building calculation list. Called separately
        for mutant and normal systems
        """
        # Set up a dictionary of external programs to use based one external progs
        progs = {'gb': self.external_progs['mmpbsa_py_energy'].full_path,
                 'sa': self.external_progs['cpptraj'].full_path,
                 'pb': self.external_progs['mmpbsa_py_energy'].full_path,
                 'rism': self.external_progs['rism3d.snglpnt'].full_path,
                 'qh': self.external_progs['cpptraj'].full_path,
                 'nmode': self.external_progs['mmpbsa_py_nabnmode'].full_path
                 }
        if self.INPUT['use_sander'] or self.INPUT['decomprun']:
            progs['gb'] = progs['pb'] = self.external_progs['sander'].full_path
        if self.INPUT['sander_apbs']:
            progs['pb'] = self.external_progs['sander.APBS'].full_path
        if self.INPUT['ifqnt']:
            progs['gb'] = self.external_progs['sander'].full_path

        # NetCDF or ASCII intermediate trajectories?
        if self.INPUT['netcdf']:
            trj_sfx = 'nc'
        else:
            trj_sfx = 'mdcrd'

        # Determine if we just copy the receptor files. This only happens if we
        # are doing mutant calculations, we're not only doing the mutant, and the
        # receptor/mutant receptor topologies are equal. Same for the ligand
        copy_receptor = (mutant and not self.INPUT['mutant_only'] and
                         self.FILES.receptor_prmtop == self.FILES.mutant_receptor_prmtop)
        copy_ligand = (mutant and not self.INPUT['mutant_only'] and
                       self.FILES.ligand_prmtop == self.FILES.mutant_ligand_prmtop)

        # First load the GB calculations
        if self.INPUT['gbrun']:
            # See if we need a PDB or restart file for the inpcrd
            if 'mmpbsa_py_energy' in progs['gb']:
                incrd = '%s%%s.pdb' % prefix
            else:
                incrd = '%sdummy%%s.inpcrd' % prefix

            # See whether we are doing molsurf or LCPO. Reduce # of arguments
            # needed to 3, filling in the others here
            if self.INPUT['molsurf']:
                SAClass = lambda a1, a2, a3: MolsurfCalc(progs['sa'], a1, a2, a3,
                                                         self.INPUT['probe'], self.INPUT['msoffset'])
            else:
                SAClass = lambda a1, a2, a3: LcpoCalc(progs['sa'], a1, a2, a3,
                                                      self.INPUT['probe'], self.INPUT['msoffset'])

            # Mdin depends on decomp or not
            if self.INPUT['decomprun']:
                mdin_template = self.pre + 'gb_decomp_%s.mdin'
            elif self.INPUT['ifqnt']:
                mdin_template = self.pre + 'gb_qmmm_%s.mdin'
            else:
                mdin_template = self.pre + 'gb.mdin'

            # Now do complex-specific stuff
            try:
                mdin = mdin_template % 'com'
            except TypeError:
                mdin = mdin_template

            self.calc_list.append(PrintCalc('\nBeginning GB calculations with %s' %
                                            progs['gb']), timer_key='gb')

            c = EnergyCalculation(progs['gb'], parm_system.complex_prmtop,
                                  incrd % 'complex',
                                  '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                  mdin, '%scomplex_gb.mdout.%%d' % (prefix),
                                  self.pre + 'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...',
                                  timer_key='gb')
            c = SAClass(parm_system.complex_prmtop,
                        '%scomplex.%s.%%d' % (prefix, trj_sfx),
                        '%scomplex_gb_surf.dat.%%d' % prefix)
            self.calc_list.append(c, '', timer_key='gb')

            if not self.stability:
                try:
                    mdin = mdin_template % 'rec'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing receptor if the mutation is in the ligand
                # or perform a receptor calculation
                if copy_receptor:
                    c = CopyCalc('%sreceptor_gb.mdout.%%d' % self.pre,
                                 '%sreceptor_gb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='gb')
                    c = CopyCalc('%sreceptor_gb_surf.dat.%%d' % self.pre,
                                 '%sreceptor_gb_surf.dat.%%d' % prefix)
                    self.calc_list.append(c, '', timer_key='gb')
                else:
                    c = EnergyCalculation(progs['gb'], parm_system.receptor_prmtop,
                                          incrd % 'receptor',
                                          '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                                          mdin, '%sreceptor_gb.mdout.%%d' % (prefix),
                                          self.pre + 'restrt.%d')
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='gb')
                c = SAClass(parm_system.receptor_prmtop,
                            '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                            '%sreceptor_gb_surf.dat.%%d' % prefix)
                self.calc_list.append(c, '', timer_key='gb')

                try:
                    mdin = mdin_template % 'lig'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing ligand if the mutation is in the receptor
                # or perform a ligand calculation
                if copy_ligand:
                    c = CopyCalc('%sligand_gb.mdout.%%d' % self.pre,
                                 '%sligand_gb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='gb')
                    c = CopyCalc('%sligand_gb_surf.dat.%%d' % self.pre,
                                 '%sligand_gb_surf.dat.%%d' % prefix)
                    self.calc_list.append(c, '', timer_key='gb')
                else:
                    c = EnergyCalculation(progs['gb'], parm_system.ligand_prmtop,
                                          incrd % 'ligand',
                                          '%sligand.%s.%%d' % (prefix, trj_sfx),
                                          mdin, '%sligand_gb.mdout.%%d' % (prefix),
                                          self.pre + 'restrt.%d')
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='gb')
                c = SAClass(parm_system.ligand_prmtop,
                            '%sligand.%s.%%d' % (prefix, trj_sfx),
                            '%sligand_gb_surf.dat.%%d' % prefix)
                self.calc_list.append(c, '', timer_key='gb')

        # end if self.INPUT['gbrun']

        # Next load the PB calculations
        if self.INPUT['pbrun']:
            # See if we need a PDB or restart file for the inpcrd
            if 'mmpbsa_py_energy' in progs['pb']:
                incrd = '%s%%s.pdb' % prefix
            else:
                incrd = '%sdummy%%s.inpcrd' % prefix

            # Mdin depends on decomp or not
            if self.INPUT['decomprun']:
                mdin_template = self.pre + 'pb_decomp_%s.mdin'
                mdin_template2 = mdin_template
            else:
                mdin_template = self.pre + 'pb.mdin'
                mdin_template2 = self.pre + 'pb.mdin2'

            # Now do complex-specific stuff
            try:
                mdin = mdin_template % 'com'
            except TypeError:
                mdin = mdin_template

            self.calc_list.append(PrintCalc('\nBeginning PB calculations with %s' %
                                            progs['pb']), timer_key='pb')

            c = PBEnergyCalculation(progs['pb'], parm_system.complex_prmtop,
                                    incrd % 'complex',
                                    '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                    mdin, '%scomplex_pb.mdout.%%d' % (prefix),
                                    self.pre + 'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...',
                                  timer_key='pb')
            if not self.stability:
                try:
                    mdin = mdin_template % 'rec'
                except TypeError:
                    mdin = mdin_template

                # Either copy the existing receptor if the mutation is in the ligand
                # or perform a receptor calculation
                if copy_receptor:
                    c = CopyCalc('%sreceptor_pb.mdout.%%d' % self.pre,
                                 '%sreceptor_pb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = PBEnergyCalculation(progs['pb'], parm_system.receptor_prmtop,
                                            incrd % 'receptor',
                                            '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                                            mdin, '%sreceptor_pb.mdout.%%d' % (prefix),
                                            self.pre + 'restrt.%d')
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='pb')

                try:
                    mdin2 = mdin_template2 % 'lig'
                except TypeError:
                    mdin2 = mdin_template2

                # Either copy the existing ligand if the mutation is in the receptor
                # or perform a ligand calculation
                if copy_ligand:
                    c = CopyCalc('%sligand_pb.mdout.%%d' % self.pre,
                                 '%sligand_pb.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = PBEnergyCalculation(progs['pb'], parm_system.ligand_prmtop,
                                            incrd % 'ligand',
                                            '%sligand.%s.%%d' % (prefix, trj_sfx),
                                            mdin2, '%sligand_pb.mdout.%%d' % (prefix),
                                            self.pre + 'restrt.%d')
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='pb')
        # end if self.INPUT['pbrun']

        if self.INPUT['rismrun']:
            self.calc_list.append(
                PrintCalc('\nBeginning 3D-RISM calculations with %s' %
                          progs['rism']), timer_key='rism')

            c = RISMCalculation(progs['rism'], parm_system.complex_prmtop,
                                '%scomplex.pdb' % prefix, '%scomplex.%s.%%d' %
                                (prefix, trj_sfx), self.FILES.xvvfile,
                                '%scomplex_rism.mdout.%%d' % prefix, self.INPUT)
            self.calc_list.append(c, '  calculating complex contribution...',
                                  timer_key='rism')

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('%sreceptor_rism.mdout.%%d' % self.pre,
                                 '%sreceptor_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = RISMCalculation(progs['rism'], parm_system.receptor_prmtop,
                                        '%sreceptor.pdb' % prefix,
                                        '%sreceptor.%s.%%d' % (prefix, trj_sfx),
                                        self.FILES.xvvfile,
                                        '%sreceptor_rism.mdout.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='rism')

                if copy_ligand:
                    c = CopyCalc('%sligand_rism.mdout.%%d' % self.pre,
                                 '%sligand_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = RISMCalculation(progs['rism'], parm_system.ligand_prmtop,
                                        '%sligand.pdb' % prefix,
                                        '%sligand.%s.%%d' % (prefix, trj_sfx),
                                        self.FILES.xvvfile,
                                        '%sligand_rism.mdout.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='rism')

        # end if self.INPUT['rismrun']

        if self.INPUT['nmoderun']:
            self.calc_list.append(
                PrintCalc('\nBeginning nmode calculations with %s' %
                          progs['nmode']), timer_key='nmode')

            c = NmodeCalc(progs['nmode'], parm_system.complex_prmtop,
                          '%scomplex.pdb' % prefix,
                          '%scomplex_nm.%s.%%d' % (prefix, trj_sfx),
                          '%scomplex_nm.out.%%d' % prefix, self.INPUT)
            self.calc_list.append(c, '  calculating complex contribution...',
                                  timer_key='nmode')

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('%sreceptor_nm.out.%%d' % self.pre,
                                 '%sreceptor_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.receptor_prmtop,
                                  '%sreceptor.pdb' % prefix,
                                  '%sreceptor_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sreceptor_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='rism')

                if copy_ligand:
                    c = CopyCalc('%sligand_nm.out.%%d' % self.pre,
                                 '%sligand_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.ligand_prmtop,
                                  '%sligand.pdb' % prefix,
                                  '%sligand_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sligand_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='rism')

        # end if self.INPUT['nmoderun']

        # Only master does entropy calculations
        if self.INPUT['entropy'] == 1:
            self.calc_list.append(
                PrintCalc('\nBeginning quasi-harmonic calculations with %s' %
                          progs['qh']), timer_key='qh')

            c = QuasiHarmCalc(progs['qh'], parm_system.complex_prmtop,
                              '%scomplex.%s' % (prefix, trj_sfx),
                              '%scpptrajentropy.in' % prefix,
                              '%scpptraj_entropy.out' % prefix,
                              self.INPUT['receptor_mask'],
                              self.INPUT['ligand_mask'], self.pre)
            self.calc_list.append(c, '', timer_key='qh')

    def loadcheck_prmtops(self):
        """ Loads the topology files and checks their consistency """



        # Start setup timer and make sure we've already set up our input
        self.timer.add_timer('setup', 'Total setup time:')
        self.timer.start_timer('setup')
        if not hasattr(self, 'FILES') or not hasattr(self, 'INPUT'):
            raise InternalError('MMPBSA_App not set up! Cannot check parms yet!')
        # create local aliases to avoid abundant selfs
        FILES, INPUT = self.FILES, self.INPUT
        # Now we're getting ready, remove existing intermediate files
        if self.master and FILES.use_mdins:
            self.remove(-1)
        elif self.master and not FILES.rewrite_output:
            self.remove(0)

        # Now load the parms and check them
        self.stdout.write('Loading and checking parameter files for '
                          'compatibility...\n')
        # Find external programs IFF we are doing a calc
        if not FILES.make_mdins:
            external_progs = {}
            if self.master:
                external_progs = find_progs(self.INPUT)
            external_progs = self.MPI.COMM_WORLD.bcast(external_progs, root=0)
            # Make external_progs an instance attribute
            self.external_progs = external_progs

        # Make amber topologies
        maketop = CheckMakeTop(FILES, INPUT, self.external_progs)
        (FILES.complex_prmtop, FILES.receptor_prmtop, FILES.ligand_prmtop, FILES.mutant_complex_prmtop,
         FILES.mutant_receptor_prmtop, FILES.mutant_ligand_prmtop) = maketop.makeToptleap()

        self.normal_system = MMPBSA_System(FILES.complex_prmtop, FILES.receptor_prmtop, FILES.ligand_prmtop)
        self.using_chamber = self.normal_system.complex_prmtop.chamber
        self.mutant_system = None
        if INPUT['alarun']:
            if (FILES.mutant_receptor_prmtop is None and FILES.mutant_ligand_prmtop is None and not self.stability):
                raise MMPBSA_Error('Alanine scanning requires either a mutated receptor or mutated ligand topology '
                                   'file!')
            if FILES.mutant_receptor_prmtop is None:
                FILES.mutant_receptor_prmtop = FILES.receptor_prmtop
            elif FILES.mutant_ligand_prmtop is None:
                FILES.mutant_ligand_prmtop = FILES.ligand_prmtop
            self.mutant_system = MMPBSA_System(FILES.mutant_complex_prmtop, FILES.mutant_receptor_prmtop,
                                               FILES.mutant_ligand_prmtop)
            if self.using_chamber is not self.mutant_system.complex_prmtop.chamber:
                raise MMPBSA_Error('CHAMBER prmtops must be used for both mutant '
                                   'and normal prmtops or neither!')
        # If we have a chamber prmtop, force using sander
        if self.using_chamber:
            INPUT['use_sander'] = True
            if INPUT['rismrun']:
                raise MMPBSA_Error('CHAMBER prmtops cannot be used with 3D-RISM')
            if INPUT['nmoderun']:
                raise MMPBSA_Error('CHAMBER prmtops cannot be used with NMODE')
            self.stdout.write('CHAMBER prmtops found. Forcing use of sander\n')

        # Print warnings if we are overwriting any masks and get default masks
        if (INPUT['ligand_mask'] is None and INPUT['receptor_mask'] is not None):
            warnings.warn('receptor_mask overwritten with default\n')
            INPUT['receptor_mask'] = None
        if (INPUT['receptor_mask'] is None and INPUT['ligand_mask'] is not None):
            warnings.warn('ligand_mask overwritten with default\n')
            INPUT['ligand_mask'] = None
        # Map the gmx_MMPBSA systems with the input masks, and get the default ones if
        # the masks were not input
        self.normal_system.Map(INPUT['receptor_mask'], INPUT['ligand_mask'])
        self.normal_system.CheckConsistency()
        if INPUT['alarun']:
            self.mutant_system.Map(INPUT['receptor_mask'], INPUT['ligand_mask'])
            self.mutant_system.CheckConsistency()
        if (INPUT['ligand_mask'] is None or INPUT['receptor_mask'] is None):
            com_mask, INPUT['receptor_mask'], INPUT['ligand_mask'] = \
                self.normal_system.Mask('all', in_complex=True)
        self.sync_mpi()
        self.timer.stop_timer('setup')

    def write_final_outputs(self):
        """ Writes the final output files for gmx_MMPBSA """
        self.timer.add_timer('output', 'Statistics calculation & output writing:')
        self.timer.start_timer('output')
        if (not hasattr(self, 'input_file_text') or not hasattr(self, 'FILES') or
                not hasattr(self, 'INPUT') or not hasattr(self, 'normal_system')):
            raise InternalError('I am not prepared to write the final output file!')
        # Only the master does this, so bail out if we are not master
        if not self.master:
            return
        # If we haven't already parsed our output files, do that now
        if not hasattr(self, 'calc_types'):
            self.parse_output_files()
        # Do the output files now
        if self.stability:
            write_stability_output(self)
        else:
            write_binding_output(self)
        if self.INPUT['decomprun']:
            if self.stability:
                write_decomp_stability_output(self.FILES, self.INPUT, self.mpi_size,
                                              self.normal_system, self.mutant_system, self.mut_str, self.pre)
            else:
                write_decomp_binding_output(self.FILES, self.INPUT, self.mpi_size,
                                            self.normal_system, self.mutant_system, self.mut_str, self.pre)
        self.timer.stop_timer('output')

    def finalize(self):
        """ We are done. Finish up timers and print out timing info """
        self.timer.done()
        if not self.master:
            self.MPI.Finalize()
            sys.exit(0)
        self.stdout.write('\nTiming:\n')
        self.timer.print_('setup', self.stdout)

        if not self.FILES.rewrite_output:
            self.timer.print_('cpptraj', self.stdout)

            if self.INPUT['alarun']:
                self.timer.print_('muttraj', self.stdout)

            self.timer.print_('calc', self.stdout)
            self.stdout.write('\n')

            if self.INPUT['gbrun']:
                self.timer.print_('gb', self.stdout)

            if self.INPUT['pbrun']:
                self.timer.print_('pb', self.stdout)

            if self.INPUT['nmoderun']:
                self.timer.print_('nmode', self.stdout)

            if self.INPUT['entropy'] == 1:
                self.timer.print_('qh', self.stdout)

            self.stdout.write('\n')

        self.timer.print_('output', self.stdout)
        self.timer.print_('global', self.stdout)

        self.remove(self.INPUT['keep_files'])

        self.stdout.write('\n\ngmx_MMPBSA Finished! Thank you for using. Please '
                          'cite us if you publish this work with this paper:\n   '
                          'Coming soon\n   '
                          ' and \n'
                          'Miller III, B. R., McGee Jr., T. D., Swails, J. M. '
                          'Homeyer, N. Gohlke, H. and Roitberg, A. E.\n   '
                          'J. Chem. Theory Comput., 2012, 8 (9) pp 3314-3321\n')
        self.MPI.Finalize()

        if self.FILES.gui and not self.FILES.stability:
            self.stdout.write('Opening GUI to analyze results...')
            GUI_run(self.FILES.prefix + 'info')
        else:
            sys.exit(0)

    def get_cl_args(self, args=None):
        """
        Gets the command-line arguments to load the INPUT array. Also determines
        if we are doing a stability calculation or not
        """
        if args is None:
            args = sys.argv
        if self.master:
            self.FILES = self.clparser.parse_args(args)
        else:
            self.FILES = object()
        # Broadcast the FILES
        self.FILES = self.MPI.COMM_WORLD.bcast(self.FILES)
        # Hand over the file prefix to the App instance
        self.pre = self.FILES.prefix
        if self.FILES.receptor_trajs or self.FILES.ligand_trajs:
            self.traj_protocol = 'MTP'  # multiple traj protocol
        else:
            self.traj_protocol = 'STP'  # single traj protocol
        # change by explicity argument
        self.stability = self.FILES.stability

    def read_input_file(self, infile=None):
        """ Reads the input file, pull it from FILES if not provided here """
        global _debug_printlevel
        if infile is None:
            if not hasattr(self, 'FILES'):
                raise InternalError('FILES not present and no input file given!')
            infile = self.FILES.input_file
        self.INPUT = self.input_file.Parse(infile)
        _debug_printlevel = self.INPUT['debug_printlevel']
        self.input_file_text = str(self.input_file)

    def process_input(self):
        """
        This handles processing of the INPUT dict if necessary if this is a 'new'
        calculation (i.e., not re-writing output). This does the following prep:
           - invert scale
           - determine trajectory file suffix
           - set decomp-dependent GBSA default
           - adjust verbose for stability calcs
           - 3D-RISM setup
           - Set temperature. Don't put it in namelist, because temp change
             for entropy requires changes to nmode and cpptraj calcs, meaning it
             is not as easily changed here.
        """
        # Invert scale
        self.INPUT['scale'] = 1 / self.INPUT['scale']

        # Set up netcdf variables and decide trajectory suffix
        if self.INPUT['netcdf'] == 0:
            self.INPUT['netcdf'] = ''
            self.trjsuffix = 'mdcrd'
        else:
            self.INPUT['netcdf'] = 'netcdf'
            self.trjsuffix = 'nc'

        # Set default GBSA for Decomp
        if self.INPUT['decomprun']:
            self.INPUT['gbsa'] = 2

        # Force to use Sander when intdiel is defined
        if self.INPUT['intdiel'] > 1.0:
            self.INPUT['use_sander'] = 1

        # Stability: no terms cancel, so print them all
        if self.stability:
            self.INPUT['verbose'] = 2

        # 3D-RISM stuff (keywords are case-insensitive)
        self.INPUT['thermo'] = self.INPUT['thermo'].lower()
        if self.INPUT['solvcut'] is None:
            self.INPUT['solvcut'] = self.INPUT['buffer']
        self.INPUT['rismrun_std'] = (self.INPUT['rismrun'] and
                                     self.INPUT['thermo'] in ['std', 'both'])
        self.INPUT['rismrun_gf'] = (self.INPUT['rismrun'] and
                                    self.INPUT['thermo'] in ['gf', 'both'])

        # Default temperature
        self.INPUT['temp'] = 298.15

    def check_for_bad_input(self, INPUT=None):
        """ Checks for bad user input """
        if INPUT is None:
            INPUT = self.INPUT

        if not INPUT['igb'] in [1, 2, 5, 7, 8]:
            raise InputError('Invalid value for IGB (%s)! ' % INPUT['igb'] +
                             'It must be 1, 2, 5, 7, or 8.')
        if INPUT['saltcon'] < 0:
            raise InputError('SALTCON must be non-negative!')
        if INPUT['surften'] < 0:
            raise InputError('SURFTEN must be non-negative!')
        if INPUT['indi'] < 0:
            raise InputError('INDI must be non-negative!')
        if INPUT['exdi'] < 0:
            raise InputError('EXDI must be non-negative!')
        if INPUT['scale'] < 0:
            raise InputError('SCALE must be non-negative!')
        if INPUT['linit'] < 0:
            raise InputError('LINIT must be a positive integer!')
        if not INPUT['prbrad'] in [1.4, 1.6]:
            raise InputError('PRBRAD (%s) must be 1.4 and 1.6!' % INPUT['prbrad'])
        if INPUT['istrng'] < 0:
            raise InputError('ISTRNG must be non-negative!')
        if not INPUT['inp'] in [0, 1, 2]:
            raise InputError('INP/NPOPT (%s) must be 0, 1, or 2!' % INPUT['inp'])
        if INPUT['cavity_surften'] < 0:
            raise InputError('CAVITY_SURFTEN must be non-negative!')
        if INPUT['fillratio'] <= 0:
            raise InputError('FILL_RATIO must be positive!')
        if not INPUT['radiopt'] in [0, 1]:
            raise InputError('RADIOPT (%s) must be 0 or 1!' % INPUT['radiopt'])
        if INPUT['dielc'] <= 0:
            raise InputError('DIELC must be positive!')
        if INPUT['maxcyc'] < 1:
            raise InputError('MAXCYC must be a positive integer!')
        if not INPUT['idecomp'] in [0, 1, 2, 3, 4]:
            raise InputError('IDECOMP (%s) must be 1, 2, 3, or 4!' %
                             INPUT['idecomp'])
        if INPUT['idecomp'] != 0 and INPUT['sander_apbs'] == 1:
            raise InputError('IDECOMP cannot be used with sander.APBS!')
        if not INPUT['entropy'] in [0, 1, 2]:
            raise InputError('ENTROPY (%s) must be 0, 1 or 2!' % INPUT['entropy'])
        if INPUT['entropy_seg'] not in range(1, 101):
            raise InputError('Entropy Segment (%s) must be in 1-100!' % INPUT['entropy_seg'])
        if not INPUT['sander_apbs'] in [0, 1]:
            raise InputError('SANDER_APBS must be 0 or 1!')
        if INPUT['alarun'] and INPUT['netcdf'] != '':
            raise InputError('Alanine scanning is incompatible with NETCDF != 0!')
        if INPUT['decomprun'] and INPUT['idecomp'] == 0:
            raise InputError('IDECOMP cannot be 0 for Decomposition analysis!')
        if INPUT['ions_parameters'] not in range(1,13):
            raise InputError('Ions parameters file name must be in %s!' % range(1,13))
        if INPUT['PBRadii'] not in [1, 2, 3, 4]:
            raise InputError('PBRadii must be 1, 2, 3 or 4!')
        if INPUT['solvated_trajectory'] not in [0, 1]:
            raise InputError('Ligand force field must be 0 or 1!')
        if not INPUT['use_sander'] in [0, 1]:
            raise InputError('USE_SANDER must be set to 0 or 1!')
        if not INPUT['ifqnt'] in [0, 1]:
            raise InputError('QMMM must be 0 or 1!')
        if INPUT['ifqnt'] == 1:
            if not INPUT['qm_theory'] in ['PM3', 'AM1', 'MNDO', 'PDDG-PM3', 'PM3PDDG',
                                          'PDDG-MNDO', 'PDDGMNDO', 'PM3-CARB1',
                                          'PM3CARB1', 'DFTB', 'SCC-DFTB', 'RM1', 'PM6',
                                          'PM3-ZnB', 'PM3-MAIS', 'PM6-D', 'PM6-DH+',
                                          'AM1-DH+', 'AM1-D*', 'PM3ZNB', 'MNDO/D',
                                          'MNDOD']:
                raise InputError('Invalid QM_THEORY (%s)! ' % INPUT['qm_theory'] +
                                 'This variable must be set to allowable options.\n' +
                                 '       See the Amber manual for allowable options.')
            if INPUT['qm_residues'] == '':
                raise InputError('QM_RESIDUES must be specified for IFQNT = 1!')
            if INPUT['decomprun']:
                raise InputError('QM/MM and decomposition are incompatible!')
            if (INPUT['qmcharge_lig'] + INPUT['qmcharge_rec'] !=
                    INPUT['qmcharge_com'] and not self.stability):
                raise InputError('The total charge of the ligand and receptor ' +
                                 'does not equal the charge of the complex!')
        if INPUT['rismrun']:
            if INPUT['rism_verbose'] > 2 or INPUT['rism_verbose'] < 0:
                raise InputError('RISM_VERBOSE must be 0, 1, or 2!')
            if INPUT['buffer'] < 0 and INPUT['solvcut'] < 0:
                raise InputError('If BUFFER < 0, SOLVCUT must be > 0!')
            if INPUT['tolerance'] < 0:
                raise InputError('TOLERANCE must be positive!')
            if INPUT['buffer'] < 0 and INPUT['ng'] == '':
                raise InputError('You must specify NG if BUFFER < 0!')
            if INPUT['closure'] == 'pse' and INPUT['closureorder'] < 1:
                raise InputError('You must specify CLOSUREORDER if CLOSURE=pse!')
            if not INPUT['polardecomp'] in [0, 1]:
                raise InputError('POLARDECOMP must be either 0 or 1!')
            if not INPUT['thermo'] in ['std', 'gf', 'both']:
                raise InputError('THERMO must be "std", "gf", or "both"!')
        if not (INPUT['gbrun'] or INPUT['pbrun'] or INPUT['rismrun'] or
                INPUT['nmoderun'] or INPUT['entropy']):
            raise InputError('You did not specify any type of calculation!')

        if INPUT['decomprun'] and not (INPUT['gbrun'] or INPUT['pbrun']):
            raise InputError('DECOMP must be run with either GB or PB!')

        if not INPUT['molsurf'] and (INPUT['msoffset'] != 0 or
                                     INPUT['probe'] != 1.4):
            warnings.warn('offset and probe are molsurf-only options',
                          InputWarning)

        # User warning when intdiel > 10
        if self.INPUT['intdiel'] > 10:
            warnings.warn('Intdiel should be less than 10, but it is {}'.format(self.INPUT['intdiel']), InputWarning)
        # check mutant definition
        if not self.INPUT['mutant'].lower() in ['rec', 'receptor', 'lig', 'ligand']:
            raise InputError('The mutant most be receptor (or rec) or ligand (or lig)')

    def remove(self, flag):
        """ Removes temporary files """
        if not self.master:
            return
        utils.remove(flag, mpi_size=self.mpi_size, fnpre=self.pre)

    def sync_mpi(self):
        """ Throws up a barrier """
        self.MPI.COMM_WORLD.Barrier()

    def calculate_interaction_entropy(self, key, mutant=False):
        """
        Calculate the interaction entropy described FIXME: article
        :param key:
        :return:
        """
        # gases constant in kcal/mol
        k = 0.001987
        if mutant:
            calc_types = self.calc_types['mutant']
        else:
            calc_types = self.calc_types
        energy_int = np.array([], dtype=np.float)
        a_energy_int = np.array([], dtype=np.float)
        d_energy_int = np.array([], dtype=np.float)
        exp_energy_int = np.array([], dtype=np.float)
        ts = np.array([], dtype=np.float)

        ggas = calc_types[key]['delta'].data['DELTA G gas']

        for eint in ggas:
            energy_int = np.append(energy_int, eint)
            aeint = energy_int.mean()
            a_energy_int = np.append(a_energy_int, aeint)
            deint = eint - aeint
            d_energy_int = np.append(d_energy_int, deint)
            eceint = exp(deint / (k * self.INPUT['entropy_temp']))
            exp_energy_int = np.append(exp_energy_int, eceint)
            aeceint = exp_energy_int.mean()
            cts = k * self.INPUT['entropy_temp'] * log(aeceint)
            ts = np.append(ts, cts)
            calc_types[key]['delta'].data['-TDS'] = ts

    def parse_output_files(self):
        """
        This parses the output files and loads them into dicts for easy access
        """
        # Only the master does this
        if not self.master:
            return
        self.calc_types = {}
        INPUT, FILES = self.INPUT, self.FILES
        # Mutant will also be a dict
        if INPUT['alarun']:
            self.calc_types['mutant'] = {}
        # Quasi-harmonic analysis is a special-case, so handle that separately
        if INPUT['entropy'] == 1:
            if not INPUT['mutant_only']:
                self.calc_types['qh'] = QHout(self.pre + 'cpptraj_entropy.out',
                                              INPUT['temp'])
            if INPUT['alarun']:
                self.calc_types['mutant']['qh'] = QHout(self.pre +
                                                        'mutant_cpptraj_entropy.out', INPUT['temp'])
        # Set BindingClass based on whether it's a single or multiple trajectory
        # analysis
        if self.traj_protocol == 'STP':
            BindClass = SingleTrajBinding
        else:
            BindClass = MultiTrajBinding
        # Determine if our GB is QM/MM or not
        if INPUT['ifqnt']:
            GBClass = QMMMout
        else:
            GBClass = GBout
        # Determine which kind of RISM output class we are based on std/gf and
        # polardecomp
        if INPUT['polardecomp']:
            RISM_GF = PolarRISM_gf_Out
            RISM_Std = PolarRISM_std_Out
        else:
            RISM_GF = RISM_gf_Out
            RISM_Std = RISM_std_Out
        # Now we make a list of the other calculation types, their INPUT triggers,
        # their key in the calc_types dict, the base name of their output files
        # without the prefix (with %s-substitution for complex, receptor, or
        # ligand), and the class for their output
        triggers = ('nmoderun', 'gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf')
        outclass = (NMODEout, GBClass, PBout, RISM_Std, RISM_GF)
        outkey = ('nmode', 'gb', 'pb', 'rism std', 'rism gf')
        basename = ('%s_nm.out', '%s_gb.mdout', '%s_pb.mdout', '%s_rism.mdout',
                    '%s_rism.mdout')
        for i, key in enumerate(outkey):
            if not INPUT[triggers[i]]:
                continue
            # Non-mutant
            if not INPUT['mutant_only']:
                self.calc_types[key] = {'complex': outclass[i](self.pre +
                                                               basename[i] % 'complex', self.INPUT, self.mpi_size,
                                                               self.using_chamber)}
                if not self.stability:
                    self.calc_types[key]['receptor'] = outclass[i](self.pre +
                                                                   basename[i] % 'receptor', self.INPUT, self.mpi_size,
                                                                   self.using_chamber)
                    self.calc_types[key]['ligand'] = outclass[i](self.pre +
                                                                 basename[i] % 'ligand', self.INPUT, self.mpi_size,
                                                                 self.using_chamber)
                    self.calc_types[key]['delta'] = BindClass(
                        self.calc_types[key]['complex'],
                        self.calc_types[key]['receptor'],
                        self.calc_types[key]['ligand'],
                        self.INPUT['verbose'], self.using_chamber)

                    if self.INPUT['entropy'] == 2:
                        self.calculate_interaction_entropy(key)
                else:
                    self.calc_types[key]['complex'].fill_composite_terms()
            # Time for mutant
            if INPUT['alarun']:
                self.calc_types['mutant'][key] = {'complex':
                                                      outclass[i](self.pre + 'mutant_' + basename[i] % 'complex',
                                                                  self.INPUT, self.mpi_size, self.using_chamber)}
                if not self.stability:
                    self.calc_types['mutant'][key]['receptor'] = outclass[i](
                        self.pre + 'mutant_' + basename[i] % 'receptor',
                        self.INPUT, self.mpi_size, self.using_chamber)
                    self.calc_types['mutant'][key]['ligand'] = outclass[i](
                        self.pre + 'mutant_' + basename[i] % 'ligand',
                        self.INPUT, self.mpi_size, self.using_chamber)
                    self.calc_types['mutant'][key]['delta'] = BindClass(
                        self.calc_types['mutant'][key]['complex'],
                        self.calc_types['mutant'][key]['receptor'],
                        self.calc_types['mutant'][key]['ligand'],
                        self.INPUT['verbose'], self.using_chamber)
                    if self.INPUT['entropy'] == 2:
                        self.calculate_interaction_entropy(key, mutant=True)
                else:
                    self.calc_types['mutant'][key]['complex'].fill_composite_terms()


# Local methods

def excepthook(exception_type, exception_value, tb):
    """
    Replaces sys.excepthook so fatal exceptions kill all MPI threads and we can
    control the printing of tracebacks. Those are helpful for debugging purposes,
    but may be unsightly to users. debug_printlevel set above controls this
    behavior
    """
    import traceback
    global _debug_printlevel, _stderr, _mpi_size, _rank
    if _debug_printlevel > 1 or not isinstance(exception_type, MMPBSA_Error):
        traceback.print_tb(tb)
    _stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
    if _mpi_size > 1:
        _stderr.write('Error occured on rank %d.' % _rank + os.linesep)
    _stderr.write('Exiting. All files have been retained.' + os.linesep)
    _MPI.COMM_WORLD.Abort(1)


def interrupt_handler(signal, frame):
    """ Handles interrupt signals for a clean exit """
    global _MPI, _stderr
    _stderr.write('\n%s interrupted! Program terminated. All files are kept.\n' %
                  os.path.split(sys.argv[0])[1])
    _MPI.COMM_WORLD.Abort(1)


def setup_run():
    """
    Replace the uncaught exception handler to control traceback printing. Also
    add a signal handler for a SIGINT (Ctrl-C). However, we only want to do this
    if we're running gmx_MMPBSA -- for the API, we don't want to clobber the
    users' python environments like this.
    """
    sys.excepthook = excepthook
    signal.signal(signal.SIGINT, interrupt_handler)
