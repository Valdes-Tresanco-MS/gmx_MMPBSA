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
import logging
# Import gmx_MMPBSA modules
from GMXMMPBSA import utils, __version__
from GMXMMPBSA.amber_outputs import (QHout, NMODEout, QMMMout, GBout, PBout, PolarRISM_std_Out, RISM_std_Out,
                                     PolarRISM_gf_Out, RISM_gf_Out, PolarRISM_pcplus_Out, RISM_pcplus_Out,
                                     BindingStatistics, IEout, C2out, DeltaDeltaStatistics, DeltaIEC2Statistic)
from GMXMMPBSA.calculation import (CalculationList, EnergyCalculation, PBEnergyCalculation,
                                   NmodeCalc, QuasiHarmCalc, CopyCalc, PrintCalc, LcpoCalc, MolsurfCalc,
                                   InteractionEntropyCalc, C2EntropyCalc)
from GMXMMPBSA.commandlineparser import parser
from GMXMMPBSA.createinput import create_inputs, SanderRISMInput
from GMXMMPBSA.exceptions import (MMPBSA_Error, InternalError, InputError, GMXMMPBSA_ERROR)
from GMXMMPBSA.infofile import InfoFile
from GMXMMPBSA.fake_mpi import MPI as FakeMPI
from GMXMMPBSA.input_parser import input_file as _input_file
from GMXMMPBSA.make_trajs import make_trajectories, make_mutant_trajectories
from GMXMMPBSA.output_file import (write_outputs, write_decomp_output, data2pkl)
from GMXMMPBSA.parm_setup import MMPBSA_System
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.timer import Timer

# Global variables for the excepthook replacement at the bottom. Override these
# in the MMPBSA_App constructor and input file reading
_unbuf_stdout = utils.Unbuffered(sys.stdout)  # unbuffered stdout
_unbuf_stderr = utils.Unbuffered(sys.stderr)  # unbuffered stderr
_stdout = sys.stdout
_stderr = sys.stderr
_mpi_size = 1
_rank = 0
_MPI = FakeMPI()


# Main class

class MMPBSA_App(object):
    """ Main MM/PBSA application for driving the entire calculation """
    # The command line parser and input file objects are class attributes here
    clparser = parser
    input_file = _input_file

    def __init__(self, MPI, stdout=None, stderr=None, size=None):
        """
        Sets up the main gmx_MMPBSA driver class. All we set up here is the output
        and error streams (unbuffered by default) and the prefix for the
        intermediate files. Also set up empty INPUT dict
        """
        global _rank, _stdout, _stderr, _mpi_size, _MPI
        _MPI = self.MPI = MPI
        self.pre = '_GMXMMPBSA_'
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
        if self.master:
            logging.info(f'Starting gmx_MMPBSA {__version__}')
            utils.get_sys_info()

        # Set up timers
        timers = [Timer() for _ in range(self.mpi_size)]
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
            GMXMMPBSA_ERROR('MMPBSA_App not set up and parms not checked!', InternalError)
        # Set up some local refs for convenience
        FILES, INPUT, master = self.FILES, self.INPUT, self.master

        # # Now we're getting ready, remove existing intermediate files
        # elif master and not FILES.rewrite_output:
        #     self.remove(0)

        # Create input files based on INPUT dict
        if master:
            create_inputs(INPUT, self.normal_system, self.pre)
        self.timer.stop_timer('setup')

        # Now create our trajectory files

        self.timer.add_timer('cpptraj', 'Creating trajectories with cpptraj:')
        self.timer.start_timer('cpptraj')

        if master:
            logging.info('Preparing trajectories for simulation...\n')
            (self.numframes, rec_frames,
             lig_frames, self.numframes_nmode) = make_trajectories(INPUT, FILES, self.mpi_size,
                                                                   self.external_progs['cpptraj'],
                                                                   self.pre)
            if self.traj_protocol == 'MTP' and not self.numframes == rec_frames == lig_frames:
                GMXMMPBSA_ERROR('The complex, receptor, and ligand trajectories must be the same length. Since v1.5.0 '
                                'we have simplified a few things to make the code easier to maintain. Please check the '
                                'documentation')

        self.MPI.COMM_WORLD.Barrier()

        self.timer.stop_timer('cpptraj')

        self.timer.add_timer('muttraj', 'Mutating trajectories:')
        self.timer.start_timer('muttraj')

        if INPUT['alarun'] and self.master:
            logging.info('Mutating trajectories...')
        _, mutant_residue = make_mutant_trajectories(INPUT, FILES, self.mpi_rank, self.external_progs['cpptraj'],
                                                     self.normal_system, self.mutant_system, self.pre)

        self.MPI.COMM_WORLD.Barrier()

        if master:
            logging.info('%d frames were processed by cpptraj for use in calculation.' % self.numframes)
            if INPUT['nmoderun']:
                logging.info('%d frames were processed by cpptraj for nmode calculations.' % self.numframes_nmode)

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
        if INPUT['qh_entropy']:
            self.timer.add_timer('qh', 'Total quasi-harmonic calculation time:')

        self.sync_mpi()

    def run_mmpbsa(self, rank=None):
        """
        Runs the MM/PBSA analysis. This assumes FILES and INPUT are already set.
        """

        if not hasattr(self, 'external_progs'):
            GMXMMPBSA_ERROR('external_progs not declared in run_mmpbsa!', InternalError)

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

        if master:
            self.timer.stop_timer('calc')
            # Write out the info file now
            info = InfoFile(self)
            info.write_info(f'{self.pre}info')

    def load_calc_list(self):
        """
        Sets up all of the calculations to be run. When adding a new
        calculation type, add a class to calculation.py, import it at the top of
        the file here, then append it to the calc list appropriately
        """
        self.calc_list = CalculationList(self.timer)

        if not self.INPUT['mutant_only']:
            self.calc_list.append(PrintCalc('Running calculations on normal system...'), timer_key=None)
            self._load_calc_list(self.pre, False, self.normal_system)
        if self.INPUT['alarun']:
            self.calc_list.append(PrintCalc('Running calculations on mutant system...'), timer_key=None)
            self._load_calc_list(f'{self.pre}mutant_', True, self.mutant_system)

    def _load_calc_list(self, prefix, mutant, parm_system):
        """
        Internal routine to handle building calculation list. Called separately
        for mutant and normal systems
        """
        # Set up a dictionary of external programs to use based one external progs
        progs = {'gb': self.external_progs['sander'],
                 'sa': self.external_progs['cpptraj'],
                 'pb': self.external_progs['sander'],
                 'rism': self.external_progs['sander'],
                 'qh': self.external_progs['cpptraj'],
                 'nmode': self.external_progs['mmpbsa_py_nabnmode']
                 }
        if self.INPUT['sander_apbs']:
            progs['pb'] = self.external_progs['sander.APBS']

        # NetCDF or ASCII intermediate trajectories?
        trj_sfx = 'nc' if self.INPUT['netcdf'] else 'mdcrd'

        # Determine if we just copy the receptor files. This only happens if we
        # are doing mutant calculations, we're not only doing the mutant, and the
        # receptor/mutant receptor topologies are equal. Same for the ligand
        copy_receptor = (mutant and not self.INPUT['mutant_only'] and
                         self.FILES.receptor_prmtop == self.FILES.mutant_receptor_prmtop)
        copy_ligand = (mutant and not self.INPUT['mutant_only'] and
                       self.FILES.ligand_prmtop == self.FILES.mutant_ligand_prmtop)

        # First load the GB calculations
        if self.INPUT['gbrun']:
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
            elif self.INPUT['alpb']:
                mdin_template = self.pre + 'gb_%s.mdin'
            else:
                mdin_template = self.pre + 'gb.mdin'

            # Now do complex-specific stuff
            try:
                mdin = mdin_template % 'com'
            except TypeError:
                mdin = mdin_template

            self.calc_list.append(PrintCalc(f"Beginning GB calculations with {progs['gb']}"), timer_key='gb')

            c = EnergyCalculation(progs['gb'], parm_system.complex_prmtop,
                                  incrd % 'complex',
                                  '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                  mdin, '%scomplex_gb.mdout.%%d' % (prefix),
                                  self.pre + 'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='gb')
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

            self.calc_list.append(PrintCalc(f"Beginning PB calculations with {progs['pb']}"), timer_key='pb')

            c = PBEnergyCalculation(progs['pb'], parm_system.complex_prmtop,
                                    incrd % 'complex',
                                    '%scomplex.%s.%%d' % (prefix, trj_sfx),
                                    mdin, '%scomplex_pb.mdout.%%d' % prefix,
                                    self.pre + 'restrt.%d')
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='pb')
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
                                            mdin, '%sreceptor_pb.mdout.%%d' % prefix,
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
            mdin = self.pre + 'rism.mdin'
            self.calc_list.append(
                PrintCalc('Beginning 3D-RISM calculations with %s' % progs['rism']), timer_key='rism')

            c = EnergyCalculation(progs['rism'], parm_system.complex_prmtop,
                                  '%sdummycomplex.inpcrd' % prefix,
                                  '%scomplex.%s.%%d' % (prefix, trj_sfx), mdin,
                                  '%scomplex_rism.mdout.%%d' % prefix,
                                  self.pre + 'restrt.%d', self.FILES.xvvfile)
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='rism')

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('%sreceptor_rism.mdout.%%d' % self.pre,
                                 '%sreceptor_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = EnergyCalculation(progs['rism'], parm_system.receptor_prmtop,
                                          '%sdummyreceptor.inpcrd' % prefix,
                                          '%sreceptor.%s.%%d' % (prefix, trj_sfx), mdin,
                                          '%sreceptor_rism.mdout.%%d' % prefix,
                                          self.pre + 'restrt.%d', self.FILES.xvvfile)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='rism')

                if copy_ligand:
                    c = CopyCalc('%sligand_rism.mdout.%%d' % self.pre,
                                 '%sligand_rism.mdout.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='pb')
                else:
                    c = EnergyCalculation(progs['rism'], parm_system.ligand_prmtop,
                                          '%sdummyligand.inpcrd' % prefix,
                                          '%sligand.%s.%%d' % (prefix, trj_sfx), mdin,
                                          '%sligand_rism.mdout.%%d' % prefix,
                                          self.pre + 'restrt.%d', self.FILES.xvvfile)
                    self.calc_list.append(c, '  calculating ligand contribution...', timer_key='rism')

        # end if self.INPUT['rismrun']

        if self.INPUT['nmoderun']:
            self.calc_list.append(
                PrintCalc('Beginning nmode calculations with %s' % progs['nmode']), timer_key='nmode')

            c = NmodeCalc(progs['nmode'], parm_system.complex_prmtop,
                          '%scomplex.pdb' % prefix,
                          '%scomplex_nm.%s.%%d' % (prefix, trj_sfx),
                          '%scomplex_nm.out.%%d' % prefix, self.INPUT)
            self.calc_list.append(c, '  calculating complex contribution...', timer_key='nmode')

            if not self.stability:
                if copy_receptor:
                    c = CopyCalc('%sreceptor_nm.out.%%d' % self.pre,
                                 '%sreceptor_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in receptor; '
                                             'using unmutated files', timer_key='nmode')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.receptor_prmtop,
                                  '%sreceptor.pdb' % prefix,
                                  '%sreceptor_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sreceptor_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating receptor contribution...',
                                          timer_key='nmode')

                if copy_ligand:
                    c = CopyCalc('%sligand_nm.out.%%d' % self.pre,
                                 '%sligand_nm.out.%%d' % prefix)
                    self.calc_list.append(c, '  no mutation found in ligand; '
                                             'using unmutated files', timer_key='nmode')
                else:
                    c = NmodeCalc(progs['nmode'], parm_system.ligand_prmtop,
                                  '%sligand.pdb' % prefix,
                                  '%sligand_nm.%s.%%d' % (prefix, trj_sfx),
                                  '%sligand_nm.out.%%d' % prefix, self.INPUT)
                    self.calc_list.append(c, '  calculating ligand contribution...',
                                          timer_key='nmode')

        # end if self.INPUT['nmoderun']

        # Only master does entropy calculations
        if self.INPUT['qh_entropy']:
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

    def make_prmtops(self):
        self.timer.add_timer('setup_gmx', 'Total GROMACS setup time:')
        self.timer.start_timer('setup_gmx')
        if not self.FILES.rewrite_output and self.master:
            self.remove(-1)

        # Find external programs IFF we are doing a calc
        external_progs = utils.find_progs(self.INPUT, self.mpi_size) if self.master else {}
        external_progs = self.MPI.COMM_WORLD.bcast(external_progs, root=0)
        # Make external_progs an instance attribute
        self.external_progs = external_progs
        if self.master:
            # Make amber topologies
            logging.info('Building AMBER topologies from GROMACS files...')
            maketop = CheckMakeTop(self.FILES, self.INPUT, self.external_progs)
            (self.FILES.complex_prmtop, self.FILES.receptor_prmtop, self.FILES.ligand_prmtop,
             self.FILES.mutant_complex_prmtop,
             self.FILES.mutant_receptor_prmtop, self.FILES.mutant_ligand_prmtop) = maketop.buildTopology()
            logging.info('Building AMBER topologies from GROMACS files... Done.\n')
            self.INPUT['receptor_mask'], self.INPUT['ligand_mask'], self.resl = maketop.get_masks()
            self.mutant_index = maketop.com_mut_index
            self.mut_str = self.resl[maketop.com_mut_index].mutant_label if self.mutant_index else ''
            self.FILES.complex_fixed = f'{self.FILES.prefix}COM_FIXED.pdb'
        self.FILES = self.MPI.COMM_WORLD.bcast(self.FILES, root=0)
        self.INPUT = self.MPI.COMM_WORLD.bcast(self.INPUT, root=0)
        self.sync_mpi()
        self.timer.stop_timer('setup_gmx')

    def loadcheck_prmtops(self):
        """ Loads the topology files and checks their consistency """
        # Start setup timer and make sure we've already set up our input
        self.timer.add_timer('setup', 'Total AMBER setup time:')
        self.timer.start_timer('setup')
        if not hasattr(self, 'FILES') or not hasattr(self, 'INPUT'):
            GMXMMPBSA_ERROR('MMPBSA_App not set up! Cannot check parms yet!', InternalError)
        # create local aliases to avoid abundant selfs
        FILES, INPUT = self.FILES, self.INPUT
        if self.master:
            # Now load the parms and check them
            logging.info('Loading and checking parameter files for compatibility...')
        self.normal_system = MMPBSA_System(FILES.complex_prmtop, FILES.receptor_prmtop, FILES.ligand_prmtop)
        self.using_chamber = self.normal_system.complex_prmtop.chamber
        self.mutant_system = None
        if INPUT['alarun']:
            if (FILES.mutant_receptor_prmtop is None and FILES.mutant_ligand_prmtop is None and not self.stability):
                GMXMMPBSA_ERROR('Alanine scanning requires either a mutated receptor or mutated ligand topology '
                                'file!')
            if FILES.mutant_receptor_prmtop is None:
                FILES.mutant_receptor_prmtop = FILES.receptor_prmtop
            elif FILES.mutant_ligand_prmtop is None:
                FILES.mutant_ligand_prmtop = FILES.ligand_prmtop
            self.mutant_system = MMPBSA_System(FILES.mutant_complex_prmtop, FILES.mutant_receptor_prmtop,
                                               FILES.mutant_ligand_prmtop)
        # If we have a chamber prmtop, force using sander
        if self.using_chamber:
            if INPUT['rismrun']:
                GMXMMPBSA_ERROR('CHAMBER prmtops cannot be used with 3D-RISM')
            if INPUT['nmoderun']:
                GMXMMPBSA_ERROR('CHAMBER prmtops cannot be used with NMODE')

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
            GMXMMPBSA_ERROR('I am not prepared to write the final output file!', InternalError)
        # Only the master does this, so bail out if we are not master
        if not self.master:
            return
        # If we haven't already parsed our output files, do that now
        # FIXME: does this make sense?
        if not hasattr(self, 'calc_types'):
            self.parse_output_files()
        # Do the output files now
        write_outputs(self)
        if self.INPUT['decomprun']:
            write_decomp_output(self)
        if self.INPUT['keep_files'] in [0, 2]:
            data2pkl(self)

        info = InfoFile(self)
        info.write_info(f'{self.pre}info')

        self.timer.stop_timer('output')

    def finalize(self):
        """ We are done. Finish up timers and print out timing info """
        self.timer.done()
        if not self.master:
            self.MPI.Finalize()
            sys.exit(0)
        logging.info('Timing:')
        if not self.FILES.rewrite_output:
            self.timer.print_('setup_gmx')
        self.timer.print_('setup')

        if not self.FILES.rewrite_output:
            self._finalize_timers()
        self.timer.print_('output')
        self.timer.print_('global', True)

        self.remove(self.INPUT['keep_files'])

        exe_info = utils.get_warnings()
        logging.info(f"\n   Finalizing gmx_MMPBSA: [ERROR  ] = {exe_info['error']}; [WARNING] = {exe_info['warning']}\n"
                     f"   Check the gmx_MMPBSA.log file for more details...\n")


        logging.info(
            '\n Thank you for using gmx_MMPBSA. Please consider supporting gmx_MMPBSA by citing our publication:'
            '\n    Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. '
            '\n    gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. '
            '\n    J Chem Theory Comput., 2021, 17 (10):6281-6291. Epub 2021 Sep 29. PMID: 34586825.'
            '\n    https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645'
            '\n\nAlso consider citing MMPBSA.py:'
            '\n    Miller III, B. R., McGee Jr., T. D., Swails, J. M. Homeyer, N. Gohlke, H. and Roitberg, A. E.'
            '\n    MMPBSA.py: An Efficient Program for End-State Free Energy Calculations.'
            '\n    J. Chem. Theory Comput., 2012, 8 (9) pp 3314-3321\n')
        self.MPI.Finalize()

        end = 0
        if self.FILES.gui and not self.FILES.rewrite_output:
            import subprocess
            from pathlib import Path
            logging.info('Opening gmx_MMPBSA_ana to analyze results...\n')
            ifile = Path(f'{self.FILES.prefix}info')
            if not ifile.exists():
                ifile = Path('COMPACT_MMXSA_RESULTS.mmxsa')

            g = subprocess.Popen(['gmx_MMPBSA_ana', '-f', ifile.as_posix()])
            if g.wait():
                end = 1
        if end:
            logging.error('Unable to start gmx_MMPBSA_ana...')
        logging.info('Finalized...')
        sys.exit(end)

    def _finalize_timers(self):
        self.timer.print_('cpptraj')

        if self.INPUT['alarun']:
            self.timer.print_('muttraj', True)

        # self.stdout.write('\n')
        self.timer.print_('calc')

        if self.INPUT['gbrun']:
            self.timer.print_('gb')

        if self.INPUT['pbrun']:
            self.timer.print_('pb')

        if self.INPUT['nmoderun']:
            self.timer.print_('nmode')

        if self.INPUT['qh_entropy']:
            self.timer.print_('qh', True)

        # self.stdout.write('\n')

    def get_cl_args(self, args=None):
        """
        Gets the command-line arguments to load the INPUT array. Also determines
        if we are doing a stability calculation or not
        """
        if args is None:
            args = sys.argv
        if self.master:
            text_args = ' '.join(args)
            _mpi = True
            # remove mpi arg before passed it to app
            if 'mpi' in args:
                args.remove('mpi')
            elif 'MPI' in args:
                args.remove('MPI')
            else:
                _mpi = False
            mpi_cl = f'  mpirun -np {self.mpi_size} ' if _mpi else '  '
            # save args in gmx_MMPBSA.log
            logging.info('Command-line\n' + mpi_cl +
                         'gmx_MMPBSA ' + text_args + '\n')
            # check if any arg is duplicated
            utils._get_dup_args(args)
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
        # change by explicit argument
        self.stability = self.FILES.stability

    def read_input_file(self, infile=None):
        """ Reads the input file, pull it from FILES if not provided here """
        if infile is None:
            if not hasattr(self, 'FILES'):
                GMXMMPBSA_ERROR('FILES not present and no input file given!', InternalError)
            infile = self.FILES.input_file
        self.INPUT = self.input_file.Parse(infile)
        self.input_file_text = str(self.input_file)
        if self.master:
            for line in self.input_file_text.split('\n'):
                logging.debug(line)

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

        # Stability: no terms cancel, so print them all
        # if self.stability:
        #     self.INPUT['verbose'] = 2

        # 3D-RISM stuff (keywords are case-insensitive)
        if self.INPUT['solvcut'] is None:
            self.INPUT['solvcut'] = self.INPUT['buffer']

        self.INPUT['rismrun_std'] = bool(self.INPUT['rismrun'])
        self.INPUT['rismrun_gf'] = self.INPUT['rismrun'] and self.INPUT['gfcorrection']
        self.INPUT['rismrun_pcplus'] = self.INPUT['rismrun'] and self.INPUT['pcpluscorrection']

        # Default temperature
        # self.INPUT['temp'] = 298.15

    def check_for_bad_input(self, INPUT=None):
        """ Checks for bad user input """
        if INPUT is None:
            INPUT = self.INPUT
        if not self.master:
            return
        # Check deprecated variables
        # check force fields

        logging.info(f'Checking {self.FILES.input_file} input file...')

        if self.FILES.ligand_mol2:
            if 'leaprc.gaff' in self.INPUT['forcefields'] or 'leaprc.gaff2' in self.INPUT['forcefields']:
                pass
            else:
                logging.error(
                    "When using -lm flag, leaprc.gaff or leaprc.gaff2 should be included in the forcefields "
                    "variable. Check this tutorial for "
                    "more details https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand/ST/")

        if INPUT['igb'] not in [1, 2, 5, 7, 8]:
            GMXMMPBSA_ERROR('Invalid value for IGB (%s)! ' % INPUT['igb'] + 'IGB must be 1, 2, 5, 7, or 8.', InputError)
        if INPUT['intdiel'] < 0:
            GMXMMPBSA_ERROR('INDI must be non-negative!', InputError)
        if INPUT['extdiel'] < 0:
            GMXMMPBSA_ERROR('EXDI must be non-negative!', InputError)
        if INPUT['saltcon'] < 0:
            GMXMMPBSA_ERROR('SALTCON must be non-negative!', InputError)
        if INPUT['surften'] < 0:
            GMXMMPBSA_ERROR('SURFTEN must be non-negative!', InputError)
        if INPUT['alpb'] == 1 and INPUT['igb'] == 8:
            GMXMMPBSA_ERROR('IGB=8 is incompatible with ALPB=1! IGB must be 1, 2, 5, or 7 if ALPB=1.', InputError)
        if INPUT['arad_method'] not in [1, 2, 3]:
            GMXMMPBSA_ERROR('ARAD_METHOD must be 1, 2, or 3!', InputError)
        if INPUT['indi'] < 0:
            GMXMMPBSA_ERROR('INDI must be non-negative!', InputError)
        if INPUT['exdi'] < 0:
            GMXMMPBSA_ERROR('EXDI must be non-negative!', InputError)
        if INPUT['memopt'] > 0 and (INPUT['emem'] < INPUT['indi'] or INPUT['emem'] > INPUT['exdi']):
            logging.warning(
                "Membrane dielectric constant (emem) should be between indi and exdi's or there may be errors."
            )
        if INPUT['scale'] < 0:
            GMXMMPBSA_ERROR('SCALE must be non-negative!', InputError)
        if INPUT['linit'] < 0:
            GMXMMPBSA_ERROR('LINIT must be a positive integer!', InputError)
        if INPUT['prbrad'] not in [1.4, 1.6]:
            GMXMMPBSA_ERROR('PRBRAD (%s) must be 1.4 and 1.6!' % INPUT['prbrad'], InputError)
        if INPUT['istrng'] < 0:
            GMXMMPBSA_ERROR('ISTRNG must be non-negative!', InputError)
        if INPUT['inp'] not in [1, 2]:
            GMXMMPBSA_ERROR('INP/NPOPT (%s) must be 1, or 2!' % INPUT['inp'], InputError)
        if INPUT['cavity_surften'] < 0:
            GMXMMPBSA_ERROR('CAVITY_SURFTEN must be non-negative!', InputError)
        if INPUT['fillratio'] <= 0:
            GMXMMPBSA_ERROR('FILL_RATIO must be positive!', InputError)
        if INPUT['radiopt'] not in [0, 1]:
            GMXMMPBSA_ERROR('RADIOPT (%s) must be 0 or 1!' % INPUT['radiopt'], InputError)
        if INPUT['dielc'] <= 0:
            GMXMMPBSA_ERROR('DIELC must be positive!', InputError)
        if INPUT['maxcyc'] < 1:
            GMXMMPBSA_ERROR('MAXCYC must be a positive integer!', InputError)
        if INPUT['idecomp'] not in [0, 1, 2, 3, 4]:
            GMXMMPBSA_ERROR('IDECOMP (%s) must be 1, 2, 3, or 4!' % INPUT['idecomp'], InputError)
        if INPUT['idecomp'] != 0 and INPUT['sander_apbs'] == 1:
            GMXMMPBSA_ERROR('IDECOMP cannot be used with sander.APBS!', InputError)
        if INPUT['sander_apbs'] not in [0, 1]:
            GMXMMPBSA_ERROR('SANDER_APBS must be 0 or 1!', InputError)
        if INPUT['alarun'] and INPUT['netcdf'] != '':
            GMXMMPBSA_ERROR('Alanine scanning is incompatible with NETCDF != 0!', InputError)
        if INPUT['decomprun'] and INPUT['idecomp'] == 0:
            GMXMMPBSA_ERROR('IDECOMP cannot be 0 for Decomposition analysis!', InputError)
        if INPUT['ions_parameters'] not in range(1, 17):
            GMXMMPBSA_ERROR('Ions parameters file name must be in %s!' % range(1, 17), InputError)
        if INPUT['PBRadii'] not in range(1, 8):
            GMXMMPBSA_ERROR('PBRadii must be 1, 2, 3, 4, 5, 6, or 7!', InputError)
        if INPUT['solvated_trajectory'] not in [0, 1]:
            GMXMMPBSA_ERROR('SOLVATED_TRAJECTORY must be 0 or 1!', InputError)
        if INPUT['ifqnt'] not in [0, 1]:
            GMXMMPBSA_ERROR('QMMM must be 0 or 1!', InputError)
        if INPUT['ifqnt'] == 0 and (INPUT['qm_theory'] or INPUT['qm_residues']):
            logging.warning('qm_theory/qm_residues variable has been defined, however the potential function is '
                            'strictly classical (ifqnt=0). Please, set ifqnt=1 if you want to use Use QM/MM')
        if INPUT['ifqnt'] == 1:
            if INPUT['qm_theory'] not in ['PM3', 'AM1', 'MNDO', 'PDDG-PM3', 'PM3PDDG',
                                          'PDDG-MNDO', 'PDDGMNDO', 'PM3-CARB1',
                                          'PM3CARB1', 'DFTB', 'SCC-DFTB', 'RM1', 'PM6',
                                          'PM3-ZnB', 'PM3-MAIS', 'PM6-D', 'PM6-DH+',
                                          'AM1-DH+', 'AM1-D*', 'PM3ZNB', 'MNDO/D',
                                          'MNDOD']:
                GMXMMPBSA_ERROR('Invalid QM_THEORY (%s)! ' % INPUT['qm_theory'] +
                                'This variable must be set to allowable options.\n' +
                                'PM3, AM1, MNDO, PDDG-PM3, PM3PDDG, PDDG-MNDO, PDDGMNDO, \n'
                                'PM3-CARB1, PM3CARB1, DFTB, SCC-DFTB, RM1, PM6, PM3-ZnB, \n'
                                'PM3-MAIS, PM6-D, PM6-DH+, AM1-DH+, AM1-D*, PM3ZNB, MNDO/D, MNDOD', InputError)
            if INPUT['qm_residues'] == '':
                GMXMMPBSA_ERROR('QM_RESIDUES must be specified for IFQNT = 1!', InputError)
            if INPUT['decomprun']:
                GMXMMPBSA_ERROR('QM/MM and decomposition are incompatible!', InputError)
            if (INPUT['qmcharge_lig'] + INPUT['qmcharge_rec'] !=
                    INPUT['qmcharge_com'] and not self.stability):
                GMXMMPBSA_ERROR('The total charge of the ligand and receptor ' +
                                'does not equal the charge of the complex!', InputError)
            if INPUT['scfconv'] < 1.0e-12:
                logging.warning('There is a risk of convergence problems when the requested convergence is less than '
                                '1.0e-12 kcal/mol')
            if INPUT['writepdb']:
                logging.info('Writing qmmm_region.pdb PDB file of the selected QM region...')
            if INPUT['verbosity'] not in [0, 1, 2, 3, 4, 5]:
                GMXMMPBSA_ERROR('VERBOSITY must be 0, 1, 2, 3, 4 or 5!', InputError)
            if INPUT['verbosity'] >= 2:
                logging.warning('VERBOSITY values of 2 or higher will produce a lot of output')

        if INPUT['rismrun']:
            if INPUT['rism_verbose'] not in [0, 1, 2]:
                GMXMMPBSA_ERROR('RISM_VERBOSE must be 0, 1, or 2!', InputError)
            if INPUT['buffer'] < 0 and INPUT['solvcut'] < 0:
                GMXMMPBSA_ERROR('If BUFFER < 0, SOLVCUT must be > 0!', InputError)
            for tol in INPUT['tolerance']:
                if tol <= 0:
                    GMXMMPBSA_ERROR('TOLERANCE must be positive!', InputError)
            if INPUT['tolerance'][-1] > 0.00001:
                logging.warning(f"Default TOLERANCE value is 0.00001! However {INPUT['tolerance'][-1]} is been used. "
                                f"Check documentation for more details...")
            if INPUT['buffer'] < 0 and INPUT['ng'] == '':
                GMXMMPBSA_ERROR('You must specify NG if BUFFER < 0!', InputError)
            if INPUT['polardecomp'] not in [0, 1]:
                GMXMMPBSA_ERROR('POLARDECOMP must be either 0 or 1!', InputError)
            if INPUT['entropicdecomp'] not in [0, 1]:
                GMXMMPBSA_ERROR('ENTROPICDECOMP must be either 0 or 1!', InputError)
            for i in zip(['treeDCF', 'treeTCF', 'treeCoulomb'], [INPUT['treeDCF'], INPUT['treeTCF'],
                                                                 INPUT['treeCoulomb']]):
                if i[1] not in [0, 1]:
                    GMXMMPBSA_ERROR(f'{i[0]} must be either 0 or 1!', InputError)
            # if INPUT['thermo'] not in ['std', 'gf', 'both']:
            #     GMXMMPBSA_ERROR('THERMO must be "std", "gf", "both"!', InputError)
            # TODO: include other corrections? pc+?
            # if INPUT['thermo'] not in ['std', 'gf', 'pc+', 'all']:
            #     GMXMMPBSA_ERROR('THERMO must be "std", "gf", "pc+" or "all"!', InputError)
        if (
                not INPUT['gbrun']
                and not INPUT['pbrun']
                and not INPUT['rismrun']
                and not INPUT['nmoderun']
                and not INPUT['qh_entropy']
        ):
            GMXMMPBSA_ERROR('You did not specify any type of calculation!', InputError)

        if INPUT['gbrun'] and INPUT['PBRadii'] == 7:
            GMXMMPBSA_ERROR('PBRadii = 7 (charmm_radii) is compatible only with &pb!', InputError)

        if INPUT['decomprun'] and not INPUT['gbrun'] and not INPUT['pbrun']:
            GMXMMPBSA_ERROR('DECOMP must be run with either GB or PB!', InputError)

        if '-deo' in sys.argv and not INPUT['decomprun']:
            logging.warning("&decomp namelist has not been defined in the input file. Ignoring '-deo' flag... ")

        if (
                not INPUT['molsurf']
                and (INPUT['msoffset'] != 0 or INPUT['probe'] != 1.4)
                and self.master
        ):
            logging.warning('offset and probe are molsurf-only options')
        if INPUT['cas_intdiel'] not in [0, 1]:
            GMXMMPBSA_ERROR('cas_intdiel must be set to 0 or 1!', InputError)

        # User warning when intdiel > 10
        if self.INPUT['intdiel'] > 10:
            logging.warning('Intdiel is greater than 10...')
        # check mutant definition
        if self.INPUT['mutant'].upper() not in ['ALA', 'A', 'GLY', 'G']:
            GMXMMPBSA_ERROR('The mutant most be ALA (or A) or GLY (or G)', InputError)

        # fixed the error when try to open gmx_MMPBSA_ana in the issue
        # https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/33
        if self.INPUT['startframe'] < 1:
            # GMXMMPBSA_ERROR('The startframe variable must be >= 1')
            logging.warning(f"The startframe variable must be >= 1. Changing startframe from"
                            f" {self.INPUT['startframe']} to 1")
            self.INPUT['startframe'] = 1
        if INPUT['nmoderun']:
            if self.INPUT['nmstartframe'] < 1:
                logging.warning(f"The nmstartframe variable must be >= 1. Changing nmstartframe from"
                                f" {self.INPUT['nmstartframe']} to 1")
                self.INPUT['nmstartframe'] = 1
            if INPUT['drms'] > 0.001:
                logging.warning(f"Default DRMS value is 0.001! However {INPUT['drms']} is been used. Check "
                                f'documentation for more details...')
            if INPUT['maxcyc'] < 10000:
                logging.warning(f"Default MAXCYC value is 10000! However {INPUT['maxcyc']} is been used. Check "
                                f'documentation for more details...')

        # set the pbtemp = temperature
        self.INPUT['pbtemp'] = self.INPUT['temperature']

        logging.info(f'Checking {self.FILES.input_file} input file...Done.\n')

    def remove(self, flag):
        """ Removes temporary files """
        if not self.master:
            return
        utils.remove(flag, fnpre=self.pre)

    def sync_mpi(self):
        """ Throws up a barrier """
        self.MPI.COMM_WORLD.Barrier()

    def parse_output_files(self, from_calc=True):
        """
        This parses the output files and loads them into dicts for easy access
        """
        # Only the master does this
        from types import SimpleNamespace
        if not self.master:
            return
        logging.info('Parsing results to output files...\n')
        self.calc_types = SimpleNamespace(normal={}, mutant={}, mut_norm={}, decomp_normal={}, decomp_mutant={})
        INPUT, FILES = self.INPUT, self.FILES
        # Quasi-harmonic analysis is a special-case, so handle that separately
        if INPUT['qh_entropy']:
            if not INPUT['mutant_only']:
                self.calc_types.normal['qh'] = QHout(f'{self.pre}cpptraj_entropy.out', INPUT['temperature'])

            if INPUT['alarun']:
                self.calc_types.mutant['qh'] = QHout(f'{self.pre}mutant_cpptraj_entropy.out', INPUT['temperature'])

        # Determine if our GB is QM/MM or not
        GBClass = QMMMout if INPUT['ifqnt'] else GBout
        # Determine which kind of RISM output class we are based on std/gf and
        # polardecomp
        if INPUT['polardecomp']:
            RISM_GF = PolarRISM_gf_Out
            RISM_Std = PolarRISM_std_Out
            RISM_PCplus = PolarRISM_pcplus_Out
        else:
            RISM_GF = RISM_gf_Out
            RISM_Std = RISM_std_Out
            RISM_PCplus = RISM_pcplus_Out
        # Now we make a list of the other calculation types, their INPUT triggers,
        # their key in the calc_types dict, the base name of their output files
        # without the prefix (with %s-substitution for complex, receptor, or
        # ligand), and the class for their output
        triggers = ('nmoderun', 'gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf', 'rismrun_pcplus')
        outclass = (NMODEout, GBClass, PBout, RISM_Std, RISM_GF, RISM_PCplus)
        outkey = ('nmode', 'gb', 'pb', 'rism std', 'rism gf', 'rism pcplus')
        basename = ('%s_nm.out', '%s_gb.mdout', '%s_pb.mdout', '%s_rism.mdout', '%s_rism.mdout', '%s_rism.mdout')

        for i, key in enumerate(outkey):
            if triggers[i] not in INPUT or not INPUT[triggers[i]]:
                continue
            numframes = self.numframes_nmode if key == 'nmode' else self.numframes
            # Non-mutant
            if not INPUT['mutant_only']:
                self.calc_types.normal[key] = {'complex': outclass[i]('complex', self.INPUT, self.using_chamber)}
                self.calc_types.normal[key]['complex'].parse_from_file(self.pre + basename[i] % 'complex',
                                                                       self.mpi_size, numframes)
                # check if the nmode output is valid
                if self.calc_types.normal[key]['complex'].no_nmode_convergence:
                    self.INPUT['nmoderun'] = False
                    del self.calc_types.normal[key]
                    continue

                if not self.stability:
                    self.calc_types.normal[key]['receptor'] = outclass[i]('receptor', self.INPUT, self.using_chamber)
                    self.calc_types.normal[key]['receptor'].parse_from_file(self.pre + basename[i] % 'receptor',
                                                                            self.mpi_size, numframes)
                    self.calc_types.normal[key]['ligand'] = outclass[i]('ligand', self.INPUT, self.using_chamber)
                    self.calc_types.normal[key]['ligand'].parse_from_file(self.pre + basename[i] % 'ligand',
                                                                          self.mpi_size, numframes)
                    self.calc_types.normal[key]['delta'] = BindingStatistics(self.calc_types.normal[key]['complex'],
                                                                             self.calc_types.normal[key]['receptor'],
                                                                             self.calc_types.normal[key]['ligand'],
                                                                             self.using_chamber, self.traj_protocol)
            # Time for mutant
            if INPUT['alarun']:
                self.calc_types.mutant[key] = {'complex': outclass[i]('Mutant-Complex', self.INPUT, self.using_chamber)}
                self.calc_types.mutant[key]['complex'].parse_from_file(self.pre + 'mutant_' + basename[i] % 'complex',
                                                                       self.mpi_size, numframes)
                if not self.stability:
                    self.calc_types.mutant[key]['receptor'] = outclass[i]('Mutant-Receptor', self.INPUT,
                                                                          self.using_chamber)
                    self.calc_types.mutant[key]['receptor'].parse_from_file(self.pre + 'mutant_' + basename[i] %
                                                                            'receptor', self.mpi_size, numframes)
                    self.calc_types.mutant[key]['ligand'] = outclass[i]('Mutant-Ligand', self.INPUT,
                                                                        self.using_chamber)
                    self.calc_types.mutant[key]['ligand'].parse_from_file(self.pre + 'mutant_' + basename[i] % 'ligand',
                                                                          self.mpi_size, numframes)
                    self.calc_types.mutant[key]['delta'] = BindingStatistics(self.calc_types.mutant[key]['complex'],
                                                                             self.calc_types.mutant[key]['receptor'],
                                                                             self.calc_types.mutant[key]['ligand'],
                                                                             self.using_chamber, self.traj_protocol)
            if INPUT['alarun'] and not INPUT['mutant_only']:
                self.calc_types.mut_norm[key] = {'complex': DeltaDeltaStatistics(
                    self.calc_types.mutant[key]['complex'], self.calc_types.normal[key]['complex'])}
                if not self.stability:
                    if self.FILES.receptor_prmtop != self.FILES.mutant_receptor_prmtop:
                        self.calc_types.mut_norm[key]['receptor'] = DeltaDeltaStatistics(
                            self.calc_types.mutant[key]['receptor'], self.calc_types.normal[key]['receptor'])
                    else:
                        self.calc_types.mut_norm[key]['ligand'] = DeltaDeltaStatistics(
                            self.calc_types.mutant[key]['ligand'], self.calc_types.normal[key]['ligand'])
                    self.calc_types.mut_norm[key]['delta'] = DeltaDeltaStatistics(
                        self.calc_types.mutant[key]['delta'], self.calc_types.normal[key]['delta'])

            self.get_iec2entropy(from_calc)

        if not hasattr(self, 'resl'):
            from GMXMMPBSA.utils import mask2list
            self.resl = mask2list(FILES.complex_fixed, INPUT['receptor_mask'], INPUT['ligand_mask'])
            if INPUT['alarun']:
                self.resl[self.mutant_index].set_mut(INPUT['mutant'])

        if INPUT['decomprun']:
            self._get_decomp()

    def get_iec2entropy(self, from_calc):
        allowed_met = ['gb', 'pb', 'rism std', 'rism gf', 'rism pcplus', 'gbnsr6']
        calculated = False
        for key in allowed_met:
            if self.INPUT['interaction_entropy']:
                if not self.INPUT['mutant_only'] and key in self.calc_types.normal:
                    if from_calc:
                        edata = self.calc_types.normal[key]['delta']['GGAS']
                        ie = InteractionEntropyCalc(edata, self.INPUT)
                        ie.save_output(f'{self.pre}normal_interaction_entropy.dat')

                    self.calc_types.normal['ie'] = IEout(self.INPUT)
                    self.calc_types.normal['ie'].parse_from_file(f'{self.pre}normal_interaction_entropy.dat',
                                                                 self.numframes)
                    calculated = True
                if key in self.calc_types.mutant:
                    if from_calc:
                        edata = self.calc_types.mutant[key]['delta']['GGAS']
                        mie = InteractionEntropyCalc(edata, self.INPUT)
                        mie.save_output(f'{self.pre}mutant_interaction_entropy.dat')

                    self.calc_types.mutant['ie'] = IEout(self.INPUT)
                    self.calc_types.mutant['ie'].parse_from_file(f'{self.pre}mutant_interaction_entropy.dat',
                                                                 self.numframes)
                    calculated = True

                if self.INPUT['alarun'] and not self.INPUT['mutant_only']:
                    self.calc_types.mut_norm['ie'] = DeltaIEC2Statistic(
                        self.calc_types.mutant['ie'], self.calc_types.normal['ie'])

            if self.INPUT['c2_entropy']:
                if not self.INPUT['mutant_only'] and key in self.calc_types.normal:
                    if from_calc:
                        edata = self.calc_types.normal[key]['delta']['GGAS']
                        c2 = C2EntropyCalc(edata, self.INPUT)
                        c2.save_output(f'{self.pre}normal_c2_entropy.dat')

                    self.calc_types.normal['c2'] = C2out()
                    self.calc_types.normal['c2'].parse_from_file(f'{self.pre}normal_c2_entropy.dat')
                    calculated = True
                if key in self.calc_types.mutant:
                    if from_calc:
                        edata = self.calc_types.mutant[key]['delta']['GGAS']
                        c2 = C2EntropyCalc(edata, self.INPUT)
                        c2.save_output(f'{self.pre}mutant_c2_entropy.dat')

                    self.calc_types.mutant['c2'] = C2out()
                    self.calc_types.mutant['c2'].parse_from_file(f'{self.pre}mutant_c2_entropy.dat')
                    calculated = True

                if self.INPUT['alarun'] and not self.INPUT['mutant_only']:
                    self.calc_types.mut_norm['c2'] = DeltaIEC2Statistic(
                        self.calc_types.mutant['c2'], self.calc_types.normal['c2'])
            if calculated:
                break

    def _res2print(self):
        """
        Get residues list from print_res variable
        Returns: residues to print list
        """
        print_res = []
        for x in self.INPUT['print_res'].split(','):
            r = list(map(int, x.split('-')))
            s = r[0]
            e = r[-1] + 1
            print_res.extend(range(s, e))
        return print_res

    def _get_decomp(self):
        from GMXMMPBSA.amber_outputs import (DecompOut, PairDecompOut, DecompBinding, PairDecompBinding)
        outkey = ('gb', 'pb')
        triggers = ('gbrun', 'pbrun')
        basename = ('%s_gb.mdout', '%s_pb.mdout')
        INPUT, FILES = self.INPUT, self.FILES
        headers = {'gb': 'Generalized Born', 'pb': 'Poisson Boltzmann'}
        if INPUT['idecomp'] in [1, 2]:
            DecompBindingClass = DecompBinding
            DecompClass = DecompOut
        # Pairwise
        else:
            DecompBindingClass = PairDecompBinding
            DecompClass = PairDecompOut

        # get residues list from print_res variable
        print_res = self._res2print()
        com_list = {}
        rec_list = {}
        lig_list = {}
        for x in self.resl:
            if x.index in print_res:
                com_list[x.index] = x
                if x.is_receptor():
                    rec_list[x.id_index] = x
                else:
                    lig_list[x.id_index] = x

        for i, key in enumerate(outkey):
            if triggers[i] not in INPUT or not INPUT[triggers[i]]:
                continue
            surften = INPUT['surften'] if key == 'gb' else INPUT['cavity_surften']

            if not self.INPUT['mutant_only']:
                self.calc_types.decomp_normal[key] = {'complex': DecompClass('complex')}
                self.calc_types.decomp_normal[key]['complex'].parse_from_file(self.pre + basename[i] % 'complex',
                                                                              com_list, INPUT, surften,
                                                                              self.mpi_size, self.numframes)
                if not self.stability:
                    self.calc_types.decomp_normal[key]['receptor'] = DecompClass('receptor')
                    self.calc_types.decomp_normal[key]['receptor'].parse_from_file(self.pre + basename[i] % 'receptor',
                                                                                   rec_list, INPUT, surften,
                                                                                   self.mpi_size, self.numframes)
                    self.calc_types.decomp_normal[key]['ligand'] = DecompClass('ligand')
                    self.calc_types.decomp_normal[key]['ligand'].parse_from_file(self.pre + basename[i] % 'ligand',
                                                                                 lig_list, INPUT, surften,
                                                                                 self.mpi_size, self.numframes)
                    self.calc_types.decomp_normal[key]['delta'] = DecompBindingClass(
                        self.calc_types.decomp_normal[key]['complex'], self.calc_types.decomp_normal[key]['receptor'],
                        self.calc_types.decomp_normal[key]['ligand'], INPUT,
                        f'Energy Decomposition Analysis (All units kcal/mol): {headers[key]} model')

            if INPUT['alarun']:
                # Do mutant
                self.calc_types.decomp_mutant[key] = {'complex': DecompClass('Mutant-Complex')}
                self.calc_types.decomp_mutant[key]['complex'].parse_from_file(
                    (f'{self.pre}mutant_' + basename[i] % 'complex'),
                    com_list,
                    INPUT,
                    surften,
                    self.mpi_size,
                    self.numframes,
                    True
                )

                if not self.stability:
                    self.calc_types.decomp_mutant[key]['receptor'] = DecompClass('Mutant-Receptor')
                    self.calc_types.decomp_mutant[key]['receptor'].parse_from_file(
                        (f'{self.pre}mutant_' + basename[i] % 'receptor'),
                        rec_list,
                        INPUT,
                        surften,
                        self.mpi_size,
                        self.numframes,
                        True
                    )

                    self.calc_types.decomp_mutant[key]['ligand'] = DecompClass('Mutant-Ligand')
                    self.calc_types.decomp_mutant[key]['ligand'].parse_from_file(
                        (f'{self.pre}mutant_' + basename[i] % 'ligand'),
                        lig_list,
                        INPUT,
                        surften,
                        self.mpi_size,
                        self.numframes,
                        True
                    )

                    self.calc_types.decomp_mutant[key]['delta'] = DecompBindingClass(
                        self.calc_types.decomp_mutant[key]['complex'], self.calc_types.decomp_mutant[key]['receptor'],
                        self.calc_types.decomp_mutant[key]['ligand'], INPUT,
                        f'Energy Decomposition Analysis (All units kcal/mol): {headers[key]} model ({self.mut_str})')


# Local methods

def excepthook(exception_type, exception_value, tb):
    """
    Replaces sys.excepthook so fatal exceptions kill all MPI threads and we can
    control the printing of tracebacks. Those are helpful for debugging purposes,
    but may be unsightly to users. debug_printlevel set above controls this
    behavior
    """
    import traceback
    global _stderr, _mpi_size, _rank
    if not isinstance(exception_type, MMPBSA_Error):
        traceback.print_tb(tb)
    _stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
    if _mpi_size > 1:
        _stderr.write('Error occurred on rank %d.' % _rank + os.linesep)
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
