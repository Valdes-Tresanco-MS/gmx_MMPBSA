"""
This module contains classes for creating input files for various MM/PBSA
calculations.

Methods:
   create_inputs -- determines which input files need to be written then 
                    handles it

Classes: 
    SanderInput -- Base class for sander input files
    SanderGBInput -- writes input files for sander GB (derived from SanderInput)
    SanderPBSAInput- writes input files for sander PB (derived from SanderInput)
    SanderAPBSInput -- writes input files for sander.APBS PB 
                       (derived from SanderInput)
    SanderGBDecomp -- writes GB decomp input files (derived from SanderInput)
    SanderPBDecomp -- writes PB decomp input files (derived from SanderInput)
    QuasiHarmonicInput -- writes a cpptraj input file for quasi-harmonic calcs
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

from __future__ import division
import logging

import parmed
from parmed.amber.mdin import Mdin
from parmed.exceptions import AmberError
from copy import deepcopy

ROH = {1: 0.586, 2: 0.699, 3: 0.734, 4: 0.183}


def create_inputs(INPUT, prmtop_system, pre):
    """ Creates the input files for all necessary calculations """
    stability = prmtop_system.stability

    # since gbnsr6 is independent of sander MM calculation is the same input when decomp or not
    if INPUT['gbnsr6']['gbnsr6run']:
        temp_input = deepcopy(INPUT)
        temp_input['gbnsr6']['roh'] = ROH[INPUT['gbnsr6']['roh']]
        if INPUT['decomp']['decomprun']:
            temp_input['gbnsr6']['dgij'] = 1
        gbnsr6_mdin = GBNSR6Input(temp_input)
        gbnsr6_mdin.make_mdin()
        gbnsr6_mdin.write_input(f'{pre}gbnsr6.mdin')

    # First check if we are running decomp
    if INPUT['decomp']['decomprun']:
        # Get the cards that will go in the complex mdin file
        com_card, rc_card, lc_card = prmtop_system.Group(INPUT['decomp']['print_res'], True)
        junk, rec_card, lig_card = prmtop_system.Group(INPUT['decomp']['print_res'], False,
                                                       rec_group_type='RES', lig_group_type='RES')
        full_com, full_rc, full_lc = prmtop_system.Group('all', True)
        junk, full_rec, full_lig = prmtop_system.Group('all', False)

        # Convert the strings into card objects
        # Now create the mdin objects
        if INPUT['gb']['gbrun']:
            com_arad = None
            rec_arad = None
            lig_arad = None
            if INPUT['gb']['alpb']:
                import subprocess
                # FIXME: save pqr file from prmtop_system
                com_top = parmed.load_file('COM.prmtop', xyz=f"{pre}COM.inpcrd")
                com_top.save('COM.pqr', format='pqr', overwrite=True)
                stdoutdata = subprocess.getoutput("elsize COM.pqr -hea")
                com_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]
                if not stability:
                    rec_top = parmed.load_file('REC.prmtop', xyz=f"{pre}REC.inpcrd")
                    rec_top.save('REC.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize REC.pqr -hea")
                    rec_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]

                    lig_top = parmed.load_file('LIG.prmtop', xyz=f"{pre}LIG.inpcrd")
                    lig_top.save('LIG.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize LIG.pqr -hea")
                    lig_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]

            rec_res = ['Residues considered as REC', full_rc]
            if stability:
                pri_res = ['Residues to print', com_card]
                com_input = deepcopy(INPUT)
                if com_arad:
                    com_input['gb']['arad'] = com_arad
                    logging.warning(f'Setting complex arad = {com_arad}')
                com_mdin = SanderGBDecomp(com_input, rec_res, pri_res)
                com_mdin.write_input(f"{pre}gb_decomp_com.mdin")
            else:
                lig_res = ['Residues considered as LIG', full_lc]
                pri_res = ['Residues to print', com_card]
                com_input = deepcopy(INPUT)
                if com_arad:
                    com_input['gb']['arad'] = com_arad
                    logging.warning(f'Setting complex arad = {com_arad}')
                com_mdin = SanderGBDecomp(com_input, rec_res, lig_res, pri_res)
                rec_res = ['Residues considered as REC', full_rec]
                pri_res = ['Residues to print', rec_card]
                rec_input = deepcopy(INPUT)
                if rec_arad:
                    rec_input['gb']['arad'] = rec_arad
                    logging.warning(f'Setting receptor arad = {rec_arad}')
                rec_mdin = SanderGBDecomp(rec_input, rec_res, pri_res)
                lig_res = ['Residues considered as LIG', full_lig]
                pri_res = ['Residues to print', lig_card]
                lig_input = deepcopy(INPUT)
                if lig_arad:
                    lig_input['gb']['arad'] = lig_arad
                    logging.warning(f'Setting ligand arad = {lig_arad}')
                lig_mdin = SanderGBDecomp(lig_input, lig_res, pri_res)
                com_mdin.write_input(f"{pre}gb_decomp_com.mdin")
                rec_mdin.write_input(f"{pre}gb_decomp_rec.mdin")
                lig_mdin.write_input(f"{pre}gb_decomp_lig.mdin")

        if INPUT['pb']['pbrun']:
            rec_res = ['Residues considered as REC', full_rc]
            if stability:
                pri_res = ['Residues to print', com_card]
                com_mdin = SanderPBDecomp(INPUT, rec_res, pri_res)
                com_mdin.write_input(f"{pre}pb_decomp_com.mdin")
            else:
                lig_res = ['Residues considered as LIG', full_lc]
                pri_res = ['Residues to print', com_card]
                com_mdin = SanderPBDecomp(INPUT, rec_res, lig_res, pri_res)
                rec_res = ['Residues considered as REC', full_rec]
                pri_res = ['Residues to print', rec_card]
                rec_mdin = SanderPBDecomp(INPUT, rec_res, pri_res)
                lig_res = ['Residues considered as LIG', full_lig]
                pri_res = ['Residues to print', lig_card]
                lig_mdin = SanderPBDecomp(INPUT, lig_res, pri_res)
                com_mdin.write_input(f"{pre}pb_decomp_com.mdin")
                rec_mdin.write_input(f"{pre}pb_decomp_rec.mdin")
                lig_mdin.write_input(f"{pre}pb_decomp_lig.mdin")

        # Require one file for GBNSR6, pbsa.cuda and APBS calculations.
        # TODO: we need to define the intdiel for decomp because the eel term is computed by sander and not by GBNSR6
        if INPUT['gbnsr6']['gbnsr6run']:
            rec_res = ['Residues considered as REC', full_rc]
            if stability:
                pri_res = ['Residues to print', com_card]
                com_mdin = SanderMMDecomp(INPUT, 'gbnsr6', rec_res, pri_res)
                com_mdin.write_input(f"{pre}mm_gbnsr6_decomp_com.mdin")
            else:
                lig_res = ['Residues considered as LIG', full_lc]
                pri_res = ['Residues to print', com_card]
                com_mdin = SanderMMDecomp(INPUT, 'gbnsr6', rec_res, lig_res, pri_res)
                rec_res = ['Residues considered as REC', full_rec]
                pri_res = ['Residues to print', rec_card]
                rec_mdin = SanderMMDecomp(INPUT, 'gbnsr6', rec_res, pri_res)
                lig_res = ['Residues considered as LIG', full_lig]
                pri_res = ['Residues to print', lig_card]
                lig_mdin = SanderMMDecomp(INPUT, 'gbnsr6', lig_res, pri_res)
                com_mdin.write_input(f"{pre}mm_gbnsr6_decomp_com.mdin")
                rec_mdin.write_input(f"{pre}mm_gbnsr6_decomp_rec.mdin")
                lig_mdin.write_input(f"{pre}mm_gbnsr6_decomp_lig.mdin")
    else:  # not decomp

        if INPUT['gb']['gbrun']:
            # We need separate input files for QM/gmx_MMPBSA
            com_arad = None
            rec_arad = None
            lig_arad = None
            if INPUT['gb']['alpb']:
                import subprocess
                com_top = parmed.load_file('COM.prmtop', xyz=f"{pre}COM.inpcrd")
                com_top.save('COM.pqr', format='pqr', overwrite=True)
                stdoutdata = subprocess.getoutput("elsize COM.pqr -hea")
                com_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]
                if not stability:
                    rec_top = parmed.load_file('REC.prmtop', xyz=f"{pre}REC.inpcrd")
                    rec_top.save('REC.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize REC.pqr -hea")
                    rec_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]

                    lig_top = parmed.load_file('LIG.prmtop', xyz=f"{pre}LIG.inpcrd")
                    lig_top.save('LIG.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize LIG.pqr -hea")
                    lig_arad = stdoutdata.split()[INPUT['gb']['arad_method'] * -1]

            if INPUT['gb']['ifqnt']:
                com_input = deepcopy(INPUT)
                rec_input = deepcopy(INPUT)
                lig_input = deepcopy(INPUT)
                (com_input['gb']['qmmask'], rec_input['gb']['qmmask'],
                 lig_input['gb']['qmmask']) = prmtop_system.Mask(INPUT['gb']['qm_residues'], in_complex=False)
                if not com_input['gb']['qmmask']:
                    raise AmberError('No valid QM residues chosen!')
                com_input['gb']['qm_theory'] = "'%s'" % com_input['gb']['qm_theory']
                com_input['gb']['qmmask'] = "'%s'" % com_input['gb']['qmmask']
                com_input['gb']['qmcharge'] = com_input['gb']['qmcharge_com']
                # check if alpb
                if com_arad:
                    com_input['gb']['arad'] = com_arad
                    logging.warning(f'Setting complex arad = {com_arad}')
                gb_mdin = SanderGBInput(com_input)
                gb_mdin.make_mdin()
                gb_mdin.write_input(f'{pre}gb_qmmm_com.mdin')
                if not stability:
                    if not rec_input['gb']['qmmask']:
                        rec_input['gb']['ifqnt'] = 0
                    else:
                        rec_input['gb']['qmmask'] = "'%s'" % rec_input['gb']['qmmask']
                    rec_input['gb']['qm_theory'] = "'%s'" % rec_input['gb']['qm_theory']
                    rec_input['gb']['qmcharge'] = rec_input['gb']['qmcharge_rec']
                    # check if alpb
                    if rec_arad:
                        rec_input['gb']['arad'] = rec_arad
                        logging.warning(f'Setting receptor arad = {rec_arad}')
                    gb_mdin = SanderGBInput(rec_input)
                    gb_mdin.make_mdin()
                    gb_mdin.write_input(f'{pre}gb_qmmm_rec.mdin')
                    if not lig_input['gb']['qmmask']:
                        lig_input['gb']['ifqnt'] = 0
                    else:
                        lig_input['gb']['qmmask'] = "'%s'" % lig_input['gb']['qmmask']
                    lig_input['gb']['qm_theory'] = "'%s'" % lig_input['gb']['qm_theory']
                    lig_input['gb']['qmcharge'] = lig_input['gb']['qmcharge_lig']
                    # check if alpb
                    if lig_arad:
                        lig_input['gb']['arad'] = lig_arad
                        logging.warning(f'Setting ligand arad = {lig_arad}')
                    gb_mdin = SanderGBInput(lig_input)
                    gb_mdin.make_mdin()
                    gb_mdin.write_input(f'{pre}gb_qmmm_lig.mdin')

            elif INPUT['gb']['alpb']:
                com_input = deepcopy(INPUT)
                rec_input = deepcopy(INPUT)
                lig_input = deepcopy(INPUT)
                com_input['gb']['arad'] = com_arad
                logging.warning(f'Setting complex arad = {com_arad}')
                gb_mdin = SanderGBInput(com_input)
                gb_mdin.make_mdin()
                gb_mdin.write_input(f'{pre}gb_com.mdin')
                rec_input['gb']['arad'] = rec_arad
                logging.warning(f'Setting receptor arad = {rec_arad}')
                gb_mdin = SanderGBInput(rec_input)
                gb_mdin.make_mdin()
                gb_mdin.write_input(f'{pre}gb_rec.mdin')
                lig_input['gb']['arad'] = lig_arad
                logging.warning(f'Setting ligand arad = {lig_arad}')
                gb_mdin = SanderGBInput(lig_input)
                gb_mdin.make_mdin()
                gb_mdin.write_input(f'{pre}gb_lig.mdin')
            else:
                gb_mdin = SanderGBInput(INPUT)
                gb_mdin.make_mdin()
                gb_mdin.write_input(f'{pre}gb.mdin')

        # We only need to run it once for the GBNSR6, pbsa.cuda and APBS calculations.
        if INPUT['gbnsr6']['gbnsr6run']:
            mm_mdin = SanderMMInput(INPUT)
            mm_mdin.set_gbnsr6_param()
            mm_mdin.make_mdin()
            mm_mdin.write_input(f'{pre}mm.mdin')

        if INPUT['pb']['pbrun']:
            pb_prog = 'sander.APBS' if INPUT['pb']['sander_apbs'] else 'sander'
            if pb_prog == 'sander.APBS':
                pb_mdin = SanderAPBSInput(INPUT)
                pb_mdin.make_mdin()
            else:
                pb_mdin = SanderPBSAInput(INPUT)
                pb_mdin.make_mdin()
            pb_mdin.write_input(f'{pre}pb.mdin')

        if INPUT['rism']['rismrun']:
            rism_mdin = SanderRISMInput(INPUT)
            rism_mdin.make_mdin()
            rism_mdin.write_input(f'{pre}rism.mdin')

    # end if decomprun

    if INPUT['general']['qh_entropy']:  # quasi-harmonic approximation input file
        trj_suffix = 'nc' if INPUT['general']['netcdf'] else 'mdcrd'
        com_mask, rec_mask, lig_mask = prmtop_system.Mask('all', True)
        if not INPUT['ala']['mutant_only']:
            qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask, temperature=INPUT['temperature'],
                                       stability=stability, prefix=pre, trj_suffix=trj_suffix)
            qh_in.write_input(f'{pre}cpptrajentropy.in')
        if INPUT['ala']['alarun']:
            qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask, temperature=INPUT['temperature'],
                                       stability=stability, prefix=pre + 'mutant_', trj_suffix=trj_suffix)
            qh_in.write_input(f'{pre}mutant_cpptrajentropy.in')


class SanderInput(object):
    """ Base class sander input file """

    def __init__(self, INPUT):
        self.program = 'sander'  # This runs with sander
        self.input_items = {'foo': 'bar'}  # replace this in derived classes
        self.name_map = {'foo': 'orig'}  # replace this in derived classes
        self.parent_namelist = {'foo': 'foo_namelist'}  # replace this in derived classes
        self.namelist = 'foo'
        self.INPUT = INPUT

    def make_mdin(self):
        self.mdin = Mdin(self.program)
        self.mdin.title = 'File generated by gmx_MMPBSA'
        for key, value in self.input_items.items():
            # Skip ioutfm since it is handled explicitly later
            if key == 'ioutfm':
                self.mdin.change('cntrl', 'ioutfm', int(bool(self.INPUT['general']['netcdf'])))
                continue
            if self.INPUT.get(self.namelist) and self.name_map[key] in self.INPUT[self.namelist]:
                self.mdin.change(self.parent_namelist[key], key, self.INPUT[self.namelist][self.name_map[key]])
            elif key in ['dec_verbose', 'idecomp']:
                self.mdin.change(self.parent_namelist[key], key, self.INPUT['decomp'][self.name_map[key]])
            else:
                self.mdin.change(self.parent_namelist[key], key, value)
        if self.namelist in {'gb', 'pb', 'rism'}:
            self.mdin.change('cntrl', 'ioutfm', int(bool(self.INPUT['general']['netcdf'])))

        if self.namelist in ['pb', 'gbnsr6']:
            # in parmed.amber.mdin gbnsr6 namelist not exists, in this case, is gb, we need to change this namelist
            # variable in mdin, but use gbnsr6 as namelist to get parameters from the the INPUT
            nl = 'pb' if self.namelist == 'pb' else 'gb'
            self.mdin.change(nl, 'istrng', self.INPUT[self.namelist]['istrng'] * 1000)

    def write_input(self, filename):
        """ Write the mdin file """
        self.mdin.write(filename)


class SanderGBInput(SanderInput):
    """ GB sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'idecomp': 0, 'offset': -999999.0,
                   'imin': 5, 'maxcyc': 1, 'ncyc': 0, 'gbsa': 0, 'ioutfm': 0, 'dec_verbose': 0,
                   # Basic options
                   'igb': 5, 'intdiel': 1.0, 'extdiel': 78.5, 'saltcon': 0.0, 'surften': 0.0072,
                   # QM options
                   'ifqnt': 0, 'qmmask': '', 'qm_theory': '', 'qmcharge': 0,
                   'qmgb': 2, 'qmcut': 999.0,
                   'scfconv': 1.0e-8, 'peptide_corr': 0, 'writepdb': 1, 'verbosity': 0, 'alpb': 0, 'arad': 15}

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'idecomp': 'cntrl', 'offset': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ncyc': 'cntrl', 'gbsa': 'cntrl', 'ioutfm': 'cntrl',
                       'dec_verbose': 'cntrl',
                       # Basic options
                       'igb': 'cntrl', 'intdiel': 'cntrl', 'extdiel': 'cntrl', 'saltcon': 'cntrl',
                       'surften': 'cntrl', 'alpb': 'cntrl', 'arad': 'cntrl',
                       # QM options
                       'ifqnt': 'cntrl', 'qmmask': 'qmmm', 'qm_theory': 'qmmm', 'qmcharge': 'qmmm',
                       'qmgb': 'qmmm', 'qmcut': 'qmmm',
                       'scfconv': 'qmmm', 'peptide_corr': 'qmmm', 'writepdb': 'qmmm', 'verbosity': 'qmmm'}

        self.name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'idecomp': 'idecomp', 'offset': 'offset',
                'imin': 'imin', 'gbsa': 'gbsa', 'ioutfm': 'netcdf', 'dec_verbose': 'dec_verbose',
                'maxcyc': 'gb_maxcyc', 'ncyc': 'ncyc',
                # Basic options
                'igb': 'igb', 'intdiel': 'intdiel', 'extdiel': 'extdiel', 'saltcon': 'saltcon',
                'surften': 'surften', 'alpb': 'alpb', 'arad': 'arad',
                # QM options
                'ifqnt': 'ifqnt',  'qmmask': 'qmmask', 'qm_theory': 'qm_theory',
                'qmcharge': 'qmcharge', 'qmgb': 'qmgb', 'qmcut': 'qmcut',
                'scfconv': 'scfconv', 'peptide_corr': 'peptide_corr', 'writepdb': 'writepdb', 'verbosity': 'verbosity'}
        self.namelist = 'gb'


class GBNSR6Input(SanderInput):
    """ GB sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.program = 'gbnsr6'
        self.input_items = {'inp': 1,
                            'b': 0.028, 'alpb': 1, 'epsin': 1.0, 'epsout': 78.5, 'istrng': 0.0, 'rs': 0.52,
                            'dprob': 1.4, 'space': 0.5, 'arcres': 0.2, 'rbornstat': 0, 'dgij': 0, 'radiopt': 0,
                            'chagb': 0, 'roh': 1, 'tau': 1.47, 'cavity_surften': 0.005}

        self.parent_namelist = {'inp': 'cntrl',
                                'b': 'gb', 'alpb': 'gb', 'epsin': 'gb', 'epsout': 'gb', 'istrng': 'gb', 'rs': 'gb',
                                'dprob': 'gb', 'space': 'gb', 'arcres': 'gb', 'rbornstat': 'gb', 'dgij': 'gb',
                                'radiopt': 'gb', 'chagb': 'gb', 'roh': 'gb', 'tau': 'gb', 'cavity_surften': 'gb'}

        self.name_map = {'inp': 'inp',
                         'b': 'b', 'alpb': 'alpb', 'epsin': 'epsin', 'epsout': 'epsout', 'istrng': 'istrng',
                         'rs': 'rs', 'dprob': 'dprob', 'space': 'space', 'arcres': 'arcres', 'rbornstat': 'rbornstat',
                         'dgij': 'dgij', 'radiopt': 'radiopt', 'chagb': 'chagb', 'roh': 'roh', 'tau': 'tau',
                         'cavity_surften': 'cavity_surften'}
        self.namelist = 'gbnsr6'


class SanderMMInput(SanderInput):
    """ GB sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'idecomp': 0, 'offset': -999999.0,
                   'imin': 5, 'maxcyc': 1, 'ncyc': 0, 'gbsa': 0, 'ioutfm': 0, 'dec_verbose': 0,
                   # Basic options
                   'igb': 5, 'intdiel': 1.0, 'extdiel': 78.5, 'saltcon': 0.0, 'surften': 0.0072,
                   }

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'idecomp': 'cntrl', 'offset': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ncyc': 'cntrl', 'gbsa': 'cntrl', 'ioutfm': 'cntrl',
                       'dec_verbose': 'cntrl',
                       # Basic options
                       'igb': 'cntrl', 'intdiel': 'cntrl', 'extdiel': 'cntrl', 'saltcon': 'cntrl', 'surften': 'cntrl'
                       }

        self.name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'idecomp': 'idecomp', 'offset': 'offset',
                'imin': 'imin', 'gbsa': 'gbsa', 'ioutfm': 'netcdf', 'dec_verbose': 'dec_verbose',
                'maxcyc': 'gb_maxcyc', 'ncyc': 'ncyc',
                # Basic options
                'igb': 'igb', 'intdiel': 'intdiel', 'extdiel': 'extdiel', 'saltcon': 'saltcon', 'surften': 'surften'
                }

        self.namelist = 'mm'

    def set_gbnsr6_param(self):
        transferable = {'epsin': 'intdiel', 'epsout': 'extdiel', 'istrng': 'saltcon', 'cavity_surften': 'surften'}
        for gbnsr6k, gbk in transferable.items():
            self.input_items[gbk] = self.INPUT['gbnsr6'][gbnsr6k]


class SanderMMDecomp(SanderMMInput):
    """ MM decomposition input file for sander. In addition to the INPUT dictionary,
       this class also needs several GROUP cards to define the
   """

    def __init__(self, INPUT, method, *cards):
        super().__init__(INPUT)
        # SanderMMInput.__init__(self, INPUT)
        self.INPUT = INPUT
        # update parameters before make mdin and adding cards
        if method == 'gbnsr6':
            self.set_gbnsr6_param()

        self.input_items['gbsa'] = 2
        self.input_items['dec_verbose'] = self.INPUT['decomp']['dec_verbose']

        self.make_mdin()
        for card in cards:
            self.mdin.AddCard(card[0], card[1])

    def set_gbnsr6_param(self):
        transferable = {'epsin': 'intdiel', 'epsout': 'extdiel', 'istrng':'saltcon', 'cavity_surften': 'surften'}
        for gbnsr6k, gbk in transferable.items():
            self.input_items[gbk] = self.INPUT['gbnsr6'][gbnsr6k]


class SanderRISMInput(SanderInput):
    """ RISM sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)

        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'idecomp': 0, 'offset': -999999.0,
                   'imin': 5, 'maxcyc': 1, 'ncyc': 0, 'gbsa': 0, 'ioutfm': 0, 'dec_verbose': 0, 'irism': 1,
                   # Closure approximations
                   'closure': "kh",
                   # Solvation free energy corrections
                   'gfcorrection': 0, 'pcpluscorrection': 0,
                   # Long-range asymptotics
                   'noasympcorr': 1, 'treedcf': 1, 'treetcf': 1, 'treecoulomb': 0, 'treedcfmac': 0.1,
                   'treetcfmac': 0.1, 'treecoulombmac': 0.1, 'treedcforder': 2, 'treetcforder': 2,
                   'treecoulomborder': 2, 'treedcfn0': 500, 'treetcfn0': 500, 'treecoulombn0': 500,
                   # Solvation box
                   'buffer': 14.0, 'grdspc': 0.5, 'ng3': '-1,-1,-1', 'solvbox': '-1,-1,-1', 'solvcut': -1.0,
                   # Solution convergence
                   'tolerance': 1e-05, 'ljtolerance': -1.0, 'asympkspacetolerance': -1.0, 'mdiis_del': 0.7,
                   'mdiis_nvec': 5, 'mdiis_restart': 10.0, 'maxstep': 10000, 'npropagate': 5,
                   # Output
                   'polardecomp': 0, 'entropicdecomp': 0, 'verbose': 0}

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'idecomp': 'cntrl', 'offset': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ncyc': 'cntrl', 'gbsa': 'cntrl', 'ioutfm': 'cntrl',
                       'dec_verbose': 'cntrl', 'irism': 'cntrl',
                       # Closure approximations
                       'closure': 'rism',
                       # Solvation free energy corrections
                       'gfcorrection': 'rism', 'pcpluscorrection': 'rism',
                       # Long-range asymptotics
                       'noasympcorr': 'rism', 'treedcf': 'rism', 'treetcf': 'rism', 'treecoulomb': 'rism',
                       'treedcfmac': 'rism', 'treetcfmac': 'rism', 'treecoulombmac': 'rism', 'treedcforder': 'rism',
                       'treetcforder': 'rism', 'treecoulomborder': 'rism', 'treedcfn0': 'rism', 'treetcfn0': 'rism',
                       'treecoulombn0': 'rism',
                       # Solvation box
                       'buffer': 'rism', 'grdspc': 'rism', 'ng3': 'rism', 'solvbox': 'rism', 'solvcut': 'rism',
                       # Solution convergence
                       'tolerance': 'rism', 'ljtolerance': 'rism', 'asympkspacetolerance': 'rism',  'mdiis_del':
                       'rism', 'mdiis_nvec': 'rism', 'mdiis_restart': 'rism', 'maxstep': 'rism',
                       'npropagate': 'rism',
                       # Output
                       'polardecomp': 'rism', 'entropicdecomp': 'rism', 'verbose': 'rism'}

        self.name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'idecomp': 'idecomp', 'offset': 'offset',
                'imin': 'imin', 'gbsa': 'gbsa', 'ioutfm': 'netcdf', 'dec_verbose': 'dec_verbose',
                'maxcyc': 'rism_maxcyc', 'ncyc': 'ncyc', 'irism':'irism',
                # Closure approximations
                'closure': 'closure',
                # Solvation free energy corrections
                'gfcorrection': 'gfcorrection', 'pcpluscorrection': 'pcpluscorrection',
                # Long-range asymptotics
                'noasympcorr': 'noasympcorr', 'treedcf': 'treedcf',
                'treetcf': 'treetcf', 'treecoulomb': 'treecoulomb', 'treedcfmac': 'treedcfmac',
                'treetcfmac': 'treetcfmac', 'treecoulombmac': 'treecoulombmac', 'treedcforder': 'treedcforder',
                'treetcforder': 'treetcforder', 'treecoulomborder': 'treecoulomborder', 'treedcfn0': 'treedcfn0',
                'treetcfn0': 'treetcfn0', 'treecoulombn0': 'treecoulombn0',
                # Solvation box
                'buffer': 'buffer', 'grdspc': 'grdspc', 'ng3': 'ng', 'solvbox': 'solvbox', 'solvcut': 'solvcut',
                # Solution convergence
                'tolerance': 'tolerance', 'ljtolerance': 'ljtolerance', 'asympkspacetolerance': 'asympkspacetolerance',
                'mdiis_del': 'mdiis_del', 'mdiis_nvec': 'mdiis_nvec', 'mdiis_restart': 'mdiis_restart',
                'maxstep': 'maxstep', 'npropagate': 'npropagate',
                # Output
                'polardecomp': 'polardecomp', 'entropicdecomp': 'entropicdecomp', 'verbose': 'rism_verbose'}

        self.namelist = 'rism'


class SanderPBSADECOMPInput(SanderInput):
    """ PB sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'ioutfm': 0, 'idecomp': 0, 'dec_verbose': 0,
                   'imin': 5, 'maxcyc': 1, 'ntx': 1, 'pbtemp': 300,
                   # Basic input options
                   'ipb': 2, 'inp': 2,
                   # Options to define the physical constants
                   'epsin': 1.0, 'epsout': 80.0, 'epsmem': 1.0, 'smoothopt': 1, 'istrng': 0.0,
                   'radiopt': 1, 'dprob': 1.4, 'iprob': 2.0, 'sasopt': 0, 'arcres': 0.25,
                   # Options for Implicit Membranes
                   'membraneopt': 0, 'mprob': 2.70, 'mthick': 40, 'mctrdz': 0.0,
                   'poretype': 1,
                   # Options to select numerical procedures
                   'npbopt': 0, 'solvopt': 1, 'accept': 0.001, 'maxitn': 1000, 'fillratio': 4.0,
                   'space': 0.5, 'nbuffer': 0, 'nfocus': 2, 'fscale': 8, 'npbgrid': 1,
                   # Options to compute energy and forces
                   'bcopt': 5, 'eneopt': 2, 'frcopt': 0, 'scalec': 0, 'cutfd': 5.0, 'cutnb': 0.0,
                   'nsnba': 1,
                   # Options to select a non - polar solvation treatment
                   'decompopt': 2, 'use_rmin': 1, 'sprob': 0.557, 'vprob': 1.300,
                   'rhow_effect': 1.129, 'use_sav': 1, 'cavity_surften': 0.0378,
                   'cavity_offset': -0.5692, 'maxsph': 400, 'maxarcdot': 1500,
                   # Options for output
                   'npbverb': 0}

        self.name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'ioutfm': 'netcdf',
                'idecomp': 'idecomp', 'dec_verbose': 'dec_verbose',
                'imin': 'imin', 'maxcyc': 'pb_maxcyc', 'ntx': 'ntx', 'pbtemp': 'temperature',
                # Basic input options
                'ipb': 'ipb', 'inp': 'inp',
                # Options to define the physical constants
                'epsin': 'indi', 'epsout': 'exdi', 'epsmem': 'emem', 'smoothopt': 'smoothopt',
                'istrng': 'istrng', 'radiopt': 'radiopt', 'dprob': 'prbrad', 'iprob': 'iprob',
                'sasopt': 'sasopt', 'arcres': 'arcres',
                # Options for Implicit Membranes
                'membraneopt': 'memopt', 'mprob': 'mprob', 'mthick': 'mthick', 'mctrdz': 'mctrdz',
                'poretype': 'poretype',
                # Options to select numerical procedures
                'npbopt': 'npbopt', 'solvopt': 'solvopt', 'accept': 'accept', 'maxitn': 'linit',
                'fillratio': 'fillratio', 'space': 'scale', 'nbuffer': 'nbuffer', 'nfocus': 'nfocus',
                'fscale': 'fscale', 'npbgrid': 'npbgrid',
                # Options to compute energy and forces
                'bcopt': 'bcopt', 'eneopt': 'eneopt', 'frcopt': 'frcopt', 'scalec': 'scalec',
                'cutfd': 'cutfd', 'cutnb': 'cutnb', 'nsnba': 'nsnba',
                # Options to select a non-polar solvation treatment
                'decompopt': 'decompopt', 'use_rmin': 'use_rmin', 'sprob': 'sprob',
                'vprob': 'vprob', 'rhow_effect': 'rhow_effect', 'use_sav': 'use_sav',
                'cavity_surften': 'cavity_surften', 'cavity_offset': 'cavity_offset',
                'maxsph': 'maxsph', 'maxarcdot': 'maxarcdot',
                # Options for output
                'npbverb': 'npbverb'}

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'ioutfm': 'cntrl',
                       'idecomp': 'cntrl', 'dec_verbose': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ntx': 'cntrl', 'pbtemp': 'pb',
                       # Basic input options
                       'ipb': 'cntrl', 'inp': 'cntrl',
                       # Options to define the physical constants
                       'epsin': 'pb', 'epsout': 'pb', 'epsmem': 'pb', 'smoothopt': 'pb',
                       'istrng': 'pb', 'radiopt': 'pb', 'dprob': 'pb', 'iprob': 'pb',
                       'sasopt': 'pb', 'arcres': 'pb',
                       # Options for Implicit Membranes
                       'membraneopt': 'pb', 'mprob': 'pb', 'mthick': 'pb', 'mctrdz': 'pb',
                       'poretype': 'pb',
                       # Options to select numerical procedures
                       'npbopt': 'pb', 'solvopt': 'pb', 'accept': 'pb', 'maxitn': 'pb',
                       'fillratio': 'pb', 'space': 'pb', 'nbuffer': 'pb', 'nfocus': 'pb',
                       'fscale': 'pb', 'npbgrid': 'pb',
                       # Options to compute energy and forces
                       'bcopt': 'pb', 'eneopt': 'pb', 'frcopt': 'pb', 'scalec': 'pb',
                       'cutfd': 'pb', 'cutnb': 'pb', 'nsnba': 'pb',
                       # Options to select a non-polar solvation treatment
                       'decompopt': 'pb', 'use_rmin': 'pb', 'sprob': 'pb', 'vprob': 'pb',
                       'rhow_effect': 'pb', 'use_sav': 'pb', 'cavity_surften': 'pb',
                       'cavity_offset': 'pb', 'maxsph': 'pb', 'maxarcdot': 'pb',
                       # Options for output
                       'npbverb': 'pb'}
        self.namelist = 'pb'

        # self.mdin.change('pb', 'istrng', INPUT['pb']['istrng'] * 1000)
    # def __init__(self, INPUT):
    #     # We need to change istrng to mM (from M).
    #     SanderInput.__init__(self, INPUT)
    #     self.mdin.change('pb', 'istrng', INPUT['pb']['istrng'] * 1000)


class SanderPBSAInput(SanderInput):
    """ PB sander input file """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'ioutfm': 0, 'idecomp': 0, 'dec_verbose': 0,
                   'imin': 5, 'maxcyc': 1, 'ntx': 1, 'pbtemp': 300,
                   # Basic input options
                   'ipb': 2, 'inp': 2,
                   # Options to define the physical constants
                   'epsin': 1.0, 'epsout': 80.0, 'epsmem': 1.0, 'smoothopt': 1, 'istrng': 0.0,
                   'radiopt': 1, 'dprob': 1.4, 'iprob': 2.0, 'sasopt': 0, 'arcres': 0.25,
                   # Options for Implicit Membranes
                   'membraneopt': 0, 'mprob': 2.70, 'mthick': 40, 'mctrdz': 0.0,
                   'poretype': 1,
                   # Options to select numerical procedures
                   'npbopt': 0, 'solvopt': 1, 'accept': 0.001, 'maxitn': 1000, 'fillratio': 4.0,
                   'space': 0.5, 'nbuffer': 0, 'nfocus': 2, 'fscale': 8, 'npbgrid': 1,
                   # Options to compute energy and forces
                   'bcopt': 5, 'eneopt': 2, 'frcopt': 0, 'scalec': 0, 'cutfd': 5.0, 'cutnb': 0.0,
                   'nsnba': 1,
                   # Options to select a non - polar solvation treatment
                   'decompopt': 2, 'use_rmin': 1, 'sprob': 0.557, 'vprob': 1.300,
                   'rhow_effect': 1.129, 'use_sav': 1, 'cavity_surften': 0.0378,
                   'cavity_offset': -0.5692, 'maxsph': 400, 'maxarcdot': 1500,
                   # Options for output
                   'npbverb': 0}

        self.name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'ioutfm': 'netcdf',
                'idecomp': 'idecomp', 'dec_verbose': 'dec_verbose',
                'imin': 'imin', 'maxcyc': 'pb_maxcyc', 'ntx': 'ntx', 'pbtemp': 'temperature',
                # Basic input options
                'ipb': 'ipb', 'inp': 'inp',
                # Options to define the physical constants
                'epsin': 'indi', 'epsout': 'exdi', 'epsmem': 'emem', 'smoothopt': 'smoothopt',
                'istrng': 'istrng', 'radiopt': 'radiopt', 'dprob': 'prbrad', 'iprob': 'iprob',
                'sasopt': 'sasopt', 'arcres': 'arcres',
                # Options for Implicit Membranes
                'membraneopt': 'memopt', 'mprob': 'mprob', 'mthick': 'mthick', 'mctrdz': 'mctrdz',
                'poretype': 'poretype',
                # Options to select numerical procedures
                'npbopt': 'npbopt', 'solvopt': 'solvopt', 'accept': 'accept', 'maxitn': 'linit',
                'fillratio': 'fillratio', 'space': 'scale', 'nbuffer': 'nbuffer', 'nfocus': 'nfocus',
                'fscale': 'fscale', 'npbgrid': 'npbgrid',
                # Options to compute energy and forces
                'bcopt': 'bcopt', 'eneopt': 'eneopt', 'frcopt': 'frcopt', 'scalec': 'scalec',
                'cutfd': 'cutfd', 'cutnb': 'cutnb', 'nsnba': 'nsnba',
                # Options to select a non-polar solvation treatment
                'decompopt': 'decompopt', 'use_rmin': 'use_rmin', 'sprob': 'sprob',
                'vprob': 'vprob', 'rhow_effect': 'rhow_effect', 'use_sav': 'use_sav',
                'cavity_surften': 'cavity_surften', 'cavity_offset': 'cavity_offset',
                'maxsph': 'maxsph', 'maxarcdot': 'maxarcdot',
                # Options for output
                'npbverb': 'npbverb'}

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'ioutfm': 'cntrl',
                       'idecomp': 'cntrl', 'dec_verbose': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ntx': 'cntrl', 'pbtemp': 'pb',
                       # Basic input options
                       'ipb': 'cntrl', 'inp': 'cntrl',
                       # Options to define the physical constants
                       'epsin': 'pb', 'epsout': 'pb', 'epsmem': 'pb', 'smoothopt': 'pb',
                       'istrng': 'pb', 'radiopt': 'pb', 'dprob': 'pb', 'iprob': 'pb',
                       'sasopt': 'pb', 'arcres': 'pb',
                       # Options for Implicit Membranes
                       'membraneopt': 'pb', 'mprob': 'pb', 'mthick': 'pb', 'mctrdz': 'pb',
                       'poretype': 'pb',
                       # Options to select numerical procedures
                       'npbopt': 'pb', 'solvopt': 'pb', 'accept': 'pb', 'maxitn': 'pb',
                       'fillratio': 'pb', 'space': 'pb', 'nbuffer': 'pb', 'nfocus': 'pb',
                       'fscale': 'pb', 'npbgrid': 'pb',
                       # Options to compute energy and forces
                       'bcopt': 'pb', 'eneopt': 'pb', 'frcopt': 'pb', 'scalec': 'pb',
                       'cutfd': 'pb', 'cutnb': 'pb', 'nsnba': 'pb',
                       # Options to select a non-polar solvation treatment
                       'decompopt': 'pb', 'use_rmin': 'pb', 'sprob': 'pb', 'vprob': 'pb',
                       'rhow_effect': 'pb', 'use_sav': 'pb', 'cavity_surften': 'pb',
                       'cavity_offset': 'pb', 'maxsph': 'pb', 'maxarcdot': 'pb',
                       # Options for output
                       'npbverb': 'pb'}
        self.namelist = 'pb'


    # def __init__(self, INPUT):
    #     # We need to change istrng to mM (from M).
    #     SanderInput.__init__(self, INPUT)
    #     self.mdin.change('pb', 'istrng', INPUT['pb']['istrng'] * 1000)


class SanderAPBSInput(SanderInput):
    """ PB sander input file using APBS as the PB solver """
    def __init__(self, INPUT):
        super().__init__(INPUT)
        self.program = 'sander.APBS'
        self.input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999,
                       'imin': 5, 'maxcyc': 1, 'igb': 6,
                       'apbs_print': 1, 'calc_type': 0, 'cmeth': 1,
                       'bcfl': 2, 'srfm': 2, 'chgm': 1, 'pdie': 1.0,
                       'sdie': 80.0, 'srad': 1.4, 'nion': 2, 'ionq': '1.0,-1.0',
                       'ionc': '0.0,0.0', 'ionrr': '2.0,2.0', 'radiopt': 0,
                       'calcforce': 0, 'calcnpenergy': 1, 'grid': '0.5,0.5,0.5',
                       'gamma': 0.00542, 'ioutfm': 0}

        self.name_map = {'inp': 'inp', 'idecomp': 'idecomp', 'pdie': 'indi',
                    'sdie': 'exdi', 'ionc': 'istrng', 'radiopt': 'radiopt',
                    'srad': 'prbrad', 'grid': 'scale', 'gamma': 'cavity_surften',
                    'ioutfm': 'netcdf', 'maxcyc': 'pb_maxcyc'}

        self.parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl',
                           'imin': 'cntrl', 'maxcyc': 'cntrl', 'igb': 'cntrl',
                           'ipb': 'cntrl', 'inp': 'cntrl', 'idecomp': 'cntrl',
                           'ioutfm': 'cntrl', 'apbs_print': 'apbs', 'calc_type': 'apbs',
                           'cmeth': 'apbs', 'bcfl': 'apbs', 'srfm': 'apbs',
                           'chgm': 'apbs', 'pdie': 'apbs', 'sdie': 'apbs',
                           'srad': 'apbs', 'nion': 'apbs', 'ionq': 'apbs',
                           'ionc': 'apbs', 'ionrr': 'apbs', 'radiopt': 'apbs',
                           'calcforce': 'apbs', 'grid': 'apbs', 'gamma': 'apbs',
                           'calcnpenergy': 'apbs'}

        self.mdin.change('apbs', 'grid', f"{INPUT['pb']['scale']},{INPUT['pb']['scale']},{INPUT['pb']['scale']}")
        self.mdin.change('apbs', 'ionc', f"{INPUT['pb']['istrng']},{INPUT['pb']['istrng']}")
        self.mdin.change('apbs', 'gamma', INPUT['pb']['cavity_surften'] * 4.184)
    # def __init__(self, INPUT):
    #     SanderInput.__init__(self, INPUT)
    #     # We also have to make some modifications specific to sander.APBS
    #     self.mdin.change('apbs', 'grid', '%s,%s,%s' % (INPUT['pb']['scale'],
    #                                                    INPUT['pb']['scale'], INPUT['pb']['scale']))
    #     self.mdin.change('apbs', 'ionc', '%s,%s' % (INPUT['pb']['istrng'],
    #                                                 INPUT['pb']['istrng']))
    #     self.mdin.change('apbs', 'gamma', INPUT['pb']['cavity_surften'] * 4.184)


class SanderGBDecomp(SanderGBInput):
    """ GB decomposition input file for sander. In addition to the INPUT dictionary,
       this class also needs several GROUP cards to define the 
   """

    def __init__(self, INPUT, *cards):
        SanderGBInput.__init__(self, INPUT)
        self.make_mdin()
        for card in cards:
            self.mdin.AddCard(card[0], card[1])


class SanderPBDecomp(SanderPBSAInput):
    """ PB decomposition input file for sander """

    def __init__(self, INPUT, *cards):
        SanderPBSAInput.__init__(self, INPUT)
        self.make_mdin()
        for card in cards:
            self.mdin.AddCard(card[0], card[1])


class QuasiHarmonicInput(object):
    """ Write a cpptraj input file to do a quasi-harmonic entropy calculation """

    def __init__(self, com_mask, rec_mask, lig_mask, temperature=298.15, stability=False,
                 prefix='_GMXMMPBSA_', trj_suffix='mdcrd'):
        self.file_string = ("trajin %scomplex.%s\n" % (prefix, trj_suffix) +
                            "reference %savgcomplex.pdb\n" % (prefix) +
                            "rms mass reference %s\n" % (com_mask) +
                            "matrix mwcovar name comp.matrix %s\n" % (com_mask) +
                            "analyze matrix comp.matrix out " +
                            "%scomplex_entropy.out thermo reduce temp %s\n" % (prefix, temperature))
        if not stability:
            self.file_string = (self.file_string +
                                "rms mass reference %s\n" % (rec_mask) +
                                "matrix mwcovar name rec.matrix %s\n" % (rec_mask) +
                                "analyze matrix rec.matrix out " +
                                "%sreceptor_entropy.out thermo reduce temp %s\n" % (prefix, temperature) +
                                "rms mass reference %s\n" % (lig_mask) +
                                "matrix mwcovar name lig.matrix %s\n" % (lig_mask) +
                                "analyze matrix lig.matrix out " +
                                "%sligand_entropy.out thermo reduce temp %s\n" % (prefix, temperature))

    def write_input(self, filename):
        """ Writes the input file """
        open(filename, 'w').write(self.file_string)

