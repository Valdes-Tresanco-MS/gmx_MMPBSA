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

import parmed
from parmed.amber.mdin import Mdin
from parmed.exceptions import AmberError
from copy import deepcopy


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def create_inputs(INPUT, prmtop_system, pre):
    """ Creates the input files for all necessary calculations """
    stability = prmtop_system.stability

    # First check if we are running decomp
    if INPUT['decomprun']:
        # Get the cards that will go in the complex mdin file
        com_card, rc_card, lc_card = prmtop_system.Group(INPUT['print_res'], True)
        junk, rec_card, lig_card = prmtop_system.Group(INPUT['print_res'], False,
                                                       rec_group_type='RES', lig_group_type='RES')
        full_com, full_rc, full_lc = prmtop_system.Group('all', True)
        junk, full_rec, full_lig = prmtop_system.Group('all', False)

        # Convert the strings into card objects
        # Now create the mdin objects
        if INPUT['gbrun']:
            com_arad = None
            rec_arad = None
            lig_arad = None
            if INPUT['alpb']:
                import subprocess
                com_top = parmed.load_file('COM.prmtop', xyz=f"{pre}COM.inpcrd")
                com_top.save('COM.pqr', format='pqr', overwrite=True)
                stdoutdata = subprocess.getoutput("elsize COM.pqr -hea")
                com_arad = stdoutdata.split()[INPUT['arad_method'] * -1]
                if not stability:
                    rec_top = parmed.load_file('REC.prmtop', xyz=f"{pre}REC.inpcrd")
                    rec_top.save('REC.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize REC.pqr -hea")
                    rec_arad = stdoutdata.split()[INPUT['arad_method'] * -1]

                    lig_top = parmed.load_file('LIG.prmtop', xyz=f"{pre}LIG.inpcrd")
                    lig_top.save('LIG.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize LIG.pqr -hea")
                    lig_arad = stdoutdata.split()[INPUT['arad_method'] * -1]

            rec_res = ['Residues considered as REC', full_rc]
            if stability:
                pri_res = ['Residues to print', com_card]
                com_input = deepcopy(INPUT)
                if com_arad:
                    com_input['arad'] = com_arad
                com_mdin = SanderGBDecomp(com_input, rec_res, pri_res)
                com_mdin.write_input(f"{pre}gb_decomp_com.mdin")
            else:
                lig_res = ['Residues considered as LIG', full_lc]
                pri_res = ['Residues to print', com_card]
                com_input = deepcopy(INPUT)
                if com_arad:
                    com_input['arad'] = com_arad
                com_mdin = SanderGBDecomp(com_input, rec_res, lig_res, pri_res)
                rec_res = ['Residues considered as REC', full_rec]
                pri_res = ['Residues to print', rec_card]
                rec_input = deepcopy(INPUT)
                if rec_arad:
                    rec_input['arad'] = rec_arad
                rec_mdin = SanderGBDecomp(rec_input, rec_res, pri_res)
                lig_res = ['Residues considered as LIG', full_lig]
                pri_res = ['Residues to print', lig_card]
                lig_input = deepcopy(INPUT)
                if lig_arad:
                    lig_input['arad'] = lig_arad
                lig_mdin = SanderGBDecomp(lig_input, lig_res, pri_res)
                com_mdin.write_input(f"{pre}gb_decomp_com.mdin")
                rec_mdin.write_input(f"{pre}gb_decomp_rec.mdin")
                lig_mdin.write_input(f"{pre}gb_decomp_lig.mdin")

        if INPUT['pbrun']:
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

    else:  # not decomp

        if INPUT['gbrun']:
            # We need separate input files for QM/gmx_MMPBSA
            com_arad = None
            rec_arad = None
            lig_arad = None
            if INPUT['alpb']:
                import subprocess
                com_top = parmed.load_file('COM.prmtop', xyz=f"{pre}COM.inpcrd")
                com_top.save('COM.pqr', format='pqr', overwrite=True)
                stdoutdata = subprocess.getoutput("elsize COM.pqr -hea")
                com_arad = stdoutdata.split()[INPUT['arad_method'] * -1]
                if not stability:
                    rec_top = parmed.load_file('REC.prmtop', xyz=f"{pre}REC.inpcrd")
                    rec_top.save('REC.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize REC.pqr -hea")
                    rec_arad = stdoutdata.split()[INPUT['arad_method'] * -1]

                    lig_top = parmed.load_file('LIG.prmtop', xyz=f"{pre}LIG.inpcrd")
                    lig_top.save('LIG.pqr', format='pqr', overwrite=True)
                    stdoutdata = subprocess.getoutput("elsize LIG.pqr -hea")
                    lig_arad = stdoutdata.split()[INPUT['arad_method'] * -1]

            if INPUT['ifqnt']:
                com_input = deepcopy(INPUT)
                rec_input = deepcopy(INPUT)
                lig_input = deepcopy(INPUT)
                (com_input['qmmask'], rec_input['qmmask'],
                 lig_input['qmmask']) = prmtop_system.Mask(INPUT['qm_residues'], in_complex=False)
                if not com_input['qmmask']:
                    raise AmberError('No valid QM residues chosen!')
                com_input['qm_theory'] = "'%s'" % com_input['qm_theory']
                com_input['qmmask'] = "'%s'" % com_input['qmmask']
                com_input['qmcharge'] = com_input['qmcharge_com']
                # check if alpb
                if com_arad:
                    com_input['arad'] = com_arad
                gb_mdin = SanderGBInput(com_input)
                gb_mdin.write_input(f'{pre}gb_qmmm_com.mdin')
                if not stability:
                    if not rec_input['qmmask']:
                        rec_input['ifqnt'] = 0
                    else:
                        rec_input['qmmask'] = "'%s'" % rec_input['qmmask']
                    rec_input['qm_theory'] = "'%s'" % rec_input['qm_theory']
                    rec_input['qmcharge'] = rec_input['qmcharge_rec']
                    # check if alpb
                    if rec_arad:
                        rec_input['arad'] = rec_arad
                    gb_mdin = SanderGBInput(rec_input)
                    gb_mdin.write_input(f'{pre}gb_qmmm_rec.mdin')
                    if not lig_input['qmmask']:
                        lig_input['ifqnt'] = 0
                    else:
                        lig_input['qmmask'] = "'%s'" % lig_input['qmmask']
                    lig_input['qm_theory'] = "'%s'" % lig_input['qm_theory']
                    lig_input['qmcharge'] = lig_input['qmcharge_lig']
                    # check if alpb
                    if lig_arad:
                        lig_input['arad'] = lig_arad
                    gb_mdin = SanderGBInput(lig_input)
                    gb_mdin.write_input(f'{pre}gb_qmmm_lig.mdin')

            elif INPUT['alpb']:
                com_input = deepcopy(INPUT)
                rec_input = deepcopy(INPUT)
                lig_input = deepcopy(INPUT)
                com_input['arad'] = com_arad
                gb_mdin = SanderGBInput(com_input)
                gb_mdin.write_input(f'{pre}gb_com.mdin')
                rec_input['arad'] = rec_arad
                gb_mdin = SanderGBInput(rec_input)
                gb_mdin.write_input(f'{pre}gb_rec.mdin')
                lig_input['arad'] = lig_arad
                gb_mdin = SanderGBInput(lig_input)
                gb_mdin.write_input(f'{pre}gb_lig.mdin')
            else:
                gb_mdin = SanderGBInput(INPUT)
                gb_mdin.write_input(f'{pre}gb.mdin')

        if INPUT['pbrun']:
            pb_prog = 'sander.APBS' if INPUT['sander_apbs'] else 'sander'
            if pb_prog == 'sander.APBS':
                pb_mdin = SanderAPBSInput(INPUT)
                pb_mdin2 = SanderAPBSInput(INPUT)
            else:
                pb_mdin = SanderPBSAInput(INPUT)
                pb_mdin2 = SanderPBSA2Input(INPUT)

            pb_mdin.write_input(f'{pre}pb.mdin')
            pb_mdin2.write_input(f'{pre}pb.mdin2')

        if INPUT['rismrun']:
            rism_mdin = SanderRISMInput(INPUT)
            rism_mdin.write_input(pre + 'rism.mdin')

    # end if decomprun

    if INPUT['qh_entropy']:  # quasi-harmonic approximation input file
        trj_suffix = 'nc' if INPUT['netcdf'] else 'mdcrd'
        com_mask, rec_mask, lig_mask = prmtop_system.Mask('all', True)
        if not INPUT['mutant_only']:
            qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask, temperature=INPUT['temperature'],
                                       stability=stability, prefix=pre, trj_suffix=trj_suffix)
            qh_in.write_input(f'{pre}cpptrajentropy.in')
        if INPUT['alarun']:
            qh_in = QuasiHarmonicInput(com_mask, rec_mask, lig_mask, temperature=INPUT['temperature'],
                                       stability=stability, prefix=pre + 'mutant_', trj_suffix=trj_suffix)
            qh_in.write_input(f'{pre}mutant_cpptrajentropy.in')


class SanderInput(object):
    """ Base class sander input file """
    program = 'sander'  # This runs with sander
    input_items = {'foo': 'bar'}  # replace this in derived classes
    name_map = {'foo': 'orig'}  # replace this in derived classes
    parent_namelist = {'foo': 'foo_namelist'}  # replace this in derived classes

    def __init__(self, INPUT):
        self.mdin = Mdin(self.program)
        self.mdin.title = 'File generated by gmx_MMPBSA'
        for key in self.input_items:
            # Skip ioutfm since it is handled explicitly later
            if key == 'ioutfm':
                continue
            try:
                self.mdin.change(self.parent_namelist[key], key, INPUT[self.name_map[key]])
            except KeyError:
                self.mdin.change(self.parent_namelist[key], key, self.input_items[key])

        self.mdin.change('cntrl', 'ioutfm', int(bool(INPUT['netcdf'])))

    def write_input(self, filename):
        """ Write the mdin file """
        self.mdin.write(filename)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderGBInput(SanderInput):
    """ GB sander input file """
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'idecomp': 0, 'offset': -999999.0,
                   'imin': 5, 'maxcyc': 1, 'ncyc': 0, 'gbsa': 0, 'ioutfm': 0, 'dec_verbose': 0,
                   # Basic options
                   'igb': 5, 'intdiel': 1.0, 'extdiel': 78.5, 'saltcon': 0.0, 'surften': 0.0072,
                   # QM options
                   'ifqnt': 0, 'qmmask': '', 'qm_theory': '', 'qmcharge': 0,
                   'qmgb': 2, 'qmcut': 999.0,
                   'scfconv': 1.0e-8, 'peptide_corr': 0, 'writepdb': 1, 'verbosity': 0, 'alpb': 0, 'arad': 15}

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'idecomp': 'cntrl', 'offset': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'ncyc': 'cntrl', 'gbsa': 'cntrl', 'ioutfm': 'cntrl',
                       'dec_verbose': 'cntrl',
                       # Basic options
                       'igb': 'cntrl', 'intdiel': 'cntrl', 'extdiel': 'cntrl', 'saltcon': 'cntrl',
                       'surften': 'cntrl', 'alpb': 'cntrl', 'arad': 'cntrl',
                       # QM options
                       'ifqnt': 'cntrl', 'qmmask': 'qmmm', 'qm_theory': 'qmmm', 'qmcharge': 'qmmm',
                       'qmgb': 'qmmm', 'qmcut': 'qmmm',
                       'scfconv': 'qmmm', 'peptide_corr': 'qmmm', 'writepdb': 'qmmm', 'verbosity': 'qmmm'}

    name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'idecomp': 'idecomp', 'offset': 'offset',
                'imin': 'imin', 'gbsa': 'gbsa', 'ioutfm': 'netcdf', 'dec_verbose': 'dec_verbose',
                'maxcyc': 'gb_maxcyc', 'ncyc': 'ncyc',
                # Basic options
                'igb': 'igb', 'intdiel': 'intdiel', 'extdiel': 'extdiel', 'saltcon': 'saltcon',
                'surften': 'surften', 'alpb': 'alpb', 'arad': 'arad',
                # QM options
                'ifqnt': 'ifqnt',  'qmmask': 'qmmask', 'qm_theory': 'qm_theory',
                'qmcharge': 'qmcharge', 'qmgb': 'qmgb', 'qmcut': 'qmcut',
                'scfconv': 'scfconv', 'peptide_corr': 'peptide_corr', 'writepdb': 'writepdb', 'verbosity': 'verbosity'}


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderRISMInput(SanderInput):
    """ RISM sander input file """
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'idecomp': 0, 'offset': -999999.0,
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

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'idecomp': 'cntrl', 'offset': 'cntrl',
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

    name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'idecomp': 'idecomp', 'offset': 'offset',
                'imin': 'imin', 'gbsa': 'gbsa', 'ioutfm': 'netcdf', 'dec_verbose': 'dec_verbose',
                'maxcyc': 'rism_maxcyc', 'ncyc': 'ncyc',
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


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBSADECOMPInput(SanderInput):
    """ PB sander input file """
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'ioutfm': 0, 'idecomp': 0, 'dec_verbose': 0,
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

    name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'ioutfm': 'netcdf',
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

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'ioutfm': 'cntrl',
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

    def __init__(self, INPUT):
        # We need to change istrng to mM (from M).
        SanderInput.__init__(self, INPUT)
        self.mdin.change('pb', 'istrng', INPUT['istrng'] * 1000)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBSAInput(SanderInput):
    """ PB sander input file """
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'ioutfm': 0,
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

    name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'ioutfm': 'netcdf',
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

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'ioutfm': 'cntrl',
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

    def __init__(self, INPUT):
        # We need to change istrng to mM (from M).
        SanderInput.__init__(self, INPUT)
        self.mdin.change('pb', 'istrng', INPUT['istrng'] * 1000)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBSA2Input(SanderInput):
    """ PB sander input file """
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999, 'ioutfm': 0,
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

    name_map = {'ntb': 'ntb', 'cut': 'cut', 'nsnb': 'nsnb', 'ioutfm': 'netcdf',
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

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl', 'ioutfm': 'cntrl',
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

    def __init__(self, INPUT):
        # We need to change istrng to mM (from M).
        SanderInput.__init__(self, INPUT)
        self.mdin.change('pb', 'istrng', INPUT['istrng'] * 1000)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderAPBSInput(SanderInput):
    """ PB sander input file using APBS as the PB solver """
    program = 'sander.APBS'
    input_items = {'ntb': 0, 'cut': 999.0, 'nsnb': 99999,
                   'imin': 5, 'maxcyc': 1, 'igb': 6,
                   'apbs_print': 1, 'calc_type': 0, 'cmeth': 1,
                   'bcfl': 2, 'srfm': 2, 'chgm': 1, 'pdie': 1.0,
                   'sdie': 80.0, 'srad': 1.4, 'nion': 2, 'ionq': '1.0,-1.0',
                   'ionc': '0.0,0.0', 'ionrr': '2.0,2.0', 'radiopt': 0,
                   'calcforce': 0, 'calcnpenergy': 1, 'grid': '0.5,0.5,0.5',
                   'gamma': 0.00542, 'ioutfm': 0}

    name_map = {'inp': 'inp', 'idecomp': 'idecomp', 'pdie': 'indi',
                'sdie': 'exdi', 'ionc': 'istrng', 'radiopt': 'radiopt',
                'srad': 'prbrad', 'grid': 'scale', 'gamma': 'cavity_surften',
                'ioutfm': 'netcdf', 'maxcyc': 'pb_maxcyc'}

    parent_namelist = {'ntb': 'cntrl', 'cut': 'cntrl', 'nsnb': 'cntrl',
                       'imin': 'cntrl', 'maxcyc': 'cntrl', 'igb': 'cntrl',
                       'ipb': 'cntrl', 'inp': 'cntrl', 'idecomp': 'cntrl',
                       'ioutfm': 'cntrl', 'apbs_print': 'apbs', 'calc_type': 'apbs',
                       'cmeth': 'apbs', 'bcfl': 'apbs', 'srfm': 'apbs',
                       'chgm': 'apbs', 'pdie': 'apbs', 'sdie': 'apbs',
                       'srad': 'apbs', 'nion': 'apbs', 'ionq': 'apbs',
                       'ionc': 'apbs', 'ionrr': 'apbs', 'radiopt': 'apbs',
                       'calcforce': 'apbs', 'grid': 'apbs', 'gamma': 'apbs',
                       'calcnpenergy': 'apbs'}

    def __init__(self, INPUT):
        SanderInput.__init__(self, INPUT)
        # We also have to make some modifications specific to sander.APBS
        self.mdin.change('apbs', 'grid', '%s,%s,%s' % (INPUT['scale'],
                                                       INPUT['scale'], INPUT['scale']))
        self.mdin.change('apbs', 'ionc', '%s,%s' % (INPUT['istrng'],
                                                    INPUT['istrng']))
        self.mdin.change('apbs', 'gamma', INPUT['cavity_surften'] * 4.184)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderGBDecomp(SanderGBInput):
    """ GB decomposition input file for sander. In addition to the INPUT dictionary,
       this class also needs several GROUP cards to define the 
   """

    def __init__(self, INPUT, *cards):
        SanderGBInput.__init__(self, INPUT)
        for card in cards:
            self.mdin.AddCard(card[0], card[1])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SanderPBDecomp(SanderPBSADECOMPInput):
    """ PB decomposition input file for sander """

    def __init__(self, INPUT, *cards):
        SanderPBSADECOMPInput.__init__(self, INPUT)
        for card in cards:
            self.mdin.AddCard(card[0], card[1])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

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

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
