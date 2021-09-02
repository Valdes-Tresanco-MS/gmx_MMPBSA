"""
Generate Amber topology files from GROMACS files
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

import os
import parmed
from GMXMMPBSA.exceptions import *
from GMXMMPBSA.utils import checkff, selector, get_dist, list2range, Residue
from GMXMMPBSA.alamdcrd import _scaledistance
import subprocess
from pathlib import Path
import logging
import string
from copy import deepcopy
from parmed.tools.changeradii import ChRad

chains_letters = list(string.ascii_uppercase)
his = ['HIS', 'HIE', 'HID', 'HIP']
cys_name = ['CYS', 'CYX', 'CYM']
lys = ['LYS', 'LYN']
asp = ['ASP', 'ASH']
glu = ['GLU', 'GLH']
positive_aa = ['LYS', 'ARG', 'HIP']
negative_aa = ['GLU', 'ASP']
nonpolar_aa = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
polar_aa = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']

PBRadii = {1: 'bondi', 2: 'mbondi', 3: 'mbondi2', 4: 'mbondi3'}

ions_para_files = {1: 'frcmod.ions234lm_126_tip3p', 2: 'frcmod.ions234lm_iod_tip4pew', 3: 'frcmod.ions234lm_iod_spce',
                   4: 'frcmod.ions234lm_hfe_spce', 5: 'frcmod.ions234lm_126_tip4pew', 6: 'frcmod.ions234lm_126_spce',
                   7: 'frcmod.ions234lm_1264_tip4pew', 8: 'frcmod.ions234lm_1264_tip3p',
                   9: 'frcmod.ions234lm_1264_spce', 10: 'frcmod.ions234lm_iod_tip3p',
                   11: 'frcmod.ions234lm_hfe_tip4pew', 12: 'frcmod.ions234lm_hfe_tip3p}'}

ions = ["AG", "AL", "Ag", "BA", "BR", "Be", "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "CU1", "Ce", "Cl-", "Cr",
        "Dy", "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", "HE+", "HG", "HZ+", "Hf", "IN", "IOD", "K", "K+",
        "LA", "LI", "LU", "MG", "MN", "NA", "NH4", "NI", "Na+", "Nd", "PB", "PD", "PR", "PT", "Pu", "RB", "Ra", "SM",
        "SR", "Sm", "Sn", "TB", "TL", "Th", "Tl", "Tm", "U4+", "V2+", "Y", "YB2", "ZN", "Zr"]


class CheckMakeTop:
    def __init__(self, FILES, INPUT, external_programs):
        self.FILES = FILES
        self.INPUT = INPUT
        self.external_progs = external_programs
        self.use_temp = False
        self.log = open('gmx_MMPBSA.log', 'a')

        self.mut_label = None

        # Define Gromacs executable
        self.make_ndx = self.external_progs['make_ndx']
        self.trjconv = self.external_progs['trjconv']
        self.editconf = self.external_progs['editconf']

        self.ref_str = None

        self.ligand_tpr = None
        self.ligand_mol2 = None

        self.rec_str_ions = False
        self.lig_str_ions = False

        # create the * prmtop variables for compatibility with the original code
        self.complex_pmrtop = 'COM.prmtop'
        self.receptor_pmrtop = 'REC.prmtop'
        self.ligand_pmrtop = 'LIG.prmtop'

        self.mutant_complex_pmrtop = 'MUT_COM.prmtop'
        self.mutant_receptor_pmrtop = 'MUT_REC.prmtop'
        self.mutant_ligand_pmrtop = 'MUT_LIG.prmtop'

        self.complex_temp_top = self.FILES.prefix + 'COM.top'
        self.receptor_temp_top = self.FILES.prefix + 'COM.top'
        self.ligand_temp_top = self.FILES.prefix + 'COM.top'

        self.complex_str_file = self.FILES.prefix + 'COM.pdb'
        self.receptor_str_file = self.FILES.prefix + 'REC.pdb'
        self.ligand_str_file = self.FILES.prefix + 'LIG.pdb'

        if self.FILES.reference_structure:
            self.ref_str = parmed.read_PDB(self.FILES.reference_structure)

        checkff()
        self.checkFiles()

    def checkFiles(self):
        if (not self.FILES.complex_tpr or not self.FILES.complex_index or not self.FILES.complex_trajs or not
        self.FILES.complex_groups):
            GMXMMPBSA_ERROR('You must define the structure, topology and index files, as well as the groups!')

    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """
        self.gmx2pdb()
        if self.FILES.complex_top:
            tops = self.gmxtop2prmtop()
        else:
            self.pdb2prmtop()
            tops = self.makeToptleap()

        if self.INPUT['decomprun']:
            self.INPUT['print_res'] = ','.join(str(x) for x in self.get_selected_residues(self.INPUT['print_res']))
        if self.INPUT['qm_residues']:
            self.INPUT['qm_residues'] = ','.join(str(x) for x in self.get_selected_residues(self.INPUT['qm_residues']))

        self.cleanup_trajs()
        return tops

    def gmx2pdb(self):
        """
        Generate PDB file to generate topology
        :return:
        """

        logging.info('Get PDB files from GROMACS structures files...')

        # wt complex
        # make index for extract pdb structure
        rec_group, lig_group = self.FILES.complex_groups

        logging.info('Making gmx_MMPBSA index for complex...')
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        make_ndx_echo_args = ['echo', 'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | '
                                      '{l}\n q\n'.format(r=rec_group, l=lig_group)]
        c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

        com_ndx = self.FILES.prefix + 'COM_index.ndx'
        make_ndx_args = self.make_ndx + ['-n', self.FILES.complex_index, '-o', com_ndx]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + (' '.join(make_ndx_echo_args).replace('\n', '\\n')) + ' | ' +
                         ' '.join(make_ndx_args))
        c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=self.log, stderr=self.log)
        if c2.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.complex_index))
        self.FILES.complex_index = com_ndx

        logging.info('Normal Complex: Saving group {}_{} in {} file as '
                     '{}'.format(rec_group, lig_group, self.FILES.complex_index, self.complex_str_file))
        editconf_echo_args = ['echo', 'GMXMMPBSA_REC_GMXMMPBSA_LIG']
        c3 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
        # we extract a pdb from structure file to make amber topology
        editconf_args = self.editconf + ['-f', self.FILES.complex_tpr, '-o', self.complex_str_file, '-n',
                                         self.FILES.complex_index]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
        c4 = subprocess.Popen(editconf_args, stdin=c3.stdout, stdout=self.log, stderr=self.log)
        if c4.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.complex_tpr))

        # Put receptor and ligand (explicitly defined) to avoid overwrite them
        # check if ligand is not protein. In any case, non-protein ligand always most be processed
        if self.FILES.ligand_mol2:
            logging.info(f'Generating ligand parameters from {self.FILES.ligand_mol2} file...')
            lig_name = os.path.splitext(os.path.split(self.FILES.ligand_mol2)[1])[0]
            self.ligand_frcmod = self.FILES.prefix + lig_name + '.frcmod'
            # run parmchk2
            parmchk2 = self.external_progs['parmchk2']
            lig_ff = '1'
            if "gaff2" in self.INPUT['forcefields']:
                lig_ff = '2'
            parmchk2_args = [parmchk2, '-i', self.FILES.ligand_mol2, '-f', 'mol2', '-o', self.ligand_frcmod, '-s',
                             lig_ff]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + ' '.join(parmchk2_args))
            l3 = subprocess.Popen(parmchk2_args, stdout=self.log, stderr=self.log)
            if l3.wait():
                GMXMMPBSA_ERROR('%s failed when querying %s' % (parmchk2, self.FILES.ligand_mol2))

        # check if the ligand force field is gaff or gaff2 and get if the ligand mol2 was defined
        else:
            if ("leaprc.gaff" in self.INPUT['forcefields'] or "leaprc.gaff2" in self.INPUT['forcefields'] and not
            self.FILES.complex_top):
                GMXMMPBSA_WARNING('You must define the ligand mol2 file (-lm) if the ligand forcefield is '
                                  '"leaprc.gaff" or "leaprc.gaff2". If the ligand is parametrized in Amber force '
                                  'fields ignore this warning')

        # make a temp receptor pdb (even when stability) if decomp to get correct receptor residues from complex. This
        # avoid get multiples molecules from complex.split()
        if self.INPUT['decomprun'] and self.FILES.stability:
            self.use_temp = True
            logging.warning('When decomp is defined, we generate a receptor file in order to extract interface '
                            'residues')
            rec_echo_args = ['echo', '{}'.format(rec_group)]
            cp1 = subprocess.Popen(rec_echo_args, stdout=subprocess.PIPE)
            # we extract a pdb from structure file to make amber topology
            editconf_args = self.editconf + ['-f', self.FILES.complex_tpr, '-o', 'rec_temp.pdb', '-n',
                                        self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(rec_echo_args)) + ' | ' + ' '.join(editconf_args))
            cp2 = subprocess.Popen(editconf_args, stdin=cp1.stdout, stdout=self.log, stderr=self.log)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.complex_tpr))

        # check if stability
        if self.FILES.stability:
            if (self.FILES.receptor_tpr or self.FILES.ligand_tpr):
                logging.warning('When Stability calculation mode is selected, receptor and ligand files are not '
                                'needed...')
            # if self.INPUT['alarun'] and (self.FILES.mutant_receptor_tpr or self.FILES.mutant_ligand_tpr):
            #     logging.warning('When Stability calculation mode is selected, mutant receptor/mutant ligand files '
            #                     'are not needed...')

        # wt receptor
        if self.FILES.receptor_tpr:
            logging.info('A receptor structure file was defined. Using MT approach...')
            logging.info('Normal receptor: Saving group {} in {} file as {}'.format(
                self.FILES.receptor_group, self.FILES.receptor_index, self.receptor_str_file))
            editconf_echo_args = ['echo', '{}'.format(self.FILES.receptor_group)]
            p1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we extract a pdb from structure file to make amber topology
            editconf_args = self.editconf + ['-f', self.FILES.receptor_tpr, '-o', self.receptor_str_file, '-n',
                                             self.FILES.receptor_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cp2 = subprocess.Popen(editconf_args, stdin=p1.stdout, stdout=self.log, stderr=self.log)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.receptor_tpr))

        else:
            logging.info('No receptor structure file was defined. Using ST approach...')
            logging.info('Using receptor structure from complex to generate AMBER topology')
            # wt complex receptor
            logging.info('Normal Complex: Saving group {} in {} file as {}'.format(
                rec_group, self.FILES.complex_index, self.receptor_str_file))
            editconf_echo_args = ['echo', '{}'.format(rec_group)]
            cp1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we extract a pdb from structure file to make amber topology
            editconf_args = self.editconf + ['-f', self.FILES.complex_tpr, '-o', self.receptor_str_file, '-n',
                                             self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cp2 = subprocess.Popen(editconf_args, stdin=cp1.stdout, stdout=self.log, stderr=self.log)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.complex_tpr))

        # ligand
        # # check consistence
        if self.FILES.ligand_tpr:  # ligand is protein
            # FIXME: if ligand is a zwitterionic aa fail
            logging.info('A ligand structure file was defined. Using MT approach...')
            logging.info('Normal Ligand: Saving group {} in {} file as {}'.format(
                self.FILES.ligand_group, self.FILES.ligand_index, self.ligand_str_file))
            # wt ligand
            editconf_echo_args = ['echo', '{}'.format(self.FILES.ligand_group)]
            l1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we extract a pdb from structure file to make amber topology
            editconf_args = self.editconf + ['-f', self.FILES.ligand_tpr, '-o', self.ligand_str_file, '-n',
                                             self.FILES.ligand_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            l2 = subprocess.Popen(editconf_args, stdin=l1.stdout, stdout=self.log, stderr=self.log)
            if l2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.ligand_tpr))
        else:
            # wt complex ligand
            logging.info('No ligand structure file was defined. Using ST approach...')
            logging.info('Using ligand structure from complex to generate AMBER topology')
            logging.info('Normal ligand: Saving group {} in {} file as {}'.format(lig_group, self.FILES.complex_index,
                                                                                  self.ligand_str_file))
            editconf_echo_args = ['echo', '{}'.format(lig_group)]
            cl1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we extract a pdb from structure file to make amber topology
            editconf_args = self.editconf + ['-f', self.FILES.complex_tpr, '-o', self.ligand_str_file, '-n',
                                             self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cl2 = subprocess.Popen(editconf_args, stdin=cl1.stdout, stdout=self.log, stderr=self.log)
            if cl2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.editconf), self.FILES.complex_tpr))

        # check for IE variable
        if self.FILES.receptor_tpr or self.FILES.ligand_tpr:
            if self.INPUT['interaction_entropy'] or self.INPUT['c2_entropy']:
                GMXMMPBSA_WARNING("The IE or C2 entropy method don't support the MTP approach...")
                self.INPUT['interaction_entropy'] = self.INPUT['c2_entropy'] = 0

        # initialize receptor and ligand structures. Needed to get residues map
        self.complex_str = self.molstr(self.complex_str_file)
        self.receptor_str = self.molstr(self.receptor_str_file)
        self.ligand_str = self.molstr(self.ligand_str_file)
        self.check4water()
        self.resi, self.resl, self.orderl, self.indexes = self.res2map()
        self.fix_chains_IDs(self.complex_str, self.receptor_str, self.ligand_str, self.ref_str)

    def check4water(self):
        counter = sum(
            res.name
            in [
                'SOD',
                'CLA',
                'TIP3P',
                'TIP4P',
                'TIPS3P',
                'TIP5P',
                'SPC',
                'SPC/E',
                'SPCE',
                'TIP3o',
                'WAT',
            ]
            for res in self.complex_str
        )

        if counter:
            GMXMMPBSA_ERROR(f'gmx_MMPBSA does not support water molecules in any structure, but we found {counter} '
                            f'molecules in the complex.')

    def gmxtop2prmtop(self):
        logging.info('Building Normal Complex Amber Topology...')
        com_top = self.cleantop(self.FILES.complex_top, self.indexes[0] + self.indexes[1])
        # parmed.gromacs.GromacsTopologyFile(self.complex_temp_top,xyz=self.complex_str_file)
        com_top.coordinates = self.complex_str.coordinates
        # try:
        if com_top.impropers or com_top.urey_bradleys or com_top.cmaps:
            com_amb_prm = parmed.amber.ChamberParm.from_structure(com_top)
            com_top_parm = 'chamber'
        else:
            com_amb_prm = parmed.amber.AmberParm.from_structure(com_top)
            com_top_parm = 'amber'

        # IMPORTANT: make_trajs ends in error if the box is defined
        com_amb_prm.box = None

        self.fixparm2amber(com_amb_prm)

        logging.info('Writing Normal Complex Amber Topology...')
        # change de PBRadii
        action = ChRad(com_amb_prm, PBRadii[self.INPUT['PBRadii']])
        com_amb_prm.write_parm(self.complex_pmrtop)

        rec_indexes_string = ','.join(self.resi['REC']['string'])

        rec_hastop = True
        if self.FILES.receptor_top:
            logging.info('A Receptor topology file was defined. Using MT approach...')
            logging.info('Building AMBER Receptor Topology from GROMACS Receptor Topology...')
            rec_top = self.cleantop(self.FILES.receptor_top, self.indexes[0])
            rec_top.coordinates = self.receptor_str.coordinates
            if rec_top.impropers or rec_top.urey_bradleys or rec_top.cmaps:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is AMBER type while the '
                                    'Receptor is CHAMBER type!')
                rec_amb_prm = parmed.amber.ChamberParm.from_structure(rec_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Receptor is AMBER type!')
                rec_amb_prm = parmed.amber.AmberParm.from_structure(rec_top)
            self.fixparm2amber(rec_amb_prm)
        else:
            logging.info('No Receptor topology files was defined. Using ST approach...')
            logging.info('Building AMBER Receptor Topology from Complex...')
            # we make a copy for receptor topology
            rec_amb_prm = self.molstr(com_amb_prm)
            rec_amb_prm.strip(f'!:{rec_indexes_string}')
            rec_hastop = False
        # change de PBRadii
        action = ChRad(rec_amb_prm, PBRadii[self.INPUT['PBRadii']])
        rec_amb_prm.write_parm(self.receptor_pmrtop)

        lig_hastop = True
        if self.FILES.ligand_top:
            logging.info('A Ligand Topology file was defined. Using MT approach...')
            logging.info('Building AMBER Ligand Topology from GROMACS Ligand Topology...')
            lig_top = self.cleantop(self.FILES.ligand_top, self.indexes[1])
            lig_top.coordinates = self.ligand_str.coordinates
            if lig_top.impropers or lig_top.urey_bradleys or lig_top.cmaps:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is AMBER type while the '
                                    'Ligand is CHAMBER type!')
                lig_amb_prm = parmed.amber.ChamberParm.from_structure(lig_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Ligand is AMBER type!')
                lig_amb_prm = parmed.amber.AmberParm.from_structure(lig_top)
            self.fixparm2amber(lig_amb_prm)
        else:
            logging.info('No Ligand Topology files was defined. Using ST approach...')
            logging.info('Building AMBER Ligand Topology from Complex...')
            # we make a copy for ligand topology
            lig_amb_prm = self.molstr(com_amb_prm)
            lig_amb_prm.strip(f':{rec_indexes_string}')
            lig_hastop = False

        # change de PBRadii
        action = ChRad(lig_amb_prm, PBRadii[self.INPUT['PBRadii']])
        lig_amb_prm.write_parm(self.ligand_pmrtop)

        if self.INPUT['alarun']:
            logging.info('Building Mutant Complex Topology...')
            # get mutation index in complex
            self.com_mut_index, self.part_mut, self.part_index, self.mut_label = self.getMutationInfo()
            mut_com_amb_prm = self.makeMutTop(com_amb_prm, self.com_mut_index)
            # change de PBRadii
            action = ChRad(mut_com_amb_prm, PBRadii[self.INPUT['PBRadii']])
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)

            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor Topology...')
                out_prmtop = self.mutant_receptor_pmrtop
                self.mutant_ligand_pmrtop = None
                if rec_hastop:
                    mtop = self.makeMutTop(rec_amb_prm, self.part_index)
                else:
                    mut_com_amb_prm.strip(f'!:{rec_indexes_string}')
                    mtop = mut_com_amb_prm
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand Topology...')
                out_prmtop = self.mutant_ligand_pmrtop
                self.mutant_receptor_pmrtop = None
                if lig_hastop:
                    mtop = self.makeMutTop(lig_amb_prm, self.part_index)
                else:
                    mut_com_amb_prm.strip(f':{rec_indexes_string}')
                    mtop = mut_com_amb_prm

            if com_top_parm == 'charmm':
                mut_prot_amb_prm = parmed.amber.ChamberParm.from_structure(mtop)
            else:
                mut_prot_amb_prm = parmed.amber.AmberParm.from_structure(mtop)
            # change de PBRadii
            action = ChRad(mut_prot_amb_prm, PBRadii[self.INPUT['PBRadii']])
            mut_prot_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

    def _split_str(self, start, r, c, basename, struct, mut_index=0):
        end = start + (r[1] - r[0])
        mask = f'!:{start}-{end}'
        # start += end
        str_ = self.molstr(struct)
        if mut_index:
            str_ = self.makeMutTop(str_, mut_index, True)
        str_.strip(mask)
        str_file = self.FILES.prefix + f'{basename}_F{c}.pdb'
        str_.save(str_file, 'pdb', True, renumber=False)
        return end, str_file

    def pdb2prmtop(self):
        """
        Generate parmed structure object for complex, receptor and ligand ( if it is protein-like)
        :return:
        """
        logging.info('Generating AMBER Compatible PDB Files...')

        # fix receptor and structures
        self.fixparm2amber(self.complex_str, removeH=True)
        self.fixparm2amber(self.receptor_str, removeH=True)
        self.fixparm2amber(self.ligand_str, removeH=True)

        self.receptor_list = {}
        start = 1
        for c, r in enumerate(self.resi['REC']['num'], start=1):
            end, sfile = self._split_str(start, r, c, 'REC', self.receptor_str)
            self.receptor_list[f'REC{c}'] = sfile
            start += end

        self.ligand_list = {}
        start = 1
        for c, r in enumerate(self.resi['LIG']['num'], start=1):
            end, sfile = self._split_str(start, r, c, 'LIG', self.ligand_str)
            self.ligand_list[f'LIG{c}'] = sfile
            start += end

        self.mut_receptor_list = {}
        self.mut_ligand_list = {}
        if self.INPUT['alarun']:
            self.com_mut_index, self.part_mut, self.part_index, self.mut_label = self.getMutationInfo()
            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor Structure...')
                self.mutant_ligand_pmrtop = None
                start = 1
                for c, r in enumerate(self.resi['REC']['num']):
                    end, sfile = self._split_str(start, r, c, f'MUT_REC', self.receptor_str, self.part_index)
                    self.mut_receptor_list[f'MREC{c}'] = sfile
                    start += end
            else:
                logging.info('Detecting mutation in Ligand.Building Mutant Ligand Structure...')
                self.mutant_receptor_pmrtop = None
                start = 1
                for c, r in enumerate(self.resi['LIG']['num']):
                    end, sfile = self._split_str(start, r, c, f'MUT_LIG', self.ligand_str, self.part_index)
                    self.mut_ligand_list[f'MLIG{c}'] =  sfile
                    start += end

    def reswithin(self):
        # Get residue form receptor-ligand interface
        if self.INPUT['decomprun']:
            if self.INPUT['print_res'] == 'all':
                return
            else:
                dist, res_selection = selector(self.INPUT['print_res'])
                res_list = []

                if dist:
                    for i in self.resl['REC']:
                        for j in self.resl['LIG']:
                            for rat in self.complex_str.residues[i - 1].atoms:
                                rat_coor = [rat.xx, rat.xy, rat.xz]
                                for lat in self.complex_str.residues[j - 1].atoms:
                                    lat_coor = [lat.xx, lat.xy, lat.xz]
                                    if get_dist(rat_coor, lat_coor) <= dist:
                                        if i not in res_list:
                                            res_list.append(i)
                                        if j not in res_list:
                                            res_list.append(j)
                                        break
                elif res_selection:
                    for i in self.resl['REC']:
                        rres = self.complex_str.residues[i - 1]
                        if [rres.chain, rres.number, rres.insertion_code] in res_selection:
                            res_list.append(i)
                            res_selection.remove([rres.chain, rres.number, rres.insertion_code])
                    for j in self.resl['LIG']:
                        lres = self.complex_str.residues[j - 1]
                        if [lres.chain, lres.number, lres.insertion_code] in res_selection:
                            res_list.append(j)
                            res_selection.remove([lres.chain, lres.number, lres.insertion_code])
                res_list.sort()
                self.INPUT['print_res'] = ','.join([str(x) for x in res_list])
                if res_selection:
                    for res in res_selection:
                        GMXMMPBSA_WARNING("We couldn't find this residue CHAIN:{} RES_NUM:{} ICODE: {}".format(*res))

    @staticmethod
    def cleantop(top_file, ndx):
        """
        Create a new top file with selected groups and without SOL and IONS
        :param top_file: User-defined topology file
        :param ndx: atoms index
        :return: new and clean top instance
        """
        top_file = Path(top_file)
        molsect = False

        ttp_file = top_file.parent.joinpath('_temp_top.top')
        temp_top = ttp_file.open(mode='w')
        # temp_top.write('; Modified by gmx_MMPBSA\n')

        with open(top_file) as topf:
            for line in topf:
                if line.startswith('#include') and 'forcefield.itp' in line:
                    if not ('charmm' in line.lower() or 'toppar' in line.lower() or 'amber' in line.lower()):
                        GMXMMPBSA_ERROR(f'Unknown force field in GROMACS topology in line:\n {line}')
                elif '[ molecules ]' in line:
                    molsect = True
                if molsect:
                    # not copy ions and solvent
                    sol_ion = [
                        # standard gmx form
                        'NA', 'CL', 'SOL',
                        # charmm-GUI form ??
                        'SOD', 'CLA', 'TIP3P', 'TIP4P', 'TIPS3P', 'TIP5P', 'SPC', 'SPC/E', 'SPCE', 'TIP3o', 'WAT']
                    if not line.split():
                        continue
                    if line.split()[0].strip().upper() in sol_ion:
                        continue
                temp_top.write(line)
        temp_top.close()

        # read the temp topology with parmed
        rtemp_top = parmed.gromacs.GromacsTopologyFile(ttp_file.as_posix())
        # get the residues in the top from the com_ndx
        res_list = []
        for i in ndx:
            if rtemp_top.atoms[i - 1].residue.number + 1 not in res_list:
                res_list.append(rtemp_top.atoms[i - 1].residue.number + 1)

        ranges = list2range(res_list)
        rtemp_top.strip(f"!:{','.join(ranges['string'])}")
        ttp_file.unlink()
        return rtemp_top

    def get_masks(self):
        rec_mask = ':' + ','.join(self.resi['REC']['string'])
        lig_mask = ':' + ','.join(self.resi['LIG']['string'])

        # to change the self.resl to get the mutant label in decomp analysis
        self.resl['MUT_COM'] = deepcopy(self.resl['COM'])
        self.resl['MUT_COM'][self.com_mut_index].name = self.INPUT['mutant']
        # self.com_mut_index, self.part_mut, self.part_index,
        self.resl['MUT_REC'] = deepcopy(self.resl['REC'])
        self.resl['MUT_LIG'] = deepcopy(self.resl['LIG'])
        if self.part_mut == 'REC':
            self.resl['MUT_REC'][self.part_index].name = self.INPUT['mutant']
        else:
            self.resl['MUT_LIG'][self.part_index].name = self.INPUT['mutant']

        return rec_mask, lig_mask, self.resl

    def get_selected_residues(self, select):
        """
        Convert string selection format to amber index list
        """
        dist, res_selection = selector(select)
        sele_res = []
        if dist:
            for i in self.resl['REC']:
                for j in self.resl['LIG']:
                    for rat in self.complex_str.residues[i - 1].atoms:
                        rat_coor = [rat.xx, rat.xy, rat.xz]
                        for lat in self.complex_str.residues[j - 1].atoms:
                            lat_coor = [lat.xx, lat.xy, lat.xz]
                            if get_dist(rat_coor, lat_coor) <= dist:
                                if i not in sele_res:
                                    sele_res.append(i)
                                if j not in sele_res:
                                    sele_res.append(j)
                                break
        elif res_selection:
            for i in self.resl['REC']:
                rres = self.complex_str.residues[i - 1]
                if [rres.chain, rres.number, rres.insertion_code] in res_selection:
                    sele_res.append(i)
                    res_selection.remove([rres.chain, rres.number, rres.insertion_code])
            for j in self.resl['LIG']:
                lres = self.complex_str.residues[j - 1]
                if [lres.chain, lres.number, lres.insertion_code] in res_selection:
                    sele_res.append(j)
                    res_selection.remove([lres.chain, lres.number, lres.insertion_code])
        sele_res.sort()
        if res_selection:
            for res in res_selection:
                GMXMMPBSA_WARNING("We couldn't find this residue CHAIN:{} RES_NUM:{} ICODE: "
                                  "{}".format(*res))
        return sele_res

    def res2map(self):
        """

        :param com_str:
        :return:
        """
        # read the index file
        ndx = {}
        with open(self.FILES.complex_index) as indexf:
            header = None
            for line in indexf:
                if line.startswith('['):
                    header = line.strip('\n[] ')
                    ndx[header] = []
                else:
                    ndx[header].extend(map(int, line.split()))
        com_str = self.complex_str
        order_list = []

        masks = {'REC': [], 'LIG': []}
        res_list = {'REC': [], 'LIG': [], 'COM': []}
        com_ndx = ndx['GMXMMPBSA_REC_GMXMMPBSA_LIG']
        com_len = len(ndx['GMXMMPBSA_REC_GMXMMPBSA_LIG'])
        resindex = 1
        proc_res = None
        for i in range(com_len):
            res = [com_str.atoms[i].residue.chain, com_str.atoms[i].residue.number, com_str.atoms[
                i].residue.insertion_code]
            # We check who owns the residue corresponding to this atom
            if com_ndx[i] in ndx['GMXMMPBSA_REC']:
                # save residue number in the rec list
                if res != proc_res and resindex not in res_list['REC']:
                    res_list['REC'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                                   com_str.atoms[i].residue.chain, 'REC', com_str.atoms[i].residue.name,
                                                   com_str.atoms[i].residue.insertion_code))
                    res_list['COM'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                                   com_str.atoms[i].residue.chain, 'COM', com_str.atoms[i].residue.name,
                                                   com_str.atoms[i].residue.insertion_code))
                    resindex += 1
                    proc_res = res
            # save residue number in the lig list
            elif res != proc_res and resindex not in res_list['LIG']:
                    res_list['LIG'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                                   com_str.atoms[i].residue.chain, 'LIG', com_str.atoms[i].residue.name,
                                                   com_str.atoms[i].residue.insertion_code))
                    res_list['COM'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                                   com_str.atoms[i].residue.chain, 'COM', com_str.atoms[i].residue.name,
                                                   com_str.atoms[i].residue.insertion_code))
                    resindex += 1
                    proc_res = res

        masks['REC'] = list2range(res_list['REC'])
        masks['LIG'] = list2range(res_list['LIG'])

        temp = []
        for m, value in masks.items():
            for e in value['num']:
                v = e[0] if isinstance(e, list) else e
                temp.append([v, m])
        temp.sort(key=lambda x: x[0])
        order_list = [c[1] for c in temp]

        return masks, res_list, order_list, [ndx['GMXMMPBSA_REC'], ndx['GMXMMPBSA_LIG']]

    def fixparm2amber(self, structure, removeH=False):

        for residue in structure.residues:
            # change atoms name from GROMACS to AMBER
            for atom in residue.atoms:
                if atom.name == 'OC1':
                    atom.name = 'O'
                elif atom.name == 'OC2':
                    atom.name = 'OXT'
                    residue.ter = True  # parmed terminal
            # change residues name according to AMBER
            if residue.name == 'ILE':
                for atom in residue.atoms:
                    if atom.name == 'CD':
                        atom.name = 'CD1'
                        break
            elif residue.name == 'LYS':
                atoms = [atom.name for atom in residue.atoms]
                if 'HZ3' not in atoms:
                    residue.name = 'LYN'
            elif residue.name == 'ASP':
                atoms = [atom.name for atom in residue.atoms]
                if 'HD2' in atoms:
                    residue.name = 'ASH'
            elif residue.name == 'GLU':
                atoms = [atom.name for atom in residue.atoms]
                if 'HE2' in atoms:
                    residue.name = 'GLH'
            elif residue.name in his:
                atoms = [atom.name for atom in residue.atoms if atom.atomic_number == 1]
                if 'HD1' in atoms and 'HE2' in atoms:
                    residue.name = 'HIP'
                elif 'HD1' in atoms:
                    residue.name = 'HID'
                elif 'HE2' in atoms:
                    residue.name = 'HIE'
            elif residue.name in cys_name:
                if residue.name == 'CYX':
                    continue
                for atom in residue.atoms:
                    if 'SG' in atom.name:
                        for bondedatm in atom.bond_partners:
                            if bondedatm.name == 'SG':
                                residue.name = 'CYX'
                                bondedatm.residue.name = 'CYX'
                        break
            # GROMACS 4.x save the pdb without atom element column, so parmed does not recognize some H atoms.
            # Parmed assigns 0 to the atomic number of these atoms. In order to correctly eliminate hydrogens,
            # it is necessary to assign the atomic number.
            if len(self.make_ndx) == 2:
                for atom in residue.atoms:
                    if 'H' in atom.name and atom.atomic_number == 0:
                        atom.atomic_number = 1
            # Remove H atoms. Only when using the pdb files with tleap to build the topologies
        if removeH:
            structure.strip('@/H')

    def getMutationInfo(self):
        if not self.INPUT['mutant_res']:
            GMXMMPBSA_ERROR("No residue for mutation was defined")
        # dict = { resind: [chain, resnum, icode]
        sele_res_dict = self.get_selected_residues(self.INPUT['mutant_res'])
        if not len(sele_res_dict) == 1:
            GMXMMPBSA_ERROR('Only ONE mutant residue is allowed.')
        r = sele_res_dict[0]
        res = self.complex_str.residues[r - 1]
        icode = ':' + res.insertion_code if res.insertion_code else ''
        if not parmed.residue.AminoAcidResidue.has(res.name):
            GMXMMPBSA_WARNING(f"Selecting residue {res.chain}:{res.name}:{res.number}{icode} can't be mutated and "
                              f"will be ignored...")
        label = f"{res.name}[{res.chain}:{res.number}]{self.INPUT['mutant']}"
        if icode:
            label = f"{res.name}[{res.chain}:{res.number}:{res.insertion_code}]{self.INPUT['mutant']}"

        if r in self.resl['REC']:
            part_index = self.resl['REC'].index(r)
            part_mut = 'REC'
        elif r in self.resl['LIG']:
            part_index = self.resl['LIG'].index(r)
            part_mut = 'LIG'
        else:
            part_index = None
            part_mut = None
            if icode:
                GMXMMPBSA_ERROR(f'Residue {res.chain}:{res.number}:{res.insertion_code} not found')
            else:
                GMXMMPBSA_ERROR(f'Residue {res.chain}:{res.number} not found')

        # return r - 1 since r is the complex mutant index from amber selection format. Needed for top mutation only
        return r - 1, part_mut, part_index, label

    def makeMutTop(self, wt_top, mut_index, pdb=False):
        """

        :param wt_top: Amber parm from GROMACS topology
        :param mut_index: index of mutation in structure
        :return: Mutant AmberParm
        """
        mut_top = self.molstr(wt_top)
        mut_aa = self.INPUT['mutant']

        bb_atoms = 'N,H,CA,HA,C,O,'
        nterm_atoms = 'H1,H2,H3,'
        cterm_atoms = 'OXT'
        sc_cb_atom = 'CB,'
        sc_ala_atoms = ('HB,' +  # VAL, ILE, THR
                        'HB1,HB2,' +
                        'CG1,CG2,OG1,' +  # VAL, ILE, THR
                        'OG,' +  # SER
                        'SG,' +  # CYS
                        'CG,')

        if mut_aa in ['GLY', 'G']:
            # FIXME: allow terminal residues to mutate?
            strip_mask = f':{mut_index + 1} &!@' + bb_atoms + nterm_atoms + cterm_atoms
            if not pdb:
                strip_mask += sc_cb_atom
        else:
            # FIXME: allow terminal residues to mutate?
            strip_mask = f':{mut_index + 1} &!@' + bb_atoms + sc_cb_atom + nterm_atoms + cterm_atoms
            if not pdb:
                strip_mask += sc_ala_atoms
        mut_top.strip(strip_mask)

        h_atoms_prop = {}
        # get an example HB atom if not PDB
        if not pdb:
            for res in mut_top.residues:
                if res.name == mut_aa:
                    for at in res.atoms:
                        if mut_aa == 'GLY':
                            if at.name in ['HA2']:
                                h_atoms_prop['mass'] = at.mass
                                h_atoms_prop['element'] = at.element
                                h_atoms_prop['atomic_number'] = at.atomic_number
                                h_atoms_prop['charge'] = at.charge
                                h_atoms_prop['atom_type'] = at.atom_type
                                h_atoms_prop['type'] = at.type
                                break
                        else:
                            if at.name in ['HB2']:
                                h_atoms_prop['mass'] = at.mass
                                h_atoms_prop['element'] = at.element
                                h_atoms_prop['atomic_number'] = at.atomic_number
                                h_atoms_prop['charge'] = at.charge
                                h_atoms_prop['atom_type'] = at.atom_type
                                h_atoms_prop['type'] = at.type
                                break
                    break
        cb_atom = None
        ca_atom = None
        logging.info(f"Mutating {mut_top.residues[mut_index].name} by {mut_aa}")

        # change intdiel if cas_intdiel was define
        if self.INPUT['cas_intdiel']:
            if self.INPUT['gbrun']:
                if self.INPUT['intdiel'] != 1.0:
                    GMXMMPBSA_WARNING('Both cas_intdiel and intdiel were defined. The dielectric constants associated '
                                      'with cas_intdiel will be ignored and intdiel will be used instead')
                else:
                    if mut_top.residues[mut_index].name in polar_aa:
                        self.INPUT['intdiel'] = self.INPUT['intdiel_polar']
                        self.INPUT['indi'] = self.INPUT['intdiel_polar']
                        logging.info(f"Setting intdiel = indi = intdiel_polar = {self.INPUT['intdiel_polar']} for "
                                     f"Alanine scanning")
                    elif mut_top.residues[mut_index].name in nonpolar_aa:
                        self.INPUT['intdiel'] = self.INPUT['intdiel_nonpolar']
                        self.INPUT['indi'] = self.INPUT['intdiel_nonpolar']
                        logging.info(f"Setting intdiel = indi = intdiel_nonpolar = {self.INPUT['intdiel_nonpolar']} "
                                     f"for Alanine scanning")
                    elif mut_top.residues[mut_index].name in positive_aa:
                        self.INPUT['intdiel'] = self.INPUT['intdiel_positive']
                        self.INPUT['indi'] = self.INPUT['intdiel_positive']
                        logging.info(f"Setting intdiel = indi = intdiel_positive = {self.INPUT['intdiel_positive']} "
                                     f"for Alanine scanning")
                    elif mut_top.residues[mut_index].name in negative_aa:
                        self.INPUT['intdiel'] = self.INPUT['intdiel_negative']
                        self.INPUT['indi'] = self.INPUT['intdiel_negative']
                        logging.info(f"Setting intdiel = indi = intdiel_negative = {self.INPUT['intdiel_negative']} "
                                     f"for Alanine scanning")
                    else:
                        GMXMMPBSA_WARNING(f"Unclassified mutant residue {mut_top.residues[mut_index].name}. The "
                                          f"default intdiel will be used")
            if self.INPUT['pbrun']:
                if self.INPUT['indi'] != 1.0:
                    GMXMMPBSA_WARNING('Both cas_intdiel and indi were defined. The dielectric constants associated '
                                      'with cas_intdiel will be ignored and indi will be used instead')
                else:
                    if mut_top.residues[mut_index].name in polar_aa:
                        self.INPUT['indi'] = self.INPUT['intdiel_polar']
                        logging.info(f"Setting indi = intdiel_polar = {self.INPUT['intdiel_polar']} for Alanine "
                                     f"scanning")
                    elif mut_top.residues[mut_index].name in nonpolar_aa:
                        self.INPUT['indi'] = self.INPUT['intdiel_nonpolar']
                        logging.info(f"Setting indi = intdiel_nonpolar = {self.INPUT['intdiel_nonpolar']} for "
                                     f"Alanine scanning")
                    elif mut_top.residues[mut_index].name in positive_aa:
                        self.INPUT['indi'] = self.INPUT['intdiel_positive']
                        logging.info(f"Setting intdiel = indi = intdiel_positive = {self.INPUT['intdiel_positive']} "
                                     f"for Alanine scanning")
                    elif mut_top.residues[mut_index].name in negative_aa:
                        self.INPUT['indi'] = self.INPUT['intdiel_negative']
                        logging.info(f"Setting indi = intdiel_negative = {self.INPUT['intdiel_negative']} for "
                                     f"Alanine scanning")
                    else:
                        GMXMMPBSA_WARNING(f"Unclassified mutant residue {mut_top.residues[mut_index].name}. The "
                                          f"default indi will be used")
        mut_top.residues[mut_index].name = mut_aa

        for at in mut_top.residues[mut_index].atoms:
            if mut_aa == 'GlY':
                if at.name == 'CA':
                    ca_atom = at
                if at.name in ['CB']:
                    at.name = 'HA2'
                    ca_atom.xx, ca_atom.xy, ca_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [ca_atom.xx, ca_atom.xy,
                         ca_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for ref_at in h_atoms_prop:
                        setattr(at, ref_at, h_atoms_prop[ref_at])
                elif at.name in ['HA']:
                    at.name = 'HA1'
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
            else:
                if at.name == 'CB':
                    cb_atom = at
                if at.name in ['CG1', 'CG2', 'OG1', 'OG', 'SG', 'CG']:
                    at.name = 'HB3'
                    cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [cb_atom.xx, cb_atom.xy,
                         cb_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for ref_at in h_atoms_prop:
                        setattr(at, ref_at, h_atoms_prop[ref_at])
                elif at.name in ['HB']:
                    at.name = 'HB1'
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
                elif at.name in ['HB1', 'HB2', 'HB3']:
                    at.type = h_atoms_prop['type']
                    at.atom_type = h_atoms_prop['atom_type']
        return mut_top

    def cleanup_trajs(self):
        # clear trajectory
        if not self.INPUT['solvated_trajectory']:
            return
        logging.info('Cleaning normal complex trajectories...')
        new_trajs = []
        for i in range(len(self.FILES.complex_trajs)):
            trjconv_echo_args = ['echo', 'GMXMMPBSA_REC_GMXMMPBSA_LIG']
            c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file and make amber topology for complex
            trjconv_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr,
                                      '-o', 'COM_traj_{}.xtc'.format(i), '-n', self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(trjconv_args))
            c6 = subprocess.Popen(trjconv_args,  # FIXME: start and end frames???
                                  stdin=c5.stdout, stdout=self.log, stderr=self.log)
            if c6.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.complex_tpr))
            new_trajs.append('COM_traj_{}.xtc'.format(i))
        self.FILES.complex_trajs = new_trajs

        # clear trajectory
        if self.FILES.receptor_tpr:
            logging.info('Cleaning normal receptor trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trjconv_echo_args = ['echo', '{}'.format(self.FILES.receptor_group)]
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.receptor_trajs[0], '-s', self.FILES.receptor_tpr,
                                          '-o', 'REC_traj_{}.xtc'.format(i), '-n',
                                          self.FILES.receptor_index]
                if self.INPUT['debug_printlevel']:
                    logging.info('Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(
                        trjconv_args))
                c6 = subprocess.Popen(trjconv_args,  # FIXME: start and end frames???
                                      stdin=c5.stdout, stdout=self.log, stderr=self.log)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR(
                        '%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.receptor_tpr))
                new_trajs.append('REC_traj_{}.xtc'.format(i))
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_tpr:
            logging.info('Cleanig normal ligand trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.ligand_trajs)):
                trjconv_echo_args = ['echo', '{}'.format(self.FILES.ligand_group)]
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.ligand_trajs[0], '-s', self.FILES.ligand_tpr,
                                          '-o', 'LIG_traj_{}.xtc'.format(i), '-n', self.FILES.ligand_index]
                if self.INPUT['debug_printlevel']:
                    logging.info(
                        'Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=self.log, stderr=self.log)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.ligand_tpr))
                new_trajs.append('LIG_traj_{}.xtc'.format(i))
            self.FILES.ligand_trajs = new_trajs

    def fix_chains_IDs(self, com_str, rec_str=None, lig_str=None, ref_str=None):
        if ref_str:
            if len(ref_str.residues) != len(com_str.residues):
                GMXMMPBSA_ERROR('The number of amino acids in the reference structure is different from that of the '
                                'complex...')
            for c, res in enumerate(ref_str.residues, start=1):
                # TODO: check if the residues are the same
                com_str.residues[c-1].chain = res.chain
                if c in self.resl['REC']:
                    i = self.resl['REC'].index(c)
                    rec_str.residues[i].chain = res.chain
                else:
                    i = self.resl['LIG'].index(c)
                    lig_str.residues[i].chain = res.chain
        else:
            assign = False
            if self.INPUT['assign_chainID'] == 1:
                assign = not com_str.residues[0].chain  # pretty simple
                if assign:
                    logging.info('Chains ID not found. Assigning chains IDs...')
                else:
                    logging.info('Chains ID found. Ignoring chains ID assignation...')
            elif self.INPUT['assign_chainID'] == 2:
                assign = True
                if com_str.residues[0].chain:
                    logging.warning('Assigning chains ID...')
                else:
                    logging.warning('Already have chain ID. Re-assigning ID...')
            elif self.INPUT['assign_chainID'] == 0 and self.FILES.complex_tpr[-3:] == 'gro':
                assign = True
                logging.warning('No reference structure was found and a gro file was used for the complex '
                                'structure. Assigning chains ID...')

            if assign:
                self._assign_chains_IDs(com_str, rec_str, lig_str)
        # Save fixed complex structure for analysis and set it in FILES to save in info file
        com_str.save(self.FILES.prefix + 'COM_FIXED.pdb', 'pdb', True, renumber=False)

    def _assign_chains_IDs(self, com_str, rec_str, lig_str):
        chains_ids = []
        chain_by_num = False
        chain_by_ter = False
        previous_res_number = 0
        curr_chain_id = 'A'
        has_nucl = 0
        for c, res in enumerate(com_str.residues, start=1):
            if res.chain:
                if res.chain != curr_chain_id:
                    res.chain = curr_chain_id
                    if c in self.resl['REC']:
                        i = self.resl['REC'].index(c)
                        rec_str.residues[i].chain = res.chain
                    else:
                        i = self.resl['LIG'].index(c)
                        lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            else:
                res.chain = curr_chain_id

                if c in self.resl['REC']:
                    i = self.resl['REC'].index(c)
                    rec_str.residues[i].chain = res.chain
                else:
                    i = self.resl['LIG'].index(c)
                    lig_str.residues[i].chain = res.chain
                if curr_chain_id not in chains_ids:
                    chains_ids.append(curr_chain_id)
                    # see if it is the end of chain
            if res.number != previous_res_number + 1 and previous_res_number != 0:
                chain_by_num = True
            if chain_by_num and chain_by_ter:
                chain_by_num = False
                chain_by_ter = False
                curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                res.chain = curr_chain_id
                if c in self.resl['REC']:
                    i = self.resl['REC'].index(c)
                    rec_str.residues[i].chain = res.chain
                else:
                    i = self.resl['LIG'].index(c)
                    lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            elif chain_by_ter:
                chain_by_ter = False
            elif chain_by_num:
                chain_by_num = False
                curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                res.chain = curr_chain_id

                if c + 1 in self.resl['REC']:
                    i = self.resl['REC'].index(c)
                    rec_str.residues[i].chain = res.chain
                else:
                    i = self.resl['LIG'].index(c)
                    lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            for atm in res.atoms:
                if atm.name == 'OXT':  # only for protein
                    res.ter = True
                    chain_by_ter = True
            if parmed.residue.RNAResidue.has(res.name) or parmed.residue.DNAResidue.has(res.name):
                has_nucl += 1

            previous_res_number = res.number
        if has_nucl == 1:
            logging.warning('This structure contains nucleotides. We recommend that you use the reference '
                            'structure')
    @staticmethod
    def molstr(data):

        if type(data) == str:
            # data is a pdb file
            pdb_file = data
            try:
                with open(pdb_file) as fo:
                    fo = fo.readlines()
                    for line in fo:
                        if 'MODEL' in line or 'ENDMDL' in line:
                            fo.remove(line)
                with open(pdb_file, 'w') as fw:
                    for x in fo:
                        fw.write(x)
            except IOError as e:
                GMXMMPBSA_ERROR('', str(e))

            structure = parmed.read_PDB(pdb_file)
        else:
            # data is Structure, AmberParm, ChamberParm or GromacsTopologyFile. This make a copy
            structure = data.__copy__()
            for c, at in enumerate(structure.atoms, start=1):
                at.number = c
        return structure

    def _write_ff(self, ofile):
        for ff in self.INPUT['forcefields']:
            ofile.write(f'source {ff}\n')
        ofile.write('loadOff atomic_ions.lib\n')
        ofile.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['ions_parameters']]))
        ofile.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['PBRadii']]))


    def makeToptleap(self):
        logging.info('Building tleap input files...')
        with open(self.FILES.prefix + 'leap.in', 'w') as tif:
            self._write_ff(tif)
            REC = []
            LIG = []
            for rec in self.receptor_list:
                REC.append(f'{rec}')
                tif.write(f'{rec} = loadpdb {self.receptor_list[rec]}\n')
            rec_out = ' '.join(REC)

            # check if ligand is not protein and always load
            if self.FILES.ligand_mol2:
                tif.write('LIG1 = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                tif.write('loadamberparams {}\n'.format(self.ligand_frcmod))
                tif.write('check LIG1\n')

                if self.FILES.stability:
                    self.ligand_pmrtop = None
                else:
                    tif.write(f'saveamberparm LIG1 {self.ligand_pmrtop} {self.FILES.prefix}LIG.inpcrd\n')
                for lig in self.ligand_list:
                    LIG.append(f'{lig}')
            else:
                for lig in self.ligand_list:
                    LIG.append(f'{lig}')
                    tif.write(f'{lig} = loadpdb {self.ligand_list[lig]}\n')
                lig_out = ' '.join(LIG)
                if self.FILES.stability:
                    self.ligand_pmrtop = None
                else:
                    tif.write(f'LIG_OUT = combine {{ {lig_out} }}\n')
                    tif.write(f'saveamberparm LIG_OUT {self.ligand_pmrtop} {self.FILES.prefix}LIG.inpcrd\n')
            COM = self._set_com_order(REC, LIG)
            if self.FILES.stability:
                self.receptor_pmrtop = None
            else:
                tif.write(f'REC_OUT = combine {{ { rec_out } }}\n')
                tif.write(f'saveamberparm REC_OUT {self.receptor_pmrtop} {self.FILES.prefix}REC.inpcrd\n')
            com_out = ' '.join(COM)
            tif.write(f'COM_OUT = combine {{ {com_out} }}\n')
            tif.write('saveamberparm COM_OUT {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            tif.write('quit')
        # changed in v1.4.3. We source the gmxMMPBSA ff directly from the data folder instead of copy to the Amber/dat
        data_path = Path(__file__).parent.joinpath('data')
        tleap = self.external_progs['tleap']
        self._run_tleap(tleap, 'leap.in', data_path)

        if self.INPUT['alarun']:
            with open(self.FILES.prefix + 'mut_leap.in', 'w') as mtif:
                self._write_ff(mtif)

                if self.mutant_receptor_pmrtop:
                    REC = []
                    for mrec in self.mut_receptor_list:
                        REC.append(f'{mrec}')
                        mtif.write(f'{mrec} = loadpdb {self.mut_receptor_list[mrec]}\n')
                    mrec_out = ' '.join(REC)

                    if not self.FILES.stability:
                        mtif.write(f'MREC_OUT = combine {{ {mrec_out} }}\n')
                        mtif.write(
                            'saveamberparm MREC_OUT {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
                                                                                    p=self.FILES.prefix))
                    else:
                        self.mutant_receptor_pmrtop = None
                    # check if ligand is not protein and always load
                    if self.FILES.ligand_mol2:
                        mtif.write('LIG1 = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                        self.mutant_ligand_pmrtop = None
                        if not self.FILES.stability:
                            mtif.write('check LIG1\n')
                            mtif.write('loadamberparams {}\n'.format(self.ligand_frcmod))
                        else:
                            self.mutant_ligand_pmrtop = None
                    else:
                        for lig in self.ligand_list:
                            mtif.write(f'{lig} = loadpdb {self.ligand_list[lig]}\n')
                else:
                    LIG = []
                    for mlig in self.mut_ligand_list:
                        LIG.append(f'{mlig}')
                        mtif.write(f'{mlig} = loadpdb {self.mut_ligand_list[mlig]}\n')
                    mlig_out = ' '.join(LIG)

                    if not self.FILES.stability:
                        mtif.write(f'MLIG_OUT = combine {{ {mlig_out} }}\n')
                        mtif.write('saveamberparm MLIG_OUT {t} {p}MUT_LIG.inpcrd\n'.format(
                            t=self.mutant_ligand_pmrtop, p=self.FILES.prefix))
                    else:
                        self.mutant_ligand_pmrtop = None
                    for rec in self.receptor_list:
                        mtif.write(f'{rec} = loadpdb {self.receptor_list[rec]}\n')

                MCOM = self._set_com_order(REC, LIG)
                mcom_out = ' '.join(MCOM)
                mtif.write(f'MCOM_OUT = combine {{ {mcom_out} }}\n')
                mtif.write('saveamberparm MCOM_OUT {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                                   p=self.FILES.prefix))
                mtif.write('quit')

            self._run_tleap(tleap, 'mut_leap.in', data_path)

        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

    def _run_tleap(self, tleap, arg1, data_path):
        tleap_args = [
            tleap,
            '-f',
            '{}'.format(self.FILES.prefix + arg1),
            '-I',
            data_path.as_posix(),
        ]

        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + ' '.join(tleap_args))
        p1 = subprocess.Popen(tleap_args, stdout=self.log, stderr=self.log)
        if p1.wait():
            GMXMMPBSA_ERROR(
                '%s failed when querying %s' % (tleap, self.FILES.prefix + arg1)
            )

    def _set_com_order(self, REC, LIG):
        result = []
        l = 0
        r = 0
        for i, e in enumerate(self.orderl):
            if e in ['R', 'REC']:
                result.append(REC[r])
                r += 1
            else:
                result.append(LIG[l])
                l += 1
        return result
