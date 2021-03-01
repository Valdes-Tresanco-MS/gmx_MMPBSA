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
from GMXMMPBSA.utils import checkff
from GMXMMPBSA.alamdcrd import _scaledistance
import subprocess
from math import sqrt
from pathlib import Path
import logging
import string
from parmed.tools.changeradii import ChRad

chains_letters = list(string.ascii_uppercase)
his = ['HIS', 'HIE', 'HID', 'HIP']
cys_name = ['CYS', 'CYX', 'CYM']
lys = ['LYS', 'LYN']
asp = ['ASP', 'ASH']
glu = ['GLU', 'GLH']
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


def dist(coor1, coor2):
    return sqrt((coor2[0] - coor1[0]) ** 2 + (coor2[1] - coor1[1]) ** 2 + (coor2[2] - coor1[2]) ** 2)


class CheckMakeTop:
    def __init__(self, FILES, INPUT, external_programs):
        self.FILES = FILES
        self.INPUT = INPUT
        self.external_progs = external_programs
        self.use_temp = False
        self.log = open('gmx_MMPBSA.log', 'a')
        self.print_residues = 'within' in self.INPUT['print_res']  # FIXME: this is pretty ugly
        self.within = 4
        if self.print_residues:
            self.within = float(self.INPUT['print_res'].split()[1])

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

        checkff(self.INPUT['overwrite_data'])
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
            logging.info('Cleaning GROMACS topologies...')
            self.cleantop(self.FILES.complex_top, self.complex_temp_top)
            if self.FILES.receptor_top:
                self.cleantop(self.FILES.receptor_top, self.receptor_temp_top)
            if self.FILES.ligand_top:
                self.cleantop(self.FILES.ligand_top, self.ligand_temp_top)
            tops = self.gmxtop2prmtop()
        else:
            self.pdb2prmtop()
            tops = self.makeToptleap()

        self.reswithin()
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
            parmchk2_args = [parmchk2, '-i', self.FILES.ligand_mol2, '-f', 'mol2', '-o', self.ligand_frcmod]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + ' '.join(parmchk2_args))
            l3 = subprocess.Popen(parmchk2_args, stdout=self.log, stderr=self.log)
            if l3.wait():
                GMXMMPBSA_ERROR('%s failed when querying %s' % (parmchk2, self.FILES.ligand_mol2))

        # make a temp receptor pdb (even when stability) if decomp to get correct receptor residues from complex. This
        # avoid get multiples molecules from complex.split()
        if self.INPUT['decomprun'] and self.print_residues:
            if self.FILES.stability:
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
            if self.INPUT['alarun'] and (self.FILES.mutant_receptor_tpr or self.FILES.mutant_ligand_tpr):
                logging.warning('When Stability calculation mode is selected, mutant receptor/mutant ligand files '
                                'are not needed...')

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

        # initialize receptor and ligand structures. Needed to get residues map
        self.complex_str = self.molstr(self.complex_str_file)
        self.receptor_str = self.molstr(self.receptor_str_file)
        self.ligand_str = self.molstr(self.ligand_str_file)
        self.resi, self.resl, self.orderl = self.res2map()
        self.fix_chains_IDs(self.complex_str, self.receptor_str, self.ligand_str, self.ref_str)

    def gmxtop2prmtop(self):
        logging.info('Building Normal Complex Amber Topology...')
        com_top = parmed.gromacs.GromacsTopologyFile(self.complex_temp_top, xyz=self.complex_str_file)
        # try:
        if com_top.impropers or com_top.urey_bradleys or com_top.cmaps:
            com_amb_prm = parmed.amber.ChamberParm.from_structure(com_top)
            com_top_parm = 'chamber'
        else:
            com_amb_prm = parmed.amber.AmberParm.from_structure(com_top)
            com_top_parm = 'amber'

        # IMPORTANT: make_trajs ends in error if the box is defined
        com_amb_prm.box = None
        # except TypeError as err:
        #     GMXMMPBSA_ERROR(str(err))

        self.fixparm2amber(com_amb_prm)

        logging.info('Writing Normal Complex Amber Topology...')
        # change de PBRadii
        action = ChRad(com_amb_prm, PBRadii[self.INPUT['PBRadii']])
        com_amb_prm.write_parm(self.complex_pmrtop)

        text_list = []
        for r in self.resi['REC']:
            if r[0] == r[1]:
                text_list.append(f'{r[0]}')
            else:
                text_list.append(f'{r[0]}-{r[1]}')
        rec_indexes_string = ','.join(text_list)

        rec_hastop = True
        if self.FILES.receptor_top:
            logging.info('A Receptor topology file was defined. Using MT approach...')
            logging.info('Building AMBER Receptor Topology from GROMACS Receptor Topology...')
            rec_top = parmed.gromacs.GromacsTopologyFile(self.receptor_temp_top, xyz=self.receptor_str_file)
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
            lig_top = parmed.gromacs.GromacsTopologyFile(self.ligand_temp_top, xyz=self.ligand_str_file)
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
            com_mut_index, part_mut, part_index, mut_label = self.getMutationIndex()
            mut_com_amb_prm = self.makeMutTop(com_amb_prm, com_mut_index)
            # change de PBRadii
            action = ChRad(mut_com_amb_prm, PBRadii[self.INPUT['PBRadii']])
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)

            if part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor Topology...')
                mut_com_amb_prm.strip(f'!:{rec_indexes_string}')
                mdata = [self.mutant_receptor_pmrtop, 'REC']
                self.mutant_ligand_pmrtop = None
                if rec_hastop:
                    mtop = self.makeMutTop(rec_amb_prm, part_index)
                else:
                    mtop = mut_com_amb_prm
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand Topology...')
                mut_com_amb_prm.strip(f':{rec_indexes_string}')
                mdata = [self.mutant_ligand_pmrtop, 'LIG']
                self.mutant_receptor_pmrtop = None
                if lig_hastop:
                    mtop = self.makeMutTop(lig_amb_prm, part_index)
                else:
                    mtop = mut_com_amb_prm

            if self.INPUT['protein_forcefield'] == 'charmm':
                mut_prot_amb_prm = parmed.amber.ChamberParm.from_structure(mtop)
            else:
                mut_prot_amb_prm = parmed.amber.AmberParm.from_structure(mtop)
            # change de PBRadii
            action = ChRad(mut_prot_amb_prm, PBRadii[self.INPUT['PBRadii']])
            mut_prot_amb_prm.write_parm(mdata[0])
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

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
        c = 1
        for r in self.resi['REC']:
            end = start + (r[1]- r[0])
            mask = f'!:{start}-{end}'
            start += end
            rec = self.molstr(self.receptor_str)
            rec.strip(mask)
            rec_file = self.FILES.prefix + f'REC_F{c}.pdb'
            rec.save(rec_file, 'pdb', True, renumber=False)
            self.receptor_list[f'REC{c}'] = rec_file
            c += 1

        self.ligand_list = {}
        start = 1
        c = 1
        for r in self.resi['LIG']:
            end = start + (r[1] - r[0])
            mask = f'!:{start}-{end}'
            start += end
            lig = self.molstr(self.ligand_str)
            lig.strip(mask)
            lig_file = self.FILES.prefix + f'LIG_F{c}.pdb'
            lig.save(lig_file, 'pdb', True, renumber=False)
            self.ligand_list[f'LIG{c}'] = lig_file
            c += 1

        self.mut_receptor_list = {}
        self.mut_ligand_list = {}

        if self.INPUT['alarun']:
            com_mut_index, part_mut, part_index, mut_label = self.getMutationIndex()
            if part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor Structure...')
                self.mutant_ligand_pmrtop = None
                start = 0
                c = 1
                for r in self.resi['REC']:
                    mask = f'!:{start}-{(r[1] - r[0]) + 1}'
                    rec = self.molstr(self.receptor_str)
                    mut_rec = self.makeMutTop(rec, part_index)
                    mut_rec.strip(mask)
                    mut_rec_file = self.FILES.prefix + f'MUT_REC_F{c}.pdb'
                    mut_rec.save(mut_rec_file, 'pdb', True, renumber=False)
                    self.mut_receptor_list[f'MREC{c}'] = mut_rec_file
                    c += 1
            else:
                logging.info('Detecting mutation in Ligand.Building Mutant Ligand Structure...')
                self.mutant_receptor_pmrtop = None
                start = 0
                c = 1
                for r in self.resi['LIG']:
                    mask = f'!:{start}-{(r[1] - r[0]) + 1}'
                    lig = self.molstr(self.ligand_str)
                    mut_lig = self.makeMutTop(lig, part_index)
                    mut_lig.strip(mask)
                    mut_lig_file = self.FILES.prefix + f'MUT_LIG_F{c}.pdb'
                    mut_lig.save(mut_lig_file, 'pdb', True, renumber=False)
                    self.mut_ligand_list[f'MLIG{c}'] = mut_lig_file
                    c += 1

    def reswithin(self):
        # Get residue form receptor-ligand interface
        if self.print_residues and self.INPUT['decomprun']:
            res_list = []
            for i in self.resl['REC']:
                if i in res_list:
                    continue
                for j in self.resl['LIG']:
                    if j in res_list:
                        continue
                    for rat in self.complex_str.residues[i - 1].atoms:
                        rat_coor = [rat.xx, rat.xy, rat.xz]
                        if i in res_list:
                            break
                        for lat in self.complex_str.residues[j - 1].atoms:
                            lat_coor = [lat.xx, lat.xy, lat.xz]
                            if dist(rat_coor, lat_coor) <= self.within:
                                if i not in res_list:
                                    res_list.append(i)
                                if j not in res_list:
                                    res_list.append(j)
                                break
            res_list.sort()
            self.INPUT['print_res'] = ','.join([str(x + 1) for x in res_list])

    def cleantop(self, top_file, temp_top_file):
        """
        Create a new top file without SOL and IONS
        :param top_file: User-defined topology file
        :param temp_top_file: temporary top file
        :return: detected ff
        """
        top_file = Path(top_file)
        molsect = False

        temp_top = open(top_file.parent.joinpath(temp_top_file), 'w')
        temp_top.write('; Modified by gmx_MMPBSA\n')

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
                    if line.split()[0].strip().upper() in sol_ion:
                        break
                temp_top.write(line)
        temp_top.close()

    def get_masks(self):
        rt = []
        for r in self.resi['REC']:
            if r[0] == r[1]: # only as precaution
                rt.append(f'{r[0]}')
            else:
                rt.append(f'{r[0]}-{r[1]}')
        rec_mask = ':' + ','.join(rt)
        lt = []
        for l in self.resi['LIG']:
            if l[0] == l[1]:
                lt.append(f'{l[0]}')
            else:
                lt.append(f'{l[0]}-{l[1]}')
        lig_mask = ':' + ','.join(lt)
        return rec_mask, lig_mask

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

        start = 1
        end = None
        order_list = []
        previous = None

        masks = {'REC': [], 'LIG': []}
        res_list = {'REC': [], 'LIG': []}
        com_ndx = ndx['GMXMMPBSA_REC_GMXMMPBSA_LIG']
        com_len = len(ndx['GMXMMPBSA_REC_GMXMMPBSA_LIG'])
        resnum = 1
        current_res = None
        for i in range(com_len):
            # We check who owns the residue corresponding to this atom
            if com_ndx[i] in ndx['GMXMMPBSA_REC']:
                current = 'R'
                # save residue number in the rec list
                if com_str.atoms[i].residue.number != current_res and not resnum in res_list['REC']:
                    res_list['REC'].append(resnum)
                    resnum += 1
                    current_res = com_str.atoms[i].residue.number
            else:
                current = 'L'
                # save residue number in the lig list
                if com_str.atoms[i].residue.number != current_res and not resnum in res_list['LIG']:
                    res_list['LIG'].append(resnum)
                    resnum += 1
                    current_res = com_str.atoms[i].residue.number
            # check for end
            if previous and current != previous:
                end = resnum - 2

            # when i is the last index
            if i == com_len - 1:
                end = resnum - 1
            if end:
                if previous == 'R':
                    masks['REC'].append([start, end])
                else:
                    masks['LIG'].append([start, end])
                # add current range identifier
                order_list.append(previous)
                # set the new start and reset end
                start = end + 1
                end = None
            # we change previous to current once it is processed
            previous = current
        return masks, res_list, order_list

    def fixparm2amber(self, structure, removeH=False):

        for residue in structure.residues:
            # change atoms name from GROMACS to AMBER
            if residue.name == 'ILE':
                for atom in residue.atoms:
                    if atom.name == 'CD':
                        atom.name = 'CD1'
            for atom in residue.atoms:
                if atom.name == 'OC1':
                    atom.name = 'O'
                elif atom.name == 'OC2':
                    atom.name = 'OXT'
                    residue.ter = True  # parmed terminal
            # change residues name according to AMBER
            if residue.name == 'LYS':
                atoms = [atom.name for atom in residue.atoms]
                if not 'HZ3' in atoms:
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

    def getMutationIndex(self):
        label = ''
        if not self.INPUT['mutant_res']:
            GMXMMPBSA_ERROR("No residue for mutation was defined")
        chain_, resnum_ = self.INPUT['mutant_res'].split(':')
        chain = str(chain_).strip().upper()
        resnum = int(str(resnum_).strip())

        if not chain or not resnum:
            GMXMMPBSA_ERROR("Wrong notation... You most define the residue to mutate as follow: CHAIN:RES_NUMBER")
        idx = 0
        for res in self.complex_str.residues:
            if res.number == int(resnum) and res.chain == chain:
                try:
                    parmed.residue.AminoAcidResidue.get(res.name, True)
                except KeyError as e:
                    GMXMMPBSA_ERROR('The mutation must be an amino acid residue ...')

                label = f'{res.chain}:{res.name}:{res.number}'
                break
            idx += 1

        if idx in self.resl['REC']:
            part_index = self.resl['REC'].index(idx)
            part_mut = 'REC'
        elif idx in self.resl['LIG']:
            part_index = self.resl['LIG'].index(idx)
            part_mut = 'LIG'
        else:
            part_index = None
            part_mut = None
            GMXMMPBSA_ERROR('Residue {}:{} not found'.format(chain, resnum))

        return (idx, part_mut, part_index, label)

    def makeMutTop(self, wt_top, mut_index):
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
            strip_mask = f':{mut_index + 1} &!@' + bb_atoms + sc_cb_atom + nterm_atoms + cterm_atoms
        else:
            # FIXME: allow terminal residues to mutate?
            strip_mask = f':{mut_index + 1} &!@' + bb_atoms + sc_cb_atom + sc_ala_atoms + nterm_atoms + cterm_atoms
        mut_top.strip(strip_mask)

        h_atoms_prop = {}
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
        mut_top.residues[mut_index].name = mut_aa
        ind = 0
        for res in mut_top.residues:
            if ind == mut_index:
                for at in res.atoms:
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
                        if at.name in ['CG', 'OG', 'CG2', 'SG']:
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
            ind += 1
        return mut_top

    def cleanup_trajs(self):
        # clear trajectory
        if self.INPUT['solvated_trajectory']:
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
            c = 0
            for res in ref_str.residues:
                res.chain = com_str.residues[c].chain
                if c + 1 in self.resl['REC']:
                    i = self.resl['REC'].index(c)
                    rec_str.residues[i].chain = res.chain
                else:
                    i = self.resl['LIG'].index(c)
                    lig_str.residues[i].chain = res.chain
                c += 1
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
                if not com_str.residues[0].chain:
                    logging.warning('Already have chain ID. Re-assigning ID...')
                else:
                    logging.warning('Assigning chains ID...')
            elif self.INPUT['assign_chainID'] == 0 and self.FILES.complex_tpr[-3] == 'gro':
                assign = True
                logging.warning('No reference structure was found and a gro file was used for the complex '
                                'structure. Assigning chains ID...')

            if assign:
                chains_ids = []
                chain_by_num = False
                chain_by_ter = False
                previous_res_number = 0
                curr_chain_id = 'A'
                has_nucl = 0
                c = 0
                for res in com_str.residues:
                    if not res.chain:
                        res.chain = curr_chain_id

                        if c + 1 in self.resl['REC']:
                            i = self.resl['REC'].index(c)
                            rec_str.residues[i].chain = res.chain
                        else:
                            i = self.resl['LIG'].index(c)
                            lig_str.residues[i].chain = res.chain
                        if curr_chain_id not in chains_ids:
                            chains_ids.append(curr_chain_id)
                    else:
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
                    # see if it is the end of chain
                    if res.number != previous_res_number + 1:
                        if previous_res_number != 0:
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
                    c += 1
                if has_nucl == 1:
                    logging.warning('This structure contains nucleotides. We recommend that you use the reference '
                                    'structure')

    def molstr(self, data):

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
            c = 1
            for at in structure.atoms:
                at.number = c
                c += 1
        return structure

    def makeToptleap(self):
        logging.info('Building Tleap input files...')
        with open(self.FILES.prefix + 'leap.in', 'w') as tif:
            tif.write('source {}\n'.format(self.INPUT['protein_forcefield']))
            tif.write('source {}\n'.format(self.INPUT['ligand_forcefield']))
            tif.write('loadOff atomic_ions.lib\n')
            tif.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['ions_parameters']]))
            tif.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['PBRadii']]))

            REC = []
            LIG = []
            for rec in self.receptor_list:
                REC.append(f'{rec}')
                tif.write(f'{rec} = loadpdb {self.receptor_list[rec]}\n')
            rec_out = ' '.join(REC)

            # check if ligand is not protein and always load
            if self.FILES.ligand_mol2:
                tif.write('LIG1 = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                tif.write('check LIG1\n')
                tif.write('loadamberparams {}\n'.format(self.ligand_frcmod))
                if not self.FILES.stability:
                    tif.write('saveamberparm LIG1 {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop,
                                                                              p=self.FILES.prefix))
                else:
                    self.ligand_pmrtop = None
                for lig in self.ligand_list:
                    LIG.append(f'{lig}')
            else:
                for lig in self.ligand_list:
                    LIG.append(f'{lig}')
                    tif.write(f'{lig} = loadpdb {self.ligand_list[lig]}\n')
                lig_out = ' '.join(LIG)
                if not self.FILES.stability:
                    tif.write(f'LIG_OUT = combine {{ {lig_out} }}\n')
                    tif.write('saveamberparm LIG_OUT {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop,
                                                                                 p=self.FILES.prefix))
                else:
                    self.ligand_pmrtop = None
            COM = []
            l = 0
            r = 0
            i = 0
            for e in self.orderl:
                if e == 'R':
                    COM.append(REC[r])
                    r += 1
                else:
                    COM.append(LIG[l])
                    l += 1
                i += 1

            if not self.FILES.stability:
                tif.write(f'REC_OUT = combine {{ { rec_out } }}\n')
                tif.write('saveamberparm REC_OUT {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop, p=self.FILES.prefix))
            else:
                self.receptor_pmrtop = None
            com_out = ' '.join(COM)
            tif.write(f'COM_OUT = combine {{ {com_out} }}\n')
            tif.write('saveamberparm COM_OUT {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            tif.write('quit')

        tleap = self.external_progs['tleap']
        tleap_args = [tleap, '-f', '{}'.format(self.FILES.prefix + 'leap.in')]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + ' '.join(tleap_args))
        p = subprocess.Popen(tleap_args, stdout=self.log, stderr=self.log)
        if p.wait():
            GMXMMPBSA_ERROR('%s failed when querying %s' % (tleap, self.FILES.prefix + 'leap.in'))

        if self.INPUT['alarun']:
            with open(self.FILES.prefix + 'mut_leap.in', 'w') as mtif:
                mtif.write('source {}\n'.format(self.INPUT['protein_forcefield']))
                mtif.write('source {}\n'.format(self.INPUT['ligand_forcefield']))
                mtif.write('loadOff atomic_ions.lib\n')
                mtif.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['ions_parameters']]))
                mtif.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['PBRadii']]))


                if self.mutant_receptor_pmrtop:
                    REC = []
                    for mrec in self.mut_receptor_list:
                        REC.append(f'{mrec}')
                        mtif.write(f'{mrec} = loadpdb {self.mut_receptor_list[mrec]}\n')
                    mrec_out = ' '.join(REC)

                    if not self.FILES.stability:
                        mtif.write(f'MREC_OUT = combine {{ {mrec_out} }}\n')
                        mtif.write('saveamberparm MREC_OUT {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
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
                        mtif.write(f'{mlig} = loadpdb {self.mut_ligand_list[lig]}\n')
                    mlig_out = ' '.join(LIG)

                    if not self.FILES.stability:
                        mtif.write(f'MLIG_OUT = combine {{ {mlig_out} }}\n')
                        mtif.write('saveamberparm MLIG_OUT {t} {p}MUT_LIG.inpcrd\n'.format(
                                t=self.mutant_ligand_pmrtop, p=self.FILES.prefix))
                    else:
                        self.mutant_ligand_pmrtop = None
                    for rec in self.receptor_list:
                        mtif.write(f'{rec} = loadpdb {self.receptor_list[rec]}\n')

                MCOM = []
                l = r = i = 0
                for e in self.orderl:
                    if e == 'R':
                        MCOM.append(REC[r])
                        r += 1
                    else:
                        MCOM.append(LIG[l])
                        l += 1
                    i += 1
                mcom_out = ' '.join(MCOM)
                mtif.write(f'MCOM_OUT = combine {{ {mcom_out} }}\n')
                mtif.write('saveamberparm MCOM_OUT {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                               p=self.FILES.prefix))
                mtif.write('quit')

            tleap_args = [tleap, '-f', '{}'.format(self.FILES.prefix + 'mut_leap.in')]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + ' '.join(tleap_args))
            p1 = subprocess.Popen(tleap_args, stdout=self.log, stderr=self.log)
            if p1.wait():
                GMXMMPBSA_ERROR('%s failed when querying %s' % (tleap, self.FILES.prefix + 'mut_leap.in'))
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)
