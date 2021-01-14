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
import subprocess
from math import sqrt
import logging
import string

chains_letters = list(string.ascii_uppercase)

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
        self.print_residues = self.INPUT['print_res'].split()[0] == 'within'  # FIXME: this is pretty ugly
        self.within = 4
        if self.print_residues:
            self.within = float(self.INPUT['print_res'].split()[1])

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

        self.complex_pdb = self.FILES.prefix + 'COM.pdb'
        self.receptor_pdb = self.FILES.prefix + 'REC.pdb'
        self.ligand_pdb = self.FILES.prefix + 'LIG.pdb'
        self.complex_pdb_fixed = self.FILES.prefix + 'COM_FIXED.pdb'
        self.receptor_pdb_fixed = self.FILES.prefix + 'REC_FIXED.pdb'
        self.ligand_pdb_fixed = self.FILES.prefix + 'LIG_FIXED.pdb'

        self.rec_ions_pdb = self.FILES.prefix + 'REC_IONS.pdb'
        self.lig_ions_pdb = self.FILES.prefix + 'LIG_IONS.pdb'

        self.mutant_complex_pdb = self.FILES.prefix + 'MUT_COM.pdb'
        self.mutant_receptor_pdb = self.FILES.prefix + 'MUT_REC.pdb'
        self.mutant_ligand_pdb = self.FILES.prefix + 'MUT_LIG.pdb'
        self.mutant_complex_pdb_fixed = self.FILES.prefix + 'MUT_COM_FIXED.pdb'
        self.mutant_receptor_pdb_fixed = self.FILES.prefix + 'MUT_REC_FIXED.pdb'
        self.mutant_ligand_pdb_fixed = self.FILES.prefix + 'MUT_LIG_FIXED.pdb'
        # self.default_ff = 'leaprc.protein.ff14SB'

        checkff()
        self.getPDBfromTpr()
        self.checkPDB()

    def getPDBfromTpr(self):
        """
        Generate PDB file to generate topology
        :return:
        """
        logging.info('Get PDB files from structures files...')
        gmx = self.external_progs['gmx'].full_path
        # check if GROMACS 4.x exists
        make_ndx = [self.external_progs['make_ndx'].full_path]
        trjconv = [self.external_progs['trjconv'].full_path]
        editconf = [self.external_progs['editconf'].full_path]
        if gmx:
            make_ndx = [gmx, 'make_ndx']
            trjconv = [gmx, 'trjconv']
            editconf = [gmx, 'editconf']

        # wt complex
        # make index for extract pdb structure
        rec_group, lig_group = self.FILES.complex_groups

        logging.info('Making gmx_MMPBSA index for complex...')
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        make_ndx_echo_args = ['echo', 'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | {l}\n q\n'.format(
            r=rec_group, l=lig_group)]
        c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)
        # FIXME: overwrite the user index file???
        com_ndx = self.FILES.prefix + 'COM_index.ndx'
        make_ndx_args = make_ndx + ['-n', self.FILES.complex_index, '-o', com_ndx]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + (' '.join(make_ndx_echo_args).replace('\n', '\\n')) + ' | ' + ' '.join(
                make_ndx_args))
        c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=self.log, stderr=self.log)
        if c2.wait():  # if it quits with return code != 0
            raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(make_ndx), self.FILES.complex_index))
        self.FILES.complex_index = com_ndx

        logging.info('Normal Complex: Saving group {}_{} in {} file as '
                     '{}'.format(rec_group, lig_group, self.FILES.complex_index, self.complex_pdb))
        editconf_echo_args = ['echo', 'GMXMMPBSA_REC_GMXMMPBSA_LIG']
        c3 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
        # we get only first trajectory to extract a pdb file and make amber topology for complex
        editconf_args = editconf + ['-f', self.FILES.complex_tpr, '-o', self.complex_pdb, '-n',
                                    self.FILES.complex_index]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
        c4 = subprocess.Popen(editconf_args, stdin=c3.stdout, stdout=self.log, stderr=self.log)
        if c4.wait():  # if it quits with return code != 0
            raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.complex_tpr))

        # clear trajectory
        if self.INPUT['solvated_trajectory']:
            logging.info('Cleaning normal complex trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.complex_trajs)):
                trjconv_echo_args = ['echo', 'GMXMMPBSA_REC_GMXMMPBSA_LIG']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr,
                                          '-o', 'COM_traj_{}.xtc'.format(i), '-n', self.FILES.complex_index]
                if self.INPUT['debug_printlevel']:
                    logging.info('Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args,  # FIXME: start and end frames???
                                      stdin=c5.stdout, stdout=self.log, stderr=self.log)
                if c6.wait():  # if it quits with return code != 0
                    raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(trjconv), self.FILES.complex_tpr))
                new_trajs.append('COM_traj_{}.xtc'.format(i))
            self.FILES.complex_trajs = new_trajs

        # Put receptor and ligand (explicitly defined) to avoid overwrite them
        # check if ligand is not protein. In any case, non-protein ligand always most be processed
        if self.FILES.ligand_mol2:
            logging.info(f'Generating ligand parameters from {self.FILES.ligand_mol2} file...')
            lig_name = os.path.splitext(os.path.split(self.FILES.ligand_mol2)[1])[0]
            self.ligand_frcmod = self.FILES.prefix + lig_name + '.frcmod'
            # run parmchk2
            parmchk2 = self.external_progs['parmchk2'].full_path
            parmchk2_args = [parmchk2, '-i', self.FILES.ligand_mol2, '-f', 'mol2', '-o', self.ligand_frcmod]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + ' '.join(parmchk2_args))
            l3 = subprocess.Popen(parmchk2_args, stdout=self.log, stderr=self.log)
            if l3.wait():
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (parmchk2, self.FILES.ligand_mol2))

        # make a temp receptor pdb (even when stability) if decomp to get correct receptor residues from complex. This
        # avoid get multiples molecules from complex.split()
        if self.INPUT['decomprun'] and self.print_residues:
            if self.FILES.stability:
                self.use_temp = True
                logging.warning('When decomp is defined, we generate a receptor file in order to extract interface '
                                'residues')
                rec_echo_args = ['echo', '{}'.format(rec_group)]
                cp1 = subprocess.Popen(rec_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file to generate amber topology
                editconf_args = editconf + ['-f', self.FILES.complex_tpr, '-o', 'rec_temp.pdb', '-n',
                                            self.FILES.complex_index]
                if self.INPUT['debug_printlevel']:
                    logging.info('Running command: ' + (' '.join(rec_echo_args)) + ' | ' + ' '.join(editconf_args))
                cp2 = subprocess.Popen(editconf_args, stdin=cp1.stdout, stdout=self.log, stderr=self.log)
                if cp2.wait():  # if it quits with return code != 0
                    raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.complex_tpr))

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
                self.FILES.receptor_group, self.FILES.receptor_index, self.receptor_pdb))
            editconf_echo_args = ['echo', '{}'.format(self.FILES.receptor_group)]
            p1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file for make amber topology
            editconf_args = editconf + ['-f', self.FILES.receptor_tpr, '-o', self.receptor_pdb, '-n',
                                        self.FILES.receptor_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cp2 = subprocess.Popen(editconf_args, stdin=p1.stdout, stdout=self.log, stderr=self.log)
            if cp2.wait():  # if it quits with return code != 0
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.receptor_tpr))
            # clear trajectory
            if self.INPUT['solvated_trajectory']:
                logging.info('Clear normal receptor trajectories...')
                new_trajs = []
                for i in range(len(self.FILES.receptor_trajs)):
                    trjconv_echo_args = ['echo', '{}'.format(self.FILES.receptor_group)]
                    c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                    # we get only first trajectory to extract a pdb file and make amber topology for complex
                    trjconv_args = trjconv + ['-f', self.FILES.receptor_trajs[0], '-s', self.FILES.receptor_tpr,
                                              '-o', 'REC_traj_{}.xtc'.format(i), '-n',
                                              self.FILES.receptor_index]
                    if self.INPUT['debug_printlevel']:
                        logging.info('Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(
                            trjconv_args))
                    c6 = subprocess.Popen(trjconv_args,  # FIXME: start and end frames???
                                          stdin=c5.stdout, stdout=self.log, stderr=self.log)
                    if c6.wait():  # if it quits with return code != 0
                        raise GMXMMPBSA_ERROR(
                            '%s failed when querying %s' % (' '.join(trjconv), self.FILES.receptor_tpr))
                    new_trajs.append('REC_traj_{}.xtc'.format(i))
                self.FILES.receptor_trajs = new_trajs
        else:
            logging.info('No receptor structure file was defined. Using ST approach...')
            logging.info('Using receptor structure from complex to generate AMBER topology')
            # wt complex receptor
            logging.info('Normal Complex: Saving group {} in {} file as {}'.format(
                rec_group, self.FILES.complex_index, self.receptor_pdb))
            editconf_echo_args = ['echo', '{}'.format(rec_group)]
            cp1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file for make amber topology
            editconf_args = editconf + ['-f', self.FILES.complex_tpr, '-o', self.receptor_pdb, '-n',
                                        self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cp2 = subprocess.Popen(editconf_args, stdin=cp1.stdout, stdout=self.log, stderr=self.log)
            if cp2.wait():  # if it quits with return code != 0
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.complex_tpr))

        # ligand
        # # check consistence
        if self.FILES.ligand_tpr:  # ligand is protein
            # FIXME: if ligand is a zwitterionic aa fail
            logging.info('A ligand structure file was defined. Using MT approach...')
            logging.info('Normal Ligand: Saving group {} in {} file as {}'.format(
                self.FILES.ligand_group, self.FILES.ligand_index, self.ligand_pdb))
            # wt ligand
            editconf_echo_args = ['echo', '{}'.format(self.FILES.ligand_group)]
            l1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory for extract a pdb file for make amber topology
            editconf_args = editconf + ['-f', self.FILES.ligand_tpr, '-o', self.ligand_pdb, '-n',
                                        self.FILES.ligand_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            l2 = subprocess.Popen(editconf_args, stdin=l1.stdout, stdout=self.log, stderr=self.log)
            if l2.wait():  # if it quits with return code != 0
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.ligand_tpr))

            # clear trajectory
            if self.INPUT['solvated_trajectory']:
                logging.info('Clear normal ligand trajectories...')
                new_trajs = []
                for i in range(len(self.FILES.ligand_trajs)):
                    trjconv_echo_args = ['echo', '{}'.format(self.FILES.ligand_group)]
                    c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                    # we get only first trajectory to extract a pdb file and make amber topology for complex
                    trjconv_args = trjconv + ['-f', self.FILES.ligand_trajs[0], '-s', self.FILES.ligand_tpr,
                                              '-o', 'LIG_traj_{}.xtc'.format(i), '-n', self.FILES.ligand_index]
                    if self.INPUT['debug_printlevel']:
                        logging.info(
                            'Running command: ' + (' '.join(trjconv_echo_args)) + ' | ' + ' '.join(trjconv_args))
                    c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=self.log, stderr=self.log)
                    if c6.wait():  # if it quits with return code != 0
                        raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(trjconv), self.FILES.ligand_tpr))
                    new_trajs.append('LIG_traj_{}.xtc'.format(i))
                self.FILES.ligand_trajs = new_trajs
        else:
            # wt complex ligand
            logging.info('No ligand structure file was defined. Using ST approach...')
            logging.info('Using ligand structure from complex to generate AMBER topology')
            logging.info('Normal ligand: Saving group {} in {} file as {}'.format(lig_group, self.FILES.complex_index,
                                                                                                  self.ligand_pdb))
            editconf_echo_args = ['echo', '{}'.format(lig_group)]
            cl1 = subprocess.Popen(editconf_echo_args, stdout=subprocess.PIPE)
            # we get only  first trajectory to extract a pdb file for make amber topology
            editconf_args = editconf + ['-f', self.FILES.complex_tpr, '-o', self.ligand_pdb, '-n',
                                        self.FILES.complex_index]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + (' '.join(editconf_echo_args)) + ' | ' + ' '.join(editconf_args))
            cl2 = subprocess.Popen(editconf_args, stdin=cl1.stdout, stdout=self.log, stderr=self.log)
            if cl2.wait():  # if it quits with return code != 0
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(editconf), self.FILES.complex_tpr))

    def checkPDB(self):
        """
        Generate parmed structure object for complex, receptor and ligand ( if it is protein-like)

        1 - Rename HIS
        2 - Rename CYS
        3 - Delete H
        4 - Rename oxygen in termini from GROMACS to AMBER name
          - Rename CD in ILE from GROMACS to AMBER name
        5 - Save
        :return:
        """
        logging.info('Generating AMBER Compatible PDB Files...')

        self.remove_MODEL(self.complex_pdb)
        self.remove_MODEL(self.receptor_pdb)
        self.remove_MODEL(self.ligand_pdb)
        # try:
        #     self.complex_str = parmed.read_PDB(self.complex_pdb)  # can always be initialized
        # except: # when structure file has no chains ids
        #     remove_MODEL(self.complex_pdb)
        self.complex_str = parmed.read_PDB(self.complex_pdb)  # can always be initialized
        self.receptor_str = parmed.read_PDB(self.receptor_pdb)
        self.ligand_str = parmed.read_PDB(self.ligand_pdb)
        if self.FILES.reference_structure:
            self.ref_str = parmed.read_PDB(self.FILES.reference_structure)

        self.fix_chains_IDs(self.complex_str, self.receptor_str, self.ligand_str, self.ref_str)

        # fix receptor structure
        self.properHIS(self.receptor_str)
        self.properCYS(self.receptor_str)
        self.properAspGluLys(self.receptor_str)
        self.fix_H_ATOMS(self.receptor_str)
        # For some reason removing the hydrogens returns the hydrogen-bound atoms to their original names. This is
        # problematic with ILE switching from CD to CD1. parmed bug?
        self.receptor_str.strip('@/H')
        self.properATOMS(self.receptor_str)

        # check if rec contain ions (metals)
        self.rec_ions = self.receptor_str[:, ions, :]
        self.rec_ions_after = False
        if self.rec_ions.atoms:
            # fix atom number, avoid core dump in tleap
            i = 1
            for at in self.rec_ions.atoms:
                at.number = i
                i += 1
            # check ions location
            count = 0
            for res in self.receptor_str.residues:
                if res.number != self.complex_str.residues[count].number:
                    self.rec_ions_after = True
                    break
                count += 1
            if self.rec_ions_after:
                self.rec_ions.save(self.rec_ions_pdb, 'pdb', True, renumber=False)
                # if exists any ions then strip them
                self.receptor_str.strip(f':{",".join(ions)}')
                self.rec_str_ions = True
        self.receptor_str.save(self.receptor_pdb_fixed, 'pdb', True, renumber=False)

        # fix ligand structure if is protein
        self.properHIS(self.ligand_str)
        self.properCYS(self.ligand_str)
        self.properAspGluLys(self.ligand_str)
        self.fix_H_ATOMS(self.ligand_str)
        self.ligand_str.strip('@/H')
        self.properATOMS(self.ligand_str)

        # check if lig contain ions (metals)
        self.lig_ions = self.ligand_str[:, ions, :]
        if self.lig_ions.atoms:
            # fix atom number, avoid core dump in tleap
            i = 1
            for at in self.lig_ions.atoms:
                at.number = i
                i += 1
            self.lig_ions.save(self.lig_ions_pdb, 'pdb', True, renumber=False)
            # if exists any ions then strip them
            self.ligand_str.strip(f':{",".join(ions)}')
            self.lig_str_ions = True
        self.ligand_str.save(self.ligand_pdb_fixed, 'pdb', True, renumber=False)

        if self.INPUT['alarun']:
            logging.info('Building Mutant receptor...')
            if self.INPUT['mutant'].lower() in ['rec', 'receptor']:
                self.mutant_receptor_str = parmed.read_PDB(self.receptor_pdb_fixed)
                # fix mutant receptor structure
                self.mutatexala(self.mutant_receptor_str)
                self.mutant_receptor_str.save(self.mutant_receptor_pdb_fixed, 'pdb', True, renumber=False)
            else:
                logging.info('Building Mutant ligand...')
                if self.FILES.ligand_mol2:
                    raise GMXMMPBSA_ERROR('Mutation is only possible if the ligand is protein-like')
                self.mutant_ligand_str = parmed.read_PDB(self.ligand_pdb_fixed)
                self.mutatexala(self.mutant_ligand_str)
                self.mutant_ligand_str.save(self.mutant_ligand_pdb_fixed, 'pdb', True, renumber=False)

        # Get residue form receptor-ligand interface
        if self.print_residues:
            if self.use_temp:
                temp_str = parmed.read_PDB('rec_temp.pdb')
                rec_resnum = len(temp_str.residues)
            else:
                rec_resnum = len(self.receptor_str.residues)
            res_list = []
            res_ndx = 1
            for rres in self.complex_str.residues[:rec_resnum]:  # iterate over receptor residues
                lres_ndx = rec_resnum + 1
                for lres in self.complex_str.residues[rec_resnum:]:  # iterate over ligand residues
                    for rat in rres.atoms:
                        rat_coor = [rat.xx, rat.xy, rat.xz]
                        for lat in lres.atoms:
                            lat_coor = [lat.xx, lat.xy, lat.xz]
                            if dist(rat_coor, lat_coor) <= self.within:
                                if res_ndx not in res_list:
                                    res_list.append(res_ndx)
                                if lres_ndx not in res_list:
                                    res_list.append(lres_ndx)
                                break
                    lres_ndx += 1
                res_ndx += 1
            res_list.sort()
            self.INPUT['print_res'] = ','.join([str(x) for x in res_list])

    def mutatexala(self, structure):
        idx = 0
        found = False
        if not self.INPUT['mutant_res']:
            raise GMXMMPBSA_ERROR("No residue for mutation was defined")
        chain, resnum = self.INPUT['mutant_res'].split(':')

        if not chain or not resnum:
            raise GMXMMPBSA_ERROR("No residue was defined")
        for res in structure.residues:
            if res.number == int(resnum) and res.chain == chain:
                found = True
                break
            idx += 1
        if found:
            structure.residues[idx].name = 'ALA'
            excluded_mask = ':{} &!@CB,C,CA,N,O'.format(idx + 1)
            structure.strip(excluded_mask)
        else:
            raise GMXMMPBSA_ERROR('Residue {}:{} not found'.format(chain, resnum))

    def fix_chains_IDs(self, com_str, rec_str, lig_str, ref_str=None):
        if ref_str:
            if len(ref_str.residues) < len(rec_str.residues):
                raise ('error')
            chains = []
            for res in ref_str.residues:
                if res.chain not in chains:
                    chains.append(res.chain)
            c = 0
            for res in com_str.residues:
                if c >= len(ref_str.residues):
                    l = c - len(rec_str.residues)
                    logging.warning('Differences in the number of residues between the reference and complex '
                                    'structure. Arbitrarily assigning chain ID to the ligand in the complex')
                    res.chain = chains_letters[chains_letters.index(chains[-1]) + 1]
                    lig_str.residues[l].chain = chains_letters[chains_letters.index(chains[-1]) + 1]
                elif c < len(rec_str.residues):
                    res.chain = ref_str.residues[c].chain
                    rec_str.residues[c].chain = ref_str.residues[c].chain
                else:
                    res.chain = ref_str.residues[c].chain
                    lig_str.residues[c - len(rec_str.residues)].chain = ref_str.residues[c].chain
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
                logging.warning('No reference structure was found and the gro format was used in the complex '
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
                        if c < len(rec_str.residues):
                            rec_str.residues[c].chain = curr_chain_id
                        else:
                            lig_str.residues[c - len(rec_str.residues)].chain = curr_chain_id
                        if curr_chain_id not in chains_ids:
                            chains_ids.append(curr_chain_id)
                    else:
                        if res.chain != curr_chain_id:
                            res.chain = curr_chain_id
                            if c < len(rec_str.residues):
                                rec_str.residues[c].chain = curr_chain_id
                            else:
                                lig_str.residues[c - len(rec_str.residues)].chain = curr_chain_id
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
                        if c < len(rec_str.residues):
                            rec_str.residues[c].chain = curr_chain_id
                        else:
                            lig_str.residues[c - len(rec_str.residues)].chain = curr_chain_id
                        if res.chain not in chains_ids:
                            chains_ids.append(res.chain)
                    elif chain_by_ter:
                        chain_by_ter = False
                    elif chain_by_num:
                        chain_by_num = False
                        curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                        res.chain = curr_chain_id
                        if c < len(rec_str.residues):
                            rec_str.residues[c].chain = curr_chain_id
                        else:
                            lig_str.residues[c - len(rec_str.residues)].chain = curr_chain_id
                        if res.chain not in chains_ids:
                            chains_ids.append(res.chain)
                    for atm in res.atoms:
                        if atm.name == 'OC2':  # only for protein
                            res.ter = True
                            chain_by_ter = True
                    if parmed.residue.RNAResidue.has(res.name) or parmed.residue.DNAResidue.has(res.name):
                        has_nucl += 1
                    if has_nucl == 1:
                        logging.warning('This structure contains nucleotides. We recommend that you use the reference '
                                        'structure')
                    previous_res_number = res.number
                    c += 1

    @staticmethod
    def remove_MODEL(pdb_file):
        try:
            with open(pdb_file) as fo:
                fo = fo.readlines()
                for line in fo:
                    if 'MODEL' in line or 'ENDMDL' in line:
                        fo.remove(line)
            with open(pdb_file, 'w') as fw:
                for x in fo:
                    fw.write(x)
            return True
        except IOError as e:
            GMXMMPBSA_ERROR('', str(e))

    @staticmethod
    def fix_H_ATOMS(structure):
        """
        GROMACS 4.x save the pdb without atom element column, so parmed does not recognize some H atoms. Parmed assigns
        0 to the atomic number of these atoms. In order to correctly eliminate hydrogens, it is necessary to assign the
        atomic number.
        """
        for residue in structure.residues:
            for atom in residue.atoms:
                if 'H' in atom.name and atom.atomic_number == 0:
                    atom.atomic_number = 1

    @staticmethod
    def properATOMS(structure):
        """
        Rename oxygen in termini from GROMACS to AMBER name
        OC1 -> 'O  '
        OC2 -> OXT -> set this res as terminal in parmed

        Rename CD in ILE from GROMACS to AMBER name
        CD   ILE -> CD1 ILE
        :return:
        """
        for residue in structure.residues:
            if residue.name == 'ILE':
                for atom in residue.atoms:
                    if atom.name == 'CD':
                        atom.name = 'CD1'

            for atom in residue.atoms:
                if atom.name == 'OC1':
                    atom.name = 'O  '
                elif atom.name == 'OC2':
                    atom.name = 'OXT'
                    residue.ter = True  # parmed terminal

    @staticmethod
    def properAspGluLys(structure):
        """
        Proper name for residues Asp, Glu and Lys
        """
        for residue in structure.residues:
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

    @staticmethod
    def properHIS(structure):
        """
        Compatible amber name for Histidines from protonation state
        """
        his = ['HIS', 'HIE', 'HID', 'HIP']

        for residue in structure.residues:
            if residue.name in his:
                atoms = [atom.name for atom in residue.atoms if atom.atomic_number == 1]
                if 'HD1' in atoms and 'HE2' in atoms:
                    residue.name = 'HIP'
                elif 'HD1' in atoms:
                    residue.name = 'HID'
                elif 'HE2' in atoms:
                    residue.name = 'HIE'

    @staticmethod
    def properCYS(structure):
        """
        Rename the cys in disulfide bond
        :return:
        """
        cys_name = ['CYS', 'CYX', 'CYM']
        allcys = [residue for residue in structure.residues if residue.name in cys_name]
        # print(llcys
        xcys = []
        for residue in allcys:
            for atom in residue.atoms:
                if 'SG' in atom.name:
                    for bondedatm in atom.bond_partners:
                        # exclude CB
                        if bondedatm.residue == residue:
                            continue
                        else:
                            # check if is bonded to cys residue
                            # TODO: Check if bonded atom is SG. Is really necessary?
                            if bondedatm.residue.name in cys_name:
                                if residue not in xcys:
                                    xcys.append(residue)
                                if bondedatm.residue not in xcys:
                                    xcys.append(bondedatm.residue)
        for cys in xcys:
            cys.name = 'CYX'

    def makeToptleap(self):
        logging.info('Building Tleap input files...')
        with open(self.FILES.prefix + 'leap.in', 'w') as tif:
            tif.write('source {}\n'.format(self.INPUT['protein_forcefield']))
            tif.write('source {}\n'.format(self.INPUT['ligand_forcefield']))
            if self.rec_ions.atoms or self.lig_ions.atoms:
                tif.write('loadOff atomic_ions.lib\n')
                tif.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['ions_parameters']]))
            tif.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['PBRadii']]))

            tif.write('REC = loadpdb {}\n'.format(self.receptor_pdb_fixed))
            if self.rec_str_ions:
                tif.write('R_IONS = loadpdb {}\n'.format(self.rec_ions_pdb))
                tif.write('REC_IONS = combine { REC R_IONS}\n')
                tif.write('saveamberparm REC_IONS {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop,
                                                                              p=self.FILES.prefix))
            else:
                tif.write('saveamberparm REC {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop, p=self.FILES.prefix))

            # check if ligand is not protein and always load
            if self.FILES.ligand_mol2:
                tif.write('LIG = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                tif.write('check LIG\n')
                tif.write('loadamberparams {}\n'.format(self.ligand_frcmod))
                tif.write('saveamberparm LIG {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop, p=self.FILES.prefix))
            else:
                tif.write('LIG = loadpdb {}\n'.format(self.ligand_pdb_fixed))
                if self.lig_str_ions:
                    tif.write('L_IONS = loadpdb {}\n'.format(self.lig_ions_pdb))
                    tif.write('LIG_IONS = combine { LIG L_IONS}\n')
                    tif.write('saveamberparm LIG_IONS {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop,
                                                                                  p=self.FILES.prefix))
                else:
                    tif.write('saveamberparm LIG {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop, p=self.FILES.prefix))

            com_string = 'complex = combine { REC '
            com_string += 'LIG '
            if self.rec_str_ions and self.rec_ions_after:
                com_string += 'R_IONS '
            # elif self.rec_str_ions:
            # com_string += 'R_IONS '
            # com_string += 'LIG '
            if self.lig_str_ions:
                com_string += 'L_IONS '
            com_string += '}\n'
            tif.write(com_string)
            tif.write('saveamberparm complex {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            tif.write('quit')

        tleap = self.external_progs['tleap'].full_path
        tleap_args = [tleap, '-f', '{}'.format(self.FILES.prefix + 'leap.in')]
        if self.INPUT['debug_printlevel']:
            logging.info('Running command: ' + ' '.join(tleap_args))
        p = subprocess.Popen(tleap_args, stdout=self.log, stderr=self.log)
        if p.wait():
            raise GMXMMPBSA_ERROR('%s failed when querying %s' % (tleap, self.FILES.prefix + 'leap.in'))

        if self.INPUT['alarun']:
            with open(self.FILES.prefix + 'mut_leap.in', 'w') as mtif:
                mtif.write('source {}\n'.format(self.INPUT['protein_forcefield']))
                mtif.write('source {}\n'.format(self.INPUT['ligand_forcefield']))
                if self.rec_str_ions or self.lig_str_ions:
                    mtif.write('loadOff atomic_ions.lib\n')
                    mtif.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['ions_parameters']]))
                mtif.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['PBRadii']]))
                # check if ligand is not protein and always load
                if self.FILES.ligand_mol2:
                    mtif.write('LIG = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                    self.mutant_ligand_pmrtop = None
                    mtif.write('check LIG\n')
                    mtif.write('loadamberparams {}\n'.format(self.ligand_frcmod))
                    mtif.write('saveamberparm LIG {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop,
                                                                              p=self.FILES.prefix))
                    mtif.write('REC = loadpdb {}\n'.format(self.mutant_receptor_pdb_fixed))
                    if self.rec_str_ions:
                        mtif.write('R_IONS = loadpdb {}\n'.format(self.rec_ions_pdb))
                        mtif.write('REC_IONS = combine { REC R_IONS}\n')
                        mtif.write('saveamberparm REC_IONS {t} {p}MUT_REC.inpcrd\n'.format(
                            t=self.mutant_receptor_pmrtop, p=self.FILES.prefix))
                    else:
                        mtif.write('saveamberparm REC {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
                                                                                      p=self.FILES.prefix))
                else:
                    if self.INPUT['mutant'].lower() in ['rec', 'receptor']:
                        mtif.write('REC = loadpdb {}\n'.format(self.mutant_receptor_pdb_fixed))
                        self.mutant_ligand_pmrtop = None
                        if self.rec_str_ions:
                            mtif.write('R_IONS = loadpdb {}\n'.format(self.rec_ions_pdb))
                            mtif.write('REC_IONS = combine { REC R_IONS}\n')
                            mtif.write('saveamberparm REC_IONS {t} {p}REC.inpcrd\n'.format(
                                t=self.mutant_receptor_pmrtop, p=self.FILES.prefix))
                        else:
                            mtif.write('saveamberparm REC {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
                                                                                          p=self.FILES.prefix))
                        mtif.write('LIG = loadpdb {}\n'.format(self.ligand_pdb_fixed))
                        mtif.write('saveamberparm LIG {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop,
                                                                                  p=self.FILES.prefix))
                    else:
                        mtif.write('LIG = loadpdb {}\n'.format(self.mutant_ligand_pdb_fixed))
                        self.mutant_receptor_pmrtop = None
                        if self.lig_str_ions:
                            mtif.write('L_IONS = loadpdb {}\n'.format(self.lig_ions_pdb))
                            mtif.write('LIG_IONS = combine { LIG L_IONS}\n')
                            mtif.write('saveamberparm LIG_IONS {t} {p}LIG.inpcrd\n'.format(t=self.mutant_ligand_pmrtop,
                                                                                           p=self.FILES.prefix))
                        else:
                            mtif.write('saveamberparm LIG {t} {p}MUT_LIG.inpcrd\n'.format(t=self.mutant_ligand_pmrtop,
                                                                                          p=self.FILES.prefix))
                        mtif.write('REC = loadpdb {}\n'.format(self.receptor_pdb_fixed))
                        if self.rec_str_ions:
                            mtif.write('R_IONS = loadpdb {}\n'.format(self.rec_ions_pdb))
                            mtif.write('REC_IONS = combine { REC R_IONS}\n')
                            mtif.write('saveamberparm REC_IONS {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop,
                                                                                           p=self.FILES.prefix))
                        else:
                            mtif.write('saveamberparm REC {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop,
                                                                                      p=self.FILES.prefix))

                    mut_com_string = 'mut_com = combine { REC LIG '
                    if self.rec_str_ions:
                        mut_com_string += 'R_IONS '
                    if self.lig_str_ions:
                        mut_com_string += 'L_IONS '
                    mut_com_string += '}\n'
                    mtif.write(mut_com_string)
                    mtif.write('saveamberparm mut_com {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                                      p=self.FILES.prefix))
                mtif.write('quit')
            tleap_args = [tleap, '-f', '{}'.format(self.FILES.prefix + 'mut_leap.in')]
            if self.INPUT['debug_printlevel']:
                logging.info('Running command: ' + ' '.join(tleap_args))
            p1 = subprocess.Popen(tleap_args, stdout=self.log, stderr=self.log)
            if p1.wait():
                raise GMXMMPBSA_ERROR('%s failed when querying %s' % (tleap, self.FILES.prefix + 'mut_leap.in'))
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)
