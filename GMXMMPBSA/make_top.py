"""
Generate Amber topology files from GROMACS files
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

import os
import platform
import textwrap

import parmed
from GMXMMPBSA.exceptions import *
from GMXMMPBSA.topology_preprocess import GromacsTopologyPreprocessor
from GMXMMPBSA.utils import (selector, get_dist, list2range, res2map, get_indexes, log_subprocess_output, check_str,
                             eq_strs, get_index_groups, reconcile_qm_charges, topology_mismatch_error, Residue,
                             residue_names_match)
from GMXMMPBSA.alamdcrd import _scaledistance
from GMXMMPBSA.make_trajs import warn_concatenated_complex_trajectories
from GMXMMPBSA.radii import source_force_field_family
import subprocess
from pathlib import Path
import logging
import string
import re
from parmed.tools.changeradii import ChRad

if platform.system() == "Darwin":
    echo_command = ['echo']
else:
    echo_command = ['echo', '-e']

chains_letters = list(string.ascii_uppercase)
his = ['HIS', 'HIE', 'HID', 'HIP', 'HSP', 'HSD', 'HSE']
cys_name = ['CYS', 'CYX', 'CYM']
lys = ['LYS', 'LYN']
asp = ['ASP', 'ASH']
glu = ['GLU', 'GLH']
positive_aa = ['LYS', 'ARG', 'HIP']
negative_aa = ['GLU', 'ASP']
nonpolar_aa = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
polar_aa = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']

PBRadii = {1: 'bondi', 2: 'mbondi', 3: 'mbondi2', 4: 'mbondi3', 5: 'mbondi_pb2', 6: 'mbondi_pb3', 7: 'charmm_radii'}
GB_RECOMMENDED_RADII = {1: 'mbondi', 2: 'mbondi2', 5: 'mbondi2', 7: 'bondi', 8: 'mbondi3'}

# ions_para_files = {1: 'frcmod.ions234lm_126_tip3p', 2: 'frcmod.ions234lm_iod_tip4pew', 3: 'frcmod.ions234lm_iod_spce',
#                    4: 'frcmod.ions234lm_hfe_spce', 5: 'frcmod.ions234lm_126_tip4pew', 6: 'frcmod.ions234lm_126_spce',
#                    7: 'frcmod.ions234lm_1264_tip4pew', 8: 'frcmod.ions234lm_1264_tip3p',
#                    9: 'frcmod.ions234lm_1264_spce', 10: 'frcmod.ions234lm_iod_tip3p',
#                    11: 'frcmod.ions234lm_hfe_tip4pew', 12: 'frcmod.ions234lm_hfe_tip3p}'}

ions_para_files = {1: 'frcmod.ions234lm_126_tip3p', 2: 'frcmod.ions234lm_126_spce', 3: 'frcmod.ions234lm_126_tip4pew',
                   4: 'frcmod.ions234lm_hfe_tip3p', 5: 'frcmod.ions234lm_hfe_spce', 6: 'frcmod.ions234lm_hfe_tip4pew',
                   7: 'frcmod.ions234lm_iod_tip3p', 8: 'frcmod.ions234lm_iod_spce', 9: 'frcmod.ions234lm_iod_tip4pew',
                   10: 'frcmod.ionslm_126_opc', 11: 'frcmod.ionslm_hfe_opc', 12: 'frcmod.ionslm_iod_opc',
                   13: 'frcmod.ions1lm_126_tip3p', 14: 'frcmod.ions1lm_126_spce', 15: 'frcmod.ions1lm_126_tip4pew',
                   16: 'frcmod.ions1lm_iod'}

ions = ["AG", "AL", "Ag", "BA", "BR", "Be", "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "CU1", "Ce", "Cl-", "Cr",
        "Dy", "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", "HE+", "HG", "HZ+", "Hf", "IN", "IOD", "K", "K+",
        "LA", "LI", "LU", "MG", "MN", "NA", "NH4", "NI", "Na+", "Nd", "PB", "PD", "PR", "PT", "Pu", "RB", "Ra", "SM",
        "SR", "Sm", "Sn", "TB", "TL", "Th", "Tl", "Tm", "U4+", "V2+", "Y", "YB2", "ZN", "Zr"]

water_residues = [
    'SOL', 'WAT',
    'TIP3P', 'TIP3', 'TP3', 'TIPS3P', 'TIP3o',
    'TIP4P', 'TIP4PEW', 'T4E', 'TIP4PD',
    'TIP5P',
    'SPC', 'SPC/E', 'SPCE',
    'OPC'
]

solvent_ion_residues = [
    'NA', 'CL', 'K',
    'SOD', 'Na+', 'CLA', 'Cl-', 'POT', 'K+',
    *water_residues
]
explicit_water_ion_mask = ':NA,CL,K,SOD,Na+,CLA,Cl-,POT,K+'
explicit_water_solvent_mask = ':' + ','.join(water_residues)
explicit_water_reference_exclusion_mask = ':' + ','.join(solvent_ion_residues)
explicit_water_group_names = ['SOLV', 'SOL', 'Water', 'WAT', *water_residues]


class CheckMakeTop:
    def __init__(self, FILES, INPUT, external_programs):
        self.FILES = FILES
        self.INPUT = INPUT
        self.external_progs = external_programs
        self.use_temp = False
        self.com_mut_index = None
        self.com_mut_indices = []
        self.part_indices = []

        # Define Gromacs executable
        self.make_ndx = self.external_progs['make_ndx']
        self.trjconv = self.external_progs['trjconv']
        self.editconf = self.external_progs['editconf']

        self.cys_bonds = {'COM': [], 'REC': [], 'LIG': []}

        self.ref_str = None

        self.ligand_tpr = None
        self.ligand_mol2 = None

        self.rec_str_ions = False
        self.lig_str_ions = False
        self.explicit_waters = self.INPUT['general']['explicit_waters']
        self.explicit_waters_mask = self.INPUT['general']['explicit_waters_mask']
        self.explicit_waters_group = self.INPUT['general']['explicit_waters_group']
        self.explicit_waters_extra_points = self.INPUT['general']['explicit_waters_extra_points'].lower()
        self.explicit_water_prmtop = None
        self.explicit_water_range = ''
        self.explicit_water_extra_point_mask = ''

        # create the * prmtop variables for compatibility with the original code
        self.complex_pmrtop = 'COM.prmtop'
        self.receptor_pmrtop = 'REC.prmtop'
        self.ligand_pmrtop = 'LIG.prmtop'

        self.mutant_complex_pmrtop = 'MUT_COM.prmtop'
        self.mutant_receptor_pmrtop = 'MUT_REC.prmtop'
        self.mutant_ligand_pmrtop = 'MUT_LIG.prmtop'

        self.complex_str_file = f'{self.FILES.prefix}COM.pdb'
        self.receptor_str_file = f'{self.FILES.prefix}REC.pdb'
        self.ligand_str_file = f'{self.FILES.prefix}LIG.pdb'

        self.checkFiles()

    def checkFiles(self):
        if (not self.FILES.complex_tpr or not self.FILES.complex_index or
                not self.FILES.complex_trajs or not self.FILES.complex_groups):
            GMXMMPBSA_ERROR(
                'You must define the complex structure (-cs), index (-ci), trajectory (-ct), and groups (-cg).'
            )
        if not self.FILES.complex_top:
            GMXMMPBSA_ERROR(
                'A GROMACS complex topology (-cp) is required. Structure-only tleap rebuilds are no longer '
                'supported; convert parameters from the GROMACS topology used in the MD. For a small-molecule '
                'ligand, include it in that topology (for example via ACPYPE) rather than relying on -lm alone.'
            )
        if (self.FILES.receptor_tpr or self.FILES.receptor_trajs) and not self.FILES.receptor_top:
            GMXMMPBSA_ERROR(
                'A receptor topology (-rp) is required when unbound receptor structure (-rs) or trajectories '
                '(-rt) are defined (multiple-trajectory approach).'
            )
        if (self.FILES.ligand_tpr or self.FILES.ligand_trajs) and not self.FILES.ligand_top:
            GMXMMPBSA_ERROR(
                'A ligand topology (-lp) is required when unbound ligand structure (-ls) or trajectories '
                '(-lt) are defined (multiple-trajectory approach).'
            )

    def buildTopology(self):
        """
        :return: complex, receptor, ligand topologies and their mutants
        """
        self.gmx2pdb()
        # dASA needs the full solvated AMBER topology generated below. The
        # distance and explicit-mask selectors can be resolved immediately.
        if self.explicit_waters_mask.strip().lower() != 'dasa':
            self._resolve_explicit_waters_mask()
        # GROMACS calculations always convert the user topology (-cp). The legacy
        # structure→PDB→tleap path is removed to avoid termini/FF mismatches.
        tops = self.gmxtop2prmtop()
        if self.explicit_waters_mask.strip().lower() == 'dasa':
            self._resolve_explicit_waters_mask()

        if self.INPUT['decomp']['decomprun']:
            explicit_water_res = []
            if self.explicit_waters and self.INPUT['decomp']['print_res'] == 'all':
                explicit_water_res = self._ensure_explicit_water_residues_mapped()
            decomp_res = self.get_selected_residues(self.INPUT['decomp']['print_res'])
            if 'within' in self.INPUT['decomp']['print_res']:
                if len(decomp_res) < 2:
                    logging.info(f"Number of decomp residues to print using "
                                 f"print_res = '{self.INPUT['decomp']['print_res']}' < 2; expanding the cutoff")
                    logging.info(
                        'Increasing cutoff value by 0.1 until number of decomp residues to print >= 2'
                    )
                    cutoff = float(self.INPUT['decomp']['print_res'].split()[1])
                    it = 0
                    while len(decomp_res) < 2:
                        cutoff = round(cutoff, 1) + 0.25
                        decomp_res = self.get_selected_residues(f'within {cutoff}')
                        if it == 20:
                            # probably not needed, but...
                            GMXMMPBSA_ERROR('The maximum number of iterations to select interaction residues was '
                                            'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for decomposition analysis...")
                else:
                    logging.info(
                        f"Selecting residues by distance ({self.INPUT['decomp']['print_res'].split()[1]} Å) between "
                        f"receptor and ligand for decomposition analysis...")
                explicit_water_res = self._ensure_explicit_water_residues_mapped()
                decomp_res = self._include_explicit_waters_in_decomp(decomp_res, explicit_water_res)
            elif self.INPUT['decomp']['print_res'] == 'all':
                logging.info('Selecting all residues for decomposition analysis...')
            else:
                logging.info('User-selected residues for decomposition analysis...')

            textwraped = textwrap.wrap('\t'.join(x.string for x in decomp_res), tabsize=4, width=120)
            logging.info(f'Selected {len(decomp_res)} residues:\n' + '\n'.join(textwraped) + '\n')

            if self.INPUT['decomp']['idecomp'] in [3, 4]:
                if self.INPUT['decomp']['dec_verbose'] == 0:
                    mol_terms = 1
                elif self.INPUT['decomp']['dec_verbose'] == 1:
                    mol_terms = 3
                elif self.INPUT['decomp']['dec_verbose'] == 2:
                    mol_terms = 4
                else:
                    mol_terms = 12
                energy_terms = 6
                num_res = len(decomp_res)
                total_items = energy_terms * mol_terms * num_res ** 2
                if total_items > 250:
                    logging.warning(f"Using idecomp = {self.INPUT['decomp']['idecomp']} and dec_verbose ="
                                    f" {self.INPUT['decomp']['dec_verbose']} will generate approximately {total_items} items. "
                                    f"Large print selections demand a large amount of memory and take a "
                                    f"significant amount of time to print!")

            self.INPUT['decomp']['print_res'] = ','.join(list2range(decomp_res)['string'])
        if self.INPUT['gb']['ifqnt'] and self.INPUT['gb']['com_qmmask'] == '':
            qm_residues, (rec_charge, lig_charge) = self.get_selected_residues(self.INPUT['gb']['qm_residues'], True)

            if 'within' in self.INPUT['gb']['qm_residues']:
                if len(qm_residues) == 0:
                    logging.info(f"Number of qm_residues using print_res = '{self.INPUT['gb']['qm_residues']}' = 0; "
                                 'expanding the cutoff')
                    logging.info(
                        'Increasing cutoff value by 0.1 until number of qm_residues > 0'
                    )
                    cutoff = float(self.INPUT['gb']['qm_residues'].split()[1])
                    it = 0
                    while len(qm_residues) == 0:
                        cutoff = round(cutoff, 1) + 0.25
                        qm_residues, (rec_charge, lig_charge) = self.get_selected_residues(f'within {cutoff}', True)
                        if it == 20:
                            # probably not needed, but...
                            GMXMMPBSA_ERROR('The maximum number of iterations to select interaction residues was '
                                            'reached. Please set print_res with a valid selection.')
                        it += 1

                    logging.info(f"Selecting residues by distance ({round(cutoff, 1)} Å) between "
                                 f"receptor and ligand for QM/MM calculation...")
                else:
                    logging.info(
                        f"Selecting residues by distance ({self.INPUT['gb']['qm_residues'].split()[1]} Å) between "
                        f"receptor and ligand for QM calculation...")
            elif self.INPUT['gb']['qm_residues'] == 'all':
                logging.info('Selecting all residues for QM calculation...')
            else:
                logging.info('User-selected residues for QM calculation...')

            textwraped = textwrap.wrap('\t'.join(x.string for x in qm_residues), tabsize=4, width=120)
            logging.info(f'Selected {len(qm_residues)} residues:\n' + '\n'.join(textwraped) + '\n')
            self.INPUT['gb']['qm_residues'] = ','.join(list2range(qm_residues)['string'])

            reconcile_qm_charges(self.INPUT['gb'], rec_charge, lig_charge)

        elif self.INPUT['gb']['com_qmmask'] != '':
            logging.warning('Overriding automatic assigment of qmcharge_com, qmcharge_rec, and qmcharge_lig. Using '
                            'default or user defined qmcharge_com, qmcharge_rec, and qmcharge_lig instead...')

        self.cleanup_trajs()
        return tops

    def _warn_gmx_gb_radius_compatibility(self):
        """Warn when the selected GROMACS radius set differs from the usual GB choice."""
        if not self.INPUT.get('gb', {}).get('gbrun', False):
            return

        igb = self.INPUT['gb']['igb']
        recommended = GB_RECOMMENDED_RADII.get(igb)
        if recommended is None:
            return

        selected = PBRadii[self.INPUT['general']['PBRadii']]
        source_family = source_force_field_family(self.INPUT)
        if source_family == 'charmm' and selected.startswith(('bondi', 'mbondi')):
            logging.warning(
                "CHARMM input with AMBER %s radii for GB is a cross-parameterization protocol; "
                "charmm_radii is available for CHARMM PB but is not selected automatically.", selected
            )
        elif source_family == 'opls' and selected.startswith(('bondi', 'mbondi')):
            logging.warning(
                "OPLS input with AMBER %s radii for GB is empirically unvalidated; interpret this "
                "combination as a calibrated scoring protocol.", selected
            )
        elif source_family == 'gromos':
            logging.warning(
                "GROMOS input uses a united-atom-oriented force field; standard all-atom continuum-radius "
                "rules have strong experimental-support limitations.")
        if selected != recommended:
            logging.warning(
                f"PBRadii='{selected}' is selected for the GROMACS topology, while igb={igb} "
                f"is conventionally used with '{recommended}' radii. The selected PBRadii will "
                f"be used; change PBRadii if this combination is not intentional."
            )

    @staticmethod
    def _get_index_group_names(index_file):
        with open(index_file) as ndx_file:
            return [line.split('[', 1)[1].split(']', 1)[0].strip()
                    for line in ndx_file if line.lstrip().startswith('[')]

    def _explicit_water_group_candidates(self):
        if self.explicit_waters_group.strip():
            return [self.explicit_waters_group.strip()]

        candidates = []
        for group_name in explicit_water_group_names:
            if group_name.lower() not in [name.lower() for name in candidates]:
                candidates.append(group_name)
        return candidates

    def gmx2pdb(self):
        """
        Generate PDB file to generate topology
        :return:
        """

        logging.info('Generating PDB files from GROMACS structure files...')

        # wt complex
        # make index for extract pdb structure
        com_rec_group, com_lig_group = self.FILES.complex_groups
        if com_rec_group == com_lig_group:
            GMXMMPBSA_ERROR('The receptor and ligand groups must be different')
        num_com_rec_group, str_com_rec_group = get_index_groups(self.FILES.complex_index, com_rec_group)
        num_com_lig_group, str_com_lig_group = get_index_groups(self.FILES.complex_index, com_lig_group)
        with open(self.FILES.complex_index) as ndx_file:
            combined_group_number = sum(1 for line in ndx_file if line.startswith('['))
        num_com_wat_group = None
        str_com_wat_group = None
        if self.explicit_waters:
            groups = self._get_index_group_names(self.FILES.complex_index)
            solvent_group_candidates = self._explicit_water_group_candidates()
            for candidate in solvent_group_candidates:
                for group_index, group_name in enumerate(groups):
                    if group_name.lower() == candidate.lower():
                        num_com_wat_group = group_index
                        str_com_wat_group = group_name
                        break
                if num_com_wat_group is not None:
                    break
            if num_com_wat_group is None:
                candidate_list = ', '.join(solvent_group_candidates)
                GMXMMPBSA_ERROR(f'EXPLICIT_WATERS requires a solvent group named one of [{candidate_list}] in the '
                                'complex index file. If your solvent group has another name, set '
                                'explicit_waters_group in &general.')
            logging.info(f'Using solvent group {str_com_wat_group} ({num_com_wat_group}) for explicit waters.')

        logging.info('Making gmx_MMPBSA index for complex...')
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        if self.explicit_waters:
            make_ndx_echo_args = echo_command + [
                'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n name {w} GMXMMPBSA_WAT\n'
                ' {r} | {l} | {w}\n name {c} GMXMMPBSA_REC_GMXMMPBSA_LIG\n q\n'.format(
                    r=num_com_rec_group, l=num_com_lig_group, w=num_com_wat_group, c=combined_group_number
                )
            ]
        else:
            make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | '
                                                 '{l}\n name {c} GMXMMPBSA_REC_GMXMMPBSA_LIG\n q\n'.format(
                                                     r=num_com_rec_group, l=num_com_lig_group,
                                                     c=combined_group_number
                                                 )]
        c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

        com_ndx = self.FILES.prefix + 'COM_index.ndx'
        make_ndx_args = self.make_ndx + ['-n', self.FILES.complex_index, '-o', com_ndx, '-f', self.FILES.complex_tpr]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                      ' '.join(make_ndx_args))
        c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c2)
        if c2.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.complex_index))
        self.FILES.complex_index = com_ndx

        if self.explicit_waters:
            logging.info(f'Normal Complex: Saving group {str_com_rec_group}_{str_com_lig_group}_{str_com_wat_group} '
                         f'({num_com_rec_group}_{num_com_lig_group}_{num_com_wat_group}) in '
                         f'{self.FILES.complex_index} file as {self.complex_str_file}')
        else:
            logging.info(f'Normal Complex: Saving group {str_com_rec_group}_{str_com_lig_group} '
                         f'({num_com_rec_group}_{num_com_lig_group}) in {self.FILES.complex_index} file as '
                         f'{self.complex_str_file}')
        # avoid PBC and not chain ID problems
        pdbcom_echo_args = echo_command + ['GMXMMPBSA_REC_GMXMMPBSA_LIG']
        c3 = subprocess.Popen(pdbcom_echo_args, stdout=subprocess.PIPE)

        str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
        if str_format == 'tpr':
            comprog = self.trjconv
            # we extract the pdb from the first frame of trajs to make amber topology
            pdbcom_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                          self.complex_str_file, '-n', self.FILES.complex_index, '-dump', '0']
        else:
            comprog = self.editconf
            pdbcom_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                           self.complex_str_file]
        logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                      (' '.join(pdbcom_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                      '| ' + ' '.join(pdbcom_args))
        c4 = subprocess.Popen(pdbcom_args, stdin=c3.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(c4)
        if c4.wait():  # if it quits with return code != 0
            GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # Put receptor and ligand (explicitly defined) to avoid overwrite them
        # -lm is legacy for the removed tleap path; ligand parameters come from -cp/-lp.
        if self.FILES.ligand_mol2:
            logging.warning(
                '-lm (%s) is ignored when using GROMACS topology conversion (-cp). '
                'Ligand parameters must already be present in the complex (and, for MT, ligand) topology.',
                self.FILES.ligand_mol2,
            )

        # make a temp receptor pdb (even when stability) if decomp to get correct receptor residues from complex. This
        # avoids get multiples molecules from complex.split()
        if self.INPUT['decomp']['decomprun'] and self.FILES.stability:
            self.use_temp = True
            logging.info('Generating a receptor file internally to extract interface residues for decomposition.')
            rec_echo_args = echo_command + ['{}'.format(num_com_rec_group)]
            cp1 = subprocess.Popen(rec_echo_args, stdout=subprocess.PIPE)
            if str_format == 'tpr':
                # we extract the pdb from the first frame of trajs to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              'rec_temp.pdb', '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               'rec_temp.pdb']
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(rec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check if stability
        if self.FILES.stability and (
                (self.FILES.receptor_tpr or self.FILES.ligand_tpr)
        ):
            logging.info('Stability calculation mode does not need separate receptor or ligand files; ignoring them.')
        # wt receptor
        if self.FILES.receptor_tpr:
            logging.info('A receptor structure file was defined. Using MT approach...')
            num_rec_group, str_rec_group = get_index_groups(self.FILES.receptor_index, self.FILES.receptor_group)

            logging.info('Making gmx_MMPBSA index for receptor...')
            make_ndx_echo_args = echo_command + ['name {r} GMXMMPBSA_REC\n q\n'.format(r=num_rec_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            rec_ndx = self.FILES.prefix + 'REC_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.receptor_index, '-o', rec_ndx, '-f',
                                             self.FILES.receptor_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.receptor_index))
            self.FILES.receptor_index = rec_ndx

            logging.info(f'Normal Receptor: Saving group {str_rec_group} ({num_rec_group}) in '
                         f'{self.FILES.receptor_index} file as {self.receptor_str_file}')
            pdbrec_echo_args = echo_command + ['{}'.format(num_rec_group)]
            p1 = subprocess.Popen(pdbrec_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.receptor_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                prog = self.trjconv
                # we extract a pdb from structure file to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.receptor_trajs[0], '-s', self.FILES.receptor_tpr, '-o',
                                              self.receptor_str_file, '-n', self.FILES.receptor_index, '-dump', '0']
            else:
                prog = self.editconf
                pdbrec_args = self.editconf + ['-f', self.FILES.receptor_tpr, '-n', self.FILES.receptor_index, '-o',
                                               self.receptor_str_file]

            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdbrec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(prog), self.FILES.receptor_trajs[0]))
        else:
            logging.info('No receptor structure file was defined. Using ST approach...')
            logging.info('Using receptor structure from complex to generate AMBER topology')
            logging.info(f'Normal Receptor: Saving group {str_com_rec_group} ({num_com_rec_group}) in '
                         f'{self.FILES.complex_index} file as {self.receptor_str_file}')
            pdbrec_echo_args = echo_command + ['{}'.format(num_com_rec_group)]
            cp1 = subprocess.Popen(pdbrec_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdbrec_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              self.receptor_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdbrec_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               self.receptor_str_file]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdbrec_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdbrec_args))
            cp2 = subprocess.Popen(pdbrec_args, stdin=cp1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(cp2)
            if cp2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # ligand
        # # check consistence
        if self.FILES.ligand_tpr:  # unbound ligand structure/trajs (MT)
            logging.info('A ligand structure file was defined. Using MT approach...')
            num_lig_group, str_lig_group = get_index_groups(self.FILES.ligand_index, self.FILES.ligand_group)

            logging.info('Making gmx_MMPBSA index for ligand...')
            make_ndx_echo_args = echo_command + ['name {l} GMXMMPBSA_LIG\n q\n'.format(l=num_lig_group)]
            c1 = subprocess.Popen(make_ndx_echo_args, stdout=subprocess.PIPE)

            lig_ndx = self.FILES.prefix + 'LIG_index.ndx'
            make_ndx_args = self.make_ndx + ['-n', self.FILES.ligand_index, '-o', lig_ndx, '-f', self.FILES.ligand_tpr]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(make_ndx_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' + ' | ' +
                          ' '.join(make_ndx_args))
            c2 = subprocess.Popen(make_ndx_args, stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c2)
            if c2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.make_ndx), self.FILES.ligand_index))
            self.FILES.ligand_index = lig_ndx

            logging.info(f'Normal Ligand: Saving group {str_lig_group} ({num_lig_group}) in {self.FILES.ligand_index}'
                         f' file as {self.ligand_str_file}')
            # wt ligand
            pdblig_echo_args = echo_command + ['{}'.format(num_lig_group)]
            l1 = subprocess.Popen(pdblig_echo_args, stdout=subprocess.PIPE)
            str_format = 'tpr' if self.FILES.ligand_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                prog = self.trjconv
                # we extract a pdb from structure file to make amber topology
                pdblig_args = self.trjconv + ['-f', self.FILES.ligand_trajs[0], '-s', self.FILES.ligand_tpr, '-o',
                                              self.ligand_str_file, '-n', self.FILES.ligand_index, '-dump', '0']
            else:
                prog = self.editconf
                pdblig_args = self.editconf + ['-f', self.FILES.ligand_tpr, '-n', self.FILES.ligand_index, '-o',
                                               self.ligand_str_file]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdblig_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdblig_args))
            l2 = subprocess.Popen(pdblig_args, stdin=l1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(l2)
            if l2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(prog), self.FILES.ligand_trajs[0]))
        else:
            # wt complex ligand
            logging.info('No ligand structure file was defined. Using ST approach...')
            logging.info('Using ligand structure from complex to generate AMBER topology')
            logging.info(f'Normal Ligand: Saving group {str_com_lig_group} ({num_com_lig_group}) in '
                         f'{self.FILES.complex_index} file as {self.ligand_str_file}')
            pdblig_echo_args = echo_command + ['{}'.format(num_com_lig_group)]
            l1 = subprocess.Popen(pdblig_echo_args, stdout=subprocess.PIPE)

            str_format = 'tpr' if self.FILES.complex_tpr[-3:] == 'tpr' else 'pdb'
            if str_format == 'tpr':
                # we extract a pdb from structure file to make amber topology
                pdblig_args = self.trjconv + ['-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                              self.ligand_str_file, '-n', self.FILES.complex_index, '-dump', '0']
            else:
                pdblig_args = self.editconf + ['-f', self.FILES.complex_tpr, '-n', self.FILES.complex_index, '-o',
                                               self.ligand_str_file]

            # we extract a pdb from structure file to make amber topology
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(pdblig_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(pdblig_args))
            l2 = subprocess.Popen(pdblig_args, stdin=l1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(l2)
            if l2.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(comprog), self.FILES.complex_trajs[0]))
        # check for IE variable
        if (self.FILES.receptor_tpr or self.FILES.ligand_tpr) and (
                self.INPUT['general']['interaction_entropy'] or self.INPUT['general']['c2_entropy']
        ):
            logging.warning("The IE or C2 entropy method doesn't support the MTP approach...")
            self.INPUT['general']['interaction_entropy'] = self.INPUT['general']['c2_entropy'] = 0

        # initialize receptor and ligand structures. Needed to get residues map
        logging.info('Loading extracted complex, receptor, and ligand PDB files with ParmEd...')
        self.complex_str = self.molstr(self.complex_str_file)
        self.receptor_str = self.molstr(self.receptor_str_file)
        self.ligand_str = self.molstr(self.ligand_str_file)
        logging.info('Loaded structures: complex %d atoms/%d residues, receptor %d atoms/%d residues, '
                     'ligand %d atoms/%d residues.',
                     len(self.complex_str.atoms), len(self.complex_str.residues),
                     len(self.receptor_str.atoms), len(self.receptor_str.residues),
                     len(self.ligand_str.atoms), len(self.ligand_str.residues))
        if self.FILES.reference_structure:
            logging.info('Loading reference structure for chain/residue consistency checks...')
            self.ref_str = check_str(self.FILES.reference_structure, ref=True)
        self.check4water()
        logging.info('Reading receptor/ligand atom indexes and building residue maps...')
        self.indexes = get_indexes(com_ndx=self.FILES.complex_index,
                                   rec_ndx=self.FILES.receptor_index,
                                   lig_ndx=self.FILES.ligand_index)
        self.resi, self.resl, self.orderl = res2map(self.indexes, self.complex_str)
        logging.info('Residue map built: %d receptor residues, %d ligand residues.',
                     sum(end - start + 1 for start, end in self.resi['REC']['num']),
                     sum(end - start + 1 for start, end in self.resi['LIG']['num']))
        self.check_structures(self.complex_str, self.receptor_str, self.ligand_str)

    def check4water(self):
        if self.explicit_waters:
            return

        if counter := sum(res.name in solvent_ion_residues for res in self.complex_str.residues):
            GMXMMPBSA_ERROR(f'gmx_MMPBSA does not support water/ions molecules in any structure, but we found'
                            f' {counter} molecules in the complex.')

    def _check_periodicity(self, parm, system):
        """
        check for periodicities == 0 and change them to 1. This is required especially for nmode calculations
        """
        invalid_per = 0

        for dt in parm.dihedral_types:
            if dt.per == 0:
                invalid_per += 1
                dt.per = 1

        if invalid_per:
            logging.warning(f'{invalid_per} invalid DIHEDRAL_PERIODICITY = 0 found in {system.capitalize()} '
                            f'topology... Setting DIHEDRAL_PERIODICITY = 1')

        return parm

    def _explicit_water_residues(self, parm):
        return [res.idx + 1 for res in parm.residues if res.name in water_residues]

    def _explicit_ion_residues(self, parm):
        return [res.idx + 1 for res in parm.residues if res.name in solvent_ion_residues and res.name not in water_residues]

    def _explicit_water_range(self, parm):
        if not self.explicit_waters:
            return ''
        water_res = self._explicit_water_residues(parm)
        if len(water_res) < self.explicit_waters:
            GMXMMPBSA_ERROR(f'EXPLICIT_WATERS requested {self.explicit_waters} waters, but only '
                            f'{len(water_res)} water residues were found in the complex topology.')
        return ','.join(list2range(water_res[:self.explicit_waters])['string'])

    def _strip_extra_explicit_waters(self, parm):
        if not self.explicit_waters:
            return ''
        ion_res = self._explicit_ion_residues(parm)
        if ion_res:
            parm.strip(f":{','.join(list2range(ion_res)['string'])}")
        water_res = self._explicit_water_residues(parm)
        water_range = self._explicit_water_range(parm)
        extra_waters = water_res[self.explicit_waters:]
        if extra_waters:
            parm.strip(f":{','.join(list2range(extra_waters)['string'])}")
        return water_range

    @staticmethod
    def _extra_point_atom_indices(parm):
        return [
            getattr(atom, 'idx', index - 1) + 1 for index, atom in enumerate(parm.atoms, start=1)
            if (
                getattr(atom, 'atomic_number', None) == 0 or
                getattr(atom, 'mass', None) == 0 or
                getattr(atom, 'type', '').upper() == 'EP' or
                getattr(atom, 'name', '').upper().startswith('EP')
            )
        ]

    def _extra_point_atom_mask(self, parm):
        atom_indices = self._extra_point_atom_indices(parm)
        if not atom_indices:
            return ''
        return '@' + ','.join(list2range(atom_indices)['string'])

    def _check_explicit_waters_supported_by_energy_model(self, parm):
        if not self.explicit_waters:
            return
        extra_point_mask = self._extra_point_atom_mask(parm)
        if not extra_point_mask:
            return

        if getattr(self, 'explicit_waters_extra_points', 'error') == 'error':
            GMXMMPBSA_ERROR(
                'EXPLICIT_WATERS found extra-point water atoms from a virtual-site water model such as OPC/TIP4P. '
                'sander calculations can fail with these atoms. Set EXPLICIT_WATERS_EXTRA_POINTS="strip" to remove '
                'the extra points and continue, or use a 3-site water model such as TIP3P/SPC.'
            )

        self.explicit_water_extra_point_mask = extra_point_mask
        extra_point_count = len(self._range_string_to_list(extra_point_mask[1:]))
        parm.strip(extra_point_mask)
        logging.warning(
            'EXPLICIT_WATERS_EXTRA_POINTS="strip" removed %d extra-point atom(s) from the selected explicit waters. '
            'The selected water model is approximated after removing virtual sites; use this only for controlled '
            'relative comparisons.',
            extra_point_count
        )
        logging.warning(
            'GB/PB energies with stripped OPC/TIP4P extra points should be interpreted cautiously because the '
            'water electrostatics no longer correspond to the original virtual-site model.'
        )

    def _resolve_explicit_waters_mask(self):
        if not self.explicit_waters:
            return

        selection = self.explicit_waters_mask.strip()
        if selection.lower() == 'dasa':
            self._resolve_dasa_explicit_waters_mask()
            return

        if not selection.startswith('within'):
            return

        selected_residues = self.get_selected_residues(selection)
        cutoff = float(selection.split()[1])
        if len(selected_residues) < 2:
            logging.info(f"Number of interface residues selected using explicit_waters_mask = '{selection}' < 2; "
                         'expanding the cutoff')
            logging.info('Increasing cutoff value by 0.25 until number of interface residues selected >= 2')
            it = 0
            while len(selected_residues) < 2:
                cutoff = round(cutoff, 1) + 0.25
                selected_residues = self.get_selected_residues(f'within {cutoff}')
                if it == 20:
                    GMXMMPBSA_ERROR('The maximum number of iterations to select interaction residues was reached. '
                                    'Please set explicit_waters_mask with a valid selection.')
                it += 1

        if not selected_residues:
            GMXMMPBSA_ERROR('EXPLICIT_WATERS_MASK did not select any interface residues.')

        textwraped = textwrap.wrap('\t'.join(x.string for x in selected_residues), tabsize=4, width=120)
        logging.info(f"Selecting waters closest to interface residues defined by explicit_waters_mask = "
                     f"'within {round(cutoff, 1)} Å'...")
        logging.info(f'Selected {len(selected_residues)} interface residues:\n' + '\n'.join(textwraped) + '\n')

        resolved_mask = ':' + ','.join(list2range(selected_residues)['string'])
        self.explicit_waters_mask = resolved_mask
        self.INPUT['general']['explicit_waters_mask'] = resolved_mask
        logging.info(f'Resolved explicit water reference mask for cpptraj closest: {resolved_mask}')

    def _explicit_water_closest_reference_mask(self):
        return f'({self.explicit_waters_mask})&(!{explicit_water_reference_exclusion_mask})'

    @staticmethod
    def _dasa_residue_mask(residues):
        return ':' + ','.join(list2range(residues)['string'])

    def _resolve_dasa_explicit_waters_mask(self):
        receptor_residues = [res for res in self.resl if res.is_receptor()]
        ligand_residues = [res for res in self.resl if res.is_ligand()]
        if not receptor_residues or not ligand_residues:
            GMXMMPBSA_ERROR('EXPLICIT_WATERS_MASK="dASA" requires receptor and ligand residues.')
        if not self.explicit_water_prmtop:
            GMXMMPBSA_ERROR('EXPLICIT_WATERS_MASK="dASA" requires a full solvated AMBER topology.')
        if not self.FILES.complex_trajs:
            GMXMMPBSA_ERROR('EXPLICIT_WATERS_MASK="dASA" requires a complex trajectory.')

        receptor_mask = self._dasa_residue_mask(receptor_residues)
        ligand_mask = self._dasa_residue_mask(ligand_residues)
        solute_mask = f'({receptor_mask}|{ligand_mask})'
        all_residues = receptor_residues + ligand_residues
        dataset_names = []
        actions = [f'trajin {self.FILES.complex_trajs[0]} 1 1', 'noprogress']
        for prefix, residues, environment in (
                ('com', all_residues, solute_mask),
                ('rec', receptor_residues, receptor_mask),
                ('lig', ligand_residues, ligand_mask)):
            for residue in residues:
                name = f'{prefix}{residue.index}'
                dataset_names.append(name)
                actions.append(f'surf {name} :{residue.index} solutemask {environment}')

        output = f'{self.FILES.prefix}explicit_waters_dasa.dat'
        log = f'{self.FILES.prefix}explicit_waters_dasa.out'
        actions.extend(['run', f'writedata {output} ' + ' '.join(dataset_names)])
        cutoff = self.INPUT['general']['explicit_waters_dasa_cutoff']
        logging.info('Selecting interface residues with cpptraj dASA cutoff %.3g for explicit waters...', cutoff)
        logging.debug('Running cpptraj dASA calculation with topology %s and first frame of %s',
                      self.explicit_water_prmtop, self.FILES.complex_trajs[0])
        with open(log, 'w') as log_file:
            process = subprocess.Popen([self.external_progs['cpptraj'], self.explicit_water_prmtop],
                                       stdin=subprocess.PIPE, stdout=log_file, stderr=subprocess.STDOUT)
            process.communicate(('\n'.join(actions) + '\n').encode())
        if process.wait():
            GMXMMPBSA_ERROR(f'{self.external_progs["cpptraj"]} failed when calculating dASA. Check {log}.')

        try:
            with open(output) as data_file:
                data_lines = [line.split() for line in data_file
                              if line.strip() and not line.lstrip().startswith('#')]
        except FileNotFoundError:
            GMXMMPBSA_ERROR(f'cpptraj did not write {output} when calculating dASA. Check {log}.')
        if not data_lines or len(data_lines[0]) != len(dataset_names) + 1:
            GMXMMPBSA_ERROR(f'cpptraj returned an unexpected dASA dataset in {output}. Check {log}.')

        values = dict(zip(dataset_names, map(float, data_lines[0][1:])))
        selected_residues = []
        for residue in all_residues:
            complex_area = values[f'com{residue.index}']
            isolated_prefix = 'rec' if residue.is_receptor() else 'lig'
            isolated_area = values[f'{isolated_prefix}{residue.index}']
            if abs(isolated_area - complex_area) >= cutoff:
                selected_residues.append(residue)
        if not selected_residues:
            GMXMMPBSA_ERROR('cpptraj dASA did not select any interface residues for '
                            'EXPLICIT_WATERS_MASK="dASA".')

        textwraped = textwrap.wrap('\t'.join(x.string for x in selected_residues), tabsize=4, width=120)
        logging.info('Selected %d cpptraj dASA interface residues:\n%s\n', len(selected_residues),
                     '\n'.join(textwraped))
        resolved_mask = ':' + ','.join(list2range(selected_residues)['string'])
        self.explicit_waters_mask = resolved_mask
        self.INPUT['general']['explicit_waters_mask'] = resolved_mask
        logging.info(f'Resolved explicit water reference mask for cpptraj closest: {resolved_mask}')

    def _assign_component_radii(self, parm, component, has_topology):
        if has_topology:
            logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to {component}...")
            ChRad(parm, PBRadii[self.INPUT['general']['PBRadii']])
        else:
            logging.info(
                f"Preserving {component} GB radii inherited from Complex: "
                f"{parm.parm_data.get('RADIUS_SET', ['unknown'])[0]}"
            )

    def gmxtop2prmtop(self):
        logging.info('Using topology conversion. Setting radiopt = 0...')
        self.INPUT['pb']['radiopt'] = 0
        self._warn_gmx_gb_radius_compatibility()
        logging.info('Building Normal Complex Amber topology...')
        com_top, error_info = self._cleantop_with_retry(
            self.FILES.complex_top, self.indexes['COM']['COM'], self.complex_str
        )
        if error_info:
            topology_mismatch_error('complex', self.FILES.complex_top, self.complex_str_file, error_info)

        logging.info('Assigning complex coordinates to the selected topology...')
        com_top.coordinates = self.complex_str.coordinates
        logging.info('Writing complex restart coordinates to %sCOM.inpcrd...', self.FILES.prefix)
        com_top.save(f"{self.FILES.prefix}COM.inpcrd", format='rst7', overwrite=True)
        # try:
        if com_top.impropers or com_top.urey_bradleys:
            logging.info('Converting selected complex topology to AMBER ChamberParm...')
            com_amb_prm = parmed.amber.ChamberParm.from_structure(com_top)
            com_top_parm = 'chamber'

            title = com_amb_prm.parm_data['CTITLE']
            com_amb_prm.add_flag('TITLE', '20a4', title or '', after='CTITLE')

            logging.info('Detected CHARMM force field topology format...')
        else:
            logging.info('Converting selected complex topology to AMBER AmberParm...')
            com_amb_prm = parmed.amber.AmberParm.from_structure(com_top)
            com_top_parm = 'amber'
            logging.info('Detected Amber/OPLS force field topology format...')

        # IMPORTANT: make_trajs ends in error if the box is defined
        com_amb_prm.box = None

        # check periodicity
        com_amb_prm = self._check_periodicity(com_amb_prm, 'complex')

        self.fixparm2amber(com_amb_prm)
        explicit_water_range = ''
        if self.explicit_waters:
            self.explicit_water_prmtop = f'{self.FILES.prefix}COM_FULL_SOLVENT.prmtop'
            com_amb_prm.write_parm(self.explicit_water_prmtop)
            explicit_water_range = self._strip_extra_explicit_waters(com_amb_prm)
            self.explicit_water_range = explicit_water_range
            self._check_explicit_waters_supported_by_energy_model(com_amb_prm)
            logging.info(f'Keeping {self.explicit_waters} explicit water residues assigned to the receptor.')

        logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Complex...")
        if com_top_parm == 'amber' and self.INPUT['general']['PBRadii'] == 7:
            GMXMMPBSA_ERROR(
                f"The PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} is not compatible with Amber/OPLS "
                f"topologies...")
        action = ChRad(com_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
        logging.info('Writing Normal Complex AMBER topology...')
        com_amb_prm.write_parm(self.complex_pmrtop)

        rec_indexes_string = ','.join(self.resi['REC']['string'])

        rec_hastop = True
        if self.FILES.receptor_top:
            logging.info('A Receptor topology file was defined. Using MT approach...')
            logging.info('Building AMBER Receptor Topology from GROMACS Receptor Topology...')
            rec_top, error_info = self._cleantop_with_retry(
                self.FILES.receptor_top, self.indexes['REC'], self.receptor_str, 'receptor'
            )

            if error_info:
                topology_mismatch_error('receptor', self.FILES.receptor_top, self.receptor_str_file, error_info)

            rec_top.coordinates = self.receptor_str.coordinates
            # rec_top.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)
            if rec_top.impropers or rec_top.urey_bradleys:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is Amber/OPLS type while the '
                                    'Receptor is CHAMBER type!')
                rec_amb_prm = parmed.amber.ChamberParm.from_structure(rec_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Receptor is Amber/OPLS type!')
                rec_amb_prm = parmed.amber.AmberParm.from_structure(rec_top)
            logging.info('Converting receptor residue names from GROMACS to AMBER...')

            # check periodicity
            rec_amb_prm = self._check_periodicity(rec_amb_prm, 'receptor')

            self.fixparm2amber(rec_amb_prm)
        else:
            logging.info('No Receptor topology file was defined. Using ST approach...')
            logging.info('Building AMBER Receptor topology from Complex...')
            # we make a copy for receptor topology
            rec_amb_prm = self.molstr(com_amb_prm)
            rec_keep = rec_indexes_string
            if explicit_water_range:
                rec_keep = f'{rec_keep},{explicit_water_range}'
            rec_amb_prm.strip(f'!:{rec_keep}')
            rec_hastop = False

        self._assign_component_radii(rec_amb_prm, 'Receptor', rec_hastop)
        logging.info('Writing Normal Receptor AMBER topology...')
        rec_amb_prm.write_parm(self.receptor_pmrtop)
        rec_amb_prm.save(f"{self.FILES.prefix}REC.inpcrd", format='rst7', overwrite=True)

        lig_hastop = True
        if self.FILES.ligand_top:
            logging.info('A Ligand Topology file was defined. Using MT approach...')
            logging.info('Building AMBER Ligand Topology from GROMACS Ligand Topology...')
            lig_top, error_info = self._cleantop_with_retry(
                self.FILES.ligand_top, self.indexes['LIG'], self.ligand_str, 'ligand'
            )

            if error_info:
                topology_mismatch_error('ligand', self.FILES.ligand_top, self.ligand_str_file, error_info)

            lig_top.coordinates = self.ligand_str.coordinates
            # lig_top.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)
            if lig_top.impropers or lig_top.urey_bradleys:
                if com_top_parm == 'amber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is Amber/OPLS type while the '
                                    'Ligand is CHAMBER type!')
                lig_amb_prm = parmed.amber.ChamberParm.from_structure(lig_top)
            else:
                if com_top_parm == 'chamber':
                    GMXMMPBSA_ERROR('Inconsistent parameter format. The defined Complex is CHAMBER type while the '
                                    'Ligand is Amber/OPLS type!')
                lig_amb_prm = parmed.amber.AmberParm.from_structure(lig_top)
            logging.info('Converting ligand residue names from GROMACS to AMBER...')

            # check periodicity
            lig_amb_prm = self._check_periodicity(lig_amb_prm, 'ligand')

            self.fixparm2amber(lig_amb_prm)
        else:
            logging.info('No Ligand topology file was defined. Using ST approach...')
            logging.info('Building AMBER Ligand topology from Complex...')
            # we make a copy for ligand topology
            lig_amb_prm = self.molstr(com_amb_prm)
            lig_strip = rec_indexes_string
            if explicit_water_range:
                lig_strip = f'{lig_strip},{explicit_water_range}'
            lig_amb_prm.strip(f':{lig_strip}')
            lig_hastop = False
        self._assign_component_radii(lig_amb_prm, 'Ligand', lig_hastop)
        logging.info('Writing Normal Ligand AMBER topology...')
        lig_amb_prm.write_parm(self.ligand_pmrtop)
        lig_amb_prm.save(f"{self.FILES.prefix}LIG.inpcrd", format='rst7', overwrite=True)

        if self.INPUT['ala']['alarun']:
            logging.info('Building Mutant Complex Topology...')
            # get mutation index in complex
            self.com_mut_indices, self.part_mut, self.part_indices = self.getMutationInfo()
            self.com_mut_index = self.com_mut_indices[0] if len(self.com_mut_indices) == 1 else None
            self.part_index = self.part_indices[0] if len(self.part_indices) == 1 else None
            mut_com_amb_prm = self.makeMutTop(com_amb_prm, self.com_mut_indices)
            logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant Complex...")
            action = ChRad(mut_com_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
            logging.info('Writing Mutant Complex AMBER topology...')
            mut_com_amb_prm.write_parm(self.mutant_complex_pmrtop)

            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor topology...')
                out_prmtop = self.mutant_receptor_pmrtop
                self.mutant_ligand_pmrtop = None
                if rec_hastop:
                    mtop = self.makeMutTop(rec_amb_prm, self.part_indices)
                else:
                    mut_rec_keep = rec_indexes_string
                    if explicit_water_range:
                        mut_rec_keep = f'{mut_rec_keep},{explicit_water_range}'
                    mut_com_amb_prm.strip(f'!:{mut_rec_keep}')
                    mtop = mut_com_amb_prm
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand topology...')
                out_prmtop = self.mutant_ligand_pmrtop
                self.mutant_receptor_pmrtop = None
                if lig_hastop:
                    mtop = self.makeMutTop(lig_amb_prm, self.part_indices)
                else:
                    mut_lig_strip = rec_indexes_string
                    if explicit_water_range:
                        mut_lig_strip = f'{mut_lig_strip},{explicit_water_range}'
                    mut_com_amb_prm.strip(f':{mut_lig_strip}')
                    mtop = mut_com_amb_prm

            if com_top_parm == 'chamber':
                mut_prot_amb_prm = parmed.amber.ChamberParm.from_structure(mtop)
            else:
                mut_prot_amb_prm = parmed.amber.AmberParm.from_structure(mtop)
            logging.info(f"Assigning PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant "
                         f"{'Receptor' if self.part_mut == 'REC' else 'Ligand'}...")
            action = ChRad(mut_prot_amb_prm, PBRadii[self.INPUT['general']['PBRadii']])
            logging.info(f"Writing Mutant {'Receptor' if self.part_mut == 'REC' else 'Ligand'} AMBER topology...")
            mut_prot_amb_prm.write_parm(out_prmtop)
        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)

    def _split_str(self, start, r, c, basename, struct, mut_index=None):
        end = start + (r[1] - r[0])
        mask = f'!:{start}-{end}'
        str_ = self.molstr(struct)
        if mut_index is not None and (not isinstance(mut_index, (list, tuple)) or mut_index):
            str_ = self.makeMutTop(str_, mut_index, True)
        str_.strip(mask)
        str_file = f'{self.FILES.prefix}{basename}_F{c}.pdb'
        str_.save(str_file, 'pdb', True, renumber=False)
        return end, str_file

    def pdb2prmtop(self):
        """Legacy structure→PDB prep for the removed tleap topology path.

        Retained only for unit tests that exercise PDB splitting/mutation helpers.
        Production GROMACS runs use ``gmxtop2prmtop`` exclusively.
        :return:
        """
        self._warn_gmx_gb_radius_compatibility()
        if self.INPUT['general']['PBRadii'] == 7:
            GMXMMPBSA_ERROR(f"The PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} is not compatible with "
                            f"Amber topologies...")

        logging.info('Generating AMBER Compatible PDB Files...')
        # fix receptor and structures
        logging.info('Converting complex residue names from GROMACS to AMBER...')
        self.fixparm2amber(self.complex_str, 'COM')
        logging.info('Converting receptor residue names from GROMACS to AMBER...')
        self.fixparm2amber(self.receptor_str, 'REC')
        logging.info('Converting ligand residue names from GROMACS to AMBER...')
        self.fixparm2amber(self.ligand_str, 'LIG')

        logging.info('Splitting  receptor and ligand in PDB files..')
        self.receptor_list = {}
        start = 1
        for c, r in enumerate(self.resi['REC']['num'], start=1):
            end, sfile = self._split_str(start, r, c, 'REC', self.receptor_str)
            self.receptor_list[f'REC{c}'] = sfile
            start = end + 1

        self.ligand_list = {}
        start = 1
        for c, r in enumerate(self.resi['LIG']['num'], start=1):
            end, sfile = self._split_str(start, r, c, 'LIG', self.ligand_str)
            self.ligand_list[f'LIG{c}'] = sfile
            start = end + 1

        self.mut_receptor_list = {}
        self.mut_ligand_list = {}
        if self.INPUT['ala']['alarun']:
            self.com_mut_indices, self.part_mut, self.part_indices = self.getMutationInfo()
            self.com_mut_index = self.com_mut_indices[0] if len(self.com_mut_indices) == 1 else None
            self.part_index = self.part_indices[0] if len(self.part_indices) == 1 else None
            start = 1
            if self.part_mut == 'REC':
                logging.info('Detecting mutation in Receptor. Building Mutant Receptor structure...')
                self.mutant_ligand_pmrtop = None
                for c, r in enumerate(self.resi['REC']['num']):
                    segment_end = start + (r[1] - r[0])
                    segment_mut_indices = [index for index in self.part_indices
                                           if start - 1 <= index < segment_end]
                    end, sfile = self._split_str(
                        start, r, c, 'MUT_REC', self.receptor_str, segment_mut_indices
                    )
                    self.mut_receptor_list[f'MREC{c}'] = sfile
                    start = end + 1
            else:
                logging.info('Detecting mutation in Ligand. Building Mutant Ligand Structure...')
                self.mutant_receptor_pmrtop = None
                for c, r in enumerate(self.resi['LIG']['num']):
                    segment_end = start + (r[1] - r[0])
                    segment_mut_indices = [index for index in self.part_indices
                                           if start - 1 <= index < segment_end]
                    end, sfile = self._split_str(
                        start, r, c, 'MUT_LIG', self.ligand_str, segment_mut_indices
                    )
                    self.mut_ligand_list[f'MLIG{c}'] = sfile
                    start = end + 1

    def _cleantop_with_retry(self, top_file, ndx, structure, id='complex'):
        logging.info('Preparing %s topology from %s using %d selected atom indexes...',
                     id, top_file, len(ndx))
        # Keep solvent on the first pass when explicit receptor waters are
        # requested; otherwise water index atoms vanish and only the retry path
        # recovers them.
        keep_solvent = bool(getattr(self, 'explicit_waters', 0))
        try:
            top = self.cleantop(top_file, ndx, id, remove_solvent=not keep_solvent)
        except IndexError as err:
            logging.warning(f'{err} Retrying with the full topology before applying the index...')
            try:
                logging.info('Preparing %s topology again without removing solvent before applying the index...', id)
                top = self.cleantop(top_file, ndx, id, remove_solvent=False)
            except IndexError:
                GMXMMPBSA_ERROR(f'The atom index in the {id} index is not found in the topology file. Please check '
                                'that the files are consistent.')
            return top, eq_strs(top, structure)

        error_info = eq_strs(top, structure)
        if error_info:
            logging.warning(
                f'The {id} topology generated after removing solvent/ions before applying the index is inconsistent '
                'with the structure. Retrying with the full topology before applying the index...'
            )
            try:
                logging.info('Preparing %s topology again without removing solvent before applying the index...', id)
                top = self.cleantop(top_file, ndx, id, remove_solvent=False)
            except IndexError:
                GMXMMPBSA_ERROR(f'The atom index in the {id} index is not found in the topology file. Please check '
                                'that the files are consistent.')
            error_info = eq_strs(top, structure)

        return top, error_info

    @staticmethod
    def cleantop(top_file, ndx, id='complex', remove_solvent=True):
        """
        Create a new top file with selected groups and without SOL and IONS
        :param top_file: User-defined topology file
        :param ndx: atoms index
        :param remove_solvent: remove solvent/ions from the temporary topology before applying the index
        :return: new and clean top instance
        """
        top_file = Path(top_file)
        logging.info('Preprocessing %s include tree for %s topology conversion%s...',
                     top_file, id, ' with solvent/ion removal' if remove_solvent else '')
        preprocessor = GromacsTopologyPreprocessor()
        temp_top = preprocessor.preprocess(top_file, remove_solvent, solvent_ion_residues)
        if preprocessor.cmap_found:
            logging.warning(
                'Ignoring CMAP terms in %s include tree for GROMACS topology conversion. '
                'The converted topology omits CMAP energy terms. For STP (single-trajectory) '
                'MM/PB(GB)SA this is not an issue: CMAP contributions cancel in the C−R−L '
                'difference. Consider the omission only for MTP (multiple-trajectory) '
                'calculations, where receptor and ligand come from separate ensembles.',
                top_file)

        # read the temp topology with parmed
        logging.info('Reading preprocessed %s topology with ParmEd...', id)
        rtemp_top = parmed.gromacs.GromacsTopologyFile(temp_top.as_posix())
        # get the residues in the top from the com_ndx
        logging.info('Applying %s index selection to topology (%d atoms selected)...', id, len(ndx))
        res_list = []

        for i in ndx:
            try:
                idx = rtemp_top.atoms[i - 1].residue.idx + 1
                if idx not in res_list:
                    res_list.append(rtemp_top.atoms[i - 1].residue.number + 1)
            except IndexError:
                for temp_file in preprocessor.created_files:
                    temp_file.unlink(missing_ok=True)
                raise IndexError(
                    f'The atom {i} in the {id} index is not found in the topology generated from {top_file}'
                )

        ranges = list2range(res_list)
        rtemp_top.strip(f"!:{','.join(ranges['string'])}")
        logging.info('Prepared %s topology selection: %d atoms/%d residues retained.',
                     id, len(rtemp_top.atoms), len(rtemp_top.residues))

        # Clean temporal file
        for temp_file in preprocessor.created_files:
            temp_file.unlink(missing_ok=True)
        return rtemp_top

    def get_masks(self):
        rec_mask = ':' + ','.join(self.resi['REC']['string'])
        lig_mask = ':' + ','.join(self.resi['LIG']['string'])
        if self.explicit_waters:
            dry_residues = []
            for part in ['REC', 'LIG']:
                for start, end in self.resi[part]['num']:
                    dry_residues.extend(range(start, end + 1))
            first_water = max(dry_residues) + 1
            last_water = first_water + self.explicit_waters - 1
            rec_mask = f'{rec_mask},{first_water}-{last_water}'

        if self.INPUT['ala']['alarun']:
            for index in self.com_mut_indices:
                self.resl[index].set_mut(self.INPUT['ala']['mutant'])
        return rec_mask, lig_mask, self.resl

    @staticmethod
    def _range_string_to_list(range_string):
        residue_numbers = []
        if not range_string:
            return residue_numbers
        for item in range_string.split(','):
            item = item.strip()
            if not item:
                continue
            if '-' in item:
                start, end = map(int, item.split('-', 1))
                residue_numbers.extend(range(start, end + 1))
            else:
                residue_numbers.append(int(item))
        return residue_numbers

    def _ensure_explicit_water_residues_mapped(self):
        if not self.explicit_waters or not self.explicit_water_range:
            return []

        rec_count = sum(end - start + 1 for start, end in self.resi['REC']['num'])
        mapped_residue_indexes = {res.index for res in self.resl}
        complex_prmtop = parmed.load_file(self.complex_pmrtop)
        water_res = []
        for offset, res_index in enumerate(self._range_string_to_list(self.explicit_water_range), start=1):
            if res_index > len(complex_prmtop.residues):
                GMXMMPBSA_ERROR(f'Explicit water residue {res_index} is not present in the complex topology.')
            top_res = complex_prmtop.residues[res_index - 1]
            res = Residue(res_index, res_index, '', 'R', rec_count + offset, top_res.name)
            water_res.append(res)
            if res.index not in mapped_residue_indexes:
                self.resl.append(res)
                mapped_residue_indexes.add(res.index)
        self.resl.sort(key=lambda res: res.index)
        return water_res

    @staticmethod
    def _include_explicit_waters_in_decomp(decomp_res, explicit_water_res):
        if not explicit_water_res:
            return decomp_res

        selected = {res.index for res in decomp_res}
        added = [res for res in explicit_water_res if res.index not in selected]
        if added:
            logging.info(f'Including {len(added)} explicit water residues assigned to the receptor in '
                         'decomposition print_res.')
            decomp_res = decomp_res + added
        return sorted(decomp_res, key=lambda res: res.index)

    @staticmethod
    def _global_frame_ranges(frame_counts, startframe, endframe, interval):
        """Map global trajectory frame selection to per-file cpptraj ranges.

        ``startframe``, ``endframe`` and ``interval`` refer to the concatenated
        trajectory supplied by the user.  cpptraj applies a ``trajin`` range to
        each file independently, so a range must be calculated for every file
        before building the explicit-water preprocessing script.
        """
        ranges = []
        offset = 0
        for file_index, frame_count in enumerate(frame_counts):
            file_start = offset + 1
            file_end = offset + frame_count
            selected_start = max(startframe, file_start)
            selected_end = min(endframe, file_end)
            if selected_start <= selected_end:
                first = startframe + ((selected_start - startframe + interval - 1) // interval) * interval
                if first <= selected_end:
                    last = first + ((selected_end - first) // interval) * interval
                    ranges.append((file_index, first - offset, last - offset, interval))
            offset += frame_count
        return ranges

    def _cpptraj_frame_count(self, trajectory):
        """Return the number of frames in one trajectory file."""
        process = subprocess.Popen(
            [self.external_progs['cpptraj'], '-p', self.explicit_water_prmtop, '-y', trajectory, '-tl'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        output, _ = process.communicate(b'')
        if process.wait():
            GMXMMPBSA_ERROR(f"{self.external_progs['cpptraj']} failed when querying {trajectory}")
        if isinstance(output, bytes):
            output = output.decode()
        frame_counts = re.findall(r'Frames:\s+(\d+)', output)
        if len(frame_counts) != 1:
            GMXMMPBSA_ERROR(f'Could not determine the number of frames in {trajectory} with '
                            f'{self.external_progs["cpptraj"]}.')
        return int(frame_counts[0])

    def get_selected_residues(self, select, qm_sele=False):
        """
        Convert string selection format to amber index list
        """
        if qm_sele:
            com_top = parmed.load_file(self.complex_pmrtop)

        dist, res_selection = selector(select)
        residues_selection = {'rec': [], 'lig': []}
        rec_charge = 0
        lig_charge = 0
        if dist:
            rec_residues = [res for res in self.resl if res.is_receptor()]
            lig_residues = [res for res in self.resl if res.is_ligand()]
            logging.info('Scanning receptor/ligand residue contacts within %.3g Å (%d receptor residues x %d '
                         'ligand residues)...', dist, len(rec_residues), len(lig_residues))
            for rres in self.resl:
                if rres.is_ligand():
                    continue
                for lres in self.resl:
                    if lres.is_receptor():
                        continue
                    for rat in self.complex_str.residues[rres - 1].atoms:
                        rat_coor = [rat.xx, rat.xy, rat.xz]
                        for lat in self.complex_str.residues[lres - 1].atoms:
                            lat_coor = [lat.xx, lat.xy, lat.xz]
                            if get_dist(rat_coor, lat_coor) <= dist:
                                if rres not in residues_selection['rec']:
                                    residues_selection['rec'].append(rres)
                                    if qm_sele:
                                        rec_charge += sum(atm.charge for atm in com_top.residues[rres - 1].atoms)
                                if lres not in residues_selection['lig']:
                                    residues_selection['lig'].append(lres)
                                    if qm_sele:
                                        lig_charge += sum(atm.charge for atm in com_top.residues[lres - 1].atoms)
                                break
            logging.info('Contact scan selected %d receptor and %d ligand residues.',
                         len(residues_selection['rec']), len(residues_selection['lig']))
        elif res_selection:
            for i in self.resl:
                rres = self.complex_str.residues[i - 1]
                if [rres.chain, rres.number, rres.insertion_code] in res_selection:
                    if i.is_ligand():
                        residues_selection['lig'].append(i)
                        if qm_sele:
                            # Sum raw partial charges; round once at the end so
                            # list / within / all selection paths agree.
                            lig_charge += sum(atm.charge for atm in com_top.residues[i - 1].atoms)
                    else:
                        residues_selection['rec'].append(i)
                        if qm_sele:
                            rec_charge += sum(atm.charge for atm in com_top.residues[i - 1].atoms)
                    res_selection.remove([rres.chain, rres.number, rres.insertion_code])
            for res in res_selection:
                logging.warning("We couldn't find this residue CHAIN:{} RES_NUM:{} ICODE: {}".format(*res))
            # check if residues in receptor and ligand was defined
            if not residues_selection['rec'] or not residues_selection['lig']:
                if not self.INPUT['ala']['alarun']:
                    GMXMMPBSA_ERROR(
                        'For decomposition analysis, you must define residues for both receptor and ligand!')
        else:
            for i in self.resl:
                if i.is_ligand():
                    residues_selection['lig'].append(i)
                    if qm_sele:
                        lig_charge += sum(atm.charge for atm in com_top.residues[i - 1].atoms)
                else:
                    residues_selection['rec'].append(i)
                    if qm_sele:
                        rec_charge += sum(atm.charge for atm in com_top.residues[i - 1].atoms)
        if qm_sele:
            rec_charge = int(round(rec_charge))
            lig_charge = int(round(lig_charge))
        sele_res = sorted([r for m in residues_selection.values() for r in m], key=lambda x: x.index)
        return (sele_res, (rec_charge, lig_charge)) if qm_sele else sele_res

    def fixparm2amber(self, structure, str_name=None):

        for c, residue in enumerate(structure.residues, start=1):
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
                for atom in residue.atoms:
                    if 'SG' in atom.name:
                        for bondedatm in atom.bond_partners:
                            if bondedatm.name == 'SG':
                                if str_name:
                                    if str_name == 'COM':
                                        cys1 = c
                                        cys2 = structure.residues.index(bondedatm.residue) + 1
                                    else:
                                        cys1 = residue.number
                                        cys2 = bondedatm.residue.number
                                    if ([cys1, cys2] not in self.cys_bonds[str_name] and
                                            [cys2, cys1] not in self.cys_bonds[str_name]):
                                        self.cys_bonds[str_name].append([cys1, cys2])
                                if residue.name == 'CYX' and bondedatm.residue.name == 'CYX':
                                    continue
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
        if str_name:
            structure.strip('@/H')

    def getMutationInfo(self):
        if not self.INPUT['ala']['mutant_res']:
            GMXMMPBSA_ERROR("No residue for mutation was defined")
        # dict = { resind: [chain, resnum, icode]
        sele_res_dict = self.get_selected_residues(self.INPUT['ala']['mutant_res'])
        if not sele_res_dict:
            GMXMMPBSA_ERROR('No valid residue was found for mutation')

        parts = {'REC' if residue.is_receptor() else 'LIG' if residue.is_ligand() else None
                 for residue in sele_res_dict}
        if None in parts:
            residue = next(residue for residue in sele_res_dict
                           if not residue.is_receptor() and not residue.is_ligand())
            GMXMMPBSA_ERROR(f'Residue {residue.chain}:{residue.number} not found')
        if len(parts) > 1:
            GMXMMPBSA_ERROR('Composite alanine/glycine mutations cannot mix receptor and ligand residues.')
        if len(sele_res_dict) > 1 and self.INPUT['ala']['cas_intdiel']:
            GMXMMPBSA_ERROR('cas_intdiel=1 is ambiguous for composite mutations. Set cas_intdiel=0 or select one residue.')

        com_mut_indices = []
        part_indices = []
        for residue_ref in sele_res_dict:
            res = self.complex_str.residues[residue_ref - 1]
            icode = f':{res.insertion_code}' if res.insertion_code else ''
            if (
                not parmed.residue.AminoAcidResidue.has(res.name) and res.name not in ['HSP', 'HSE', 'HSD']
                or res.name in ['CYX', 'PRO', 'GLY']
                or res.name == 'ALA' and self.INPUT['ala']['mutant'] == 'ALA'
            ):
                GMXMMPBSA_ERROR(f"Selecting residue {res.chain}:{res.name}:{res.number}{icode} can't be mutated. Please, "
                                f"define a valid residue...")
            com_mut_indices.append(residue_ref - 1)
            part_indices.append(residue_ref.id_index - 1)

        return com_mut_indices, next(iter(parts)), part_indices

    def _assign_ter(self, structure=None):
        structure = self.complex_str if structure is None else structure
        for res in structure.residues:
            # evident terminal
            if len(res.name) == 4:
                if res.name.startswith('N'):
                    res.ter = 'N'  # ter is a boolean variable, but it is not used anyway
                elif res.name.startswith('C'):
                    res.ter = 'C'  # ter is a boolean variable, but it is not used anyway
                else:
                    res.ter = False
            else:
                atms_name = [at.name for at in res.atoms]
                if 'OXT' in atms_name or 'OT1' in atms_name or 'OT2' in atms_name:  # already used, but is better to get here anyway
                    res.ter = 'C'  # ter is a boolean variable, but it is not used anyway
                elif (('H3' in atms_name or 'HT3' in atms_name) and
                      parmed.residue.AminoAcidResidue.has(res.name)  # exclude other residues with H3, for examples ligs
                ):
                    res.ter = 'N'  # ter is a boolean variable, but it is not used anyway
                else:
                    res.ter = False

    def makeMutTop(self, wt_top, mut_index, pdb=False):
        """Apply the requested mutation to one or more residue indices."""
        mut_indices = [mut_index] if isinstance(mut_index, int) else list(mut_index)
        if not mut_indices:
            return self.molstr(wt_top)
        if len(mut_indices) == 1:
            return self._makeMutTopSingle(wt_top, mut_indices[0], pdb)
        mut_top = self.molstr(wt_top)
        for index in mut_indices:
            mut_top = self._makeMutTopSingle(mut_top, index, pdb, copy=False)
        return mut_top

    def _makeMutTopSingle(self, wt_top, mut_index, pdb=False, copy=True):
        """

        :param wt_top: Amber parm from GROMACS topology
        :param mut_index: index of mutation in structure
        :return: Mutant AmberParm
        """
        mut_top = self.molstr(wt_top) if copy else wt_top
        mut_aa = self.INPUT['ala']['mutant']

        bb_atoms = 'N,H,CA,HA,C,O,HN'
        nterm_atoms = 'H1,H2,H3,HT1,HT2,HT3'
        cterm_atoms = 'OXT,OT1,OT2'  # OXT amber, OTx charmm
        sc_cb_atom = 'CB'
        sc_ala_atoms = ('HB,' +  # VAL, ILE, THR
                        'HB1,HB2,' +  # charmm -> HB1, HB2
                        'HB3,' +  # amber -> HB2, HB3
                        'CG1,CG2,OG1,' +  # VAL, ILE, THR
                        'OG,' +  # SER
                        'SG,' +  # CYS
                        'CG')

        if mut_aa in ['GLY', 'G']:
            strip_mask = f":{mut_index + 1} &!@{','.join([bb_atoms, nterm_atoms, cterm_atoms])}"
            if not pdb:
                strip_mask += f",{sc_cb_atom}"
        else:
            strip_mask = f":{mut_index + 1} &!@{','.join([bb_atoms, sc_cb_atom, nterm_atoms, cterm_atoms])}"
            if not pdb:
                strip_mask += f",{sc_ala_atoms}"
        mut_top.strip(strip_mask)

        # add terminals only for mutation
        self._assign_ter(mut_top)

        # solution for issue #364
        # NOTE: We selected the charge from  Amber14SB because it is the same as amber99sb, amber99SB-ILDN, amber12SB,
        # etc. In any case, the error associated with this charge must be small.
        # GLY-HB atom type - amber 14SB
        hc = parmed.AtomType('HC', None, 1.008, 1)
        hc.set_lj_params(eps=0.0157, rmin=1.4870, eps14=0.0157, rmin14=1.4870)
        # GLY-HA atom type - amber 14SB
        h1 = parmed.AtomType('H1', None, 1.008, 1)
        h1.set_lj_params(eps=0.0157, rmin=1.3870, eps14=0.0157, rmin14=1.3870)
        # GLY-HB atom type - charmm
        ha3 = parmed.AtomType('HA3', None, 1.008, 1)
        ha3.set_lj_params(eps=1.3400, rmin=0.0240, eps14=1.3400, rmin14=0.0240)
        # GLY-HA atom type - charmm
        hb2 = parmed.AtomType('HB2', None, 1.008, 1)
        hb2.set_lj_params(eps=1.3400, rmin=0.0280, eps14=1.3400, rmin14=0.0280)

        # In amber the charge depend on aa position in the sequence.
        # ALA: C: 0.0764, N: 0.0300, int: 0.0603
        # GLY: C: 0.1056, N: 0.0895, int: 0.0698

        mutation_residue = mut_top.residues[mut_index]
        if mutation_residue.ter == 'C':
            h_ala_charge = 0.0764
            h_gly_charge = 0.1056
        elif mutation_residue.ter == 'N':
            h_ala_charge = 0.0300
            h_gly_charge = 0.0895
        else:
            h_ala_charge = 0.0603
            h_gly_charge = 0.0698

        h_atoms_prop = {
            'charmm': {
                'ALA': {
                    'mass': 1.008, 'element': 'H', 'atomic_number': 1, 'atom_type': ha3, 'type': 'HA3',
                    'charge': 0.09},
                'GLY': {
                    'mass': 1.008, 'element': 'H', 'atomic_number': 1, 'atom_type': hb2, 'type': 'HB2',
                    'charge': 0.09}},
            'amber': {
                'ALA': {
                    'mass': 1.008, 'element': 'H', 'atomic_number': 1, 'atom_type': hc, 'type': 'HC',
                    'charge': h_ala_charge},
                'GLY': {
                    'mass': 1.008, 'element': 'H', 'atomic_number': 1, 'atom_type': h1, 'type': 'H1',
                    'charge': h_gly_charge}}
        }
        ff_rep = 'charmm' if isinstance(mut_top, parmed.amber.ChamberParm) else 'amber'

        cb_atom = None
        ca_atom = None
        logging.info(
            f"Mutating {mutation_residue.chain}/{mutation_residue.number} "
            f"{mutation_residue.name} to {mut_aa}")

        mutant_resname = mut_top.residues[mut_index].name

        mut_top.residues[mut_index].name = mut_aa

        for at in mut_top.residues[mut_index].atoms:
            if mut_aa == 'GLY':
                if at.name == 'CA':
                    ca_atom = at
                if at.name in ['CB']:
                    at.name = 'HA2'
                    ca_atom.xx, ca_atom.xy, ca_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [ca_atom.xx, ca_atom.xy,
                         ca_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for p in h_atoms_prop[ff_rep]['GLY']:
                        setattr(at, p, h_atoms_prop[ff_rep]['GLY'][p])
                elif at.name in ['HA']:
                    at.name = 'HA1'
                    for p in h_atoms_prop[ff_rep]['GLY']:
                        setattr(at, p, h_atoms_prop[ff_rep]['GLY'][p])
            else:
                # ARG, ASN, ASP, GLN, GLU, HIS, LEU, LYS, MET, PHE, TYR, TRP
                #    |
                # HN-N
                #    |   HB1
                #    |   |
                # HA-CA--CB--CG (ARG, ASN, ASP, GLN, GLU, HIS, LEU, LYS, MET, PHE, TYR, TRP)
                #    |   |
                #    |   HB2
                #  O=C
                #    |
                #
                # CYS (no S-S), SER
                #    |
                # HN-N
                #    |   HB1
                #    |   |
                # HA-CA--CB--SG|OG (CYS | SER)
                #    |   |
                #    |   HB2
                #  O=C
                #    |
                #
                # ILE, THR, VAL
                #    |
                # HN-N
                #    |     CG2 (ILE, VAL, THR)
                #    |    /
                # HA-CA--CB-HB
                #    |    \
                #    |     CG1|OG1  (ILE, VAL | THR)
                #  O=C
                #    |
                #
                # EXCLUDE: GLY, PRO, ALA

                if at.name == 'CB':
                    cb_atom = at
                    continue
                if at.name == 'CG2':  # VAL, ILE and THR
                    at.name = 'HB2'
                    cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [cb_atom.xx, cb_atom.xy,
                         cb_atom.xz, at.xx, at.xy,
                         at.xz], 1.09)
                    for p in h_atoms_prop[ff_rep]['ALA']:
                        setattr(at, p, h_atoms_prop[ff_rep]['ALA'][p])
                elif at.name in ['HB']:  # ILE, VAL and THR
                    at.name = 'HB1'
                    for p in h_atoms_prop[ff_rep]['ALA']:
                        setattr(at, p, h_atoms_prop[ff_rep]['ALA'][p])
                elif at.name in ['CG',  # LEU, PHE, TRP, MET, TYR, ARG, LYS, ASN, GLN, ASP, GLU, HIS,
                                 'OG',  # SER
                                 'SG',  # CYS (no S-S)
                                 'CG1',  # VAL, ILE and THR
                                 'OG1'  # THR
                                 ]:
                    # check if it was assigned. In some cases, the HB can be HB2 and HB3 instead HB1 and HB2
                    if 'HB3' not in [atm.name for atm in mut_top.residues[mut_index].atoms]:
                        at.name = 'HB3'
                    else:
                        at.name = 'HB1'
                    cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz = _scaledistance(
                        [cb_atom.xx, cb_atom.xy, cb_atom.xz, at.xx, at.xy, at.xz], 1.09)
                    for p in h_atoms_prop[ff_rep]['ALA']:
                        setattr(at, p, h_atoms_prop[ff_rep]['ALA'][p])
                elif at.name in ['HB1', 'HB2', 'HB3']:
                    for p in h_atoms_prop[ff_rep]['ALA']:
                        setattr(at, p, h_atoms_prop[ff_rep]['ALA'][p])

        # change intdiel if cas_intdiel was defined before end the mutation process
        if self.INPUT['ala']['cas_intdiel']:
            if self.INPUT['gb']['gbrun']:
                if self.INPUT['gb']['intdiel'] != 1.0:
                    logging.warning('Both cas_intdiel and intdiel were defined. The dielectric constants associated '
                                    'with cas_intdiel will be ignored and intdiel will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(f"Setting intdiel = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for "
                                 f"Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(
                        f"Setting intdiel = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting intdiel = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                        f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gb']['intdiel'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting intdiel = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                        f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default "
                                    f"intdiel will be used")
            if self.INPUT['gbnsr6']['gbnsr6run']:
                if self.INPUT['gbnsr6']['epsin'] != 1.0:
                    logging.warning('Both cas_intdiel and epsin were defined. The dielectric constants associated '
                                    'with cas_intdiel will be ignored and epsin will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(f"Setting epsin = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for "
                                 f"Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(
                        f"Setting epsin = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting epsin = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for Alanine "
                        f"scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['gbnsr6']['epsin'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting epsin = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                        f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default "
                                    f"intdiel will be used")

            if self.INPUT['pb']['pbrun']:
                if self.INPUT['pb']['indi'] != 1.0:
                    logging.warning('Both cas_intdiel and indi were defined. The dielectric constants associated with '
                                    'cas_intdiel will be ignored and indi will be used instead')
                elif mutant_resname in polar_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_polar']
                    logging.info(
                        f"Setting indi = intdiel_polar = {self.INPUT['ala']['intdiel_polar']} for Alanine scanning")
                elif mutant_resname in nonpolar_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_nonpolar']
                    logging.info(
                        f"Setting indi = intdiel_nonpolar = {self.INPUT['ala']['intdiel_nonpolar']} for Alanine "
                        f"scanning")
                elif mutant_resname in positive_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_positive']
                    logging.info(
                        f"Setting intdiel = indi = intdiel_positive = {self.INPUT['ala']['intdiel_positive']} for "
                        f"Alanine scanning")
                elif mutant_resname in negative_aa:
                    self.INPUT['pb']['indi'] = self.INPUT['ala']['intdiel_negative']
                    logging.info(
                        f"Setting indi = intdiel_negative = {self.INPUT['ala']['intdiel_negative']} for Alanine "
                        f"scanning")
                else:
                    logging.warning(f"Unclassified mutant residue {mutant_resname}. The default indi will be used")
        return mut_top

    def cleanup_trajs(self):
        # clear trajectory
        if not self.INPUT['general']['solvated_trajectory']:
            return
        logging.info('Cleaning normal complex trajectories...')
        if self.explicit_waters:
            warn_concatenated_complex_trajectories(self.FILES.complex_trajs)
        new_trajs = []
        full_trajs = []
        for i in range(len(self.FILES.complex_trajs)):
            trjconv_echo_args = echo_command + ['GMXMMPBSA_REC_GMXMMPBSA_LIG']
            c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file and make amber topology for complex
            com_traj_name = (f'{self.FILES.prefix}COM_full_traj_{i}.xtc' if self.explicit_waters
                             else f'COM_traj_{i}.xtc')
            trjconv_args = self.trjconv + ['-f', self.FILES.complex_trajs[i], '-s', self.FILES.complex_tpr, '-o',
                                           com_traj_name, '-n', self.FILES.complex_index]
            logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                          (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                          '| ' + ' '.join(trjconv_args))
            c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_subprocess_output(c6)
            if c6.wait():  # if it quits with return code != 0
                GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.complex_trajs[i]))
            if self.explicit_waters:
                full_trajs.append(com_traj_name)
            else:
                new_trajs.append(f'COM_traj_{i}.xtc')

        if self.explicit_waters:
            startframe = self.INPUT['general']['startframe']
            endframe = self.INPUT['general']['endframe']
            interval = self.INPUT['general']['interval']
            if len(full_trajs) == 1:
                trajin_commands = f'trajin {full_trajs[0]} {startframe} {endframe} {interval}\n'
            else:
                frame_counts = [self._cpptraj_frame_count(trajectory) for trajectory in full_trajs]
                frame_ranges = self._global_frame_ranges(frame_counts, startframe, endframe, interval)
                if not frame_ranges:
                    GMXMMPBSA_ERROR('No frames were selected across the explicit-water trajectories.')
                logging.info('Applying global frame selection across %d concatenated trajectories.', len(full_trajs))
                trajin_commands = ''.join(
                    f'trajin {full_trajs[file_index]} {local_start} {local_end} {local_interval}\n'
                    for file_index, local_start, local_end, local_interval in frame_ranges
                )

            filtered_traj = f'{self.FILES.prefix}COM_traj_0.mdcrd'
            extra_point_strip = (f'strip {self.explicit_water_extra_point_mask}\n'
                                 if self.explicit_water_extra_point_mask else '')
            cpptraj_input = (
                f'{trajin_commands}'
                f'strip {explicit_water_ion_mask}\n'
                f'closest {self.explicit_waters} {self._explicit_water_closest_reference_mask()} '
                f'solventmask {explicit_water_solvent_mask} noimage '
                f'closestout {self.FILES.prefix}explicit_waters_closest_0.dat\n'
                f'{extra_point_strip}'
                f'trajout {filtered_traj} nobox\n'
            )
            logging.debug('Running command: %s %s', self.external_progs['cpptraj'], self.explicit_water_prmtop)
            with open(f'{self.FILES.prefix}explicit_waters_cpptraj_0.out', 'w') as cpptraj_out:
                c7 = subprocess.Popen([self.external_progs['cpptraj'], self.explicit_water_prmtop],
                                      stdin=subprocess.PIPE, stdout=cpptraj_out, stderr=subprocess.STDOUT)
                c7.communicate(cpptraj_input.encode())
            if c7.wait():
                GMXMMPBSA_ERROR(f"{self.external_progs['cpptraj']} failed when selecting explicit waters from "
                                f"{', '.join(full_trajs)}")
            new_trajs.append(filtered_traj)
            self.FILES.explicit_waters_preselected = True
        self.FILES.complex_trajs = new_trajs

        # clear trajectory
        if self.FILES.receptor_tpr:
            logging.info('Cleaning normal receptor trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.receptor_trajs)):
                trjconv_echo_args = echo_command + ['GMXMMPBSA_REC']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.receptor_trajs[i], '-s', self.FILES.receptor_tpr,
                                               '-o', 'REC_traj_{}.xtc'.format(i), '-n', self.FILES.receptor_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR(
                        '%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.receptor_trajs[i]))
                new_trajs.append('REC_traj_{}.xtc'.format(i))
            self.FILES.receptor_trajs = new_trajs

        if self.FILES.ligand_tpr:
            logging.info('Cleaning normal ligand trajectories...')
            new_trajs = []
            for i in range(len(self.FILES.ligand_trajs)):
                trjconv_echo_args = echo_command + ['GMXMMPBSA_LIG']
                c5 = subprocess.Popen(trjconv_echo_args, stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file and make amber topology for complex
                trjconv_args = self.trjconv + ['-f', self.FILES.ligand_trajs[i], '-s', self.FILES.ligand_tpr, '-o',
                                               'LIG_traj_{}.xtc'.format(i), '-n', self.FILES.ligand_index]
                logging.debug('Running command: ' + ' '.join(echo_command) + ' "' +
                              (' '.join(trjconv_echo_args[len(echo_command):]).replace('\n', '\\n')) + '"' +
                              '| ' + ' '.join(trjconv_args))
                c6 = subprocess.Popen(trjconv_args, stdin=c5.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_subprocess_output(c6)
                if c6.wait():  # if it quits with return code != 0
                    GMXMMPBSA_ERROR('%s failed when querying %s' % (' '.join(self.trjconv), self.FILES.ligand_trajs[i]))
                new_trajs.append('LIG_traj_{}.xtc'.format(i))
            self.FILES.ligand_trajs = new_trajs

    def check_structures(self, com_str, rec_str=None, lig_str=None):
        logging.info('Checking structural consistency...')
        logging.info('Validating complex structure...')
        check_str(com_str)
        logging.info('Validating receptor structure...')
        check_str(rec_str, skip=True)
        logging.info('Validating ligand structure...')
        check_str(lig_str, skip=True)

        if self.FILES.reference_structure:
            logging.info('Assigning chain IDs and insertion codes to structure files according to the reference structure...')
            ref_str = check_str(self.FILES.reference_structure)
            if len(ref_str.residues) != len(com_str.residues):
                GMXMMPBSA_ERROR(f'The number of residues of the complex ({len(com_str.residues)}) and of the '
                                f'reference structure ({len(ref_str.residues)}) are different. Please check that the '
                                f'reference structure is correct')
            for c, res in enumerate(ref_str.residues):
                if com_str.residues[c].number != res.number or not residue_names_match(com_str.residues[c].name,
                                                                                       res.name):
                    GMXMMPBSA_ERROR('There is no match between the complex and the reference structure used. An '
                                    f'attempt was made to assign the chain ID to "{com_str.residues[c].name}'
                                    f':{com_str.residues[c].number}:{com_str.residues[c].insertion_code}" in the '
                                    f'complex, but "{res.name}:{res.number}:{res.insertion_code}" was expected '
                                    'based on the reference structure. Please check that the reference structure is '
                                    'correct')
                com_str.residues[c].chain = res.chain
                com_str.residues[c].insertion_code = res.insertion_code
                # The explicit-water workflow can include solvent residues in
                # the full reference structure, while ``self.resl`` contains
                # only receptor/ligand residues.  Assign the chain to the
                # complex residue, but do not map solvent residues into the
                # receptor/ligand residue index list.
                if c >= len(self.resl):
                    continue
                # update the chain and insertion code (https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/354)
                self.resl[c].chain = res.chain
                self.resl[c].icode = res.insertion_code
                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                    rec_str.residues[i].insertion_code = res.insertion_code
                else:
                    lig_str.residues[i].chain = res.chain
                    lig_str.residues[i].insertion_code = res.insertion_code
        else:
            assign = False
            if self.INPUT['general']['assign_chainID'] == 1:
                assign = not com_str.residues[0].chain  # pretty simple
                if assign:
                    logging.info('Chain IDs not found; assigning chain IDs...')
                else:
                    logging.info('Chain IDs found; skipping chain-ID assignment...')
            elif self.INPUT['general']['assign_chainID'] == 2:
                assign = True
                if com_str.residues[0].chain:
                    logging.warning('Reassigning existing chain IDs as requested.')
                else:
                    logging.info('Assigning missing chain IDs as requested.')
            elif self.INPUT['general']['assign_chainID'] == 0 and not com_str.residues[0].chain:
                assign = True
                logging.info('No reference structure or chain IDs were provided; assigning chain IDs automatically.')
            if assign:
                self._assign_chains_IDs(com_str, rec_str, lig_str)
        # Save fixed complex structure for analysis and set it in FILES to save in info file
        logging.info('Writing fixed complex structure to %sCOM_FIXED.pdb...', self.FILES.prefix)
        com_str.save(f'{self.FILES.prefix}COM_FIXED.pdb', 'pdb', True, renumber=False)
        logging.info('Structure consistency checks complete.')

    def _assign_chains_IDs(self, com_str, rec_str, lig_str):
        chains_ids = []
        chain_by_num = False
        chain_by_ter = False
        previous_res_number = 0
        curr_chain_id = 'A'
        has_nucl = 0
        for c, res in enumerate(com_str.residues):
            if c >= len(self.resl):
                continue
            if res.chain:
                if res.chain != curr_chain_id:
                    res.chain = curr_chain_id
                    i = self.resl[c].id_index - 1
                    if self.resl[c].is_receptor():
                        rec_str.residues[i].chain = res.chain
                    else:
                        lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            else:
                res.chain = curr_chain_id

                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
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

                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
                    lig_str.residues[i].chain = res.chain
                if res.chain not in chains_ids:
                    chains_ids.append(res.chain)
            elif chain_by_ter:
                chain_by_ter = False
            elif chain_by_num:
                chain_by_num = False
                curr_chain_id = chains_letters[chains_letters.index(chains_ids[-1]) + 1]
                res.chain = curr_chain_id
                i = self.resl[c].id_index - 1
                if self.resl[c].is_receptor():
                    rec_str.residues[i].chain = res.chain
                else:
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
            logging.warning('This structure contains nucleotides. We recommend that you use the reference structure')

    @staticmethod
    def molstr(data):
        if type(data) == str:
            # data is a pdb file
            pdb_file = data
            try:
                new_str = []
                with open(pdb_file) as fo:
                    fo = fo.readlines()
                    for line in fo:
                        if 'MODEL' in line or 'ENDMDL' in line:
                            continue
                        # check new charmm-gui format for Amber ff19SB (with N- and C- terminals)
                        if 'ATOM' in line:
                            resn = line[17:21].strip()
                            if len(resn) == 4 and resn.startswith(('N', 'C')):
                                line = f'{line[:17]}{resn[1:]} {line[21:]}'
                        new_str.append(line)
                with open(pdb_file, 'w') as fw:
                    for x in new_str:
                        fw.write(x)
            except IOError as e:
                GMXMMPBSA_ERROR(str(e))

            structure = parmed.read_PDB(pdb_file)
        else:
            # data is Structure, AmberParm, ChamberParm or GromacsTopologyFile. This make a copy
            structure = data.__copy__()
            for c, at in enumerate(structure.atoms, start=1):
                at.number = c
        return structure

    def _write_ff(self, ofile):
        for ff in self.INPUT['general']['forcefields']:
            ofile.write(f'source {ff}\n')
        ofile.write('loadOff atomic_ions.lib\n')
        ofile.write('loadamberparams {}\n'.format(ions_para_files[self.INPUT['general']['ions_parameters']]))
        # check if it is a modified PBRadii
        if self.INPUT['general']['PBRadii'] in [5, 6]:
            ofile.write('set default PBRadii {}\n'.format(PBRadii[1]))
        else:
            ofile.write('set default PBRadii {}\n'.format(PBRadii[self.INPUT['general']['PBRadii']]))

    def makeToptleap(self):
        """Legacy tleap topology builder (structure→loadpdb/mol2).

        Removed from the production GROMACS path; ``-cp`` topology conversion is required.
        Kept for unit tests that still call this helper directly.
        """
        logging.info('Building tleap input files...')
        with open(f'{self.FILES.prefix}leap.in', 'w') as tif:
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
                LIG.extend(f'{lig}' for lig in self.ligand_list)
            else:
                for lig in self.ligand_list:
                    LIG.append(f'{lig}')
                    tif.write(f'{lig} = loadpdb {self.ligand_list[lig]}\n')
                lig_out = ' '.join(LIG)
                if self.FILES.stability:
                    self.ligand_pmrtop = None
                else:
                    tif.write(f'LIG_OUT = combine {{ {lig_out} }}\n')
                    for cys1, cys2 in self.cys_bonds['LIG']:
                        tif.write(f'bond LIG_OUT.{cys1}.SG LIG_OUT.{cys2}.SG\n')
                    tif.write(f'saveamberparm LIG_OUT {self.ligand_pmrtop} {self.FILES.prefix}LIG.inpcrd\n')
            COM = self._set_com_order(REC, LIG)
            if self.FILES.stability:
                self.receptor_pmrtop = None
            else:
                tif.write(f'REC_OUT = combine {{ {rec_out} }}\n')
                for cys1, cys2 in self.cys_bonds['REC']:
                    tif.write(f'bond REC_OUT.{cys1}.SG REC_OUT.{cys2}.SG\n')
                tif.write(f'saveamberparm REC_OUT {self.receptor_pmrtop} {self.FILES.prefix}REC.inpcrd\n')
            com_out = ' '.join(COM)
            tif.write(f'COM_OUT = combine {{ {com_out} }}\n')
            for cys1, cys2 in self.cys_bonds['COM']:
                tif.write(f'bond COM_OUT.{cys1}.SG COM_OUT.{cys2}.SG\n')
            tif.write('saveamberparm COM_OUT {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            tif.write('quit')
        # changed in v1.4.3. We source the gmxMMPBSA ff directly from the data folder instead of copy to the Amber/dat
        data_path = Path(__file__).parent.joinpath('data')
        tleap = self.external_progs['tleap']
        self._run_tleap(tleap, 'leap.in', data_path)

        # check if it is a modified PBRadii
        if self.INPUT['general']['PBRadii'] in [5, 6]:
            com_prmtop = parmed.load_file(self.complex_pmrtop)
            com_amb_parm = parmed.amber.AmberParm.from_structure(com_prmtop)
            action = ChRad(com_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
            logging.info(
                f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Normal Complex AMBER "
                f"topology...")
            com_amb_parm.write_parm(self.complex_pmrtop)
            if not self.FILES.stability:
                rec_prmtop = parmed.load_file(self.receptor_pmrtop)
                rec_amb_parm = parmed.amber.AmberParm.from_structure(rec_prmtop)
                action = ChRad(rec_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
                logging.info(
                    f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Normal Receptor AMBER "
                    f"topology...")
                rec_amb_parm.write_parm(self.receptor_pmrtop)

                lig_prmtop = parmed.load_file(self.ligand_pmrtop)
                lig_amb_parm = parmed.amber.AmberParm.from_structure(lig_prmtop)
                action = ChRad(lig_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
                logging.info(
                    f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Normal Ligand AMBER "
                    f"topology...")
                lig_amb_parm.write_parm(self.ligand_pmrtop)

        if self.INPUT['ala']['alarun']:
            with open(f'{self.FILES.prefix}mut_leap.in', 'w') as mtif:
                self._write_ff(mtif)

                if self.mutant_receptor_pmrtop:
                    REC = []
                    for mrec in self.mut_receptor_list:
                        REC.append(f'{mrec}')
                        mtif.write(f'{mrec} = loadpdb {self.mut_receptor_list[mrec]}\n')
                    mrec_out = ' '.join(REC)

                    if not self.FILES.stability:
                        mtif.write(f'MREC_OUT = combine {{ {mrec_out} }}\n')
                        for cys1, cys2 in self.cys_bonds['REC']:
                            mtif.write(f'bond MREC_OUT.{cys1}.SG MREC_OUT.{cys2}.SG\n')
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
                        for cys1, cys2 in self.cys_bonds['LIG']:
                            mtif.write(f'bond MLIG_OUT.{cys1}.SG MLIG_OUT.{cys2}.SG\n')
                        mtif.write('saveamberparm MLIG_OUT {t} {p}MUT_LIG.inpcrd\n'.format(
                            t=self.mutant_ligand_pmrtop, p=self.FILES.prefix))
                    else:
                        self.mutant_ligand_pmrtop = None
                    for rec in self.receptor_list:
                        mtif.write(f'{rec} = loadpdb {self.receptor_list[rec]}\n')

                MCOM = self._set_com_order(REC, LIG)
                mcom_out = ' '.join(MCOM)
                mtif.write(f'MCOM_OUT = combine {{ {mcom_out} }}\n')
                for cys1, cys2 in self.cys_bonds['COM']:
                    mtif.write(f'bond MCOM_OUT.{cys1}.SG MCOM_OUT.{cys2}.SG\n')
                mtif.write('saveamberparm MCOM_OUT {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                                   p=self.FILES.prefix))
                mtif.write('quit')

            self._run_tleap(tleap, 'mut_leap.in', data_path)

            # check if it is a modified PBRadii
            if self.INPUT['general']['PBRadii'] in [5, 6]:
                mcom_prmtop = parmed.load_file(self.mutant_complex_pmrtop)
                mcom_amb_parm = parmed.amber.AmberParm.from_structure(mcom_prmtop)
                action = ChRad(mcom_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
                logging.info(
                    f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant Complex AMBER "
                    f"topology...")
                mcom_amb_parm.write_parm(self.mutant_complex_pmrtop)
                if not self.FILES.stability:
                    mrec_prmtop = parmed.load_file(self.mutant_receptor_pmrtop)
                    mrec_amb_parm = parmed.amber.AmberParm.from_structure(mrec_prmtop)
                    action = ChRad(mrec_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
                    logging.info(
                        f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant Receptor "
                        f"AMBER topology...")
                    mrec_amb_parm.write_parm(self.mutant_receptor_pmrtop)

                    mlig_prmtop = parmed.load_file(self.mutant_ligand_pmrtop)
                    mlig_amb_parm = parmed.amber.AmberParm.from_structure(mlig_prmtop)
                    action = ChRad(mlig_amb_parm, PBRadii[self.INPUT['general']['PBRadii']])
                    logging.info(
                        f"Assigning modified PBRadii {PBRadii[self.INPUT['general']['PBRadii']]} to Mutant Ligand AMBER "
                        f"topology...")
                    mlig_amb_parm.write_parm(self.mutant_ligand_pmrtop)

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

        p1 = subprocess.Popen(tleap_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_subprocess_output(p1)
        if p1.wait():
            GMXMMPBSA_ERROR('%s failed when querying %s' % (tleap, self.FILES.prefix + arg1))

    def _set_com_order(self, REC, LIG):
        result = []
        l_idx = 0
        r_idx = 0
        for e in self.orderl:
            if e in ['R', 'REC']:
                result.append(REC[r_idx])
                r_idx += 1
            else:
                result.append(LIG[l_idx])
                l_idx += 1
        return result
