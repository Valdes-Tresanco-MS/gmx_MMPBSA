"""
Make Amber topology files from Gromacs
"""

# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/GMX-MMGBSA                  #
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

import re
import os
import parmed
import warnings
from GMXMMPBSA.exceptions import *
import subprocess
from openbabel import openbabel

ff_list = {'amber03': 'oldff/leaprc.ff03', 'amber99': 'oldff/leaprc.ff99', 'amber99sb': 'oldff/leaprc.ff99SB',
           'amber99sb-ildn': 'oldff/leaprc.ffSBildn', 'amber94': 'oldff/leaprc.ff94', 'amber96': 'oldff/leaprc.ff96',
           'amber14sb': 'leaprc.protein.ff14SB'}
for x in ff_list:
    print x
lig_ff = ['gaff', 'gaff2']
# print os.path.splitext('../COM.top')

babel_elements = {0: 'Xx', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
                  11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
                  21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
                  31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
                  41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
                  51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
                  61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
                  71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
                  81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
                  91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
                  101: 'Md', 102: 'No', 106: 'Lr', 107: 'D'}

std_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS',
          'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIS']

gmx = '/usr/local/gromacs/bin/gmx'


# FIXME: check if gromacs exists (findprog module)

class CheckMakeTop:
    def __init__(self, FILES, alarun=False):
        self.FILES = FILES
        self.ligand_isProt = True
        mutprefix = ''
        # create local variables based on mutant to avoid duplicate
        self.mutation = alarun

        self.ligand_tpr = None
        self.ligand_mol2 = None

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

        self.mutant_complex_pdb = self.FILES.prefix + 'MUT_COM.pdb'
        self.mutant_receptor_pdb = self.FILES.prefix + 'MUT_REC.pdb'
        self.mutant_ligand_pdb = self.FILES.prefix + 'MUT_LIG.pdb'
        self.mutant_complex_pdb_fixed = self.FILES.prefix + 'MUT_COM_FIXED.pdb'
        self.mutant_receptor_pdb_fixed = self.FILES.prefix + 'MUT_REC_FIXED.pdb'
        self.mutant_ligand_pdb_fixed = self.FILES.prefix + 'MUT_LIG_FIXED.pdb'
        # self.default_ff = 'leaprc.protein.ff14SB'

        self.getPDBfromTpr()
        self.checkPDB()

    def getPDBfromTpr(self):
        """
        Get PDB file to make topology
        :return:
        """
        # wt complex
        # make index for extract pdb structure
        rec_group, lig_group = self.FILES.complex_groups

        print 'WildType Complex: Save group {}_{} in {} (gromacs index) file as {}'.format(rec_group, lig_group,
                                                                                           self.FILES.complex_index,
                                                                                           self.complex_pdb)
        # merge both (rec and lig) groups into complex group, modify index and create a copy
        # 1-rename groups, 2-merge
        c1 = subprocess.Popen(['echo', 'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | {l}\n'
                                       ' q\n'.format(r=rec_group, l=lig_group)], stdout=subprocess.PIPE)
        # FIXME: overwrite the user index file???
        com_ndx = self.FILES.prefix + 'COM_index.ndx'

        c2 = subprocess.Popen([gmx, "make_ndx", '-n', self.FILES.complex_index, '-o', com_ndx],
                              stdin=c1.stdout, stdout=subprocess.PIPE)
        if c2.wait():  # if it quits with return code != 0
            raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.receptor_tpr))
        self.FILES.complex_index = com_ndx

        c3 = subprocess.Popen(['echo', 'GMXMMPBSA_REC_GMXMMPBSA_LIG'], stdout=subprocess.PIPE)
        print 'echo', '"GMXMMPBSA_REC_GMXMMPBSA_LIG"'
        # we get only first trajectory to extract a pdb file and make amber topology for complex
        c4 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr,
                               '-o', self.complex_pdb, '-n', self.FILES.complex_index, '-b', '0', '-e', '0'],
                              stdin=c3.stdout, stdout=subprocess.PIPE)
        if c4.wait():  # if it quits with return code != 0
            raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.complex_tpr))

        # check lig before extract from complex
        if os.path.splitext(self.FILES.ligand_tprOmol2)[1] == '.tpr':
            self.ligand_isProt = True
            self.ligand_tpr = self.FILES.ligand_tprOmol2
        else:
            self.ligand_isProt = False
            self.ligand_mol2 = self.FILES.ligand_tprOmol2

        # Put receptor and ligand (explicitly defined) to avoid overwrite them
        # check if ligand is not protein. In any case, non-protein ligand always most be processed
        if not self.ligand_isProt:
            lig_name = os.path.splitext(os.path.split(self.ligand_mol2)[1])[0]
            self.ligand_frcmod = self.FILES.prefix + lig_name + '.frcmod'
            # run parmchk2
            l3 = subprocess.Popen(['parmchk2', '-i', self.ligand_mol2, '-f', 'mol2', '-o', self.ligand_frcmod],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if l3.wait():
                raise MMPBSA_Error('%s failed when querying %s' % ('parmchk2', self.ligand_mol2))

        # check if stability
        if self.FILES.stability:
            if (self.FILES.receptor_tpr or self.ligand_tpr):
                warnings.warn(
                    'When Stability calculation mode is selected receptor and ligand are not needed. However, '
                    'the receptor and/or the ligand are defined, so we will ignore them.', StabilityWarning)
            if self.mutation and (self.FILES.mutant_receptor_tpr or self.FILES.mutant_ligand_tpr):
                warnings.warn(
                    'When Stability calculation mode is selected mutant receptor and/or mutant ligand are not '
                    'needed. However, the receptor or the ligand (mutant) are defined, so we will ignore them.',
                    StabilityWarning)
            return

        # wt receptor
        if self.FILES.receptor_tpr:
            print 'Save group {} in {} (gromacs index) file as {}'.format(self.FILES.receptor_group,
                                                                          self.FILES.receptor_index,
                                                                          self.receptor_pdb)
            print 'Force to use this receptor structure instead the generated from complex'
            p1 = subprocess.Popen(['echo', '{}'.format(rec_group)], stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file for make amber topology
            cp2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.receptor_trajs[0], '-s', self.FILES.receptor_tpr,
                                    '-o', self.receptor_pdb, '-n', self.FILES.receptor_index, '-b', '0', '-e', '0'],
                                   stdin=p1.stdout, stdout=subprocess.PIPE)
            if cp2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.receptor_tpr))
        else:
            print 'Using receptor structure from complex to make amber topology'
            # wt complex receptor
            print 'Complex: Save group {} in {} (gromacs index) file as {}'.format(rec_group, self.FILES.complex_index,
                                                                                   self.receptor_pdb)
            cp1 = subprocess.Popen(['echo', '{}'.format(rec_group)], stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file for make amber topology
            cp2 = subprocess.Popen(
                [gmx, "trjconv", '-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                 self.receptor_pdb, '-n', self.FILES.complex_index, '-b', '0', '-e', '0'],
                stdin=cp1.stdout, stdout=subprocess.PIPE)
            if cp2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.complex_tpr))

        # ligand
        # # check consistence
        if self.ligand_tpr:  # ligand is protein
            # wt ligand
            l1 = subprocess.Popen(['echo', '{}'.format(self.FILES.ligand_group)], stdout=subprocess.PIPE)
            # we get only first trajectory for extract a pdb file for make amber topology
            l2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.ligand_trajs[0], '-s',
                                   self.FILES.ligand_tpr, '-o', self.ligand_pdb, '-b', '0', '-e', '0'],
                                  stdin=l1.stdout, stdout=subprocess.PIPE)
            if l2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.ligand_tpr))
        elif self.ligand_mol2:
            # done above
            pass
        else:
            # wt complex ligand
            print 'Using ligand structure from complex to make amber topology'
            print 'Save group {} in {} (gromacs index) file as {}'.format(lig_group, self.FILES.complex_index,
                                                                          self.ligand_pdb)
            cl1 = subprocess.Popen(['echo', '{}'.format(lig_group)], stdout=subprocess.PIPE)
            # we get only  first trajectory to extract a pdb file for make amber topology
            cl2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr,
                                    '-o', self.ligand_pdb, '-n', self.FILES.complex_index, '-b', '0', '-e', '0'],
                                   stdin=cl1.stdout, stdout=subprocess.PIPE)
            if cl2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.complex_tpr))

    def checkPDB(self):
        """
        Get parmed structure object for complex, receptor and ligand if is protein-like

        1 - Rename HIS
        2 - Rename CYS
        3 - Delete H
        4 - Rename oxygen in termini from GROMACS to AMBER name
          - Rename CD in ILE from GROMACS to AMBER name
        5 - Save
        :return:
        """
        self.complex_str = parmed.read_PDB(self.complex_pdb)  # can always be initialized
        # fix complex structure and save
        self.properHIS(self.complex_str)
        self.properCYS(self.complex_str)
        # For some reason removing the hydrogens returns the hydrogen-bound atoms to their original names. This is
        # problematic with ILE switching from CD to CD1. parmed bug?
        self.complex_str.strip('@/H')
        self.properATOMS(self.complex_str)
        self.complex_str.save(self.complex_pdb_fixed, 'pdb', True)

        if self.mutation:
            self.mutant_complex_str = parmed.read_PDB(self.complex_pdb_fixed)
            # make mutation and save
            self.mutatexala(self.mutant_complex_str)
            self.mutant_complex_str.save(self.mutant_complex_pdb_fixed, 'pdb')

        if not self.FILES.stability:
            self.receptor_str = parmed.read_PDB(self.receptor_pdb)
            # fix receptor structure
            self.properHIS(self.receptor_str)
            self.properCYS(self.receptor_str)
            self.receptor_str.strip('@/H')
            self.properATOMS(self.receptor_str)
            self.receptor_str.save(self.receptor_pdb_fixed, 'pdb', True)

            # fix ligand structure
            if self.ligand_isProt and not self.FILES.ligand_tpr:  # ligand from complex
                # check if ligand (from complex structure) is really protein-like.
                self.ligand_str = parmed.read_PDB(self.ligand_pdb)
                for res in self.ligand_str.residues:
                    if res.name not in std_aa:
                        self.ligand_isProt = False
                        raise MMPBSA_Error(
                            'It appears that the ligand that defined based on complex is non-protein type. '
                            'This ligand type requires a structure (mol2) and a parameter (frcmod) files. '
                            'Please define these parameters to perform the calculation correctly.')

                # fix ligand structure if is protein
                self.properHIS(self.ligand_str)
                self.properCYS(self.ligand_str)
                self.ligand_str.strip('@/H')
                self.properATOMS(self.ligand_str)
                self.ligand_str.save(self.ligand_pdb_fixed, 'pdb', True)

            if self.mutation:
                if self.FILES.mutant == 'REC':
                    self.mutant_receptor_str = parmed.read_PDB(self.receptor_pdb_fixed)
                    # fix mutant receptor structure
                    self.mutatexala(self.mutant_receptor_str)
                    self.mutant_receptor_str.save(self.mutant_receptor_pdb_fixed, 'pdb', True)

                elif self.FILES.mutant == 'LIG':
                    if not self.ligand_isProt:
                        raise MMPBSA_Error('Mutation is only possible if the ligand is protein-like')
                    self.mutant_ligand_str = parmed.read_PDB(self.ligand_pdb_fixed)
                    self.mutatexala(self.mutant_ligand_str)
                    self.mutant_ligand_str.save(self.mutant_ligand_pdb_fixed, 'pdb', True)

    def mutatexala(self, structure):
        idx = 0
        found = False
        if not self.FILES.mutant_residue:
            raise MMPBSA_Error("No residue for mutation was defined")
        chain, resnum = self.FILES.mutant_residue.split(':')

        print '#@#@#@#@#', chain, resnum

        if not chain or not resnum:
            raise MMPBSA_Error("No residue was defined")
        for res in structure.residues:
            print res.name, res.chain, res.number
            if res.number == int(resnum) and res.chain == chain:
                found = True
                print 'encontrado', res.name, res.chain, res.number
                break
            idx += 1
        if found:
            structure.residues[idx].name = 'ALA'
            print structure.residues[idx].name
            excluded_mask = ':{} &!@CB,C,CA,N,O'.format(idx + 1)
            print excluded_mask
            structure.strip(excluded_mask)
            print [atm.name for atm in structure.residues[idx]]
        else:
            raise MMPBSA_Error('Residue {}:{} not found'.format(chain, resnum))

    def properATOMS(self, structure):
        """
        Rename oxygen in termini from GROMACS to AMBER name
        OC1 -> 'O  '
        OC2 -> OXT
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

    def properHIS(self, structure):
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

    def properCYS(self, structure):
        """
        Rename the cys that form a disulfide bond
        :return:
        """
        cys_name = ['CYS', 'CYX', 'CYM']
        allcys = [residue for residue in structure.residues if residue.name in cys_name]
        # print allcys
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

    def checkForceField(self):
        if self.FILES.protein_ff not in ff_list:
            raise ValueError('This forcefield {} does not match any of the allowed '
                             '({})'.format(self.FILES.protein_ff, ', '.join([x for x in ff_list])))
        if not self.ligand_isProt and self.FILES.ligand_ff not in lig_ff:
            raise ValueError('This forcefield {} does not match any of the allowed '
                             '({})'.format(self.FILES.ligand_ff, ', '.join([x for x in ff_list])))

    def makeToptleap(self):
        self.checkForceField()
        with open(self.FILES.prefix + 'leap.in', 'w') as tif:
            tif.write('source {}\n'.format(ff_list[self.FILES.protein_ff]))
            tif.write('source leaprc.DNA.bsc1\n')
            tif.write('source leaprc.RNA.OL3\n')
            tif.write('source leaprc.{}\n'.format(self.FILES.ligand_ff))
            tif.write('set default PBRadii mbondi2\n')
            # check if ligand is not protein and always load
            if not self.ligand_isProt:
                tif.write('LIG = loadmol2 {}\n'.format(self.ligand_mol2))
                tif.write('check LIG\n')
                tif.write('loadamberparams {}\n'.format(self.ligand_frcmod))

            if not self.FILES.stability:
                tif.write('REC = loadpdb {}\n'.format(self.receptor_pdb_fixed))
                tif.write('saveamberparm REC {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop, p=self.FILES.prefix))
                if self.ligand_isProt:
                    tif.write('LIG = loadpdb {}\n'.format(self.ligand_pdb_fixed))
                tif.write('saveamberparm LIG {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop, p=self.FILES.prefix))

            tif.write('complex = loadpdb {}\n'.format(self.complex_pdb_fixed))
            tif.write('saveamberparm complex {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            tif.write('quit')

        tleap = '/home/mario/programs/amber18/bin/tleap'
        p = subprocess.check_output('{t} -f {f}'.format(t=tleap, f=self.FILES.prefix + 'leap.in'),
                                    stderr=subprocess.STDOUT, shell=True)
        print p.decode()

        if self.mutation:
            with open(self.FILES.prefix + 'mut_leap.in', 'w') as mtif:
                mtif.write('source {}\n'.format(ff_list[self.FILES.protein_ff]))
                mtif.write('source leaprc.DNA.bsc1\n')
                mtif.write('source leaprc.RNA.OL3\n')
                mtif.write('source leaprc.{}\n'.format(self.FILES.ligand_ff))
                mtif.write('set default PBRadii mbondi2\n')
                # check if ligand is not protein and always load
                if not self.ligand_isProt:
                    mtif.write('LIG = loadmol2 {}\n'.format(self.ligand_mol2))
                    mtif.write('check LIG\n')
                    mtif.write('loadamberparams {}\n'.format(self.ligand_frcmod))

                if not self.FILES.stability:
                    if self.FILES.mutant == 'REC':
                        mtif.write('mut_rec = loadpdb {}\n'.format(self.mutant_receptor_pdb_fixed))
                        mtif.write('saveamberparm mut_rec {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
                                                                                          p=self.FILES.prefix))
                        self.mutant_ligand_pmrtop = None
                    else:
                        mtif.write('mut_lig = loadpdb {}\n'.format(self.mutant_ligand_pdb_fixed))
                        self.mutant_receptor_pmrtop = self.receptor_pmrtop
                        mtif.write('saveamberparm mut_lig {t} {p}MUT_LIG.inpcrd\n'.format(t=self.mutant_ligand_pmrtop,
                                                                                          p=self.FILES.prefix))
                        self.mutant_receptor_pmrtop = None
                mtif.write('mut_com = loadpdb {}\n'.format(self.mutant_complex_pdb_fixed))
                mtif.write('saveamberparm mut_com {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                                  p=self.FILES.prefix))
                mtif.write('quit')

            p1 = subprocess.check_output('{t} -f {f}'.format(t=tleap, f=self.FILES.prefix + 'mut_leap.in'),
                                         stderr=subprocess.STDOUT, shell=True)
            print p1.decode()

        else:
            self.mutant_complex_pmrtop = None

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)
