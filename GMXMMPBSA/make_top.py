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
    def __init__(self, FILES, mutant=False):
        self.FILES = FILES
        self.ligand_isProt = True
        mutprefix = ''
        # create local variables based on mutant to avoid duplicate
        self.mutant = mutant

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
        # merge both (rec and lig) groups into complex group
        # modify index and create a copy
        # 1-rename groups, 2-merge

        c1 = subprocess.Popen(['echo', 'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | {l}\n'
                                       ' q\n'.format(r=rec_group, l=lig_group)], stdout=subprocess.PIPE)
        print 'echo', 'name {r} GMXMMPBSA_REC\n name {l} GMXMMPBSA_LIG\n  {r} | {l}\n q\n'.format(r=rec_group,
                                                                                                  l=lig_group)
        # FIXME: overwrite the user index file???
        com_ndx = self.FILES.prefix + 'COM_index.ndx'

        c2 = subprocess.Popen([gmx, "make_ndx", '-n', self.FILES.complex_index, '-o', com_ndx ],
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

        # wt complex receptor
        print 'Complex: Save group {} in {} (gromacs index) file as {}'.format(rec_group, self.FILES.complex_index,
                                                                               self.receptor_pdb)
        cp1 = subprocess.Popen(['echo', '{}'.format(rec_group)], stdout=subprocess.PIPE)
        # we get only first trajectory to extract a pdb file for make amber topology
        cp2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr, '-o',
                                self.receptor_pdb, '-n', self.FILES.complex_index, '-b', '0', '-e', '0'],
                               stdin=cp1.stdout, stdout=subprocess.PIPE)
        if cp2.wait():  # if it quits with return code != 0
            raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.complex_tpr))

        # wt complex ligand
        print 'Save group {} in {} (gromacs index) file as {}'.format(lig_group, self.FILES.complex_index,
                                                                      self.ligand_pdb)
        cl1 = subprocess.Popen(['echo', '{}'.format(lig_group)], stdout=subprocess.PIPE)
        # we get only  first trajectory to extract a pdb file for make amber topology
        cl2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.complex_trajs[0], '-s', self.FILES.complex_tpr,
                                '-o', self.ligand_pdb, '-n', self.FILES.complex_index, '-b', '0', '-e', '0'],
                               stdin=cl1.stdout, stdout=subprocess.PIPE)
        if cl2.wait():  # if it quits with return code != 0
            raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.complex_tpr))

        # ============================ MUTANT =================================================================
        if self.mutant:
            # Mut complex
            mrec_group, mlig_group = self.FILES.mutant_complex_groups
            print '''Mutant Complex: Save group {}_{} in {} (gromacs index) file as 
                 {}'''.format(mrec_group, mlig_group, self.FILES.mutant_complex_index, self.mutant_complex_pdb)
            # merge both (rec and lig) groups into mutant complex group
            c1 = subprocess.Popen(['echo', '"{r}" | "{l}"\n q\n'.format(r=mrec_group, l=mlig_group)],
                                  stdout=subprocess.PIPE)
            c2 = subprocess.Popen([gmx, "make_ndx", '-f', self.FILES.mutant_complex_tpr, '-n',
                                   self.FILES.mutant_complex_index], stdin=c1.stdout, stdout=subprocess.PIPE)
            if c2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_complex_tpr))

            c3 = subprocess.Popen(['echo', '"{}_{}"'.format(mrec_group, mlig_group)], stdout=subprocess.PIPE)
            # we get only first trajectory to extract a pdb file and make amber topology for complex
            c4 = subprocess.Popen(
                [gmx, "trjconv", '-f', self.FILES.mutant_complex_trajs[0], '-s', self.FILES.mutant_complex_tpr,
                 '-o', self.mutant_complex_pdb, '-n', self.FILES.mutant_complex_index, '-b', '0', '-e',
                 '0'],
                stdin=c3.stdout, stdout=subprocess.PIPE)
            if c4.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_complex_tpr))

            # # mutant complex receptor
            # print 'Complex: Save group {} in {} (gromacs index) file as {}'.format(mrec_group,
            #                                                                        self.FILES.mutant_complex_index,
            #                                                                        self.mutant_receptor_pdb)
            # cp1 = subprocess.Popen(['echo', '{}'.format(mrec_group)], stdout=subprocess.PIPE)
            # # we get only first trajectory to extract a pdb file for make amber topology
            # cp2 = subprocess.Popen(
            #     [gmx, "trjconv", '-f', self.FILES.mutant_complex_trajs[0], '-s', self.FILES.mutant_complex_tpr, '-o',
            #      self.mutant_receptor_pdb, '-n', self.FILES.mutant_complex_index, '-b', '0', '-e', '0'],
            #     stdin=cp1.stdout, stdout=subprocess.PIPE)
            # if cp2.wait():  # if it quits with return code != 0
            #     raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_complex_tpr))
            #
            # # mutant complex ligand
            # print 'Save group {} in {} (gromacs index) file as {}'.format(mlig_group, self.FILES.mutant_complex_index,
            #                                                               self.mutant_ligand_pdb)
            # cl1 = subprocess.Popen(['echo', '{}'.format(mlig_group)], stdout=subprocess.PIPE)
            # # we get only  first trajectory to extract a pdb file for make amber topology
            # cl2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.mutant_complex_trajs[0], '-s',
            #                         self.FILES.mutant_complex_tpr, '-o', self.mutant_ligand_pdb, '-n',
            #                         self.FILES.mutant_complex_index, '-b', '0', '-e', '0'], stdin=cl1.stdout,
            #                        stdout=subprocess.PIPE)
            # if cl2.wait():  # if it quits with return code != 0
            #     raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_complex_tpr))
        # =============================================================================================

        # check if stability
        if self.FILES.stability:
            if (self.FILES.receptor_tpr or self.FILES.ligand_tpr):
                warnings.warn(
                    'When Stability calculation mode is selected receptor and ligand are not needed. However, '
                    'the receptor and/or the ligand are defined, so we will ignore them.', StabilityWarning)
            if self.mutant and (self.FILES.mutant_receptor_tpr or self.FILES.mutant_ligand_tpr):
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

        # ligand
        # check consistence
        if self.FILES.ligand_tpr and (self.FILES.ligand_mol2 or self.FILES.ligand_frcmod):
            raise MMPBSA_Error('Inconsistencies found. Both the tpr and the (mol2 or frcmod) of the ligand were '
                               'defined. Note that tpr is used for protein-like molecules, while mol2 and frcmod are '
                               'result of the parametrization of a small molecule in the antechamber. Define one of '
                               'them according to your system.')

        if self.FILES.ligand_tpr:  # ligand is protein
            l1 = subprocess.Popen(['echo', '{}'.format(self.FILES.ligand_group)], stdout=subprocess.PIPE)
            # we get only first trajectory for extract a pdb file for make amber topology
            l2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.ligand_trajs[0], '-s',
                                   self.FILES.ligand_tpr, '-o', self.ligand_pdb, '-b', '0', '-e', '0'],
                                  stdin=l1.stdout, stdout=subprocess.PIPE)
            if l2.wait():  # if it quits with return code != 0
                raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.ligand_tpr))
        elif self.FILES.ligand_mol2 and self.FILES.ligand_frcmod:
            self.ligand_isProt = False
        else:
            print 'Using ligand structure from complex to make amber topology'

        # ------------------------------------------------------------------------------------------------------
        if self.mutant:
            # mutant recept
            if self.FILES.mutant_receptor_tpr:
                print '''MUTANT RECEPTOR: Save group {} in {} (gromacs index) file as 
                {}'''.format(self.FILES.mutant_receptor_group, self.FILES.mutant_receptor_index, self.mutant_receptor_pdb)

                p1 = subprocess.Popen(['echo', '{}'.format(rec_group)], stdout=subprocess.PIPE)
                # we get only first trajectory to extract a pdb file for make amber topology
                cp2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.mutant_receptor_trajs[0], '-s',
                                        self.FILES.mutant_receptor_tpr, '-o', self.mutant_receptor_pdb, '-n',
                                        self.FILES.mutant_receptor_index, '-b', '0', '-e', '0'], stdin=p1.stdout,
                                       stdout=subprocess.PIPE)
                if cp2.wait():  # if it quits with return code != 0
                    raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_receptor_tpr))
            else:
                print 'Using receptor structure from complex to make amber topology'
                self.mutant_receptor_pdb = self.receptor_pdb
            # mutant ligand
            # check consistence
            if self.FILES.mutant_ligand_tpr and not self.ligand_isProt:
                raise MMPBSA_Error('Inconsistencies found. The mutant ligand tpr were defined but the ligand bound in '
                                   'normal complex is non-protein type. The mutation is only possible for protein type.')

            if self.FILES.mutant_ligand_tpr:  # ligand is protein
                l1 = subprocess.Popen(['echo', '{}'.format(self.FILES.mutant_ligand_group)], stdout=subprocess.PIPE)
                # we get only first trajectory for extract a pdb file for make amber topology
                l2 = subprocess.Popen([gmx, "trjconv", '-f', self.FILES.mutant_ligand_trajs[0], '-s',
                                       self.FILES.mutant_ligand_tpr, '-o', self.mutant_ligand_pdb, '-b', '0', '-e', '0'],
                                      stdin=l1.stdout, stdout=subprocess.PIPE)
                if l2.wait():  # if it quits with return code != 0
                    raise MMPBSA_Error('%s failed when querying %s' % (gmx + 'make_ndx', self.FILES.mutant_ligand_tpr))
            else:
                print 'Using ligand structure from complex to make amber topology'
                # self.mutant_ligand_pdb = self.ligand_pdb # I can't be sure what kind of type it is

    def checkPDB(self):
        """
        Get parmed structure object for complex, receptor and ligand if is protein-like

        1 - Rename HIS
        2 - Rename CYS
        3 - Rename oxygen in termini from GROMACS to AMBER name
          - Rename CD in ILE from GROMACS to AMBER name
        4 - Delete hydrogen atoms
        5 - Save
        :return:
        """
        self.complex_str = parmed.read_PDB(self.complex_pdb) # can always be initialized
        # fix complex structure and save
        self.properHIS(self.complex_str)
        self.properCYS(self.complex_str)
        # For some reason removing the hydrogens returns the hydrogen-bound atoms to their original names. This is
        # problematic with ILE switching from CD to CD1. parmed bug?
        self.complex_str.strip('@/H') #
        self.properATOMS(self.complex_str)
        self.complex_str.save(self.complex_pdb_fixed, 'pdb', True)

        if self.mutant:
            self.mutant_complex_str = parmed.read_PDB(self.mutant_complex_pdb)
            # fix mutant complex structure and save
            self.properHIS(self.mutant_complex_str)
            self.properCYS(self.mutant_complex_str)
            self.mutant_complex_str.strip('@/H')
            self.properATOMS(self.mutant_complex_str)
            self.mutant_complex_str.save(self.mutant_complex_pdb_fixed, 'pdb', True)

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

            if self.mutant:
                if self.mutant_receptor_pdb == self.receptor_pdb: # ligand mutant
                    self.mutant_receptor_pdb_fixed = self.receptor_pdb_fixed
                    self.mutant_part = 'LIG'
                    if self.ligand_isProt:
                        self.mutant_ligand_str = parmed.read_PDB(self.mutant_ligand_pdb)
                        self.properHIS(self.mutant_ligand_str)
                        self.properCYS(self.mutant_ligand_str)
                        self.mutant_ligand_str.strip('@/H')
                        self.properATOMS(self.mutant_ligand_str)
                        self.mutant_ligand_str.save(self.ligand_pdb_fixed, 'pdb', True)

                else: # receptor mutant
                    if self.ligand_isProt:
                        self.mutant_part = 'REC'
                        self.mutant_ligand_pdb_fixed = self.ligand_pdb_fixed
                        self.mutant_receptor_str = parmed.read_PDB(self.mutant_receptor_pdb)
                        # fix mutant receptor structure
                        self.properHIS(self.mutant_receptor_str)
                        self.properCYS(self.mutant_receptor_str)
                        self.mutant_receptor_str.strip('@/H')
                        self.properATOMS(self.mutant_receptor_str)
                        self.mutant_receptor_str.save(self.mutant_receptor_pdb_fixed, 'pdb', True)



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
            tif.write('com = loadpdb {}\n'.format(self.complex_pdb_fixed))
            tif.write('saveamberparm com {t} {p}COM.inpcrd\n'.format(t=self.complex_pmrtop, p=self.FILES.prefix))
            if self.mutant:
                tif.write('mut_com = loadpdb {}\n'.format(self.mutant_complex_pdb_fixed))
                tif.write('saveamberparm mut_com {t} {p}MUT_COM.inpcrd\n'.format(t=self.mutant_complex_pmrtop,
                                                                           p=self.FILES.prefix))
            else:
                self.mutant_complex_pmrtop = None # FIXME: needed?
            if not self.FILES.stability:
                tif.write('rec = loadpdb {}\n'.format(self.receptor_pdb_fixed))
                tif.write('saveamberparm rec {t} {p}REC.inpcrd\n'.format(t=self.receptor_pmrtop, p=self.FILES.prefix))
                if self.mutant and self.mutant_part == 'REC':
                    tif.write('mut_rec = loadpdb {}\n'.format(self.mutant_receptor_pdb_fixed))
                    tif.write('saveamberparm mut_rec {t} {p}MUT_REC.inpcrd\n'.format(t=self.mutant_receptor_pmrtop,
                                                                            p=self.FILES.prefix))
                    self.mutant_ligand_pmrtop = None
                if self.ligand_isProt:
                    tif.write('lig = loadpdb {}\n'.format(self.ligand_pdb_fixed))
                    if self.mutant and self.mutant_part == 'LIG':
                        tif.write('mut_lig = loadpdb {}\n'.format(self.mutant_ligand_pdb_fixed))
                        self.mutant_receptor_pmrtop = self.receptor_pmrtop
                        tif.write('saveamberparm mut_lig {t} {p}MUT_LIG.inpcrd\n'.format(t=self.mutant_ligand_pmrtop,
                                                                                p=self.FILES.prefix))
                else:
                    tif.write('lig = loadmol2 {}\n'.format(self.FILES.ligand_mol2))
                    tif.write('check lig\n')
                    tif.write('loadamberparams {}\n'.format(self.FILES.ligand_frcmod))
                    if self.mutant:
                        self.mutant_receptor_pmrtop = None
                tif.write('saveamberparm lig {t} {p}LIG.inpcrd\n'.format(t=self.ligand_pmrtop, p=self.FILES.prefix))

            tif.write('quit')
        print 'tleap -f {}'.format(self.FILES.prefix + 'leap.in')
        print os.path.exists(self.FILES.prefix + 'leap.in')
        tleap = '/home/mario/programs/amber18/bin/tleap'
        p = subprocess.check_output('{t} -f {f}'.format(t=tleap, f=self.FILES.prefix + 'leap.in'),
                                    stderr=subprocess.STDOUT, shell=True)
        print p.decode()

        return (self.complex_pmrtop, self.receptor_pmrtop, self.ligand_pmrtop, self.mutant_complex_pmrtop,
                self.mutant_receptor_pmrtop, self.mutant_ligand_pmrtop)


# ff = getForceField('../MMGBSA_test/topol_xtc_gro_tpr/test/COM.top')
#     print ff
# g = CheckTop('../MMGBSA_test/topol_xtc_gro_tpr/1X1U_md_free.tpr', None, None, 'amber14sb')
# # g.getPDB('/media/mario/Dynamics/4uwh/complex.pdb')
# g.getPDB('test_ssbond.pdb')
# g.getPDB('5ivo.pdb')
# g.properHIS()
# g.properCYS()

# s = parmed.read_PDB('/media/mario/Dynamics/4uwh/COM/com.pdb')
# s = parmed.read_PDB('/media/mario/Dynamics/4uwh/COM/com.pdb')
# s = parmed.read_PDB('/media/mario/Dynamics/4uwh/complex.pdb')
# s = parmed.read_PDB('../MMGBSA_test/complex.pdb')
# s = parmed.read_PDB('../MMGBSA_test/topol_xtc_gro_tpr/test/test_conect.pdb')
# s = parmed.read_PDB('/home/mario/Drive/scripts/MMGBSA/MMGBSA_test/topol_xtc_gro_tpr/test_lig_charged/test.pdb')
#
#
# for res in s.residues:
#     if res.name == 'ILE':
#         for atom in res.atoms:
#             if atom.name == 'CD':
#                 print atom.name
#                 atom.name = 'CD1'
#                 print atom.name
#
# s.save('kkkkk.pdb', 'pdb')
# print prot + lig

# obConversion = openbabel.OBConversion()
# obConversion.SetInAndOutFormats("pdb", "pdb")
#
# mol = openbabel.OBMol()
#
# obConversion.ReadFile(mol, '/media/mario/Dynamics/dynamics/alpha/4ykn/LIG/ligand.pdb')
#
# f = 1
# for obatom in openbabel.OBMolAtomIter(mol):
#     print f, obatom.GetExplicitValence(), obatom.GetAtomicNum(), obatom.GetTotalValence(), obatom.GetFormalCharge()
#     print [b.GetBondOrder() for b in openbabel.OBAtomBondIter(obatom)]
#     print obatom.GetFormalCharge()
#     print 'element', babel_elements[obatom.GetAtomicNum()]
#     # print
#     f += 1
# lig[0].save('kk.pdb', 'pdb')
# print s.split()[2][0].residues
# print s.split()
# for res in s.bonds:
# print res.charge
# print res.order
# for bond in res.bonds:
#     print bond
