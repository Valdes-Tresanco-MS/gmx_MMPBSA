import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--receptor', help="PDB file for receptor. Make sure that pdb file contain only one chain",
                    required=True)
parser.add_argument('-l', '--ligand', help="PDB file for ligand. Make sure that pdb file contain only one chain",
                    required=True)
parser.add_argument('-d', '--dec_file', help="Decomposition file", required=True)
parser.add_argument('-c', '--cutoff', help="Energy cutoff for representation", type=int)
parser.add_argument('-o', '--out', help="Output PDB file for receptor+ligand with Total Energy values in b-factor "
                                        "place.", default='perresidue_energy.pdb')
args = parser.parse_args()

class DECOMP:
    def __init__(self, file):
        self.file = open(file)
        self.read = False
        self.receptor_list = []
        self.ligand_list = []
        self.extract()

    def extract(self):
        for line in self.file:
            line = line.strip('\n')
            if not self.read:
                if '----------------------------------' in line:
                    self.read = True
                continue
            all_line = line.split('|')
            if len(all_line) < 8:
                continue
            ftype = all_line[1].split()[0]

            if len(all_line[0].split()) != 2:
                res = all_line[0][:3].strip()
                num = all_line[0][3:].strip()
            else:
                res, num = all_line[0].split()
            vdw = all_line[3].split('+/-')[0].strip()
            ele = all_line[4].split('+/-')[0].strip()
            psol = all_line[5].split('+/-')[0].strip()
            npsol = all_line[6].split('+/-')[0].strip()
            total = float(all_line[7].split('+/-')[0].strip())

            if total < args.cutoff and total > args.cutoff * -1 :
                total = '0.00'
            else:
                total = '%.2f' % total
            if res == 'WAT':
                res = 'HOH'
            if 'R' in ftype:
                self.receptor_list.append([res, num, vdw, ele, psol, npsol, total])
            else:
                self.ligand_list.append([res, num, vdw, ele, psol, npsol, total])


class PDBMap:
    def __init__(self, pdb_file):
        self.pdb_map = dict()
        self.pdb_file = pdb_file
        self.chains = []
        self.pdb_list = self.pdb2map()
        self.receptor = None
        self.ligand = None
        self.element = None

    def pdb2map(self):
        file = open(self.pdb_file, 'r')
        i = 1
        for line in file:
            line = line.strip('\n')
            if not (re.search("ATOM", line[0:6]) or re.search("HETATM", line[0:6])):
                continue
            if not i in self.pdb_map:
                self.pdb_map[i] = dict()
            self.pdb_map[i]["id"] = line[0:6].strip()
            self.pdb_map[i]["atom_number"] = line[6:11].strip()
            self.pdb_map[i]["atom_name"] = line[12:16]
            self.pdb_map[i]["alternate_location"] = line[16]
            self.pdb_map[i]["three_letter_code"] = line[17:21].strip()
            if self.pdb_map[i]["three_letter_code"] == 'WAT':
                self.pdb_map[i]["three_letter_code"] = 'HOH'
            if self.pdb_map[i]["three_letter_code"] == 'HOH':
                self.pdb_map[i]["id"] = 'HETATM'
            self.pdb_map[i]["chain"] = line[21].strip()
            self.pdb_map[i]["residue_number"] = line[22:26].strip()
            self.pdb_map[i]["i_code"] = line[26]
            self.pdb_map[i]["x"] = line[27:38].strip()
            self.pdb_map[i]["y"] = line[38:46].strip()
            self.pdb_map[i]["z"] = line[46:54].strip()
            self.pdb_map[i]["occupancy"] = line[54:60].strip()
            self.pdb_map[i]["b_factor"] = line[60:66].strip()
            self.element = [line[76:78].strip()][0]
            if len(self.element) == 2:
                if re.search(r'[A-Za-z]', self.element[0]) and re.search(r'[A-Za-z]', self.element[1]):
                    self.pdb_map[i]["element"] = self.element[0] + self.element[1]
                else:
                    self.pdb_map[i]["element"] = self.element[0]
            elif len(self.element) == 3:
                if re.search(r'[A-Za-z]', self.element[0]) and re.search(r'[A-Za-z]', self.element[1]):
                    self.pdb_map[i]["element"] = self.element[0] + self.element[1]
                elif re.search(r'[A-Za-z]', self.element[1]) and re.search(r'[A-Za-z]', self.element[2]):
                    self.pdb_map[i]["element"] = self.element[1] + self.element[2]
                else:
                    self.pdb_map[i]["element"] = self.element[0]
            else:
                self.pdb_map[i]["element"] = self.element
            self.pdb_map[i]["charge"] = [line[76:81].strip()][1:2]
            if not self.pdb_map[i]['chain'] in self.chains:
                self.chains.append(self.pdb_map[i]['chain'])
            i += 1
        file.close()
        if len(self.chains) > 1:
            raise IOError('Check that %s contain only one chain' % self.pdb_file)
        return self.pdb_map

    def writePDB(self, line_num):
        # Create the PDB line.
        line = (self.pdb_map[line_num]['id']).ljust(6) + \
               (self.pdb_map[line_num]['atom_number']).rjust(5) + \
               ' ' + \
               ((self.pdb_map[line_num]['atom_name']).ljust(3)).rjust(4) + \
               " " + \
               (self.pdb_map[line_num]['three_letter_code']).rjust(3) + \
               " " + \
               (self.pdb_map[line_num]["chain"]) + \
               (self.pdb_map[line_num]['residue_number']).rjust(4) + \
               (self.pdb_map[line_num]['x']).rjust(11) + \
               (self.pdb_map[line_num]['y']).rjust(8) + \
               (self.pdb_map[line_num]['z']).rjust(8) + \
               '  ' + \
               (self.pdb_map[line_num]["occupancy"]).rjust(4) + \
               (self.pdb_map[line_num]['b_factor']).rjust(6) + \
               (self.pdb_map[line_num]["element"]).rjust(12)
        return line


if args.cutoff and args.cutoff < 0 :
    raise IOError('Cutoff value must be positive integer')
dec = DECOMP(args.dec_file)
out_file = open(args.out, 'w')

rec = PDBMap(args.receptor)
lig = PDBMap(args.ligand)

match = []
f = 0
for x in rec.pdb_list:
    rec.pdb_map[x]['chain'] = 'R'
    b = False
    for d in dec.receptor_list:
        f += 1
        if d[0][:2] in rec.pdb_map[x]['three_letter_code'] and d[1] == rec.pdb_map[x]['residue_number']:
            rec.pdb_map[x]['b_factor'] = d[-1]
            if d not in match:
                match.append(d)
            b = True
            break
    if not b:
        rec.pdb_map[x]['b_factor'] = '0.00'
if len(match) != len(dec.receptor_list):
    print(match)
    print(dec.receptor_list)
    raise IOError('Mismatch between given receptor and receptor information in decomposition file ')

match = []
f = 0
for x in lig.pdb_list:
    lig.pdb_map[x]['chain'] = 'L'
    b = False
    for d in dec.ligand_list:
        f += 1
        if d[0][:2] in lig.pdb_map[x]['three_letter_code'] and d[1] == lig.pdb_map[x]['residue_number']:
            lig.pdb_map[x]['b_factor'] = d[-1]
            if d not in match:
                match.append(d)
            b = True
            break
    if not b:
        lig.pdb_map[x]['b_factor'] = '0.00'
if len(match) != len(dec.ligand_list):
    raise IOError('Mismatch between given ligand and ligand information in decomposition file ')

out_file = open('perresidue_energy.pdb', 'w')
out_file.write('REMARK   PER-RESIDUE ENERGY\n')
out_file.write('REMARK   RECEPTOR\n')
for at in rec.pdb_map:
    line = rec.writePDB(at) + '\n'
    out_file.write(line)
out_file.write('REMARK   LIGAND\n')
for at in lig.pdb_map:
    line = lig.writePDB(at) + '\n'
    out_file.write(line)
out_file.close()
