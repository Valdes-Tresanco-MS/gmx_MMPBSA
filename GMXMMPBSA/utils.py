"""
This is a module that contains functions generally useful for the
gmx_MMPBSA script. A full list of functions/subroutines is shown below.
It must be included to insure proper functioning of gmx_MMPBSA

List of functions and a brief description of their purpose
-remove: Removes temporary work files in this directory. It has a number of
    different levels to remove only a small number of files.
-concatenate: combines 2 files into a single, common file
"""

# TODO get rid of this file altogether and move these functions into the main
# app class

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
import re
import shutil
from pathlib import Path
import json
import logging
from string import ascii_letters
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
from math import sqrt
import parmed


def create_input_args(args: list):
    if not args or 'all' in args:
        return 'general', 'gb', 'pb', 'ala', 'nmode', 'decomp', 'rism'
    elif 'gb' not in args and 'pb' not in args and 'rism' not in args and 'nmode' not in 'args':
        GMXMMPBSA_ERROR('You did not specify any type of calculation!')
    elif 'gb' not in args and 'pb' not in args and 'decomp' in args:
        logging.warning('decomp calculation is only compatible with gb y pb calculations. Will be ignored!')
        args.remove('decomp')
        return ['general'] + args
    else:
        return ['general'] + args


def log_subprocess_output(process):
    while True:
        output = process.stdout.readline().decode()
        if output:
            logging.debug(output.strip('\n'))
        else:
            break


class Residue(int):
    """
    Residue class
    """
    def __init__(self, index, number, chain, mol, id, name, icode=''):
        int.__init__(index)
        self.index = index
        self.number = number
        self.chain = chain
        self.mol = mol
        self.id = id
        self.name = name
        self.icode = icode

    def __new__(cls, index, number, chain, mol, id, name, icode=''):
        i = int.__new__(cls, index)
        i.index = index
        i.number = number
        i.chain = chain
        i.mol = mol
        i.name = name
        i.icode = icode
        i.string = f"{id}:{chain}:{name}:{number}:{icode}" if icode else f"{id}:{chain}:{name}:{number}"
        return i

    def __copy__(self):
        return Residue(self.index, self.number, self.chain, self.mol, self.id, self.name, self.icode)

    def __deepcopy__(self, memo):
        cls = self.__class__
        return cls.__new__(cls, self.index, self.number, self.chain, self.mol, self.id, self.name, self.icode)

    def __repr__(self):
        text = f"{type(self).__name__}(index: {self.index}, {self.id}:{self.chain}:{self.name}:{self.number}"
        if self.icode:
            text += f":{self.icode}"
        text += ')'
        return text

    def __str__(self):
        return f"{self.index}"

    def issame(self, other):
        pass


def get_indexes(com_ndx, rec_ndx=None, rec_group=1, lig_ndx=None, lig_group=1):
    ndx_files = {'COM': com_ndx, 'REC': rec_ndx, 'LIG': lig_ndx}
    ndx = {'COM': {}, 'REC': {}, 'LIG': {}}
    for n, f in ndx_files.items():
        if f is None:
            continue
        with open(f) as indexf:
            header = None
            for line in indexf:
                if line.startswith('['):
                    header = line.strip('\n[] ')
                    ndx[n][header] = []
                else:
                    ndx[n][header].extend(map(int, line.split()))
    com_indexes = {'COM': ndx['COM']['GMXMMPBSA_REC_GMXMMPBSA_LIG'], 'REC': ndx['COM']['GMXMMPBSA_REC'],
                   'LIG': [ndx['COM']['GMXMMPBSA_LIG']]}
    rec_indexes = ndx['REC'][list(ndx['REC'].keys())[rec_group]] if rec_ndx else {}
    lig_indexes = ndx['LIG'][list(ndx['LIG'].keys())[lig_group]] if lig_ndx else {}
    return {'COM': com_indexes, 'REC': rec_indexes, 'LIG': lig_indexes}


def _get_restype(resname):
    if resname == 'LYN':
        return 'LYS'
    elif resname == 'ASH':
        return 'ASP'
    elif resname == 'GLH':
        return 'GLU'
    elif resname in ['HIP', 'HIE', 'HID']:
        return 'HIS'
    elif resname in ['CYX', 'CYM']:
        return 'CYS'
    else:
        return resname


def eq_strs(struct1, struct2):
    if len(struct1.atoms) != len(struct2.atoms):
        return 'atoms', len(struct1.atoms), len(struct2.atoms)
    elif len(struct1.residues) != len(struct2.residues):
        return 'residues', len(struct1.residues), len(struct2.residues)
    else:
        return


def check_str(structure, ref=False):
    if isinstance(structure, str):
        refstr = parmed.read_PDB(structure)
    else:
        refstr = structure

    previous = 0
    ind = 1
    res_dict = {}
    duplicates = []
    for res in refstr.residues:
        if res.chain == '':
            if ref:
                GMXMMPBSA_ERROR('The reference structure used is inconsistent. The following residue does not have a '
                                f'chain ID: {res.number}:{res.name}')
            elif not previous:
                res_dict[ind] = [[res.number, res.name, res.insertion_code]]
            elif res.number - previous in [0, 1]:
                res_dict[ind].append([res.number, res.name, res.insertion_code])
            else:
                ind += 1
                res_dict[ind] = [[res.number, res.name, res.insertion_code]]
            previous = res.number
        elif res.chain not in res_dict:
            res_dict[res.chain] = [[res.number, res.name, res.insertion_code]]
        else:
            res_dict[res.chain].append([res.number, res.name, res.insertion_code])

    for chain, resl in res_dict.items():
        res_id_list = [[x, x2] for x, x1, x2 in resl]
        for c, x in enumerate(res_id_list):
            if res_id_list.count(x) > 1:
                duplicates.append(f'{chain}:{resl[c][0]}:{resl[c][1]}:{resl[c][2]}')
    if duplicates:
        GMXMMPBSA_ERROR(f'The {"reference" if ref else "complex"} structure used is inconsistent. The following '
                        f'residues are duplicates:\n {", ".join(duplicates)}')
    return refstr


def res2map(indexes, com_file):
    """
    :param com_str:
    :return:
    """
    res_list = {'REC': [], 'LIG': [], 'COM': []}
    com_len = len(indexes['COM']['COM'])
    if isinstance(com_file, parmed.Structure):
        com_str = com_file
    else:
        com_str = parmed.load_file(com_file)

    resindex = 1
    proc_res = None
    for i in range(com_len):
        res = [com_str.atoms[i].residue.chain, com_str.atoms[i].residue.number, com_str.atoms[
            i].residue.insertion_code]
        # We check who owns the residue corresponding to this atom
        if indexes['COM']['COM'][i] in indexes['COM']['REC']:
            # save residue number in the rec list
            if res != proc_res and resindex not in res_list['REC']:
                res_list['REC'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                               com_str.atoms[i].residue.chain, 'REC', 'R',
                                               com_str.atoms[i].residue.name,
                                               com_str.atoms[i].residue.insertion_code))
                res_list['COM'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                               com_str.atoms[i].residue.chain, 'COM', 'R',
                                               com_str.atoms[i].residue.name,
                                               com_str.atoms[i].residue.insertion_code))
                resindex += 1
                proc_res = res
        # save residue number in the lig list
        elif res != proc_res and resindex not in res_list['LIG']:
            res_list['LIG'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                           com_str.atoms[i].residue.chain, 'LIG', 'L',
                                           com_str.atoms[i].residue.name,
                                           com_str.atoms[i].residue.insertion_code))
            res_list['COM'].append(Residue(resindex, com_str.atoms[i].residue.number,
                                           com_str.atoms[i].residue.chain, 'COM', 'L',
                                           com_str.atoms[i].residue.name,
                                           com_str.atoms[i].residue.insertion_code))
            resindex += 1
            proc_res = res

    masks = {'REC': list2range(res_list['REC']), 'LIG': list2range(res_list['LIG'])}

    temp = []
    for m, value in masks.items():
        for e in value['num']:
            v = e[0] if isinstance(e, list) else e
            temp.append([v, m])
    temp.sort(key=lambda x: x[0])
    order_list = [c[1] for c in temp]

    return masks, res_list, order_list


def get_dist(coor1, coor2):
    return sqrt((coor2[0] - coor1[0]) ** 2 + (coor2[1] - coor1[1]) ** 2 + (coor2[2] - coor1[2]) ** 2)


def list2range(input_list):
    """
    Convert a list in list of ranges
    :return: list of ranges, string format of the list of ranges
    """

    def _add(temp):
        if len(temp) == 1:
            ranges_str.append(f"{temp[0]}")
            ranges.append([temp[0], temp[0]])
        else:
            ranges_str.append(f"{str(temp[0])}-{str(temp[-1])}")
            ranges.append([temp[0], temp[-1]])

    ranges = []
    ranges_str = []
    if not input_list:
        return ''
    temp = []
    previous = None

    input_list.sort()

    for x in input_list:
        if not previous:
            temp.append(x)
        elif x == previous + 1:
            temp.append(x)
        else:
            _add(temp)
            temp = [x]
        if x == input_list[-1]:
            _add(temp)
        previous = x
    return {'num': ranges, 'string': ranges_str}


def selector(selection: str):
    string_list = re.split(r"\s|;\s*", selection)
    dist = None
    # exclude = None
    res_selections = []
    if selection.startswith('within'):
        try:
            dist = float(string_list[1])
        except:
            GMXMMPBSA_ERROR(f'Invalid dist, we expected a float value but we get "{string_list[1]}"')
    else:
        # try to process residue selection
        for s in string_list:
            n = re.split(r":\s*|/\s*", s)
            if len(n) != 2 or n[0] not in ascii_letters:
                GMXMMPBSA_ERROR(f'We expected something like this: A/2-10,35,41 but we get {s} instead')
            chain = n[0]
            resl = n[1].split(',')
            for r in resl:
                rr = r.split('-')
                if len(rr) == 1:
                    ci = rr[0].split(':')
                    ri = [chain, int(ci[0]), ''] if len(ci) == 1 else [chain, int(ci[0]), ci[1]]
                    if ri in res_selections:
                        logging.warning('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
                                          '{}'.format(*ri))
                        continue
                    res_selections.append(ri)
                else:
                    try:
                        start = int(rr[0])
                        end = int(rr[1]) + 1
                    except:
                        GMXMMPBSA_ERROR(f'When residues range is defined, start and end most be integer but we get'
                                        f' {rr[0]} and {rr[1]}')
                    for cr in range(start, end):
                        if [chain, cr, ''] in res_selections:
                            logging.warning('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
                                              '{}'.format(chain, cr, ''))
                            continue
                        res_selections.append([chain, cr, ''])
    return dist, res_selections

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def remove(flag, mpi_size=0, fnpre='_GMXMMPBSA_'):
    """ Removes temporary files. Allows for different levels of cleanliness """

    # A list of all input files that we keep for the flag -use-mdins
    input_files = [fnpre + 'gb.mdin', fnpre + 'pb.mdin', fnpre + 'pb.mdin2',
                   fnpre + 'gb_decomp_com.mdin', fnpre + 'gb_decomp_rec.mdin',
                   fnpre + 'gb_decomp_lig.mdin', fnpre + 'pb_decomp_com.mdin',
                   fnpre + 'pb_decomp_rec.mdin', fnpre + 'pb_decomp_lig.mdin',
                   fnpre + 'cpptrajentropy.in',
                   fnpre + 'mutant_cpptrajentropy.in',
                   fnpre + 'gb_qmmm_com.mdin', fnpre + 'gb_qmmm_rec.mdin',
                   fnpre + 'gb_qmmm_lig.mdin']

    # All the extra files we keep for keep_files = 1
    keep_files_1 = [fnpre + 'ligand.mdcrd', fnpre + 'ligand.nc',
                    fnpre + 'mutant_ligand.mdcrd', fnpre + 'mutant_ligand.nc',
                    fnpre + 'complex.mdcrd', fnpre + 'mutant_complex.mdcrd',
                    fnpre + 'complex.nc', fnpre + 'mutant_complex.nc',
                    fnpre + 'receptor.mdcrd', fnpre + 'mutant_receptor.mdcrd',
                    fnpre + 'receptor.nc', fnpre + 'mutant_receptor.nc',
                    fnpre + 'dummycomplex.inpcrd', fnpre + 'complex.pdb',
                    fnpre + 'dummyreceptor.inpcrd', fnpre + 'receptor.pdb',
                    fnpre + 'dummyligand.inpcrd', fnpre + 'ligand.pdb',
                    fnpre + 'mutant_dummycomplex.inpcrd',
                    fnpre + 'mutant_dummyreceptor.inpcrd',
                    fnpre + 'mutant_dummyligand.inpcrd',
                    fnpre + 'mutant_complex.pdb', fnpre + 'mutant_receptor.pdb',
                    fnpre + 'mutant_ligand.pdb', fnpre + 'complex_nm.mdcrd',
                    fnpre + 'complex_nm.nc', fnpre + 'mutant_complex_nm.mdcrd',
                    fnpre + 'mutant_complex_nm.nc', fnpre + 'receptor_nm.mdcrd',
                    fnpre + 'receptor_nm.nc', fnpre + 'mutant_receptor_nm.nc',
                    fnpre + 'mutant_receptor_nm.mdcrd', fnpre + 'ligand_nm.nc',
                    fnpre + 'ligand_nm.mdcrd', fnpre + 'mutant_ligand_nm.nc',
                    fnpre + 'mutant_ligand_nm.mdcrd', fnpre + 'avgcomplex.pdb',
                    fnpre + 'mutant_avgcomplex.pdb', fnpre + 'ligand_entropy.out',
                    fnpre + 'complex_entropy.out', fnpre + 'receptor_entropy.out',
                    fnpre + 'cpptraj_entropy.out',
                    fnpre + 'mutant_cpptraj_entropy.out',
                    fnpre + 'mutant_complex_entropy.out',
                    fnpre + 'mutant_receptor_entropy.out',
                    fnpre + 'mutant_ligand_entropy.out',
                    fnpre + 'complex_gb.mdout', fnpre + 'mutant_complex_gb.mdout',
                    fnpre + 'receptor_gb.mdout', fnpre + 'mutant_receptor_gb.mdout',
                    fnpre + 'ligand_gb.mdout', fnpre + 'mutant_ligand_gb.mdout',
                    fnpre + 'complex_pb.mdout', fnpre + 'mutant_complex_pb.mdout',
                    fnpre + 'receptor_pb.mdout', fnpre + 'mutant_receptor_pb.mdout',
                    fnpre + 'ligand_pb.mdout', fnpre + 'mutant_ligand_pb.mdout',
                    fnpre + 'complex_rism.mdout',
                    fnpre + 'mutant_complex_rism.mdout',
                    fnpre + 'receptor_rism.mdout',
                    fnpre + 'mutant_receptor_rism.mdout',
                    fnpre + 'ligand_rism.mdout', fnpre + 'mutant_ligand_rism.mdout',
                    fnpre + 'complex_nm.out', fnpre + 'mutant_complex_nm.out',
                    fnpre + 'receptor_nm.out', fnpre + 'mutant_receptor_nm.out',
                    fnpre + 'ligand_nm.out', fnpre + 'mutant_ligand_nm.out',
                    fnpre + 'complex_gb_surf.dat', fnpre + 'receptor_gb_surf.dat',
                    fnpre + 'ligand_gb_surf.dat',
                    fnpre + 'mutant_complex_gb_surf.dat',
                    fnpre + 'mutant_receptor_gb_surf.dat',
                    fnpre + 'mutant_ligand_gb_surf.dat', fnpre + 'info']

    # Collect all of the temporary files (those starting with _MMPBSA_)
    allfiles = os.listdir(os.getcwd())
    tempfiles = []
    for fil in allfiles:
        if fil.startswith(fnpre): tempfiles.append(fil)

    if flag == -1:  # internal -- keep all mdin files
        for fil in tempfiles:
            if not fil in input_files: os.remove(fil)
    elif flag == 0:  # remove all temporary files
        for fil in tempfiles: os.remove(fil)
    elif flag == 1:  # keep keep mdcrds, mdouts, and other relevant output files
        for fil in tempfiles:
            if fil in keep_files_1: continue  # keep this file
            # Now we have to split out this file and analyze the base. If the
            # suffix is just a number (corresponding to a thread-specific output
            # file or trajectory), then we only want to remove it if in the base
            # name is not in keep_files_1
            base, ext = os.path.splitext(fil)
            if ext.strip('.').isdigit() and base in keep_files_1: continue
            # if we've made it this far, remove the file
            os.remove(fil)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def concatenate(file1, file2):
    """ Adds contents of file2 onto beginning of file1 """
    import os
    chunksize = 1048576  # Read the file in 1 MB chunks
    # Open the 2 files, the first in append mode
    with open(file2, 'r') as fl2:
        # Add a newline (make it OS-independent) to the first file if it doesn't
        # already end in one
        file1.write(os.linesep)

        str1 = fl2.read(chunksize)
        while str1:
            file1.write(str1)
            str1 = fl2.read(chunksize)

    file1.flush()
    # Now remove the merged file (file2)
    os.remove(file2)


class Unbuffered(object):
    """ Takes a stream handle and calls flush() on it after writing """

    def __init__(self, handle):
        self._handle = handle

    def write(self, data):
        self._handle.write(data)
        self._handle.flush()

    def __getattr__(self, attr):
        return getattr(self._handle, attr)
