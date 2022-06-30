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
import platform
import re
import shutil
import sys
from pathlib import Path
import json
import logging
from string import ascii_letters

import pandas as pd

from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
from math import sqrt
import parmed
import numpy as np
from typing import Union


class EnergyVector(np.ndarray):
    def __new__(cls, values=None, com_std=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        if isinstance(values, int):
            obj = np.zeros((values,)).view(cls)
        elif isinstance(values, (list, tuple, np.ndarray)):
            obj = np.array(values).view(cls)
        else:
            obj = np.array([]).view(cls)
        obj.com_std = com_std
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.com_std = getattr(obj, 'com_std', None)

    # This fix the pickle problem. Taken from
    # https://stackoverflow.com/questions/26598109/preserve-custom-attributes-when-pickling-subclass-of-numpy-array
    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(EnergyVector, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (self.__dict__,)
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self.__dict__.update(state[-1])  # Update the internal dict from state
        # Call the parent's __setstate__ with the other tuple elements.
        super(EnergyVector, self).__setstate__(state[:-1])

    def stdev(self):
        return self.com_std or self.std()

    def sem(self):
        return float(self.std() / sqrt(len(self)))

    def semp(self):
        return float(self.stdev() / sqrt(len(self)))

    def append(self, values):
        return EnergyVector(np.append(self, values))

    def avg(self):
        return np.average(self)

    def corr_add(self, other):
        selfstd = self.com_std or float(self.std())
        comp_std = None
        if isinstance(other, EnergyVector):
            otherstd = other.com_std or float(other.std())
            comp_std = get_corrstd(selfstd, otherstd)
        return EnergyVector(np.add(self, other), comp_std)

    def corr_sub(self, other):
        self_std = self.com_std or float(np.asarray(self).std())
        comp_std = None
        if isinstance(other, EnergyVector):
            other_std = other.com_std or float(np.asarray(other).std())
            comp_std = get_corrstd(self_std, other_std)
        return EnergyVector(np.subtract(self, other), comp_std)

    def __add__(self, other):
        selfstd = self.com_std or float(self.std())
        comp_std = None
        if isinstance(other, EnergyVector):
            otherstd = other.com_std or float(other.std())
            comp_std = get_std(selfstd, otherstd)
        return EnergyVector(np.add(self, other), comp_std)

    def __sub__(self, other):
        self_std = self.com_std or float(np.asarray(self).std())
        comp_std = None
        if isinstance(other, EnergyVector):
            other_std = other.com_std or float(np.asarray(other).std())
            comp_std = get_std(self_std, other_std)
        return EnergyVector(np.subtract(self, other), comp_std)

    def __eq__(self, other):
        return np.all(np.equal(self, other))

    def __lt__(self, other):
        return np.all(np.less(self, other))

    def __le__(self, other):
        return np.all(np.less_equal(self, other))

    def __gt__(self, other):
        return np.all(np.greater(self, other))

    def __ge__(self, other):
        return np.all(np.greater_equal(self, other))

    def abs_gt(self, val):
        """ If any element's absolute value is greater than a # """
        return np.any(np.greater(self, val))


def get_std(val1, val2):
    return sqrt(val1 ** 2 + val2 ** 2)


def get_corrstd(val1, val2):
    return sqrt(val1 ** 2 + val2 ** 2 - 2 * val1 * val2)


def calc_sum(vector1, vector2, mut=False) -> (float, float):
    """
    Calculate the mean and std of the two vector/numbers sum
    Args:
        vector1: EnergyVector or float
        vector2: EnergyVector or float
        mut: If mutant, the SD is the standard deviation of the array

    Returns:
        dmean: Mean of the sum
        dstd: Standard deviation
    """
    if isinstance(vector2, EnergyVector) and isinstance(vector1, EnergyVector):
        if mut:
            d = vector2 + vector1
            dmean = float(d.mean())
            dstd = float(d.std())
        else:
            dmean = float(vector2.mean() + vector1.mean())
            dstd = float(get_std(vector2.std(), vector1.std()))
    elif isinstance(vector2, EnergyVector) and isinstance(vector1, (int, float)):
        dmean = float(vector2.mean() + vector1)
        dstd = vector2.std()
    elif isinstance(vector2, (int, float)) and isinstance(vector1, EnergyVector):
        dmean = float(vector2 + vector1.mean())
        dstd = vector1.std()
    else:
        dmean = float(vector2 + vector1)
        dstd = 0.0
    return dmean, dstd


def create_input_args(args: list):
    if not args or 'all' in args:
        return 'general', 'gb', 'pb', 'ala', 'nmode', 'decomp', 'rism'
    elif 'gb' not in args and 'pb' not in args and 'rism' not in args and 'nmode' not in 'args':
        GMXMMPBSA_ERROR('You did not specify any type of calculation!')
    elif 'gb' not in args and 'pb' not in args and 'decomp' in args:
        logging.warning('&decomp calculation is only compatible with &gb and &pb calculations. Will be ignored!')
        args.remove('decomp')
        return ['general'] + args
    else:
        return ['general'] + args


def mask2list(com_str, rec_mask, lig_mask):
    rm_list = rec_mask.strip(":").split(',')
    lm_list = lig_mask.strip(':').split(',')
    res_list = []

    for r in rm_list:
        if '-' in r:
            s, e = r.split('-')
            res_list.extend([i, 'R'] for i in range(int(s) - 1, int(e)))
        else:
            res_list.append([int(r) - 1, 'R'])
    for l in lm_list:
        if '-' in l:
            s, e = l.split('-')
            res_list.extend([i, 'L'] for i in range(int(s) - 1, int(e)))
        else:
            res_list.append([int(l) - 1, 'L'])
    res_list = sorted(res_list, key=lambda x: x[0])
    comstr = parmed.load_file(com_str)
    resl = []
    rec_index = 1
    lig_index = 1
    for res, rl in zip(comstr.residues, res_list):
        if rl[1] == 'R':
            resl.append(Residue(rl[0] + 1, res.number, res.chain, rl[1], rec_index, res.name, res.insertion_code))
            rec_index += 1
        else:
            resl.append(Residue(rl[0] + 1, res.number, res.chain, rl[1], lig_index, res.name, res.insertion_code))
            lig_index += 1
    return resl


def log_subprocess_output(process):
    while True:
        if output := process.stdout.readline().decode():
            if output.startswith(' ->  frame'):
                continue
            logging.debug(output.strip('\n'))
        else:
            break


class Residue(object):
    def __init__(self, index, number, chain, mol_id, id_index, name, icode=''):
        self.index = int(index)
        self.number = number
        self.chain = chain
        self.mol_id = mol_id
        self.id_index = id_index
        self.name = name
        self.icode = icode
        self.mutant_label = None
        self.string = f"{mol_id}:{chain}:{name}:{number}:{icode}" if icode else f"{mol_id}:{chain}:{name}:{number}"
        self.mutant_string = None

    def __repr__(self):
        text = f"{type(self).__name__}(index: {self.index}, {self.mol_id}:{self.chain}:{self.name}:{self.number}"
        if self.icode:
            text += f":{self.icode}"
        text += ')'
        return text

    def __str__(self):
        return f"{self.index}"

    def __add__(self, other):
        if isinstance(other, Residue):
            return int(self.index + other.index)
        return int(self.index + other)

    def __sub__(self, other):
        if isinstance(other, Residue):
            return int(self.index - other.index)
        return int(self.index - other)

    def __int__(self):
        return self.index

    def is_mutant(self):
        return bool(self.mutant_label)

    def is_receptor(self):
        return self.mol_id == 'R'

    def is_ligand(self):
        return self.mol_id == 'L'

    def issame(self, other):
        pass

    def set_mut(self, mut):
        self.mutant_label = f'{self.chain}/{self.number}{f":{self.icode}" if self.icode else ""} - {self.name}x{mut}'
        self.mutant_string = (f"{self.mol_id}:{self.chain}:{mut}:{self.number}:{self.icode}" if self.icode
                              else f"{self.mol_id}:{self.chain}:{mut}:{self.number}")


def multiindex2dict(p: Union[pd.MultiIndex, pd.Index, dict]) -> dict:
    """
    Converts a pandas Multiindex to a nested dict
    :parm p: As this is a recursive function, initially p is a pd.MultiIndex, but after the first iteration it takes
    the internal_dict value, so it becomes to a dictionary
    """
    internal_dict = {}
    end = False
    for x in p:
        # Since multi-indexes have a descending hierarchical structure, it is convenient to start from the last
        # element of each tuple. That is, we start by generating the lower level to the upper one. See the example
        if isinstance(p, pd.MultiIndex):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used. Only for 2 levels
            # pd.MultiIndex
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = [x[-1]]
            else:
                internal_dict[t].append(x[-1])
        elif isinstance(x, tuple):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = {x[-1]: p[x]}
            else:
                internal_dict[t][x[-1]] = p[x]
    if end:
        return internal_dict
    return multiindex2dict(internal_dict)


def flatten(dictionary, parent_key: list = False):
    """
    Turn a nested dictionary into a flattened dictionary
    :param dictionary: The dictionary to flatten
    :param parent_key:The accumulated list of keys
    :return: A flattened dictionary with key as tuples of nested keys
    """

    items = []
    for key, value in dictionary.items():
        new_key = parent_key + [key] if parent_key else [key]
        if isinstance(value, dict):
            items.extend(flatten(value, new_key).items())
        else:
            items.append((tuple(new_key), value))
    return dict(items)


def emapping(d):
    internal_dict = {}
    for k, v in d.items():
        if isinstance(v, dict):
            if v.values():
                if isinstance(list(v.values())[0], (dict, pd.DataFrame)):
                    internal_dict[k] = emapping(v)
                else:
                    internal_dict[k] = list(v.keys())
        elif isinstance(v, pd.DataFrame):
            internal_dict[k] = multiindex2dict(v.columns)
        else:
            internal_dict[k] = v
    return internal_dict


def get_indexes(com_ndx, rec_ndx=None, rec_group=1, lig_ndx=None, lig_group=1):
    ndx_files = {'COM': com_ndx, 'REC': rec_ndx, 'LIG': lig_ndx}
    ndx = {'COM': {'header': [], 'index': []}, 'REC': {'header': [], 'index': []}, 'LIG': {'header': [], 'index': []}}
    for n, f in ndx_files.items():
        if f is None:
            continue
        with open(f) as indexf:
            indexes = []
            for line in indexf:
                if line.startswith('['):
                    header = line.strip('\n[] ')
                    ndx[n]['header'].append(header)
                    if indexes:
                        ndx[n]['index'].append(indexes)
                        indexes = []
                else:
                    indexes.extend(map(int, line.split()))
            ndx[n]['index'].append(indexes)

    comind = ndx['COM']['header'].index('GMXMMPBSA_REC_GMXMMPBSA_LIG')
    crecind = ndx['COM']['header'].index('GMXMMPBSA_REC')
    cligind = ndx['COM']['header'].index('GMXMMPBSA_LIG')
    com_indexes = {'COM': ndx['COM']['index'][comind], 'REC': ndx['COM']['index'][crecind],
                   'LIG': ndx['COM']['index'][cligind]}
    rec_indexes = ndx['REC']['index'][rec_group] if rec_ndx else {}
    lig_indexes = ndx['LIG']['index'][lig_group] if lig_ndx else {}
    return {'COM': com_indexes, 'REC': rec_indexes, 'LIG': lig_indexes}


def _get_dup_args(args):
    flags = []
    flags_values = []
    cv = []
    for o in args:
        if o.startswith('-'):
            flags.append(o)
            flags_values.append(cv)
            cv = []
        else:
            cv.append(o)
    flags_values.append(cv)

    opt_duplicates = []
    args_duplicates = []
    for x in flags:
        if flags.count(x) > 1 and x not in opt_duplicates:
            opt_duplicates.append(x)
    for x in flags_values:
        if flags_values.count(x) > 1 and x and ' '.join(x) not in args_duplicates:
            args_duplicates.append(' '.join(x))

    if opt_duplicates or args_duplicates:
        text_dup = 'Duplicates were found on the command line:\n'
        if opt_duplicates:
            text_dup += f"Options: {', '.join(opt_duplicates)}\n"
        if args_duplicates:
            text_dup += f"Arguments: {', '.join(args_duplicates)}\n"
        GMXMMPBSA_ERROR(text_dup)


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


def check_str(structure, ref=False, skip=False):
    if isinstance(structure, str):
        refstr = parmed.read_PDB(structure)
    else:
        refstr = structure

    previous = 0
    ind = 1
    res_dict = {}
    duplicates = []
    for res in refstr.residues:
        if 'LP' in res.name:
            GMXMMPBSA_ERROR('The LP pseudo-atom is not supported. Please remove them following this instructions: '
                            'https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand_LPH_atoms_CHARMMff/')
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
        duplicates.extend(
            f'{chain}:{resl[c][0]}:{resl[c][1]}:{resl[c][2]}'
            for c, x in enumerate(res_id_list)
            if res_id_list.count(x) > 1
        )

    if ref:
        if duplicates:
            GMXMMPBSA_ERROR(f'The reference structure used is inconsistent. The following residues are duplicates:\n'
                            f' {", ".join(duplicates)}')
    elif skip:
        if duplicates:
            return refstr
    elif duplicates:
        logging.warning(f'The complex structure used is inconsistent. The following residues are duplicates:\n'
                        f' {", ".join(duplicates)}')
    return refstr


def res2map(indexes, com_file):
    """
    :param com_str:
    :return:
    """
    res_list = []
    rec_list = []
    lig_list = []
    com_len = len(indexes['COM']['COM'])
    if isinstance(com_file, parmed.Structure):
        com_str = com_file
    else:
        com_str = parmed.load_file(com_file)

    resindex = 1
    rec_index = 1
    lig_index = 1
    proc_res = None
    for i in range(com_len):
        res = [com_str.atoms[i].residue.chain, com_str.atoms[i].residue.number, com_str.atoms[i].residue.name,
               com_str.atoms[i].residue.insertion_code]
        # We check who owns the residue corresponding to this atom
        if indexes['COM']['COM'][i] in indexes['COM']['REC']:
            # save residue number in the rec list
            if res != proc_res and resindex not in res_list:
                rec_list.append(resindex)
                res_list.append(Residue(resindex, com_str.atoms[i].residue.number,
                                        com_str.atoms[i].residue.chain, 'R', rec_index,
                                        com_str.atoms[i].residue.name,
                                        com_str.atoms[i].residue.insertion_code))
                resindex += 1
                rec_index += 1
                proc_res = res
        # save residue number in the lig list
        elif res != proc_res and resindex not in res_list:
            lig_list.append(resindex)
            res_list.append(Residue(resindex, com_str.atoms[i].residue.number,
                                    com_str.atoms[i].residue.chain, 'L', lig_index,
                                    com_str.atoms[i].residue.name,
                                    com_str.atoms[i].residue.insertion_code))
            resindex += 1
            lig_index += 1
            proc_res = res

    masks = {'REC': list2range(rec_list), 'LIG': list2range(lig_list)}

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

    ilist = sorted(input_list, key=lambda x: x.index if isinstance(x, Residue) else x)

    for x in ilist:
        if not previous:
            temp.append(x)
        elif x == previous + 1:
            temp.append(x)
        else:
            _add(temp)
            temp = [x]
        if x == ilist[-1]:
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
                GMXMMPBSA_ERROR(f'We expected something like this: A/2-10,35,41 B/104 but we get {s} instead')
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


def remove(flag, fnpre='_GMXMMPBSA_'):
    """ Removes temporary files. Allows for different levels of cleanliness """
    # Collect all of the temporary files (those starting with _GMXMMPBSA_)
    allfiles = os.listdir(os.getcwd())

    other_files = ['COM.prmtop', 'REC.prmtop', 'LIG.prmtop', 'MUT_COM.prmtop', 'MUT_REC.prmtop', 'MUT_LIG.prmtop',
                   'leap.log']
    if flag == -1:
        result_files = ['FINAL_RESULTS_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.dat']
        for fil in allfiles:
            if (
                    fil.startswith(fnpre) or fil.startswith(f"#{fnpre}") or
                    bool(re.match('#?(COM|REC|LIG|MUT_COM|MUT_REC|MUT_LIG)_traj_(\d)\.xtc', fil)) or
                    fil == 'COMPACT_MMXSA_RESULTS.mmxsa' or
                    fil in other_files or
                    fil in result_files):
                os.remove(fil)

    elif flag == 0:  # remove all temporary files
        for fil in allfiles:

            if fil.startswith(fnpre) or bool(re.match('#?(COM|REC|LIG|MUT_COM|MUT_REC|MUT_LIG)_traj_(\d)\.xtc',
                                                      fil)) or fil in other_files:
                os.remove(fil)


def find_progs(INPUT, mpi_size=0):
    """ Find the necessary programs based in the user INPUT """
    # List all of the used programs with the conditions that they are needed
    logging.info('Checking external programs...')
    used_progs = {'cpptraj': True,
                  'tleap': True,
                  'parmchk2': True,
                  'sander': True,
                  'sander.APBS': INPUT['sander_apbs'] == 1,
                  'mmpbsa_py_nabnmode': INPUT['nmoderun'],
                  # 'rism3d.snglpnt': INPUT['rismrun']
                  'elsize': INPUT['alpb']

                  }
    gro_exe = {
        'gmx5': [
            # look for any available gromacs executable
            'gmx', 'gmx_mpi', 'gmx_d', 'gmx_mpi_d'],
        'gmx4': [
            # look for gromacs 4.x
            'make_ndx', 'trjconv', 'editconf']}

    # The returned dictionary:
    my_progs = {}

    for prog, needed in used_progs.items():
        my_progs[prog] = shutil.which(prog, path=os.environ['PATH'])
        if needed:
            if not my_progs[prog]:
                GMXMMPBSA_ERROR('Could not find necessary program [%s]' % prog)
            logging.info('%s found! Using %s' % (prog, str(my_progs[prog])))

    search_parth = INPUT['gmx_path'] or os.environ['PATH']
    g5 = False
    for gv, g_exes in gro_exe.items():
        if gv == 'gmx5':
            for prog in g_exes:
                if exe := shutil.which(prog, path=search_parth):
                    logging.info('Using GROMACS version > 5.x.x!')
                    my_progs['make_ndx'] = [exe, 'make_ndx']
                    my_progs['editconf'] = [exe, 'editconf']
                    my_progs['trjconv'] = [exe, 'trjconv']
                    g5 = True
                    if prog in ['gmx_mpi', 'gmx_mpi_d'] and mpi_size > 1:
                        GMXMMPBSA_ERROR('gmx_mpi and gmx_mpi_d are not supported when running gmx_MMPBSA in parallel '
                                        'due to incompatibility between the mpi libraries used to compile GROMACS and '
                                        'mpi4py respectively. You can still use gmx_mpi or gmx_mpi_d to run gmx_MMPBSA '
                                        'serial. For parallel calculations use gmx instead')
                    logging.info('%s found! Using %s' % (prog, exe))
                    break
            if g5:
                break
        else:
            logging.info('Using GROMACS version 4.x.x!')
            for prog in g_exes:
                if exe := shutil.which(prog, path=search_parth):
                    my_progs[prog] = [exe]
                    logging.info('%s found! Using %s' % (prog, str(my_progs[prog])))

    if 'make_ndx' not in my_progs or 'editconf' not in my_progs or 'trjconv' not in my_progs:
        GMXMMPBSA_ERROR('Could not find necessary program [ GROMACS ]')
    logging.info('Checking external programs...Done.\n')
    return my_progs


def get_sys_info():
    """
    Print relevant system info for debugging proposes in the gmx_MMPBSA.log file
    """
    logging.debug(f"WDIR          : {Path('.').absolute().as_posix()}")
    logging.debug(f"AMBERHOME     : {os.environ['AMBERHOME'] if 'AMBERHOME' in os.environ else ''}")
    logging.debug(f"PYTHON EXE    : {shutil.which('python')}")
    logging.debug("PYTHON VERSION: " + ''.join(sys.version.split('\n')))
    logging.debug(f"MPI           : {shutil.which('mpirun')}")
    logging.debug(f"ParmEd        : {parmed.__version__}")
    logging.debug(f"OS PLATFORM   : {platform.platform()}")
    logging.debug(f"OS SYSTEM     : {platform.system()}")
    logging.debug(f"OS VERSION    : {platform.version()}")
    logging.debug(f"OS PROCESSOR  : {platform.processor()}\n")


def get_warnings():
    info = {'warning': 0, 'error': 0}
    with open('gmx_MMPBSA.log') as logfile:
        for line in logfile:
            if line.startswith('[ERROR  ]'):
                info['error'] += 1
            elif line.startswith('[WARNING]'):
                info['warning'] += 1
    return info


class Unbuffered(object):
    """ Takes a stream handle and calls flush() on it after writing """

    def __init__(self, handle):
        self._handle = handle

    def write(self, data):
        self._handle.write(data)
        self._handle.flush()

    def __getattr__(self, attr):
        return getattr(self._handle, attr)
