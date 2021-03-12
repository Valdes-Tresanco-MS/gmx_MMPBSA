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
import shutil
from pathlib import Path
import json
import logging
from string import ascii_letters
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, GMXMMPBSA_WARNING
from math import sqrt

def get_dist(coor1, coor2):
    return sqrt((coor2[0] - coor1[0]) ** 2 + (coor2[1] - coor1[1]) ** 2 + (coor2[2] - coor1[2]) ** 2)


def selector(selection:str):
    string_list = selection.split()
    dist = None
    exclude = None
    res_selections = []
    if selection.startswith('within'):
        try:
            dist = float(string_list[1])
        except:
            GMXMMPBSA_ERROR(f'Invalid dist, we expected a float value but we get "{string_list[1]}"')
        if len(string_list) == 4:
            if string_list[2] != 'not':
                GMXMMPBSA_ERROR(f'We expected "not" but we get {string_list[2]} instead')
            if str(string_list[3]).lower() not in ['receptor', 'ligand', 'rec', 'lig']:
                GMXMMPBSA_ERROR(f'We expected one of this values: "receptor", "rec", "ligand" or "lig" but we get'
                      f' {(string_list[3]).lower()} instead')
            exclude = 'REC' if str(string_list[3]).lower() in ['receptor', 'rec'] else 'LIG'
    else:
        # try to process residue selection
        for s in string_list:
            n = s.split('/')
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
                        GMXMMPBSA_WARNING('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
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
                            GMXMMPBSA_WARNING('Found duplicated residue in selection: CHAIN:{} RES_NUM:{} ICODE: '
                                              '{}'.format(chain, cr, ''))
                            continue
                        res_selections.append([chain, cr, ''])
    return dist, exclude, res_selections


def checkff(overwrite):
    """

    :param overwrite:
    :param sel_ff: folder/leaprc
    Examples:
        gmxMMPBSA/leaprc.GLYCAM_06h-1
        oldff/leaprc.ff99SB
        leaprc.protein.ff14SB

    :return:
    """
    amberhome = os.getenv('AMBERHOME')
    if not amberhome:
        GMXMMPBSA_ERROR('Could not found Amber. Please make sure you have sourced %s/amber.sh (if you are using sh/ksh/'
              'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
              (amberhome, amberhome))
        return
    amberhome = Path(amberhome)

    logging.info('Checking if supported force fields exists in Amber data...')
    if overwrite:
        logging.info('The overwrite_data option is enabled. Overwriting the gmxMMPBSA data...')
    data_path = Path(__file__).parent.joinpath('data')
    info_file = data_path.joinpath('info.dat')
    with info_file.open('r') as read_file:
        data = json.load(read_file)

    leap_dat = amberhome.joinpath('dat/leap/')

    for ff in data['forcefield']:
        for p in data['forcefield'][ff]:
            gmxf = leap_dat.joinpath(p, 'gmxMMPBSA')
            if not gmxf.exists():
                gmxf.mkdir()
            if type(data['forcefield'][ff][p]) == list:
                for l in data['forcefield'][ff][p]:
                    cf = data_path.joinpath(l)
                    if overwrite:
                        shutil.copy(cf, gmxf)
                    else:
                        if cf.exists() and not gmxf.joinpath(l).exists():
                            shutil.copy(cf, gmxf)
            else:
                cf = data_path.joinpath(data['forcefield'][ff][p])
                if overwrite:
                    shutil.copy(cf, gmxf)
                else:
                    if not cf.exists():
                        shutil.copy(cf, gmxf)
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


