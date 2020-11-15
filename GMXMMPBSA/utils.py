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

def checkff(sel_ff):
    print('Checking if supported force fields exists...')
    amberhome = os.getenv('AMBERHOME')
    if not amberhome:
        print('Could not found Amber. Please make sure you have sourced %s/amber.sh (if you are using sh/ksh/'
                          'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                          (amberhome, amberhome))
        return
    amberhome = Path(amberhome)
    data_info = {'forcefield': {}}
    current_file_path = Path(__file__).parent
    with open(Path(__file__).parent.joinpath('data/info.dat')) as f:
        for line in f.readlines():
            if line.startswith('FF:'):
                ff = line.split()[1]
                folder = line.split()[2]
                files = line.split()[3:]
                if ff not in data_info['forcefield']:
                    data_info['forcefield'][ff] = {folder: files}
                else:
                    data_info['forcefield'][ff][folder] = files
    try:
        leap_dat = amberhome.joinpath('dat/leap/')
        if not leap_dat.joinpath('cmd/' + data_info['forcefield'][sel_ff]['cmd'][0]).exists():
            shutil.copy(current_file_path.joinpath('data/' + data_info['forcefield'][sel_ff]['cmd'][0]),
                        leap_dat.joinpath('cmd'))
            shutil.copy(current_file_path.joinpath('data/' + data_info['forcefield'][sel_ff]['prep'][0]),
                        leap_dat.joinpath('prep'))
            shutil.copy(current_file_path.joinpath('data/' + data_info['forcefield'][sel_ff]['parm'][0]),
                        leap_dat.joinpath('parm'))
            for f in data_info['forcefield'][sel_ff]['lib']:
                shutil.copy(current_file_path.joinpath('data/' + f), leap_dat.joinpath('lib'))
            print(f'Coping {sel_ff} to Amber data... Done')
    except:
        pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def remove(flag, mpi_size=0, fnpre='_MMPBSA_'):
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


class PDB:
    def __init__(self):

        self.allAtoms = []
        self.info = []
        self.element = None

    def parse(self, filename):
        ofile = open(filename).readlines()
        model = 0
        for line in ofile:
            if 'MODEL' in line:
                model += 1
                if model > 1:
                    raise ValueError('Only one Model is allowed')
            elif 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
                atm = {}
                atm['id'] = line[0:6].strip()
                atm["number"] = int(line[6:11].strip())
                atm["name"] = line[12:16].strip()
                atm["alternate_location"] = line[16].strip()
                atm["resname"] = line[17:21].strip()
                atm["chain"] = line[21].strip()
                atm["resnum"] = int(line[22:26].strip())
                atm["i_code"] = line[26].strip()
                atm["x"] = float(line[30:38].strip())
                atm["y"] = float(line[38:46].strip())
                atm["z"] = float(line[46:54].strip())
                atm["occupancy"] = float(line[54:60].strip())
                atm["b_factor"] = float(line[60:66].strip())
                atm['segment_id'] = line[72:76].strip()
                atm["element"] = line[76:78]
                atm["charge"] = line[78:80]
                self.allAtoms.append(atm)
            elif 'TER' in line[0:6]:
                atm = {}
                atm['id'] = 'TER'
                self.allAtoms.append(atm)
            elif 'ENDML' in line[0:6]:
                continue
            else:
                self.info.append(line)

    def getOutLine(self, atm, ter=False):
        # Create the PDB line.
        if ter:
            atomrec = '%-6s\n' % (atm['id'])
        else:
            atomrec = '%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s\n' % (
                atm['id'], atm["number"], atm["name"], atm["alternate_location"], atm["resname"], atm["chain"],
                atm["resnum"], atm["i_code"], atm["x"], atm["y"], atm["z"], atm["occupancy"], atm["b_factor"],
                atm['segment_id'], atm["element"], atm["charge"]
            )
        return atomrec
