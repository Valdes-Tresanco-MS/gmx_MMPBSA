"""
 This is a module that contains functions generally useful for the
 MMPBSA.py script. A full list of functions/subroutines is shown below.
 It must be included to insure proper functioning of MMPBSA.py

                           GPL LICENSE INFO

Copyright (C) 2009 - 2012 Jason Swails, Bill Miller III, and Dwight McGee

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

################################################################################
#  List of functions and a brief description of their purpose
#
#  remove: Removes temporary work files in this directory. It has a number of
#          different levels to remove only a small number of files.
#  concatenate: combines 2 files into a single, common file
################################################################################

# TODO get rid of this file altogether and move these functions into the main
# app class

import os

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
   keep_files_1=[fnpre + 'ligand.mdcrd', fnpre + 'ligand.nc',
                 fnpre + 'mutant_ligand.mdcrd', fnpre + 'mutant_ligand.nc',
                 fnpre + 'complex.mdcrd', fnpre + 'mutant_complex.mdcrd',
                 fnpre + 'complex.nc', fnpre + 'mutant_complex.nc',
                 fnpre + 'receptor.mdcrd', fnpre + 'mutant_receptor.mdcrd',
                 fnpre + 'receptor.nc', fnpre + 'mutant_receptor.nc',
                 fnpre + 'dummycomplex.inpcrd', fnpre + 'complex.pdb',
                 fnpre + 'dummyreceptor.inpcrd',fnpre + 'receptor.pdb',
                 fnpre + 'dummyligand.inpcrd',fnpre + 'ligand.pdb',
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
                 fnpre + 'receptor_gb.mdout',fnpre + 'mutant_receptor_gb.mdout',
                 fnpre + 'ligand_gb.mdout', fnpre + 'mutant_ligand_gb.mdout',
                 fnpre + 'complex_pb.mdout', fnpre + 'mutant_complex_pb.mdout',
                 fnpre + 'receptor_pb.mdout',fnpre + 'mutant_receptor_pb.mdout',
                 fnpre + 'ligand_pb.mdout', fnpre + 'mutant_ligand_pb.mdout',
                 fnpre + 'complex_rism.mdout',
                 fnpre + 'mutant_complex_rism.mdout',
                 fnpre + 'receptor_rism.mdout',
                 fnpre + 'mutant_receptor_rism.mdout',
                 fnpre + 'ligand_rism.mdout',fnpre + 'mutant_ligand_rism.mdout',
                 fnpre + 'complex_nm.out', fnpre + 'mutant_complex_nm.out',
                 fnpre + 'receptor_nm.out',fnpre + 'mutant_receptor_nm.out',
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

   if flag == -1: # internal -- keep all mdin files
      for fil in tempfiles:
         if not fil in input_files: os.remove(fil)
   elif flag == 0: # remove all temporary files
      for fil in tempfiles: os.remove(fil)
   elif flag == 1: # keep keep mdcrds, mdouts, and other relevant output files
      for fil in tempfiles:
         if fil in keep_files_1: continue # keep this file
         # Now we have to split out this file and analyze the base. If the
         # suffix is just a number (corresponding to a thread-specific output
         # file or trajectory), then we only want to remove it if in the base
         # name is not in keep_files_1
         base, ext = os.path.splitext(fil)
         if ext.strip('.').isdigit() and base in keep_files_1: continue
         # if we've made it this far, remove the file
         os.remove(fil)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def concatenate(file1, file2):
   """ Adds contents of file2 onto beginning of file1 """
   import os
   chunksize = 1048576 # Read the file in 1 MB chunks
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
