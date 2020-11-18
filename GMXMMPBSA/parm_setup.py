"""
This module contains classes and functions that manipulate Amber prmtop
objects with respect to MM/PBSA. It will figure out where atoms belong
in the context of a complex/receptor/ligand system. It will also be able
to translate complex atom selections into either Amber masks or Amber
group input selection strings for the complex, receptor, and ligand.
It will also check through the systems and make sure that the prmtop
files are compatible, as well. Necessary for gmx_MMPBSA functioning.
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

from parmed.amber import LoadParm
from GMXMMPBSA.exceptions import PrmtopError, SelectionError


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Residue(object):
    """ Atom class in MM/PBSA, complex/receptor/ligand context """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, number, residue):
        """ Initializes an Atom instance. NOTE: Atom numbers here index
          beginning at 1 for both number (complex) and receptor_number
          or ligand_number! Normal python arrays index from 0
      """
        self.number = number
        self.residue = residue
        self.receptor_number = None
        self.ligand_number = None
        self.selected = False

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def select(self):
        """ Label this atom as selected """
        self.selected = True

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def deselect(self):
        """ Deselect this atom """
        self.selected = False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MMPBSA_System(object):
    """
    Sets up a system of a complex, receptor, and ligand.

    Get amber topology file from GROMACS

    The prmtops get initialized to prmtop objects immediately, so they must exist.
    Use it as:

    new_sys = MMPBSA_system(complex_prmtop, receptor_prmtop, ligand_prmtop)
    new_sys.Map()

    The MMPBSA_system is now set up, and you can get a Mask or Group input
    string. You can also check charges

    The following will return an amber mask string of the selected residues
    in the comma-delimited list format assuming you want a mask that corresponds
    to residue numbers in the complex prmtop (by definition, selection MUST
    correspond to residue numbers in the complex prmtop)

    com_mask, rec_mask, lig_mask = new_sys.Mask(selection, in_complex=False)

    The following will return a GROUP input string according to the same
    criteria as above for the masks

    com_group, rec_group, lig_group = new_sys.Group(selection, in_complex=True)

    The following will run sanity checks on the gmx_MMPBSA system to make sure stuff
    is consistent.

    new_sys.CheckConsistency()
    """

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, complex_prmtop, receptor_prmtop=None, ligand_prmtop=None):
        """
         Get amber topology from gromacs tpr file and Initializes the topology file classes for all of the prmtop files
         """

        self.complex_prmtop = LoadParm(complex_prmtop)
        self.stability = not receptor_prmtop and not ligand_prmtop
        if not self.stability:
            self.receptor_prmtop = LoadParm(receptor_prmtop)
            self.ligand_prmtop = LoadParm(ligand_prmtop)
        else:
            self.receptor_prmtop = None
            self.ligand_prmtop = None
        self.mapped = False
        self.ligstart = -1

        self._validate()

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Map(self, receptor_mask=None, ligand_mask=None):
        """ This determines the receptor and ligand masks based on residue
          sequence. It will try placing the ligand inside the receptor to
          match the complex residue sequence, starting from the end.
        """
        from parmed.amber import AmberMask

        # First, if a mask is provided, we can actually map it using that
        if (receptor_mask is not None and ligand_mask is not None and
                not self.stability):
            rmask = AmberMask(self.complex_prmtop, receptor_mask).Selection()
            lmask = AmberMask(self.complex_prmtop, ligand_mask).Selection()
            # Check that every atom is selected once and only once by comparing
            # the sums of the two mask selections
            if sum(rmask) + sum(lmask) != self.complex_prmtop.ptr('natom'):
                raise PrmtopError("provided receptor/ligand masks don't select " +
                                  "every atom in the complex topology!")
            if sum(rmask) != self.receptor_prmtop.ptr('natom'):
                raise PrmtopError('mismatch in receptor mask and receptor prmtop')
            if sum(lmask) != self.ligand_prmtop.ptr('natom'):
                raise PrmtopError('mismatch in ligand mask and ligand prmtop')

            # Now check that the masks don't select an atom twice
            for i in range(len(rmask)):
                if rmask[i] == 1 and lmask[i] == 1:
                    raise PrmtopError('Atom %d selected by both receptor and ' % i +
                                      'ligand masks')
                if lmask[i] == 0 and rmask[i] == 0:
                    raise PrmtopError('Atom %d is not selected by either ' % i +
                                      'receptor or ligand masks')

            # Now that we've verified everything is OK (we could do more, but this
            # should be good enough), let's actually create the res_list mapping
            lignum = recnum = 1
            self.res_list = []
            for i in range(self.complex_prmtop.ptr('nres')):
                start_ptr = self.complex_prmtop.parm_data['RESIDUE_POINTER'][i] - 1
                in_lig = lmask[start_ptr] == 1
                new_res = Residue(i + 1, self.complex_prmtop.parm_data[
                    'RESIDUE_LABEL'][i])
                if in_lig:
                    new_res.ligand_number = lignum
                    if (self.complex_prmtop.parm_data['RESIDUE_LABEL'][i] !=
                            self.ligand_prmtop.parm_data['RESIDUE_LABEL'][lignum - 1]):
                        raise PrmtopError('Residue mismatch while mapping. ' +
                                          'Incompatible topology files or bad mask definitions')
                    lignum += 1
                else:
                    new_res.receptor_number = recnum
                    if (self.complex_prmtop.parm_data['RESIDUE_LABEL'][i] !=
                            self.receptor_prmtop.parm_data['RESIDUE_LABEL'][recnum - 1]):
                        raise PrmtopError('Residue mismatch while mapping. ' +
                                          'Incompatible topology files or bad mask definitions')
                    recnum += 1
                self.res_list.append(new_res)
            # end for i in range(self.complex_prmtop.ptr('nres'))
            self.mapped = True
            return None  # done here

        # Begin find our own masks!
        if self.stability:
            # Stability mapping -- just copy the complex residue sequence
            self.res_list = \
                [Residue(i + 1, self.complex_prmtop.parm_data['RESIDUE_LABEL'][i])
                 for i in range(self.complex_prmtop.ptr('nres'))]
            self.mapped = True
            return None

        complex_residues = self.complex_prmtop.parm_data['RESIDUE_LABEL']
        receptor_residues = self.receptor_prmtop.parm_data['RESIDUE_LABEL']
        ligand_residues = self.ligand_prmtop.parm_data['RESIDUE_LABEL']

        # The way this loop is going to work is it's going to insert the ligand
        # in each location at the receptor, starting from the end and working to
        # the beginning. At each location, we will compare the residue sequence
        # of the receptor + ligand with the complex residue. If it matches, we've
        # the start of the ligand. If it doesn't match, we move on to the next.

        for i in range(len(receptor_residues)):
            idx = len(receptor_residues) - i - 1
            found_position = True

            for j in range(len(complex_residues)):
                com_res = complex_residues[j]
                if j <= idx:
                    rec_or_lig_res = receptor_residues[j]
                elif j > idx and j - idx <= len(ligand_residues):
                    rec_or_lig_res = ligand_residues[j - idx - 1]
                else:
                    rec_or_lig_res = receptor_residues[j - len(ligand_residues)]

                if com_res != rec_or_lig_res:
                    found_position = False
                    break

            # end for j...:
            if found_position:
                self.ligstart = idx + 2

        # end for i...:

        # Note that the above will NOT have tried to put the ligand first. Try to
        # do that here if we haven't found the ligand starting location yet

        if self.ligstart == -1:

            found_position = True
            for i in range(len(complex_residues)):
                com_res = complex_residues[i]

                if i < len(ligand_residues):
                    rec_or_lig_res = ligand_residues[i]
                else:  # in the receptor
                    rec_or_lig_res = receptor_residues[i - len(ligand_residues)]

                if com_res != rec_or_lig_res:
                    found_position = False
                    break
            # end for i ...:

            if found_position:
                self.ligstart = 1

            else:
                # We've only reached here because we CAN'T find it
                raise PrmtopError("Couldn't predict mask from topology files!\n" +
                                  "Your ligand residues must be sequential in your complex.\n" +
                                  "There are likely problems with your topology files if " +
                                  "this is not the case.")

        # end if self.ligstart == -1

        # Now that we know where the ligand is, we can map out all of the atoms
        # in both their complex and receptor or ligand locations.

        self.res_list = []  # start with an empty list

        for i in range(len(complex_residues)):
            new_res = Residue(i + 1, complex_residues[i])

            # Now set the number that residue is in the receptor or ligand. First
            # we see if it's before the ligand. Then we see if it's in the ligand,
            # then we see if it's after the ligand

            if i < self.ligstart - 1:
                new_res.receptor_number = i + 1
            elif i >= self.ligstart - 1 and \
                    i < self.ligstart + len(ligand_residues) - 1:
                new_res.ligand_number = i + 2 - self.ligstart
            else:
                new_res.receptor_number = i + 1 - len(ligand_residues)

            self.res_list.append(new_res)

        # End for i...:

        self.mapped = True

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Mask(self, selection, in_complex=True):
        """ This function turns a ,- or ;-delimited list of residues into an
          Amber mask for the complex, receptor, and ligand. If in_complex is
          true, then the masks returned will correspond to the residues inside
          the complex. If it's false, then it will correspond to just the
          selection in the receptor/ligand prmtops
      """

        if selection.lower() == 'all':
            selection = '1-%d' % self.complex_prmtop.ptr('nres')

        if self.stability:
            # M.L. added check for single residue
            if self.complex_prmtop.ptr('nres') == 1:
                return (':1', None, None)
            else:
                return self._stability_mask(selection)
        else:
            return self._binding_mask(selection, in_complex)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _stability_mask(self, selection):
        """ Internal call for Mask if stability is true """

        self._deselect_residues()
        self._select_residues(selection)

        # M.L. added a first residue condition, so do not put comma before it
        first_residue = False
        com_start = 0
        mask = ':'
        # M.L. keeps i as 0-index, use resid as 1-index for printout for cpptraj
        i = 0
        resid = 1
        while i < self.complex_prmtop.ptr('nres'):

            if not self.res_list[i].selected:
                i += 1
                resid += 1
                continue

            else:
                com_start = resid
                # M.L. modified loop
                # check that the next residue j is in the selection
                while i < self.complex_prmtop.ptr('nres') - 1:  # nres-1 b/c searching for j
                    j = i + 1
                    if self.res_list[j].selected:
                        i += 1
                        resid += 1
                    else:
                        break
                if first_residue == False:
                    mask += range_string(com_start, resid)
                    first_residue = True
                else:
                    mask += ',' + range_string(com_start, resid)
                com_start = -1
                i += 1
                resid += 1

        return (mask, None, None)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _binding_mask(self, selection, in_complex):
        """ Internal call if stability is False """

        self._deselect_residues()
        self._select_residues(selection)
        com_start = -1
        rec_start = -1
        lig_start = -1
        rec_end = -1
        lig_end = -1
        complex_mask = ''
        receptor_mask = ''
        ligand_mask = ''

        # Create 0ed arrays for receptor and ligand selections. If in_complex is
        # False, then we will need these arrays to generate the masks, or we run
        # the risk of artificially fragmented masks

        receptor_selection = [0 for i in range(self.receptor_prmtop.ptr('nres'))]
        ligand_selection = [0 for i in range(self.ligand_prmtop.ptr('nres'))]

        # Now every residue in self.res_list has been selected according to
        # that selection. Now we have to loop through all of the atoms until
        # we hit one that's selected
        i = 0
        while i < len(self.res_list):

            # If we have any intermediate receptor/ligand selections, flush them
            # here. i has been added twice too many times, so we need to subtract
            # off 2 of them, but add 1 back to adjust from 0 -> 1 indexing
            if rec_start != -1:
                if rec_start == rec_end:
                    receptor_mask += ',%d' % rec_start
                else:
                    receptor_mask += ',%d-%d' % (rec_start, i - 1)
                rec_start = -1
            if lig_start != -1:
                if lig_start == lig_end:
                    ligand_mask += ',%d' % lig_start
                else:
                    ligand_mask += ',%d-%d' % (lig_start, i - 1)
                lig_start = -1

            # Skip over unselected atoms
            if not self.res_list[i].selected:
                i += 1
                continue

            # Now this atom is selected. Add it to where it belongs
            com_start = i + 1
            if self.res_list[i].receptor_number != None:  # it's in receptor
                if in_complex:
                    rec_start = i + 1
                else:
                    receptor_selection[self.res_list[i].receptor_number - 1] = 1
                in_rec = True
            elif self.res_list[i].ligand_number != None:  # it's in ligand
                if in_complex:
                    lig_start = i + 1
                else:
                    ligand_selection[self.res_list[i].ligand_number - 1] = 1
                in_rec = False

            i += 1

            while i < len(self.res_list) and self.res_list[i].selected:
                # Now we have to make sure that we didn't leave the receptor
                # or ligand, or update where we are if we did.
                if in_rec:
                    if self.res_list[i].receptor_number != None:
                        # We're still in the receptor, go on
                        receptor_selection[self.res_list[i].receptor_number - 1] = 1
                        i += 1
                        continue
                    else:
                        # Define where we end and add to mask if in_complex
                        if in_complex:
                            rec_end = i
                            receptor_mask += ',' + range_string(rec_start, rec_end)
                            rec_start = -1
                        in_rec = False
                        # Define the start of our ligand
                        ligand_selection[self.res_list[i].ligand_number - 1] = 1
                        if in_complex: lig_start = i + 1
                        i += 1
                        continue
                else:
                    if self.res_list[i].ligand_number != None:
                        # We're still in the ligand, go on
                        ligand_selection[self.res_list[i].ligand_number - 1] = 1
                        i += 1
                        continue
                    else:
                        # Define where we end and add to mask if in_complex
                        if in_complex:
                            lig_end = i
                            ligand_mask += ',' + range_string(lig_start, lig_end)
                            lig_start = -1  # so we know this one is ended
                        in_rec = True
                        # Define the start of our receptor
                        receptor_selection[self.res_list[i].receptor_number - 1] = 1
                        if in_complex: rec_start = i + 1
                        i += 1
                        continue

                raise SelectionError("About to enter an infinite loop. Oh no!")

                # end if in_rec

            # end while i ... selected

            com_end = i
            complex_mask += ',' + range_string(com_start, com_end)
            com_start = -1

            i += 1

        # end while i < len(self.res_list)

        # If we had selected the last residue, then we have to add that last
        # chunk of mask back in.
        if com_start != -1:
            complex_mask += ',' + range_string(com_start, len(self.res_list))

        # Now we have a choice in how to construct our masks. If in_complex, then
        # just use them as they've been created and replace the beginning , with :
        # If not in_complex, then we need to use receptor_selection and
        # ligand_selection to determine the masks.

        if in_complex:
            if rec_start != -1:
                if in_complex:
                    rec_end = len(self.res_list)
                else:
                    rec_end = self.res_list[len(self.res_list) - 1].receptor_number
                receptor_mask += ',' + range_string(rec_start, rec_end)

            if lig_start != -1:
                if in_complex:
                    lig_end = len(self.res_list)
                else:
                    lig_end = self.res_list[len(self.res_list) - 1].ligand_number
                ligand_mask += ',' + range_string(lig_start, lig_end)
        else:
            # Get receptor mask
            i = 0
            start = -1
            while i < len(receptor_selection):
                if receptor_selection[i] == 1:
                    start = i + 1
                    i += 1
                    while i < len(receptor_selection) and receptor_selection[i] == 1:
                        i += 1
                    if start == i:
                        receptor_mask += ',%d' % start
                    else:
                        receptor_mask += ',%d-%d' % (start, i)
                    start = -1
                i += 1

            if start != -1:
                # The end was selected
                if start == len(receptor_selection):
                    receptor_mask += ',%d' % start
                else:
                    receptor_mask += ',%d-%d' % (start, len(receptor_selection))

            # Get ligand mask
            i = 0
            start = -1
            while i < len(ligand_selection):
                if ligand_selection[i] == 1:
                    start = i + 1
                    i += 1
                    while i < len(ligand_selection) and ligand_selection[i] == 1:
                        i += 1
                    if start == i:
                        ligand_mask += ',%d' % start
                    else:
                        ligand_mask += ',%d-%d' % (start, i)
                    start = -1
                i += 1

            if start != -1:
                # The end was selected
                if start == len(ligand_selection):
                    ligand_mask += ',%d' % start
                else:
                    ligand_mask += ',%d-%d' % (start, len(ligand_selection))

        # end if in_complex

        # Now, for any existing masks, replace leading , with :

        if complex_mask != '':
            complex_mask = ':' + complex_mask[1:]
        if ligand_mask != '':
            ligand_mask = ':' + ligand_mask[1:]
        if receptor_mask != '':
            receptor_mask = ':' + receptor_mask[1:]

        return (complex_mask, receptor_mask, ligand_mask)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def Group(self, selection, in_complex=True, com_group_type='RES',
              rec_group_type='RRES', lig_group_type='LRES'):
        """ This function turns a ,- or ;-delimited list of residues into an
          Amber group input specification for the complex, receptor, and ligand.
          If in_complex is true, then the group inputs returned will correspond
          to the residues inside the complex. If it's false, then it will
          correspond to just the selection in the receptor/ligand prmtops
        """

        if selection.lower() == 'all':
            selection = '1-%d' % self.complex_prmtop.ptr('nres')

        if self.stability:
            return self._stability_group(selection, com_group_type, rec_group_type)
        else:
            return self._binding_group(selection, in_complex, com_group_type,
                                       rec_group_type, lig_group_type)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _stability_group(self, selection, com_group_type, rec_group_type):
        """ Internal function for Group for stability calcs """

        self._deselect_residues()
        self._select_residues(selection)
        MAX_GROUP = 7
        com_start = -1

        i = 0
        num_groups = 1
        com_group_string = '%s' % com_group_type
        rec_group_string = '%s' % rec_group_type
        while i < self.complex_prmtop.ptr('nres'):
            if not self.res_list[i].selected:
                i += 1
                continue

            else:
                com_start = i + 1
                i += 1
                while i < self.complex_prmtop.ptr('nres') and \
                        self.res_list[i].selected:
                    i += 1

                if num_groups == MAX_GROUP:
                    com_group_string += '\n%s' % com_group_type
                    rec_group_string += '\n%s' % rec_group_type

                com_group_string += ' %d %d' % (com_start, i)
                rec_group_string += ' %d %d' % (com_start, i)
                num_groups += 1
                i += 1

        # end while i < len(....(nres))

        return (com_group_string, rec_group_string, '')

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _binding_group(self, selection, in_complex, com_group_type,
                       rec_group_type, lig_group_type):
        """ Internal function for Group for no stability """
        import re
        numre = re.compile(r'(\d+)')

        self._deselect_residues()
        self._select_residues(selection)
        MAX_GROUP = 7
        com_start = -1
        rec_start = -1
        lig_start = -1
        complex_group = []
        receptor_group = []
        ligand_group = []

        # Create 0ed arrays for receptor and ligand selections. If in_complex is
        # False, then we will need these arrays to generate the masks, or we run
        # the risk of artificially fragmented masks

        receptor_selection = [0 for i in range(self.receptor_prmtop.ptr('nres'))]
        ligand_selection = [0 for i in range(self.ligand_prmtop.ptr('nres'))]

        # Now every residue in self.res_list has been selected according to
        # that selection. Now we have to loop through all of the atoms until
        # we hit one that's selected
        i = 0
        while i < len(self.res_list):

            # If we have any intermediate receptor/ligand selections, flush them
            # here. i has been added twice too many times, so we need to subtract
            # off 2 of them, but add 1 back to adjust from 0 -> 1 indexing
            if rec_start != -1:
                receptor_group.append('%d %d' % (rec_start, i - 1))
                rec_start = -1
            if lig_start != -1:
                ligand_group.append('%d %d' % (lig_start, i - 1))
                lig_start = -1

            # Skip over unselected atoms
            if not self.res_list[i].selected:
                i += 1
                continue

            # Now this atom is selected. Add it to where it belongs
            com_start = i + 1
            if self.res_list[i].receptor_number != None:  # it's in receptor
                if in_complex:
                    rec_start = i + 1
                else:
                    receptor_selection[self.res_list[i].receptor_number - 1] = 1
                in_rec = True
            elif self.res_list[i].ligand_number != None:  # it's in ligand
                if in_complex:
                    lig_start = i + 1
                else:
                    ligand_selection[self.res_list[i].ligand_number - 1] = 1
                in_rec = False

            i += 1

            while i < len(self.res_list) and self.res_list[i].selected:
                # Now we have to make sure that we didn't leave the receptor
                # or ligand, or update where we are if we did.
                if in_rec:
                    if self.res_list[i].receptor_number != None:
                        # We're still in the receptor, go on
                        receptor_selection[self.res_list[i].receptor_number - 1] = 1
                        i += 1
                        continue
                    else:
                        # Define where we end and add to mask if in_complex
                        if in_complex:
                            rec_end = i
                            receptor_group.append('%d %d' % (rec_start, rec_end))
                            rec_start = -1
                        in_rec = False
                        # Define the start of our ligand
                        ligand_selection[self.res_list[i].ligand_number - 1] = 1
                        if in_complex: lig_start = i + 1
                        i += 1
                        continue
                else:
                    if self.res_list[i].ligand_number != None:
                        # We're still in the ligand, go on
                        ligand_selection[self.res_list[i].ligand_number - 1] = 1
                        i += 1
                        continue
                    else:
                        # Define where we end and add to mask if in_complex
                        if in_complex:
                            lig_end = i
                            ligand_group.append('%d %d' % (lig_start, lig_end))
                            lig_start = -1  # so we know this one is ended
                        in_rec = True
                        # Define the start of our receptor
                        receptor_selection[self.res_list[i].receptor_number - 1] = 1
                        if in_complex: rec_start = i + 1
                        i += 1
                        continue

                raise SelectionError("About to enter an infinite loop. Oh no!")

                # end if in_rec

            # end while i ... selected

            com_end = i
            complex_group.append('%d %d' % (com_start, com_end))
            com_start = -1

            i += 1

        # end while i < len(self.res_list)

        # If we had selected the last residue, then we have to add that last
        # chunk of mask back in.
        if com_start != -1:
            complex_group.append('%d %d' % (com_start, len(self.res_list)))

        # Now we have a choice in how to construct our masks. If in_complex, then
        # just use them as they've been created and replace the beginning , with :
        # If not in_complex, then we need to use receptor_selection and
        # ligand_selection to determine the masks.

        if in_complex:
            if rec_start != -1:
                receptor_group.append('%d %d' % (rec_start, len(self.res_list)))

            if lig_start != -1:
                ligand_group.append('%d %d' % (lig_start, len(self.res_list)))
        else:
            # Get receptor group
            i = 0
            start = -1
            while i < len(receptor_selection):
                if receptor_selection[i] == 1:
                    start = i + 1
                    i += 1
                    while i < len(receptor_selection) and receptor_selection[i] == 1:
                        i += 1
                    receptor_group.append('%d %d' % (start, i))
                    start = -1
                i += 1

            if start != -1:
                # The end was selected
                receptor_group.append('%d %d' % (start, len(receptor_selection)))

            # Get ligand group
            i = 0
            start = -1
            while i < len(ligand_selection):
                if ligand_selection[i] == 1:
                    start = i + 1
                    i += 1
                    while i < len(ligand_selection) and ligand_selection[i] == 1:
                        i += 1
                    ligand_group.append('%d %d' % (start, i))
                    start = -1
                i += 1

            if start != -1:
                # The end was selected
                ligand_group.append('%d %d' % (start, len(ligand_selection)))

        # end if in_complex

        # Now, for any existing groups, format them for group input, making
        # sure to limit the number of groups to MAX_GROUP

        com_grp_str = com_group_type + ' '
        rec_grp_str = rec_group_type + ' '
        lig_grp_str = lig_group_type + ' '

        # First do the complex
        num_grps = 0
        for grp in complex_group:
            if num_grps % MAX_GROUP == 0 and num_grps != 0:
                com_grp_str += '\n%s ' % com_group_type
            com_grp_str += grp + ' '
            num_grps += 1

        # Next do the receptor
        num_grps = 0
        for grp in receptor_group:
            if num_grps % MAX_GROUP == 0 and num_grps != 0:
                rec_grp_str += '\n%s ' % rec_group_type
            rec_grp_str += grp + ' '
            num_grps += 1

        # Next do the ligand
        num_grps = 0
        for grp in ligand_group:
            if num_grps % MAX_GROUP == 0 and num_grps != 0:
                lig_grp_str += '\n%s ' % lig_group_type
            lig_grp_str += grp + ' '
            num_grps += 1

        # Groups are invalid if they are empty. So if either rec_grp_str or
        # lig_grp_str have no numbers, add a 1 1 to it
        if len(numre.findall(rec_grp_str)) == 0:
            rec_grp_str += "1 1"
        if len(numre.findall(lig_grp_str)) == 0:
            lig_grp_str += "1 1"
        return (com_grp_str, rec_grp_str, lig_grp_str)

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _select_residues(self, selection):
        """ This routine takes a given selection, comma- or ;-delimited,
          and selects the appropriate residues from res_list.
      """

        if not self.mapped:
            raise PrmtopError('The MMPBSA_system must be mapped before you can '
                              'request masks or groups!')

        selection = selection.replace(';', ',')  # move to one kind of delimiter
        selection_fields = selection.split(',')

        # Get rid of any blank fields

        for field in selection_fields:

            field = field.strip()

            # Skip over blank fields
            if len(field) == 0:
                continue

            # If '-' is in the field, then it's a range
            if '-' in field:
                field_items = field.split('-')

                # Make sure the fields are properly formed
                if len(field_items) != 2:
                    raise SelectionError('Selection ranges must be "# - #"!')
                if len(field_items[0].strip()) == 0:
                    raise SelectionError('Selection ranges must be "# - #"!')
                if len(field_items[1].strip()) == 0:
                    raise SelectionError('Selection ranges must be "# - #"!')

                # Make sure the fields are composed of integers
                try:
                    res1 = int(field_items[0])
                    res2 = int(field_items[1])
                except:
                    raise SelectionError('Invalid selection! Integers expected.')

                # Now make sure that they're within the legal range and not stupid
                if res2 < res1 or res1 <= 0 or \
                        res2 > self.complex_prmtop.ptr('nres'):
                    raise SelectionError('Invalid selection. Illegal residue range')

                # Now go through and select them
                for i in range(res1 - 1, res2):
                    self.res_list[i].select()

            else:
                # Here the field is just a single residue
                try:
                    res1 = int(field)
                except:
                    raise SelectionError('Invalid selection! Integers expected.')

                # select that one residue

                self.res_list[res1 - 1].select()

            # End '-' in field

        # End for field...:

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _deselect_residues(self):
        """ Deselect all of the residues in the system """

        if not self.mapped:
            raise PrmtopError('Cannot deselect residues before mapping!')

        for i in range(len(self.res_list)):
            self.res_list[i].deselect()

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _validate(self):
        """ Validates the passed prmtops. Does a preliminary compatibility
          check
      """
        if self.stability: return None  # done checking here for stability calcs

        if self.complex_prmtop.ptr('natom') != self.receptor_prmtop.ptr('natom') \
                + self.ligand_prmtop.ptr('natom'):
            raise PrmtopError('Complex natom != receptor natom + ligand natom')

        if self.complex_prmtop.ptr('nres') != self.receptor_prmtop.ptr('nres') + \
                self.ligand_prmtop.ptr('nres'):
            raise PrmtopError('Complex nres != receptor nres + ligand nres')

        if (self.complex_prmtop.chamber or self.receptor_prmtop.chamber or
                self.ligand_prmtop.chamber):
            if (not self.complex_prmtop.chamber or not self.receptor_prmtop.chamber
                    or not self.ligand_prmtop.chamber):
                raise PrmtopError('Prmtops must be either ALL chamber ' +
                                  'prmtops or NONE')

    # -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def CheckConsistency(self):
        """ Checks to make sure that the charge definitions are consistent across
          each prmtop file, but only if we know which atoms are which atoms in
          each prmtop
      """

        if self.stability: return None

        # Check that our radii sets are consistent
        try:
            com_radii = self.complex_prmtop.parm_data['RADIUS_SET'][0]
            rec_radii = self.receptor_prmtop.parm_data['RADIUS_SET'][0]
            lig_radii = self.ligand_prmtop.parm_data['RADIUS_SET'][0]
            if com_radii != rec_radii or com_radii != lig_radii:
                raise PrmtopError('Topology files have inconsistent RADIUS_SETs')
        except KeyError:
            raise PrmtopError('Topology files have no Implicit radii! You can ' +
                              'add implicit radii using xparmed.py or parmed.py')
        if not self.mapped:
            raise PrmtopError('Cannot check charges prior to mapping!')

        TINY = 0.0001

        # Loop through every residue and compare its charge value in the complex
        # to the corresponding atom in the ligand or receptor
        for i in range(self.complex_prmtop.ptr('natom')):
            resnum = self.complex_prmtop.atoms[i].residue.idx
            atom_in_res = (i + 1 -
                           self.complex_prmtop.parm_data['RESIDUE_POINTER'][resnum])
            com_chg = self.complex_prmtop.parm_data['CHARGE'][i]
            if self.res_list[resnum].ligand_number is not None:
                otherparm = self.ligand_prmtop
                resnum = self.res_list[resnum].ligand_number - 1
            else:
                otherparm = self.receptor_prmtop
                resnum = self.res_list[resnum].receptor_number - 1

            atnum = otherparm.parm_data['RESIDUE_POINTER'][resnum] - 1 + atom_in_res
            lig_or_rec_chg = otherparm.parm_data['CHARGE'][atnum]
            if abs(lig_or_rec_chg - com_chg) > TINY:
                raise PrmtopError('Inconsistent charge definition for atom %d!' %
                                  (i + 1))


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def range_string(num1, num2):
    """ This function returns a single number if they're equal, or a range if
       they're not
   """

    if num1 == num2:
        return str(num1)
    else:
        return '%d-%d' % (num1, num2)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
