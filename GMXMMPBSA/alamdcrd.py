"""
This is a module that contains functions responsible for mutating
the trajectory file for alanine scanning in gmx_MMPBSA. It must be
included with gmx_MMPBSA to insure proper functioning of alanine
scanning.
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

from GMXMMPBSA.exceptions import MutateError, MutantResError

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _getCoords(line, coordsperline, coordsize):
    """ Returns the coordinates of a line in a mdcrd file """
    holder = []
    location = 0
    for i in range(coordsperline):
        try:
            tmp = float(line[location:location+coordsize])
            holder.append(tmp)
        except:
            pass
        location += coordsize

    return holder

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _scaledistance(coords, dist):
    """ Scales the distance between 2 3-D cartesian coordinates to the specified
        distance
    """
    from math import sqrt

    if len(coords) != 6:
        raise MutateError('_scaledistance requires x,y,z coords for 2 atoms')

    coords[3] -= coords[0]  # set first 3 coordinates as origin
    coords[4] -= coords[1]
    coords[5] -= coords[2]

    actualdist = sqrt(coords[3]*coords[3] + coords[4]*coords[4] +
                      coords[5]*coords[5])

    scalefactor = dist / actualdist # determine scale factor

    coords[3] *= scalefactor  # scale original coordinates
    coords[4] *= scalefactor
    coords[5] *= scalefactor

    coords[3] += coords[0] # move back to original place
    coords[4] += coords[1]
    coords[5] += coords[2]

    return coords  # return the coordinates

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _getnumatms(resname):
    """ Returns the number of atoms in a given Amino acid residue """
    if resname in 'GLY':
        return 7
    if resname in ('ALA', 'CYM', 'CYX'):
        return 10
    if resname in ('CYS', 'SER'):
        return 11
    if resname in ('ASP'):
        return 12
    if resname in ('ASN', 'PRO', 'THR'):
        return 14
    if resname in ('GLU'):
        return 15
    if resname in ('GLH', 'VAL'):
        return 16
    if resname in ('GLN', 'HID', 'HIE', 'MET'):
        return 17
    if resname in ('HIP'):
        return 18
    if resname in ('ILE', 'LEU'):
        return 19
    if resname in ('PHE'):
        return 20
    if resname in ('LYN', 'TYR'):
        return 21
    if resname in ('LYS'):
        return 22
    if resname in ('ARG', 'TRP'):
        return 24

    raise MutateError(('Unrecognized residue! Add %s to _getnumatms(resname) ' +
                       'in alamdcrd.py and reinstall gmx_MMPBSA' % resname))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _ressymbol(resname):
    """ Return 1-letter symbol of give amino acid """

    if resname == 'ALA':
        return 'A'
    elif resname == 'ARG':
        return 'R'
    elif resname == 'ASN':
        return 'N'
    elif resname in ['ASP','ASH']:
        return 'D'
    elif resname in ['CYS','CYX','CYM']:
        return 'C'
    elif resname in ['GLU','GLH']:
        return 'E'
    elif resname == 'GLN':
        return 'Q'
    elif resname == 'GLY':
        return 'G'
    elif resname in ['HIP','HID','HIE']:
        return 'H'
    elif resname == 'ILE':
        return 'I'
    elif resname == 'LEU':
        return 'L'
    elif resname in ['LYN','LYS']:
        return 'K'
    elif resname == 'MET':
        return 'M'
    elif resname == 'PHE':
        return 'F'
    elif resname == 'PRO':
        return 'P'
    elif resname == 'SER':
        return 'S'
    elif resname == 'THR':
        return 'T'
    elif resname == 'TRP':
        return 'W'
    elif resname == 'TYR':
        return 'Y'
    elif resname == 'VAL':
        return 'V'
    else:
        return resname

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MutantMdcrd(object):
    """ Class for an alanine-mutated amber trajectory file.
        ASCII only (no netcdf)
    """

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, trajname, prm1, prm2):
        self.traj = trajname
        self.orig_prm = prm1
        self.new_prm = prm2
        self.mutres = self.FindMutantResidue()
        self.hasbox = bool(prm1.ptr('ifbox'))

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __str__(self):
        return '%s%d%s' % (_ressymbol(
            self.orig_prm.parm_data['RESIDUE_LABEL'][self.mutres-1]),
                           self.mutres, _ressymbol('ALA'))

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def FindMutantResidue(self):
        """ Finds which residue is the alanine mutant in a pair of prmtop files
        """
        origres = self.orig_prm.parm_data['RESIDUE_LABEL']
        newres = self.new_prm.parm_data['RESIDUE_LABEL']
        diffs = 0
        mutres = -1

        if len(origres) != len(newres):
            raise MutateError(('Mutant prmtop (%s) has a different number of ' +
                               'residues than the original (%s)!') %
                              (self.new_prm.prm_name, self.orig_prm.prm_name))

        for i in range(len(origres)):
            if origres[i] != newres[i]:
                diffs += 1
                if newres[i] != 'ALA':
                    raise MutantResError('Mutant residue %s is %s but must be ALA!' %
                                         (i+1, newres[i]))
                mutres = i + 1

        if diffs == 0:
            raise MutateError('Your mutant prmtop (%s) has the same sequence '
                              'as the original!' % self.new_prm.prm_name)
        elif diffs > 1:
            raise MutateError('Your mutant prmtop (%s) can only have one mutation!'
                              % self.new_prm.prm_name)

        return mutres

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def MutateTraj(self, newname):
        """ Mutates a given mdcrd file based on 2 prmtops """

        mutres = self.mutres

        orig_resname = self.orig_prm.parm_data['RESIDUE_LABEL'][mutres-1]
        resstart = self.orig_prm.parm_data['RESIDUE_POINTER'][mutres-1]
        try:
            nextresstart = self.orig_prm.parm_data['RESIDUE_POINTER'][mutres]
        except IndexError:
            nextresstart = self.orig_prm.ptr('natom')
        number_atoms_mut = self.new_prm.ptr('natom')
        if self.orig_prm.ptr('ifbox'):
            number_atoms_mut += 1

        coordsperline = 10 # number of coordinates in each line
        coordsize = 8      # how large coordinates are in characters
        counter = 0
        coords_done = 0
        coords_tomutate = []
        temp_holder = []
        new_coords = []

        # location is 0 before modified coordinates, 1 during modified
        # coordinates, and 2 after modified coordinates
        location = 0

        if orig_resname == 'GLY':
            raise MutateError('You are trying to mutate GLY to ALA! ' +
                              'Not currently supported.')

        new_mdcrd = open(newname, 'w')
        mdcrd = open(self.traj, 'r')

        for line in mdcrd:
            counter += 1
            # First line is always a comment
            if counter == 1:
                new_mdcrd.write('%-80s' % (line.strip() +
                                           ' and mutated by gmx_MMPBSA for alanine scanning'))
                continue

            if location == 0 and coords_done <= max(resstart * 3 - 4, 0) and \
                    coords_done + coordsperline > max(resstart * 3 - 4, 0):
                location = 1
                words = _getCoords(line, coordsperline, coordsize)
                for i in range(coordsperline):
                    if coords_done <= resstart * 3 - 4:
                        if coords_done % coordsperline == 0:
                            new_mdcrd.write('\n')
                        new_mdcrd.write('%8.3f' % words[i])
                        coords_done += 1
                    else:
                        coords_tomutate.append(words[i])
                continue

            elif location == 1:
                words = _getCoords(line, coordsperline, coordsize)
                if coordsperline + len(coords_tomutate) >= \
                        3 * (nextresstart - resstart):
                    location = 2
                    for i in range(coordsperline):
                        if len(coords_tomutate) < 3 * (nextresstart - resstart):
                            coords_tomutate.append(words[i])
                        else:
                            temp_holder.append(words[i])

                    new_coords = self._mutate(orig_resname, coords_tomutate)
                    for i in range(len(new_coords)):
                        if coords_done % coordsperline == 0:
                            new_mdcrd.write('\n')
                        new_mdcrd.write('%8.3f' % new_coords[i])
                        coords_done += 1
                    if len(temp_holder) != 0:
                        for i in range(len(temp_holder)):
                            if coords_done % coordsperline == 0:
                                new_mdcrd.write('\n')
                            new_mdcrd.write('%8.3f' % temp_holder[i])
                            coords_done += 1
                    coords_tomutate = []
                    temp_holder = []
                else:
                    for i in range(coordsperline):
                        coords_tomutate.append(words[i])

                continue

            elif location == 2:
                words = _getCoords(line, coordsperline, coordsize)
                for i in range(len(words)):
                    if coords_done % coordsperline == 0:
                        if not self.hasbox or coords_done < number_atoms_mut * 3 - 3:
                            new_mdcrd.write('\n')
                    new_mdcrd.write('%8.3f' % words[i])
                    coords_done += 1
                    if self.hasbox and coords_done == number_atoms_mut * 3 - 3:
                        new_mdcrd.write('\n')

                if coords_done == number_atoms_mut * 3:
                    coords_done = 0
                    location = 0
                continue
            else:
                if coords_done % coordsperline == 0:
                    new_mdcrd.write('\n')
                new_mdcrd.write(line[:len(line)-1])
                coords_done = coords_done + coordsperline
                continue

        new_mdcrd.write('\n')
        mdcrd.close()
        new_mdcrd.close()

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _mutate(self, resname, coords):
        list_one = 'ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU HID HIE HIP ' + \
                   'LEU LYN LYS MET PHE SER TRP TYR'
        list_two = 'ILE THR VAL'
        list_three = 'PRO'

        chdist = 1.09
        nhdist = 1.01

        coords_tosend = []
        new_coords  = []
        coords_received = []
        cterm = False

        if _getnumatms(resname) * 3 == len(coords):
            startindex = 0
            cterm = False
        elif (_getnumatms(resname) + 2) * 3 == len(coords):
            startindex = 2
            cterm = False
        elif (_getnumatms(resname) + 1) * 3 == len(coords):
            startindex = 0
            cterm = True
        else:
            raise MutateError(('Mismatch in atom # in residue %s. (%d in ' +
                               'alamdcrd.py and %d passed in)') % (resname,
                                                                   _getnumatms(resname), len(coords)))

        if resname in list_one:
            for i in range((7+startindex)*3):
                new_coords.append(coords[i])
            for i in range(3):
                coords_tosend.append(coords[(4+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(7+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])

        elif resname in list_two:
            for i in range((6+startindex)*3):
                new_coords.append(coords[i])

            for i in range(3):
                coords_tosend.append(coords[(4+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(6+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])

            coords_tosend = []
            coords_received = []
            for i in range(3):
                coords_tosend.append(coords[(4+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(10+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[3+i])

        elif resname in list_three:
            for i in range((1+startindex)*3):
                new_coords.append(coords[i])

            coords_tosend = coords[startindex*3:startindex*3 + 6]
            coords_received = _scaledistance(coords_tosend, nhdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])
            for i in range(6):
                new_coords.append(coords[(10+startindex)*3+i])
            for i in range(9):
                new_coords.append(coords[(7+startindex)*3+i])

            coords_tosend = []
            coords_received = []
            for i in range(3):
                coords_tosend.append(coords[(7+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(4+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])

        else:
            raise MutateError("Residue %s not recognized! Can't mutate." % resname)

        if cterm:
            for i in range(9):
                new_coords.append(coords[len(coords)-9+i])
        else:
            for i in range(6):
                new_coords.append(coords[len(coords)-6+i])

        return new_coords

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class GlyMutantMdcrd(MutantMdcrd):

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __str__(self):
        return '%s%d%s' % (_ressymbol(
            self.orig_prm.parm_data['RESIDUE_LABEL'][self.mutres-1]),
                           self.mutres, _ressymbol('GLY'))

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def _mutate(self, resname, coords):
        list_one = 'ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU HID HIE HIP ' + \
                   'LEU LYN LYS MET PHE SER TRP TYR ALA ILE THR VAL'
        list_two = 'PRO'

        chdist = 1.09
        nhdist = 1.01

        coords_tosend = []    # Coordinates to send to be scaled
        new_coords  = []      # Mutated coordinates to return
        coords_received = []  # Scaled coordinates received
        cterm = False

        if _getnumatms(resname) * 3 == len(coords):
            startindex = 0
            cterm = False
        elif (_getnumatms(resname) + 2) * 3 == len(coords):
            startindex = 2
            cterm = False
        elif (_getnumatms(resname) + 1) * 3 == len(coords):
            startindex = 0
            cterm = True
        else:
            raise MutateError(('Mismatch in atom # in residue %s. (%d in ' +
                               'alamdcrd.py and %d passed in)') % (resname,
                                                                   _getnumatms(resname), len(coords)))

        if resname in list_one:
            for i in range((4+startindex)*3):
                new_coords.append(coords[i])
            for i in range(3):
                coords_tosend.append(coords[(2+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(4+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])

        elif resname in list_two:
            for i in range((1+startindex)*3):
                new_coords.append(coords[i])

            coords_tosend = coords[startindex*3:startindex*3 + 6]
            coords_received = _scaledistance(coords_tosend, nhdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])
            for i in range(6):
                new_coords.append(coords[(10+startindex)*3+i])

            coords_tosend = []
            coords_received = []
            for i in range(3):
                coords_tosend.append(coords[(10+startindex)*3+i])
            for i in range(3):
                coords_tosend.append(coords[(7+startindex)*3+i])
            coords_received = _scaledistance(coords_tosend, chdist)

            for i in range(3):
                new_coords.append(coords_received[i+3])

        else:
            raise MutateError("Residue %s not recognized! Can't mutate." % resname)

        if cterm:
            for i in range(9):
                new_coords.append(coords[len(coords)-9+i])
        else:
            for i in range(6):
                new_coords.append(coords[len(coords)-6+i])

        return new_coords

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def FindMutantResidue(self):
        """ Finds which residue is the alanine mutant in a pair of prmtop files
        """
        origres = self.orig_prm.parm_data['RESIDUE_LABEL']
        newres = self.new_prm.parm_data['RESIDUE_LABEL']
        diffs = 0
        mutres = -1

        if len(origres) != len(newres):
            raise MutateError(('Mutant prmtop (%s) has a different number of ' +
                               'residues than the original (%s)!') %
                              (self.new_prm.prm_name, self.orig_prm.prm_name))

        for i in range(len(origres)):
            if origres[i] != newres[i]:
                diffs += 1
                if newres[i] != 'GLY':
                    raise MutantResError('Mutant residue %s is %s but must be GLY!' %
                                         (i+1, newres[i]))
                mutres = i + 1

        if diffs == 0:
            raise MutateError(('Your mutant prmtop (%s) has the same sequence ' +
                               'as the original!') % (self.new_prm.prm_name,
                                                      self.orig_prm.prm_name))
        elif diffs > 1:
            raise MutateError('Your mutant prmtop (%s) can only have one mutation!'
                              % self.new_prm.prm_name)

        return mutres
