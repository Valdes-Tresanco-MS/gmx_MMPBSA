"""
This is a module that contains functions responsible for mutating
the trajectory file for alanine scanning in gmx_MMPBSA. It must be
included with gmx_MMPBSA to insure proper functioning of alanine
scanning.
"""

# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdés-Tresanco and Mario E. Valdés-Tresanco   #
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
import contextlib
from pathlib import Path


def _getCoords(line, coordsperline, coordsize):
    """ Returns the coordinates of a line in a mdcrd file """
    holder = []
    location = 0
    for _ in range(coordsperline):
        with contextlib.suppress(Exception):
            tmp = float(line[location:location + coordsize])
            holder.append(tmp)
        location += coordsize
    return holder


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

    scalefactor = dist / actualdist  # determine scale factor

    coords[3] *= scalefactor  # scale original coordinates
    coords[4] *= scalefactor
    coords[5] *= scalefactor

    coords[3] += coords[0]  # move back to original place
    coords[4] += coords[1]
    coords[5] += coords[2]

    return coords  # return the coordinates


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
    raise MutateError(f'Unrecognized residue! Add {resname} to _getnumatms(resname) in alamdcrd.py and reinstall '
                      f'gmx_MMPBSA')


def _ressymbol(resname):
    """ Return 1-letter symbol of give amino acid """

    if resname == 'ALA':
        return 'A'
    elif resname == 'ARG':
        return 'R'
    elif resname == 'ASN':
        return 'N'
    elif resname in ['ASP', 'ASH']:
        return 'D'
    elif resname in ['CYS', 'CYX', 'CYM']:
        return 'C'
    elif resname in ['GLU', 'GLH']:
        return 'E'
    elif resname == 'GLN':
        return 'Q'
    elif resname == 'GLY':
        return 'G'
    elif resname in ['HIP', 'HID', 'HIE']:
        return 'H'
    elif resname == 'ILE':
        return 'I'
    elif resname == 'LEU':
        return 'L'
    elif resname in ['LYN', 'LYS']:
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


class MutantMdcrd(object):
    """ Class for an alanine-mutated amber trajectory file.
        ASCII only (no netcdf)
    """
    target_resname = 'ALA'

    def __init__(self, trajname, prm1, prm2):
        self.traj = trajname
        self.orig_prm = prm1
        self.new_prm = prm2
        self.mutres = self.FindMutantResidue()
        self.hasbox = bool(prm1.ptr('ifbox'))

    def __str__(self):
        mutation_indices = self.mutres if isinstance(self.mutres, list) else [self.mutres]
        return '; '.join(
            '%s%d%s' % (_ressymbol(self.orig_prm.parm_data['RESIDUE_LABEL'][index - 1]),
                        index, _ressymbol(self.target_resname))
            for index in mutation_indices
        )

    def FindMutantResidue(self):
        """ Finds which residue is the alanine mutant in a pair of prmtop files
        """
        origres = self.orig_prm.parm_data['RESIDUE_LABEL']
        newres = self.new_prm.parm_data['RESIDUE_LABEL']
        mutation_indices = []

        if len(origres) != len(newres):
            raise MutateError(('Mutant prmtop (%s) has a different number of ' +
                               'residues than the original (%s)!') %
                              (self.new_prm.prm_name, self.orig_prm.prm_name))

        for i in range(len(origres)):
            if origres[i] != newres[i]:
                if newres[i] != self.target_resname:
                    raise MutantResError(f'Mutant residue {i + 1} is {newres[i]} but must be {self.target_resname}!')
                mutation_indices.append(i + 1)
        if not mutation_indices:
            raise MutateError(f'Mutant prmtop ({self.new_prm.prm_name}) has the same sequence as the original!')
        return mutation_indices[0] if len(mutation_indices) == 1 else mutation_indices

    def MutateTraj(self, newname):
        """ Mutates a given mdcrd file based on 2 prmtops """

        mutation_indices = self.mutres if isinstance(self.mutres, list) else [self.mutres]
        original_atoms = self.orig_prm.ptr('natom')
        coordinates_per_frame = original_atoms * 3 + (3 if self.hasbox else 0)
        residue_labels = self.orig_prm.parm_data['RESIDUE_LABEL']
        residue_pointers = self.orig_prm.parm_data['RESIDUE_POINTER']

        if Path(self.traj).resolve() == Path(newname).resolve():
            raise MutateError('Original and mutated trajectory paths must differ.')
        output_started = False
        try:
            with open(self.traj) as mdcrd, open(newname, 'w') as new_mdcrd:
                output_started = True
                title = next(mdcrd, '').strip()
                new_mdcrd.write('%-80s' % f'{title} and mutated by gmx_MMPBSA for scanning')
                frames = 0
                for frame in self._iter_frames(mdcrd, coordinates_per_frame):
                    atom_coordinates = frame[:original_atoms * 3]
                    mutated_coordinates = []
                    for residue_index, residue_name in enumerate(residue_labels):
                        atom_start = (residue_pointers[residue_index] - 1) * 3
                        atom_end = ((residue_pointers[residue_index + 1] - 1) * 3
                                    if residue_index + 1 < len(residue_pointers) else original_atoms * 3)
                        residue_coordinates = atom_coordinates[atom_start:atom_end]
                        if residue_index + 1 in mutation_indices:
                            mutated_coordinates.extend(self._mutate(residue_name, residue_coordinates))
                        else:
                            mutated_coordinates.extend(residue_coordinates)
                    mutated_coordinates.extend(frame[original_atoms * 3:])
                    expected = self.new_prm.ptr('natom') * 3 + (3 if self.hasbox else 0)
                    if len(mutated_coordinates) != expected:
                        raise MutateError(
                            f'Mutated trajectory frame has {len(mutated_coordinates)} coordinates; expected {expected}.'
                        )
                    for coordinate_index, value in enumerate(mutated_coordinates):
                        if coordinate_index % 10 == 0:
                            new_mdcrd.write('\n')
                        new_mdcrd.write('%8.3f' % value)
                    frames += 1
                if not frames:
                    raise MutateError(f'Trajectory {self.traj} contains no complete frames.')
                new_mdcrd.write('\n')
        except BaseException:
            if output_started:
                Path(newname).unlink(missing_ok=True)
            raise

    def _iter_frames(self, mdcrd, coordinates_per_frame):
        """Yield complete frames while buffering at most one frame."""
        frame = []
        for line in mdcrd:
            for start in range(0, len(line.rstrip('\r\n')), 8):
                field = line[start:start + 8].strip()
                if not field:
                    continue
                try:
                    frame.append(float(field))
                except ValueError as exc:
                    raise MutateError(f'Invalid coordinate in trajectory {self.traj}: {field!r}') from exc
                if len(frame) == coordinates_per_frame:
                    yield frame
                    frame = []
        if frame:
            raise MutateError(f'Trajectory {self.traj} does not contain complete frames for the original topology.')

    def _mutate(self, resname, coords):
        list_one = 'ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU HID HIE HIP LEU LYN LYS MET PHE SER TRP TYR'

        list_two = 'ILE THR VAL'
        list_three = 'PRO'
        chdist = 1.09
        nhdist = 1.01
        coords_tosend = []
        new_coords = []
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
            raise MutateError('Mismatch in atom # in residue %s. (%d in alamdcrd.py and %d passed '
                              'in)' % (resname, _getnumatms(resname), len(coords)))

        if resname in list_one:
            new_coords.extend(coords[i] for i in range((7 + startindex) * 3))
            coords_tosend.extend(coords[(4 + startindex) * 3 + i] for i in range(3))
            coords_tosend.extend(coords[(7 + startindex) * 3 + i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)
            new_coords.extend(coords_received[i + 3] for i in range(3))
        elif resname in list_two:
            new_coords.extend(coords[i] for i in range((6 + startindex) * 3))
            coords_tosend.extend(coords[(4 + startindex) * 3 + i] for i in range(3))
            coords_tosend.extend(coords[(6 + startindex) * 3 + i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)
            new_coords.extend(coords_received[i + 3] for i in range(3))
            coords_tosend = []
            coords_tosend.extend(coords[(4 + startindex) * 3 + i] for i in range(3))
            coords_tosend.extend(coords[(10 + startindex) * 3 + i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)
            new_coords.extend(coords_received[3 + i] for i in range(3))
        elif resname in list_three:
            new_coords.extend(coords[i] for i in range((1 + startindex) * 3))
            coords_tosend = coords[startindex * 3:startindex * 3 + 6]
            coords_received = _scaledistance(coords_tosend, nhdist)
            new_coords.extend(coords_received[i + 3] for i in range(3))
            new_coords.extend(coords[(10 + startindex) * 3 + i] for i in range(6))
            new_coords.extend(coords[(7 + startindex) * 3 + i] for i in range(9))
            coords_tosend = []
            coords_tosend.extend(coords[(7 + startindex) * 3 + i] for i in range(3))
            coords_tosend.extend(coords[(4 + startindex) * 3 + i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)
            new_coords.extend(coords_received[i + 3] for i in range(3))
        else:
            raise MutateError(f"Residue {resname} not recognized! Can't mutate.")
        if cterm:
            new_coords.extend(coords[len(coords) - 9 + i] for i in range(9))
        else:
            new_coords.extend(coords[len(coords) - 6 + i] for i in range(6))
        return new_coords


class GlyMutantMdcrd(MutantMdcrd):
    target_resname = 'GLY'

    def _mutate(self, resname, coords):
        list_one = 'ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU HID HIE HIP LEU LYN LYS MET PHE SER TRP TYR ALA ILE THR VAL'
        list_two = 'PRO'

        chdist = 1.09
        nhdist = 1.01

        coords_tosend = []    # Coordinates to send to be scaled
        new_coords = []      # Mutated coordinates to return

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
                               'alamdcrd.py and %d passed in)') % (resname, _getnumatms(resname), len(coords)))

        if resname in list_one:
            new_coords.extend(coords[i] for i in range((4+startindex)*3))
            coords_tosend.extend(coords[(2+startindex)*3+i] for i in range(3))
            coords_tosend.extend(coords[(4+startindex)*3+i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)

            new_coords.extend(coords_received[i+3] for i in range(3))
        elif resname in list_two:
            new_coords.extend(coords[i] for i in range((1+startindex)*3))
            coords_tosend = coords[startindex*3:startindex*3 + 6]
            coords_received = _scaledistance(coords_tosend, nhdist)

            new_coords.extend(coords_received[i+3] for i in range(3))
            new_coords.extend(coords[(10+startindex)*3+i] for i in range(6))
            coords_tosend = []
            coords_tosend.extend(coords[(10+startindex)*3+i] for i in range(3))
            coords_tosend.extend(coords[(7+startindex)*3+i] for i in range(3))
            coords_received = _scaledistance(coords_tosend, chdist)

            new_coords.extend(coords_received[i+3] for i in range(3))
        else:
            raise MutateError("Residue %s not recognized! Can't mutate." % resname)

        if cterm:
            new_coords.extend(coords[len(coords)-9+i] for i in range(9))
        else:
            new_coords.extend(coords[len(coords)-6+i] for i in range(6))
        return new_coords
