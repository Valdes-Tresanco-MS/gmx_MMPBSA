"""
This module contains all the classes and code to collect data and calculate
statistics from the output files of various calculation types. Each calculation
type needs its own class.

All data is stored in a special class derived from the list.
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
import logging
from copy import deepcopy
from math import sqrt
from GMXMMPBSA.exceptions import (OutputError, LengthError, DecompError)
from GMXMMPBSA.utils import EnergyVector, get_std
from types import SimpleNamespace
import numpy as np
import sys
import re
from csv import writer

idecompString = ['idecomp = 0: No decomposition analysis',
                 'idecomp = 1: Per-residue decomp adding 1-4 interactions to Internal.',
                 'idecomp = 2: Per-residue decomp adding 1-4 interactions to EEL and VDW.',
                 'idecomp = 3: Pairwise decomp adding 1-4 interactions to Internal.',
                 'idecomp = 4: Pairwise decomp adding 1-4 interactions to EEL and VDW.']
sep = '-------------------------------------------------------------------------------'

data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'], 'DIHED': ['GGAS', 'TOTAL'],
                  'VDWAALS': ['GGAS', 'TOTAL'], 'EEL': ['GGAS', 'TOTAL'], '1-4 VDW': ['GGAS', 'TOTAL'],
                  '1-4 EEL': ['GGAS', 'TOTAL'],
                  # charmm
                  'UB': ['GGAS', 'TOTAL'], 'IMP': ['GGAS', 'TOTAL'], 'CMAP': ['GGAS', 'TOTAL'],
                  # non lineal PB
                  'EEL+EPB': ['TOTAL'],
                  # PB
                  'EPB': ['GSOLV', 'TOTAL'], 'ENPOLAR': ['GSOLV', 'TOTAL'], 'EDISPER': ['GSOLV', 'TOTAL'],
                  # GB
                  'EGB': ['GSOLV', 'TOTAL'], 'ESURF': ['GSOLV', 'TOTAL'],
                  # QM/GB
                  'ESCF': ['GGAS', 'TOTAL'],
                  # RISM
                  'POLAR SOLV': ['GSOLV', 'TOTAL'], 'APOLAR SOLV': ['GSOLV', 'TOTAL'],
                  'ERISM': ['GSOLV', 'TOTAL'],
                  # NMODE and QH
                  'TRANSLATIONAL': ['TOTAL'], 'ROTATIONAL': ['TOTAL'], 'VIBRATIONAL': ['TOTAL']
                  }


class AmberOutput(dict):
    """
    Base Amber output class. It takes a basename as a file name and parses
    through all of the thread-specific output files (assumed to have the suffix
    .# where # spans from 0 to num_files - 1
    """
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EPOL': 1}

    def __init__(self, mol: str, INPUT, chamber=False, **kwargs):
        super(AmberOutput, self).__init__(**kwargs)
        self.numframes = None
        self.mol = mol
        self.INPUT = INPUT
        self.chamber = chamber
        self.basename = None
        self.num_files = None
        self.frame_idx = 0
        self.extraframe_idx = 0
        self.is_read = False
        self.apbs = INPUT['sander_apbs']

        # This variable is used to get if the nmode calculation hasn't at least one frame
        self.no_nmode_convergence = False

        self.data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL']
        self.chamber_keys = ['CMAP', 'IMP', 'UB']
        self.composite_keys = ['GGAS', 'GSOLV', 'TOTAL']

        if self.chamber:
            for key in self.chamber_keys:
                if key not in self.data_keys:
                    self.data_keys.insert(3, key)

    def parse_from_file(self, basename, num_files=1, numframes=1):
        self.num_files = num_files
        self.basename = basename
        self.temperature = self.INPUT['temperature']
        self.numframes = numframes

        for key in self.data_keys:
            self[key] = EnergyVector(numframes)
        for key in self.composite_keys:
            self[key] = EnergyVector(numframes)
        AmberOutput._read(self)
        self._fill_composite_terms()

    def _print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        print_keys = list(self.data_keys)
        # Add on the composite keys
        print_keys += self.composite_keys

        # write the header
        csvwriter.writerow(['Frame #'] + print_keys)

        # write out each frame
        c = self.INPUT['nmstartframe'] if self.__class__ == NMODEout else self.INPUT['startframe']
        for i in range(self.numframes):
            csvwriter.writerow([c] + [round(self[key][i], 2) for key in print_keys])
            c += self.INPUT['nminterval'] if self.__class__ == NMODEout else self.INPUT['interval']

    def set_frame_range(self, start=None, end=None, interval=None):
        d = deepcopy(self)
        for key in d.data_keys:
            d[key] = d[key][start:end:interval]
        d._fill_composite_terms()
        return d

    def summary_output(self):
        if not self.is_read:
            raise OutputError('Cannot print summary before reading output files')
        text = [f'{self.mol.capitalize()}:']
        summary = self.summary()
        for c, row in enumerate(summary, start=1):
            key, avg, stdev, std, semp, sem = row
            if key in ['GGAS', 'TOTAL']:
                text.append('')
            if isinstance(avg, str):
                text.extend(
                    (
                        f'{key:16s} {avg:>13s} {stdev:>13s} {std:>10s} {semp:>12s} {sem:>10s}',
                        sep,
                    )
                )
            else:
                text.append(f'{key:16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}')
        return '\n'.join(text) + '\n\n'

    def summary(self):
        """ Returns a formatted string that can be printed directly to the
            output file
        """
        if not self.is_read:
            raise OutputError('Cannot print summary before reading output files')

        comp_name = 'Entropy Component' if self.__class__ in [NMODEout, QHout] else 'Energy Component'

        summary_list = [[comp_name, 'Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']]

        for key in self.data_keys:
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys:
                continue
            avg = float(self[key].mean())
            stdev = float(self[key].stdev())
            semp = float(self[key].semp())
            std = float(self[key].std())
            sem = float(self[key].sem())
            summary_list.append([key, avg, stdev, std, semp, sem])

        for key in self.composite_keys:
            # Now print out the composite terms
            avg = float(self[key].mean())
            stdev = float(self[key].stdev())
            semp = float(self[key].semp())
            std = float(self[key].std())
            sem = float(self[key].sem())
            summary_list.append([key, avg, stdev, std, semp, sem])

        return summary_list

    def _read(self):
        """        Internal reading function. This should be called at the end of __init__"""

        if self.is_read:
            return None  # don't read through them twice

        # Loop through all filenames
        for fileno in range(self.num_files):
            with open('%s.%d' % (self.basename, fileno)) as output_file:
                self._get_energies(output_file)
            self._extra_reading(fileno)
        self._fill_nmode_values()

        self.is_read = True

    def _get_energies(self, output_file):
        pass

    def _extra_reading(self, fileno):
        pass

    def _fill_nmode_values(self):
        pass

    def _fill_composite_terms(self):
        """
        Fills in the composite terms WITHOUT adding in terms we're not printing.
        This should be called after the final verbosity level has been set (based
        on whether or not certain terms need to be added in)
        """

        for key in self.composite_keys:
            self[key] = EnergyVector(self.numframes)

        for key in self.data_keys:
            for component in data_key_owner[key]:
                self[component] = self[key] + self[component]


class IEout(dict):
    """
    Interaction Entropy output
    """

    def __init__(self, INPUT, **kwargs):
        super(IEout, self).__init__(**kwargs)
        self.INPUT = INPUT

    def parse_from_dict(self, d: dict):
        self.update(d)

    def parse_from_file(self, filename, numframes=1):
        self['data'] = EnergyVector(numframes)
        with open(filename) as of:
            c = 0
            f = 0
            while line := of.readline():
                f += 1
                if line.startswith('|') or not line.split():
                    continue
                if line.startswith('IE-frames:'):
                    self['ieframes'] = int(line.strip('\n').split()[-1])
                elif line.startswith('Internal Energy SD (sigma):'):
                    self['sigma'] = float(line.strip('\n').split()[-1])
                elif line.startswith('Frame'):
                    continue
                else:
                    frame, value = line.strip('\n').split()
                    self['data'][c] = float(value)
                    c += 1
                f += 1
        self['iedata'] = self['data'][-self['ieframes']:]

    def _print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        csvwriter.writerow(['Frame #', 'Interaction Entropy'])
        f = self.INPUT['startframe']
        for d in self['data']:
            csvwriter.writerow([f] + [round(d, 2)])
            f += self.INPUT['interval']
        csvwriter.writerow([])

    def summary_output(self):
        summary = self.summary()
        text = []
        for row in summary:
            key, sigma, avg, std, sem = row
            if isinstance(avg, str):
                text.extend((f'{key:15s} {sigma:>14s} {avg:>16s} {std:>14s} {sem:>16s}', sep))
            else:
                text.append(f"{key:15s} {sigma:14.2f} {avg:16.2f} {std:14.2f} {sem:16.2f}")
        return '\n'.join(text) + '\n\n'

    def summary(self):
        """ Formatted summary of Interaction Entropy results """

        avg = float(self['data'][-self['ieframes']:].mean())
        stdev = float(self['data'][-self['ieframes']:].stdev())
        sem = float(self['data'][-self['ieframes']:].sem())

        return [
            [
                'Method',
                'σ(Int. Energy)',
                'Average',
                'SD',
                'SEM'
            ],
            ['IE', self['sigma'], avg, stdev, sem]
        ]


class C2out(dict):
    """
    C2 Entropy output
    """

    def __init__(self, **kwargs):
        super(C2out, self).__init__(**kwargs)

    def parse_from_dict(self, d):
        self.update(d)

    def parse_from_file(self, filename):
        with open(filename) as of:
            while line := of.readline():
                if line.startswith('|') or not line:
                    continue
                if line.startswith('C2 Entropy (-TΔS):'):
                    self['c2data'] = float(line.strip('\n').split()[-1])
                elif line.startswith('C2 Entropy SD:'):
                    self['c2_std'] = float(line.strip('\n').split()[-1])
                elif line.startswith('Internal Energy SD (sigma):'):
                    self['sigma'] = float(line.strip('\n').split()[-1])
                elif line.startswith('C2 Entropy CI:'):
                    self['c2_ci'] = [float(line.strip('\n').split()[-2]), float(line.strip('\n').split()[-1])]

    def summary_output(self):
        summary = self.summary()
        text = []
        for row in summary:
            key, sigma, avg, std, ci = row
            if isinstance(avg, str):
                text.extend((f'{key:15s} {sigma:>14s} {avg:>16s} {std:>14s} {ci:>16s}', sep))
            else:
                text.append(f"{key:15s} {sigma:14.2f} {avg:16.2f} {std:14.2f} {ci:>16s}")
        return '\n'.join(text) + '\n\n'

    def summary(self):
        """ Formatted summary of C2 Entropy results """

        return [
            [
                'Method',
                'σ(Int. Energy)',
                'C2 Value',
                'SD',
                'C.Inter.(95%)'
            ],
            ['C2', float(self['sigma']), float(self['c2data']), float(self['c2_std']),
             f"{self['c2_ci'][0]:.2f}-{self['c2_ci'][1]:.2f}",]
        ]


class QHout(dict):
    """ Quasi-harmonic output file class. QH output files are strange so we won't
        derive from AmberOutput
    """

    def __init__(self, filename=None, temp=298.15, **kwargs):
        super(QHout, self).__init__(**kwargs)
        self.filename = filename
        self.temperature = temp
        self.stability = False
        self._read()

    def summary_output(self):
        text = ['           Translational      Rotational      Vibrational           Total',
                'Complex   %13.4f %15.4f %16.4f %15.4f\n' % (self['complex']['TRANSLATIONAL'],
                                                             self['complex']['ROTATIONAL'],
                                                             self['complex']['VIBRATIONAL'],
                                                             self['complex']['TOTAL'],)]

        if not self.stability:
            text.extend(['Receptor  %13.4f %15.4f %16.4f %15.4f\n' % (self['receptor']['TRANSLATIONAL'],
                                                                      self['receptor']['ROTATIONAL'],
                                                                      self['receptor']['VIBRATIONAL'],
                                                                      self['receptor']['TOTAL']),
                         'Ligand    %13.4f %15.4f %16.4f %15.4f\n' % (self['ligand']['TRANSLATIONAL'],
                                                                      self['ligand']['ROTATIONAL'],
                                                                      self['ligand']['VIBRATIONAL'],
                                                                      self['ligand']['TOTAL']),
                         '',
                         '-TΔS   %13.4f %15.4f %16.4f %15.4f\n' % (self['delta']['TRANSLATIONAL'],
                                                                   self['delta']['ROTATIONAL'],
                                                                   self['delta']['VIBRATIONAL'],
                                                                   self['delta']['TOTAL'])])
        return '\n'.join(text) + '\n\n'

    def summary(self):
        """ Formatted summary of quasi-harmonic results """

        summry_list = [['', 'TRANSLATIONAL', 'ROTATIONAL', 'VIBRATIONAL', 'TOTAL'],
                       ['Complex', self['complex']['TRANSLATIONAL'], self['complex']['ROTATIONAL'],
                        self['complex']['VIBRATIONAL'], self['complex']['TOTAL']]]

        if not self.stability:
            summry_list.extend([['Receptor',
                                 self['receptor']['TRANSLATIONAL'],
                                 self['receptor']['ROTATIONAL'],
                                 self['receptor']['VIBRATIONAL'],
                                 self['receptor']['TOTAL']],
                                ['Ligand',
                                 self['ligand']['TRANSLATIONAL'],
                                 self['ligand']['ROTATIONAL'],
                                 self['ligand']['VIBRATIONAL'],
                                 self['ligand']['TOTAL']],
                                ['-TΔS',
                                 self['delta']['TRANSLATIONAL'],
                                 self['delta']['ROTATIONAL'],
                                 self['delta']['VIBRATIONAL'],
                                 self['delta']['TOTAL']]])

        return summry_list

    def _read(self):
        """ Parses the output files and fills the data arrays """
        with open(self.filename, 'r') as output:
            rawline = output.readline()
            self['complex'] = {}
            self['receptor'] = {'TOTAL': 0}
            self['ligand'] = {'TOTAL': 0}
            self['delta'] = {}
            comdone = False  # if we've done the complex yet (filled in self.com)
            recdone = False  # if we've done the receptor yet (filled in self.rec)

            # Try to fill in all found entropy values. If we can only find 1 set,
            # we're doing stability calculations
            while rawline:
                if rawline[:6] == " Total":
                    if not comdone:
                        self['complex']['TOTAL'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['TRANSLATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['ROTATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['VIBRATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        comdone = True
                    elif not recdone:
                        self['receptor']['TOTAL'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['TRANSLATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['ROTATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['VIBRATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        recdone = True
                    else:
                        self['ligand']['TOTAL'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['TRANSLATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['ROTATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['VIBRATIONAL'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        break
                rawline = output.readline()
        # end while rawline
        self.stability = not recdone

        # fill the delta if not stability
        if not self.stability:
            self['delta']['TOTAL'] = self['complex']['TOTAL'] - self['receptor']['TOTAL'] - self['ligand']['TOTAL']
            self['delta']['TRANSLATIONAL'] = (self['complex']['TRANSLATIONAL'] - self['receptor']['TRANSLATIONAL'] -
                                              self['ligand']['TRANSLATIONAL'])
            self['delta']['ROTATIONAL'] = (self['complex']['ROTATIONAL'] - self['receptor']['ROTATIONAL'] -
                                           self['ligand']['ROTATIONAL'])
            self['delta']['VIBRATIONAL'] = (self['complex']['VIBRATIONAL'] - self['receptor']['VIBRATIONAL'] -
                                            self['ligand']['VIBRATIONAL'])


class NMODEout(AmberOutput):
    """ Normal mode entropy approximation output class """
    print_levels = {'TRANSLATIONAL': 1, 'ROTATIONAL': 1, 'VIBRATIONAL': 1, 'TOTAL': 1}

    def __init__(self, mol: str, INPUT, chamber=False, **kwargs):
        super(NMODEout, self).__init__(mol, INPUT, chamber, **kwargs)

        # Ordered list of keys in the data dictionary
        self.data_keys = ['TRANSLATIONAL', 'ROTATIONAL', 'VIBRATIONAL']

        # Other aspects of AmberOutputs, which are just blank arrays
        self.composite_keys = ['TOTAL']

    def _get_energies(self, outfile):
        """ Parses the energy terms from the output file. This will parse 1 line
            at a time in order to minimize the memory requirements (we should only
            have to store a single line at a time in addition to the arrays of
            data)
        """
        while rawline := outfile.readline():
            if "|---- Entropy not Calculated---|" in rawline:
                self['TOTAL'][self.frame_idx] = np.nan
                self['TRANSLATIONAL'][self.frame_idx] = np.nan
                self['ROTATIONAL'][self.frame_idx] = np.nan
                self['VIBRATIONAL'][self.frame_idx] = np.nan
                self.frame_idx += 1

            if rawline[:6] == 'Total:':
                self['TOTAL'][self.frame_idx] = float(rawline.split()[3]) * self.temperature / 1000 * -1
                self['TRANSLATIONAL'][self.frame_idx] = (float(outfile.readline().split()[3]) * self.temperature /
                                                         1000 * -1)
                self['ROTATIONAL'][self.frame_idx] = (float(outfile.readline().split()[3]) * self.temperature / 1000
                                                      * -1)
                self['VIBRATIONAL'][self.frame_idx] = (float(outfile.readline().split()[3]) * self.temperature / 1000
                                                       * -1)
                self.frame_idx += 1


    def _fill_nmode_values(self):
        if np.isnan(self['TOTAL']).all():
            logging.warning(f'{self.mol.capitalize()}: Convergence criteria for minimized energy gradient in NMODE has not '
                            f'been\n'
                            '    satisfied in any of the frames selected. Increase the convergence criteria for\n'
                            '    minimized energy gradient (drms) or the maximum number of minimization cycles to\n '
                            '    useper snapshot in sander (maxcyc)...\n')
            self.no_nmode_convergence = True
            return

        filling = False
        for t in self.data_keys:
            if np.isnan(self[t]).any():
                filling = True
                self[t] = EnergyVector(np.nan_to_num(self[t], nan=float(np.nanmean(self[t]))))
        if filling:
            logging.warning(f'{self.mol.capitalize()}: Convergence criteria for minimized energy gradient in NMODE\n '
                            '    has not been satisfied in several frames selected. Filling "NaN" with the mean\n'
                            '    value. Please, consider to  increase the convergence criteria for minimized energy\n'
                            '    gradient (drms) or the maximum number minimization cycles to use per snapshot in\n'
                            '    sander (maxcyc)...\n')


class GBout(AmberOutput):
    """ Amber output class for normal generalized Born simulations """
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1,
                    'ESURF': 1}

    # Ordered list of keys in the data dictionary

    def __init__(self, mol, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber, **kwargs)
        self.data_keys.extend(['EGB', 'ESURF'])

    def _get_energies(self, outfile):
        """ Parses the mdout files for the GB potential terms """
        while rawline := outfile.readline():
            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'][self.frame_idx] = float(words[2])
                self['ANGLE'][self.frame_idx] = float(words[5])
                self['DIHED'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()
                if self.chamber:
                    self['UB'][self.frame_idx] = float(words[2])
                    self['IMP'][self.frame_idx] = float(words[5])
                    self['CMAP'][self.frame_idx] = float(words[8])
                    words = outfile.readline().split()
                self['VDWAALS'][self.frame_idx] = float(words[2])
                self['EEL'][self.frame_idx] = float(words[5])
                self['EGB'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()
                self['1-4 VDW'][self.frame_idx] = float(words[3])
                self['1-4 EEL'][self.frame_idx] = float(words[7])
                self.frame_idx += 1

    def _extra_reading(self, fileno):
        # Load the ESURF data from the cpptraj output
        fname = '%s.%d' % (self.basename, fileno)
        fname = fname.replace('gb.mdout', 'gb_surf.dat')
        surf_data = _get_cpptraj_surf(fname)
        for sd in surf_data:
            self['ESURF'][self.extraframe_idx] = sd * self.INPUT['surften'] + self.INPUT['surfoff']
            self.extraframe_idx += 1


class PBout(AmberOutput):

    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EPB': 1, 'ENPOLAR': 1, 'EDISPER': 1}

    def __init__(self, mol, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber, **kwargs)
        # FIXME: include Non linear PB
        self.data_keys.extend(['EPB', 'ENPOLAR', 'EDISPER'])

    def _get_energies(self, outfile):
        """ Parses the energy values from the output files """
        while rawline := outfile.readline():
            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'][self.frame_idx] = float(words[2])
                self['ANGLE'][self.frame_idx] = float(words[5])
                self['DIHED'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()
                if self.chamber:
                    self['UB'][self.frame_idx] = float(words[2])
                    self['IMP'][self.frame_idx] = float(words[5])
                    self['CMAP'][self.frame_idx] = float(words[8])
                    words = outfile.readline().split()
                self['VDWAALS'][self.frame_idx] = float(words[2])
                self['EEL'][self.frame_idx] = float(words[5])
                self['EPB'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()
                self['1-4 VDW'][self.frame_idx] = float(words[3])
                self['1-4 EEL'][self.frame_idx] = float(words[7])
                words = outfile.readline().split()
                self['ENPOLAR'][self.frame_idx] = float(words[2])
                if self.INPUT['inp'] == 2 and not self.apbs:
                    self['EDISPER'][self.frame_idx] = float(words[5])
                self.frame_idx += 1

class RISMout(AmberOutput):
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'ERISM': 1}

    def __init__(self, mol, INPUT, chamber=False, solvtype=0, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber)
        self.solvtype = solvtype
        self.data_keys.extend(['ERISM'])

    def _get_energies(self, outfile):
        """ Parses the RISM output file for energy terms """

        # Getting the RISM solvation energies requires some decision-making.
        # There are 2 possibilities (right now):
        #
        # 1. Standard free energy (solvtype==0)
        # 2. GF free energy (solvtype==1)
        # 3. PC+ GF free energy (solvtype==2)

        while rawline := outfile.readline():

            if re.match(r'(solute_epot|solutePotentialEnergy)', rawline):
                words = rawline.split()
                self['VDWAALS'][self.frame_idx] = float(words[2])
                self['EEL'][self.frame_idx] = float(words[3])
                self['BOND'][self.frame_idx] = float(words[4])
                self['ANGLE'][self.frame_idx] = float(words[5])
                self['DIHED'][self.frame_idx] = float(words[6])
                self['1-4 VDW'][self.frame_idx] = float(words[7])
                self['1-4 EEL'][self.frame_idx] = float(words[8])

            elif self.solvtype == 0 and re.match(r'(rism_exchem|rism_excessChemicalPotential)\s', rawline):
                self['ERISM'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1
            elif self.solvtype == 1 and re.match(r'(rism_exchGF|rism_excessChemicalPotentialGF)\s', rawline):
                self['ERISM'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1
            elif self.solvtype == 2 and re.match(r'(rism_exchPCPLUS|rism_excessChemicalPotentialPCPLUS)\s', rawline):
                self['ERISM'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1


class RISM_std_Out(RISMout):
    """ No polar decomp RISM output file for standard free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 0)


class RISM_gf_Out(RISMout):
    """ No polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 1)

class RISM_pcplus_Out(RISMout):
    """ No polar decomp RISM output file for PC+ free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 2)


class PolarRISMout(RISMout):
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'POLAR SOLV': 1, 'APOLAR SOLV': 1}

    def __init__(self, mol, INPUT, chamber=False, solvtype=0, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber)
        self.solvtype = solvtype
        self.data_keys.extend(['POLAR SOLV', 'APOLAR SOLV'])

    def _get_energies(self, outfile):
        """ Parses the RISM output file for energy terms """
        # Getting the RISM solvation energies requires some decision-making.
        # There are 2 possibilities (right now):
        #
        # 1. Standard free energy (solvtype==0)
        # 2. GF free energy (solvtype==1)

        while rawline := outfile.readline():

            if re.match(r'(solute_epot|solutePotentialEnergy)', rawline):
                words = rawline.split()
                self['VDWAALS'][self.frame_idx] = float(words[2])
                self['EEL'][self.frame_idx] = float(words[3])
                self['BOND'][self.frame_idx] = float(words[4])
                self['ANGLE'][self.frame_idx] = float(words[5])
                self['DIHED'][self.frame_idx] = float(words[6])
                self['1-4 VDW'][self.frame_idx] = float(words[8])
                self['1-4 EEL'][self.frame_idx] = float(words[8])

            elif self.solvtype == 0 and re.match(
                    r'(rism_polar|rism_polarExcessChemicalPotential)\s', rawline):
                self['POLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
            elif self.solvtype == 0 and re.match(
                    r'(rism_apolar|rism_apolarExcessChemicalPotential)\s', rawline):
                self['APOLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1
            elif self.solvtype == 1 and re.match(
                    r'(rism_polGF|rism_polarExcessChemicalPotentialGF)\s', rawline):
                self['POLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
            elif self.solvtype == 1 and re.match(
                    r'(rism_apolGF|rism_apolarExcessChemicalPotentialGF)\s', rawline):
                self['APOLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1
            elif self.solvtype == 2 and re.match(
                    r'(rism_polPCPLUS|rism_polarExcessChemicalPotentialPCPLUS)\s', rawline):
                self['POLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
            elif self.solvtype == 2 and re.match(
                    r'(rism_apolPCPLUS|rism_apolarExcessChemicalPotentialPCPLUS)\s', rawline):
                self['APOLAR SOLV'][self.frame_idx] = float(rawline.split()[1])
                self.frame_idx += 1


class PolarRISM_std_Out(PolarRISMout):
    """ Polar decomp RISM output file for standard free energy """

    def __init__(self, mol, INPUT, chamber=False):
        PolarRISMout.__init__(self, mol, INPUT, chamber, 0)


class PolarRISM_gf_Out(PolarRISMout):
    """ Polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, mol, INPUT, chamber=False):
        PolarRISMout.__init__(self, mol, INPUT, chamber, 1)

class PolarRISM_pcplus_Out(PolarRISMout):
    """ Polar decomp RISM output file for PC+ free energy """

    def __init__(self, mol, INPUT, chamber=False):
        PolarRISMout.__init__(self, mol, INPUT, chamber, 2)


class QMMMout(GBout):
    """ Class for QM/MM GBSA output files """
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1, 'ESURF': 1, 'ESCF': 1}

    def __init__(self, mol, INPUT, chamber=False, **kwargs):
        GBout.__init__(self, mol, INPUT, chamber, **kwargs)

        self.data_keys.extend(['ESCF'])

    def _get_energies(self, outfile):
        """ Parses the energies from a QM/MM output file. NOTE, however, that a
            QMMMout *could* just be a GBout with ESCF==0 if the QM region lies
            entirely outside this system
        """

        while rawline := outfile.readline():
            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'][self.frame_idx] = float(words[2])
                self['ANGLE'][self.frame_idx] = float(words[5])
                self['DIHED'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()

                if self.chamber:
                    self['UB'][self.frame_idx] = float(words[2])
                    self['IMP'][self.frame_idx] = float(words[5])
                    self['CMAP'][self.frame_idx] = float(words[8])
                    words = outfile.readline().split()

                self['VDWAALS'][self.frame_idx] = float(words[2])
                self['EEL'][self.frame_idx] = float(words[5])
                self['EGB'][self.frame_idx] = float(words[8])
                words = outfile.readline().split()
                self['1-4 VDW'][self.frame_idx] = float(words[3])
                self['1-4 EEL'][self.frame_idx] = float(words[7])
                words = outfile.readline().split()
                # This is where ESCF will be. Since ESCF can differ based on which
                # qmtheory was chosen, we just check to see if it's != ESURF:
                if words[0] == 'minimization':
                    continue
                elif words[0].endswith('='):
                    self['ESCF'][self.frame_idx] = float(words[1])
                else:
                    self['ESCF'][self.frame_idx] = float(words[2])
                self.frame_idx += 1


class BindingStatistics(dict):
    """ Base class for compiling the binding statistics """
    st_null = ['BOND', 'ANGLE', 'DIHED', '1-4 VDW', '1-4 EEL']

    def __init__(self, com, rec, lig, chamber=False, traj_protocol='STP', **kwargs):
        super(BindingStatistics, self).__init__(**kwargs)
        self.com = com
        self.rec = rec
        self.lig = lig
        self.mol = 'delta'
        self.numframes = self.com.numframes
        self.INPUT = self.com.INPUT
        self.chamber = chamber
        self.traj_protocol = traj_protocol
        self.inconsistent = False
        self.missing_terms = False

        self.data_keys = self.com.data_keys
        self.composite_keys = []

        try:
            self._delta()
            self.missing_terms = False
        except LengthError:
            self._delta2()
            self.missing_terms = True

    def _delta(self):
        """
        Calculates the delta statistics. Should check for any consistencies that
        would cause verbosity levels to change, and it should change them
        accordingly in the child classes
        """
        # First thing we do is check to make sure that all of the terms that
        # should *not* be printed actually cancel out (i.e. bonded terms)
        if not isinstance(self.com, NMODEout) and self.traj_protocol == 'STP':
            TINY = 0.005
            for key in self.st_null:
                diff = self.com[key] - self.rec[key] - self.lig[key]
                if diff.abs_gt(TINY):
                    self.inconsistent = True
                    logging.warning(f"{key} component is reported as inconsistent. Please, check the output file "
                                    f"for more details")
                    break

        for key in self.com.data_keys:
            if self.traj_protocol == 'STP':
                temp = self.com[key].corr_sub(self.rec[key])
                self[key] = temp.corr_sub(self.lig[key])
            else:
                self[key] = self.com[key] - self.rec[key] - self.lig[key]
        for key in self.com.composite_keys:
            self[key] = EnergyVector(self.numframes)
            self.composite_keys.append(key)
        for key in self.com.data_keys:
            if self.traj_protocol == 'STP' and key in self.st_null:
                continue
            for component in data_key_owner[key]:
                self[component] = self[key] + self[component]

    def _print_vectors(self, csvwriter):
        """ Output all of the energy terms including the differences if we're
            doing a single trajectory simulation and there are no missing terms
        """
        term_text = 'Entropy' if isinstance(self.com, NMODEout) else 'Energy'
        csvwriter.writerow([f'Complex {term_text} Terms'])
        self.com._print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow([f'Receptor {term_text} Terms'])
        self.rec._print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow([f'Ligand {term_text} Terms'])
        self.lig._print_vectors(csvwriter)
        csvwriter.writerow([])

        csvwriter.writerow([f'Delta {term_text} Terms'])
        print_keys = list(self.data_keys)
        # Add on the composite keys
        print_keys += self.composite_keys

        # write the header
        csvwriter.writerow(['Frame #'] + print_keys)

        # write out each frame
        c = self.com.INPUT['nmstartframe'] if isinstance(self.com, NMODEout) else self.com.INPUT['startframe']
        for i in range(self.numframes):
            csvwriter.writerow([c] + [round(self[key][i], 2) for key in print_keys])
            c += self.com.INPUT['nminterval'] if isinstance(self.com, NMODEout) else self.com.INPUT['interval']
        csvwriter.writerow([])

    def report_inconsistency(self, output_format: str = 'ascii'):
        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.append(['WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL TERMS AND\n'
                         'THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE!\n'
                         '\n'
                         'Some absolute differences in the internal potential terms are greater than 0.005.\n'
                         'This should not happen when using Single Trajectory Protocol!\n'
                         '\n'
                         'You can generate a detailed *.csv file with all the terms and differences as follows:\n'
                         'gmx_MMPBSA --rewrite-output -eo energy_terms_differences.csv'])
        else:
            text.append('WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL TERMS AND\n'
                         'THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE!\n'
                         '\n'
                         'Some absolute differences in the internal potential terms are greater than 0.005.\n'
                         'This should not happen when using Single Trajectory Protocol!\n'
                         '\n'
                         'You can generate a detailed *.csv file with all the terms and differences as follows:\n'
                         'gmx_MMPBSA --rewrite-output -eo energy_terms_differences.csv')
        return text if _output_format else '\n'.join(text) + '\n'

    def summary_output(self, output_format: str = 'ascii'):
        _output_format = 0 if output_format == 'ascii' else 1
        summary = self.summary()
        text = []

        if _output_format:
            text.append(['Delta (Complex - Receptor - Ligand):'])
        else:
            text.append('Delta (Complex - Receptor - Ligand):')

        for c, row in enumerate(summary, start=1):
            # Skip the composite terms, since we print those at the end
            if _output_format:
                text.append(row)
            else:
                key, avg, stdev, std, semp, sem = row
                if key in ['GGAS', 'TOTAL']:
                    text.append('')
                if isinstance(avg, str):
                    text.extend(
                        (
                            f'{key:16s} {avg:>13s} {stdev:>13s} {std:>10s} {semp:>12s} {sem:>10s}',
                            sep,
                        )
                    )

                else:
                    text.append(f'{f"Δ{key}":16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}')
        return text if _output_format else '\n'.join(text) + '\n'

    def summary(self):
        """ Returns a string printing the summary of the binding statistics """
        if isinstance(self.com, NMODEout):
            col_name = '%-16s' % 'Entropy Term'
        else:
            col_name = '%-16s' % 'Energy Component'

        summary_list = [
            [col_name] + ['Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']
        ]

        for key in self.data_keys:
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys:
                continue
            # Now print out the stats
            stdev = float(self[key].stdev())
            avg = float(self[key].mean())
            std = float(self[key].std())
            semp = float(self[key].semp())
            sem = float(self[key].sem())
            summary_list.append([key, avg, stdev, std, semp, sem])

        for key in self.composite_keys:
            # Now print out the composite terms
            stdev = float(self[key].stdev())
            avg = float(self[key].mean())
            std = float(self[key].std())
            semp = float(self[key].semp())
            sem = float(self[key].sem())
            summary_list.append([key, avg, stdev, std, semp, sem])

        return summary_list


class DeltaDeltaStatistics(dict):
    """ Base class for compiling the binding statistics. Include """
    st_null = ['BOND', 'ANGLE', 'DIHED', '1-4 VDW', '1-4 EEL']

    def __init__(self, mut, norm, **kwargs):
        super(DeltaDeltaStatistics, self).__init__(**kwargs)

        self.mut = mut
        self.norm = norm
        self.numframes = self.norm.numframes
        self.mol = self.norm.mol
        self.data_keys = self.norm.data_keys
        self.composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
        self.term_text = 'Entropy' if 'ROTATIONAL' in self.data_keys else 'Energy'
        self._delta()

    def _delta(self):
        """
        Calculates the delta statistics. Should check for any consistencies that
        would cause verbosity levels to change, and it should change them
        accordingly in the child classes
        """
        for key in self.norm:
            if key in self.composite_keys:
                continue
            self[key] = self.mut[key].corr_sub(self.norm[key])

        for key in self.composite_keys:
            self[key] = EnergyVector(self.numframes)
            # self.composite_keys.append(key)
        for key in self.data_keys:
            for component in data_key_owner[key]:
                self[component] = self[key] + self[component]

    def _print_vectors(self, csvwriter):
        """ Output all of the energy terms including the differences if we're
            doing a single trajectory simulation and there are no missing terms
        """
        csvwriter.writerow([f'Delta Delta {self.term_text} Terms (Mutant - Normal)'])

        # write the header
        csvwriter.writerow(['Frame #'] + list(self.keys()))

        # write out each frame
        c = self.norm.INPUT['nmstartframe'] if isinstance(self.norm, NMODEout) else self.norm.INPUT['startframe']
        for i in range(self.numframes):
            csvwriter.writerow([c] + [round(self[key][i], 2) for key in self])
            c += self.norm.INPUT['nminterval'] if isinstance(self.norm, NMODEout) else self.norm.INPUT['interval']
        csvwriter.writerow([])

    def summary_output(self):
        summary = self.summary()
        text = ['Delta Delta (Mutant - Normal):']
        for c, row in enumerate(summary, start=1):
            # Skip the composite terms, since we print those at the end
            key, avg, stdev, std, semp, sem = row
            if key in ['GGAS', 'TOTAL']:
                text.append('')
            if isinstance(avg, str):
                text.extend(
                    (
                        f'{key:16s} {avg:>13s} {stdev:>13s} {std:>10s} {semp:>12s} {sem:>10s}',
                        sep,
                    )
                )

            else:
                text.append(f'{f"ΔΔ{key}":16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}')
        return '\n'.join(text) + '\n'

    def summary(self):
        """ Returns a string printing the summary of the binding statistics """

        summary_list = [
            [f'{self.term_text} Component'] + ['Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']
        ]

        for key in self.norm.data_keys:
            # Now print out the stats
            stdev = float(self[key].stdev())
            avg = float(self[key].mean())
            std = float(self[key].std())
            sem = float(self[key].sem())
            semp = float(self[key].semp())
            summary_list.append([key, avg, stdev, std, semp, sem])

        for key in self.composite_keys:
            # Now print out the composite terms
            stdev = float(self[key].stdev())
            avg = float(self[key].mean())
            std = float(self[key].std())
            sem = float(self[key].sem())
            semp = float(self[key].semp())
            summary_list.append([key, avg, stdev, std, semp, sem])

        return summary_list


class DeltaIEC2Statistic(dict):
    def __init__(self, mut, norm, **kwargs):
        super(DeltaIEC2Statistic, self).__init__(**kwargs)

        self.mut = mut
        self.norm = norm
        self.meth = 'C2' if 'c2data' in self.norm else 'IE'
        self._delta()

    def _delta(self):
        """
        Calculates the delta statistics. Should check for any consistencies that
        would cause verbosity levels to change, and it should change them
        accordingly in the child classes
        """
        for key in self.norm:
            if key == 'sigma':
                self[key] = (self.mut[key] + self.norm[key]) / 2
            elif key == 'ieframes':
                self[key] = self.norm[key]
            elif key == 'c2data':
                self[key] = self.mut[key] - self.norm[key]
            elif key == 'c2_std':
                self[key] = get_std(self.mut[key], self.norm[key])
            elif key == 'c2_ci':
                continue
            else:
                self[key] = self.mut[key].corr_sub(self.norm[key])

    def _print_vectors(self, csvwriter):
        """ Output all of the energy terms including the differences if we're
            doing a single trajectory simulation and there are no missing terms
        """

        csvwriter.writerow(['Delta Delta Entropy Terms (Mutant - Normal)'])

        # write the header
        csvwriter.writerow(['Frame #'] + list(self.keys()))

        # write out each frame
        c = self.norm.INPUT['nmstartframe'] if isinstance(self.norm, NMODEout) else self.norm.INPUT['startframe']
        for i in range(self.numframes):
            csvwriter.writerow([c] + [round(self[key][i], 2) for key in self])
            c += self.norm.INPUT['nminterval'] if isinstance(self.norm, NMODEout) else self.norm.INPUT['interval']
        csvwriter.writerow([])

    def summary_output(self):
        summary = self.summary()
        text = []
        for row in summary:
            key, sigma, avg, std = row
            if isinstance(avg, str):
                text.extend((f'{key:15s} {sigma:>14s} {avg:>16s} {std:>14s}', sep))
            else:
                text.append(f"Δ{key:14s} {sigma:14.2f} {avg:16.2f} {std:14.2f}")
        return '\n'.join(text) + '\n\n'

    def summary(self):
        """ Returns a string printing the summary of the binding statistics """
        if self.meth == 'C2':
            return [
                [
                    'Method',
                    'σ(Int. Energy)',
                    'C2 Value',
                    'SD'
                ],
                ['C2', float(self['sigma']), float(self['c2data']), float(self['c2_std'])]
            ]
        avg = float(self['data'][-self['ieframes']:].mean())
        stdev = float(self['data'][-self['ieframes']:].stdev())

        return [
            [
                'Method',
                'σ(Int. Energy)',
                'Average',
                'SD'
            ],
            ['IE', self['sigma'], avg, stdev]
        ]


class DecompOut(dict):
    """ Class for decomposition output file to collect statistics and output them """
    indicator = "                    PRINT DECOMP - TOTAL ENERGIES"
    descriptions = {'TDC': 'Total Energy Decomposition:',
                    'SDC': 'Sidechain Energy Decomposition:',
                    'BDC': 'Backbone Energy Decomposition:'}

    def __init__(self, mol: str, **kwargs):
        super(DecompOut, self).__init__(**kwargs)

        self.numframes = None
        self.mut = None
        self.mol = mol
        self.decfile = None
        self.num_terms = None
        self.allowed_tokens = ('TDC',)
        self.verbose = None
        self.num_files = None
        self.resl = None
        self.basename = None
        self.csvwriter = None
        self.surften = None
        self.frame_idx = 0
        self.current_file = 0  # File counter

    def set_frame_range(self, start=0, end=None, interval=1):
        frames_updated = False
        for term in self.allowed_tokens:
            for res in self[term]:
                for et in self[term][res]:
                    self[term][res][et] = self[term][res][et][start:end:interval]
                    if not frames_updated:
                        self.numframes = len(self[term][res][et])
                        frames_updated = True
        self._fill_composite_terms()

    def parse_from_file(self, basename, resl, INPUT, surften, num_files=1, numframes=1, mut=False):
        self.basename = basename  # base name of output files
        self.resl = resl
        self.mut = mut
        self.numframes = numframes

        self.num_files = num_files  # how many MPI files we created
        self.INPUT = INPUT
        self.verbose = INPUT['dec_verbose']
        self.surften = surften  # explicitly defined since is for GB and PB models

        if self.verbose in [1, 3]:
            self.allowed_tokens = 'TDC', 'SDC', 'BDC'

        try:
            self.num_terms = int(self._get_num_terms())
        except TypeError:
            raise OutputError('DecompOut: Not a decomp output file')
        for token in self.allowed_tokens:
            self[token] = {}

        self._read()
        self._fill_composite_terms()

    def _get_num_terms(self):
        """ Gets the number of terms in the output file """
        with open('%s.%d' % (self.basename, 0), 'r') as decfile:
            lines = decfile.readlines()
        num_terms = 0
        flag = False
        for line in lines:
            if line[:3] == 'TDC':
                num_terms += 1
                flag = True
            elif flag:
                break
            # We've now gotten to the end of the Total Decomp Contribution,
            # so we know how many terms we have
        if not flag:
            raise TypeError(f"{self.basename}.0 have 0 TDC starts")
        return num_terms

    def _read(self):
        """
        Internal reading function. This should be called at the end of __init__.
        It loops through all of the output files to populate the arrays
        """
        for fileno in range(self.num_files):
            with open('%s.%d' % (self.basename, fileno)) as output_file:
                self._get_decomp_energies(output_file)

    def _get_decomp_energies(self, outfile):
        while line := outfile.readline():
            if self.frame_idx == self.numframes:
                self.frame_idx = 0
            if line[:3] in self.allowed_tokens:
                if self.mut and self.resl[int(line[4:10])].is_mutant():
                    resnum = self.resl[int(line[4:10])].mutant_string
                else:
                    resnum = self.resl[int(line[4:10])].string
                internal = float(line[11:20])
                vdw = float(line[21:30])
                eel = float(line[31:40])
                pol = float(line[41:50])
                sas = float(line[51:60]) * self.surften

                if resnum not in self[line[:3]]:
                    self[line[:3]][resnum] = {}
                    for term in ['int', 'vdw', 'eel', 'pol', 'sas']:
                        self[line[:3]][resnum][term] = EnergyVector(self.numframes)

                self[line[:3]][resnum]['int'][self.frame_idx] = internal
                self[line[:3]][resnum]['vdw'][self.frame_idx] = vdw
                self[line[:3]][resnum]['eel'][self.frame_idx] = eel
                self[line[:3]][resnum]['pol'][self.frame_idx] = pol
                self[line[:3]][resnum]['sas'][self.frame_idx] = sas

                if line[:3] == self.allowed_tokens[-1] and resnum == list(self.resl.values())[-1].string:
                    self.frame_idx += 1

    def _print_vectors(self, csvwriter):
        tokens = {'TDC': 'Total Decomposition Contribution (TDC)',
                  'SDC': 'Sidechain Decomposition Contribution (SDC)',
                  'BDC': 'Backbone Decomposition Contribution (BDC)'}
        for term in self.allowed_tokens:
            csvwriter.writerow([tokens[term]])
            csvwriter.writerow(['Frame #', 'Residue', 'Internal', 'van der Waals', 'Electrostatic', 'Polar Solvation',
                                'Non-Polar Solv.', 'TOTAL'])
            c = self.INPUT['startframe']
            for i in range(self.numframes):
                for res in self[term]:
                    csvwriter.writerow([c, res] + [round(self[term][res][key][i], 2) for key in self[term][res]])
                c += self.INPUT['interval']

    def _fill_composite_terms(self):
        for term in self:
            for res in self[term]:
                item = self[term][res][list(self[term][res].keys())[0]]
                # pair decomp scheme
                if isinstance(item, dict):
                    for res2 in self[term][res]:
                        tot = EnergyVector(self.numframes)
                        for e in self[term][res][res2]:
                            tot = tot + self[term][res][res2][e]
                        self[term][res][res2]['tot'] = tot
                else:
                    tot = EnergyVector(self.numframes)
                    for e in self[term][res]:
                        tot = tot + self[term][res][e]
                    self[term][res]['tot'] = tot

    def summary(self, output_format: str = 'ascii'):
        """ Writes the summary in ASCII format to and open output_file """

        _output_format = 0 if output_format == 'ascii' else 1

        text = []
        if _output_format:
            text.append([f'{self.mol.capitalize()}:'])
        else:
            text.append(f'{self.mol.capitalize()}:')

        for term in self:
            if _output_format:
                text.extend([[self.descriptions[term]],
                             ['Residue', 'Internal', '', '', 'van der Waals', '',
                              '', 'Electrostatic', '', '', 'Polar Solvation', '',
                              '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''],
                             [''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6])
            else:
                text.extend([self.descriptions[term],
                             'Residue        |       Internal      |    van der Waals    |    Electrostatic    |   '
                             'Polar Solvation   |   Non-Polar Solv.   |       TOTAL',
                             '-------------------------------------------------------------------------------------'
                             '-------------------------------------------------------------'])
            for res in self[term]:
                int_avg = self[term][res]['int'].mean()
                int_std = self[term][res]['int'].stdev()
                vdw_avg = self[term][res]['vdw'].mean()
                vdw_std = self[term][res]['vdw'].stdev()
                eel_avg = self[term][res]['eel'].mean()
                eel_std = self[term][res]['eel'].stdev()
                pol_avg = self[term][res]['pol'].mean()
                pol_std = self[term][res]['pol'].stdev()
                sas_avg = self[term][res]['sas'].mean()
                sas_std = self[term][res]['sas'].stdev()
                tot_avg = self[term][res]['tot'].mean()
                tot_std = self[term][res]['tot'].stdev()
                sqrt_frames = sqrt(self.numframes)

                if _output_format:
                    # FIXME: use EnergyVector.sem or EnergyVector.semp
                    text.append([res,
                                 int_avg, int_std, int_std / sqrt_frames,
                                 vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                 eel_avg, eel_std, eel_std / sqrt_frames,
                                 pol_avg, pol_std, pol_std / sqrt_frames,
                                 sas_avg, sas_std, sas_std / sqrt_frames,
                                 tot_avg, tot_std, tot_std / sqrt_frames])
                else:
                    text.append(f"{res:14s} "
                                f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            if _output_format:
                text.append([])
            else:
                text.append('')

        return text if _output_format else '\n'.join(text)


class PairDecompOut(DecompOut):
    """ Same as DecompOut, but for Pairwise decomposition """
    indicator = "                    PRINT PAIR DECOMP - TOTAL ENERGIES"

    def set_frame_range(self, start=0, end=None, interval=1):
        frames_updated = False
        for term in self.allowed_tokens:
            for res in self[term]:
                for res2 in self[term][res]:
                    for et in self[term][res][res2]:
                        self[term][res][res2][et] = self[term][res][res2][et][start:end:interval]
                        if not frames_updated:
                            self.numframes = len(self[term][res][res2][et])
                            frames_updated = True
        self._fill_composite_terms()

    def _get_decomp_energies(self, outfile):
        while line := outfile.readline():
            if line[:3] in self.allowed_tokens:
                if self.mut and self.resl[int(line[4:11])].is_mutant():
                    resnum = self.resl[int(line[4:11])].mutant_string
                else:
                    resnum = self.resl[int(line[4:11])].string
                if self.mut and self.resl[int(line[13:20])].is_mutant():
                    resnum2 = self.resl[int(line[13:20])].mutant_string
                else:
                    resnum2 = self.resl[int(line[13:20])].string

                internal = float(line[21:33])
                vdw = float(line[34:46])
                eel = float(line[47:59])
                pol = float(line[60:72])
                sas = float(line[73:85]) * self.surften

                if resnum not in self[line[:3]]:
                    self[line[:3]][resnum] = {}

                if resnum2 not in self[line[:3]][resnum]:
                    self[line[:3]][resnum][resnum2] = {}
                    for term in ['int', 'vdw', 'eel', 'pol', 'sas']:
                        self[line[:3]][resnum][resnum2][term] = EnergyVector(self.numframes)

                self[line[:3]][resnum][resnum2]['int'][self.frame_idx] = internal
                self[line[:3]][resnum][resnum2]['vdw'][self.frame_idx] = vdw
                self[line[:3]][resnum][resnum2]['eel'][self.frame_idx] = eel
                self[line[:3]][resnum][resnum2]['pol'][self.frame_idx] = pol
                self[line[:3]][resnum][resnum2]['sas'][self.frame_idx] = sas
                if line[:3] == self.allowed_tokens[-1] and resnum == list(self.resl.values())[-1].string == resnum2:
                    self.frame_idx += 1

    def _print_vectors(self, csvwriter):
        tokens = {'TDC': 'Total Decomposition Contribution (TDC)',
                  'SDC': 'Sidechain Decomposition Contribution (SDC)',
                  'BDC': 'Backbone Decomposition Contribution (BDC)'}
        for term in self.allowed_tokens:
            csvwriter.writerow([tokens[term]])
            csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal', 'van der Waals', 'Electrostatic',
                                'Polar Solvation', 'Non-Polar Solv.', 'TOTAL'])
            c = self.INPUT['startframe']
            for i in range(self.numframes):
                for res in self[term]:
                    for res2 in self[term][res]:
                        csvwriter.writerow([c, res, res2] +
                                           [round(self[term][res][res2][key][i], 2) for key in self[term][res][res2]])
                c += self.INPUT['interval']

    def summary(self, output_format: str = 'ascii'):
        """ Writes the summary in ASCII format to and open output_file """
        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[f'{self.mol.capitalize()}:']])
        else:
            text.extend([f'{self.mol.capitalize()}:'])
        for term in self:
            if _output_format:
                text.extend([[self.descriptions[term]],
                             ['Resid 1', 'Resid 2', 'Internal', '', '', 'van der Waals', '', '', 'Electrostatic',
                              '', '', 'Polar Solvation', '', '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''],
                             [''] * 2 + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6])
            else:
                text.append(self.descriptions[term] + '\n' +
                            'Resid 1        | Resid 2        |       Internal      |    van der Waals    '
                            '|    Electrostatic    |   Polar Solvation   |   Non-Polar Solv.   |       TOTAL\n' +
                            '-----------------------------------------------------------------------------'
                            '--------------------------------------------------------------------------------------')
            for res in self[term]:
                for res2 in self[term][res]:
                    int_avg = self[term][res][res2]['int'].mean()
                    int_std = self[term][res][res2]['int'].stdev()
                    vdw_avg = self[term][res][res2]['vdw'].mean()
                    vdw_std = self[term][res][res2]['vdw'].stdev()
                    eel_avg = self[term][res][res2]['eel'].mean()
                    eel_std = self[term][res][res2]['eel'].stdev()
                    pol_avg = self[term][res][res2]['pol'].mean()
                    pol_std = self[term][res][res2]['pol'].stdev()
                    sas_avg = self[term][res][res2]['sas'].mean()
                    sas_std = self[term][res][res2]['sas'].stdev()
                    tot_avg = self[term][res][res2]['tot'].mean()
                    tot_std = self[term][res][res2]['tot'].stdev()
                    sqrt_frames = sqrt(len(self[term][res][res2]['int']))
                    if _output_format:
                        # FIXME: use EnergyVector.sem or EnergyVector.semp
                        text.append([res, res2,
                                     int_avg, int_std, int_std / sqrt_frames,
                                     vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                     eel_avg, eel_std, eel_std / sqrt_frames,
                                     pol_avg, pol_std, pol_std / sqrt_frames,
                                     sas_avg, sas_std, sas_std / sqrt_frames,
                                     tot_avg, tot_std, tot_std / sqrt_frames])
                    else:
                        text.append(f"{res:14s} | {res2:14s} "
                                    f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                    f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                    f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                    f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                    f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                    f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            if _output_format:
                text.append([])
            else:
                text.append('')
        return text if _output_format else '\n'.join(text) + '\n\n'


class DecompBinding(dict):
    """ Class for decomposition binding (per-residue) """

    def __init__(self, com, rec, lig, INPUT, desc=None, **kwargs):
        """
        output should be an open file and csvfile should be a csv.writer class. If
        the output format is specified as csv, then output should be a csv.writer
        class as well.
        """
        super(DecompBinding, self).__init__(**kwargs)
        self.com, self.rec, self.lig = com, rec, lig
        self.num_terms = self.com.num_terms
        self.desc = desc  # Description
        self.INPUT = INPUT
        self.idecomp = INPUT['idecomp']
        self.verbose = INPUT['dec_verbose']
        # Set up the data for the DELTAs
        if self.verbose in [1, 3]:
            self.allowed_tokens = 'TDC', 'SDC', 'BDC'
        else:
            self.allowed_tokens = ('TDC',)
        for token in self.allowed_tokens:
            self[token] = {}

        # Parse everything
        self._parse_all_begin()

    def _print_vectors(self, csvwriter):
        tokens = {'TDC': 'Total Decomposition Contribution (TDC)',
                  'SDC': 'Sidechain Decomposition Contribution (SDC)',
                  'BDC': 'Backbone Decomposition Contribution (BDC)'}
        for term in self.allowed_tokens:
            csvwriter.writerow([tokens[term]])
            csvwriter.writerow(['Frame #', 'Residue', 'Internal', 'van der Waals', 'Electrostatic', 'Polar Solvation',
                                'Non-Polar Solv.', 'TOTAL'])
            c = self.INPUT['startframe']
            for i in range(self.com.numframes):
                for res in self[term]:
                    csvwriter.writerow([c, res] + [round(self[term][res][key][i], 2) for key in self[term][res]])
                c += self.INPUT['interval']

    def _parse_all_begin(self):
        """ Parses through all of the terms in all of the frames, but doesn't
            do any printing
        """
        # For per-residue decomp, we need terms to match up
        if self.com.num_terms != (self.rec.num_terms + self.lig.num_terms):
            raise DecompError('Mismatch in number of decomp terms!')

        for term in self.com:
            for res in self.com[term]:
                self[term][res] = {}
                other_token = self.rec[term][res] if res.startswith('R') else self.lig[term][res]
                for e in self.com[term][res]:
                    self[term][res][e] = self.com[term][res][e] - other_token[e]

    def summary(self, output_format: str = 'ascii'):

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[idecompString[self.idecomp]], [self.desc], []])
        else:
            text.extend((idecompString[self.idecomp], self.desc))
        if self.verbose > 1:
            if _output_format:
                text.extend(self.com.summary(output_format))
                text.extend(self.rec.summary(output_format))
                text.extend(self.lig.summary(output_format))
            else:
                text.extend((self.com.summary(), self.rec.summary(), self.lig.summary()))
        if _output_format:
            text.append(['DELTAS:'])
        else:
            text.append('DELTAS:')

        for term in self:
            if _output_format:
                text.extend([[DecompOut.descriptions[term]],
                             ['Residue', 'Internal', '', '', 'van der Waals', '', '', 'Electrostatic',
                              '', '', 'Polar Solvation', '', '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''],
                             [''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6])
            else:
                text.extend([DecompOut.descriptions[term],
                             'Residue        |       Internal      |    van der Waals    |    Electrostatic    '
                             '|   Polar Solvation   |    Non-Polar Solv.  |       TOTAL',
                             '----------------------------------------------------------------------------------'
                             '----------------------------------------------------------------'])
            for res in self[term]:
                int_avg = self[term][res]['int'].mean()
                int_std = self[term][res]['int'].stdev()
                vdw_avg = self[term][res]['vdw'].mean()
                vdw_std = self[term][res]['vdw'].stdev()
                eel_avg = self[term][res]['eel'].mean()
                eel_std = self[term][res]['eel'].stdev()
                pol_avg = self[term][res]['pol'].mean()
                pol_std = self[term][res]['pol'].stdev()
                sas_avg = self[term][res]['sas'].mean()
                sas_std = self[term][res]['sas'].stdev()
                tot_avg = self[term][res]['tot'].mean()
                tot_std = self[term][res]['tot'].stdev()
                sqrt_frames = sqrt(len(self[term][res]['int']))
                if _output_format:
                    # FIXME: use EnergyVector.sem or EnergyVector.semp
                    text.append([res,
                                 int_avg, int_std, int_std / sqrt_frames,
                                 vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                 eel_avg, eel_std, eel_std / sqrt_frames,
                                 pol_avg, pol_std, pol_std / sqrt_frames,
                                 sas_avg, sas_std, sas_std / sqrt_frames,
                                 tot_avg, tot_std, tot_std / sqrt_frames])
                else:
                    text.append(f"{res:14s} "
                                f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            if _output_format:
                text.append([])
            else:
                text.append('')
        return text if _output_format else '\n'.join(text)


class PairDecompBinding(DecompBinding):
    """ Class for decomposition binding (pairwise) """

    def _parse_all_begin(self):
        """ Parses through all of the terms in all of the frames, but doesn't
            do any printing
        """
        for term in self.com:
            for res in self.com[term]:
                self[term][res] = {}
                for res2 in self.com[term][res]:
                    self[term][res][res2] = {}
                    if res.startswith('R') and res2.startswith('R'):
                        # Both residues are in the receptor -- pull the next one
                        other_token = self.rec[term][res][res2]
                    elif (
                            res.startswith('R')
                            and not res2.startswith('R')
                            or not res.startswith('R')
                            and not res2.startswith('L')
                    ):
                        other_token = {}
                    else:
                        # Both residues are in the ligand -- pull the next one
                        other_token = self.lig[term][res][res2]
                    for e in self.com[term][res][res2]:
                        if other_token:
                            self[term][res][res2][e] = self.com[term][res][res2][e] - other_token[e]
                        else:
                            self[term][res][res2][e] = self.com[term][res][res2][e]

    def _print_vectors(self, csvwriter):
        tokens = {'TDC': 'Total Decomposition Contribution (TDC)',
                  'SDC': 'Sidechain Decomposition Contribution (SDC)',
                  'BDC': 'Backbone Decomposition Contribution (BDC)'}
        for term in self.allowed_tokens:
            csvwriter.writerow([tokens[term]])
            csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal', 'van der Waals', 'Electrostatic',
                                'Polar Solvation', 'Non-Polar Solv.', 'TOTAL'])
            c = self.com.INPUT['startframe']
            for i in range(self.com.numframes):
                for res in self[term]:
                    for res2 in self[term][res]:
                        csvwriter.writerow([c, res, res2] +
                                           [round(self[term][res][res2][key][i], 2) for key in self[term][res][res2]])
                c += self.com.INPUT['interval']

    def summary(self, output_format: str = 'ascii'):

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[idecompString[self.idecomp]], [self.desc], []])
        else:
            text.extend((idecompString[self.idecomp], self.desc, ''))
        if self.verbose > 1:
            if _output_format:
                text.extend(self.com.summary(output_format))
                text.extend(self.rec.summary(output_format))
                text.extend(self.lig.summary(output_format))
            else:
                text.extend((self.com.summary(), self.rec.summary(), self.lig.summary()))
        if _output_format:
            text.append(['DELTAS:'])
        else:
            text.append('DELTAS:')

        for term in self:
            if _output_format:
                text.extend([[DecompOut.descriptions[term]],
                             ['Resid 1', 'Resid 2', 'Internal', '', '', 'van der Waals', '', '', 'Electrostatic', '',
                              '', 'Polar Solvation', '', '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''],
                             ['', ''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6])
            else:
                text.extend([DecompOut.descriptions[term],
                             'Resid 1        | Resid 2        |       Internal      |    van der Waals    '
                             '|    Electrostatic    |   Polar Solvation   |    Non-Polar Solv.  |       TOTAL',
                             '-----------------------------------------------------------------------------------------'
                             '--------------------------------------------------------------------------'])

            for res in self[term]:
                for res2 in self[term][res]:
                    int_avg = self[term][res][res2]['int'].mean()
                    int_std = self[term][res][res2]['int'].stdev()
                    vdw_avg = self[term][res][res2]['vdw'].mean()
                    vdw_std = self[term][res][res2]['vdw'].stdev()
                    eel_avg = self[term][res][res2]['eel'].mean()
                    eel_std = self[term][res][res2]['eel'].stdev()
                    pol_avg = self[term][res][res2]['pol'].mean()
                    pol_std = self[term][res][res2]['pol'].stdev()
                    sas_avg = self[term][res][res2]['sas'].mean()
                    sas_std = self[term][res][res2]['sas'].stdev()
                    tot_avg = self[term][res][res2]['tot'].mean()
                    tot_std = self[term][res][res2]['tot'].stdev()
                    sqrt_frames = sqrt(len(self[term][res][res2]['int']))

                    if _output_format:
                        # FIXME: use EnergyVector.sem or EnergyVector.semp
                        text.append([res, res2,
                                     int_avg, int_std, int_std / sqrt_frames,
                                     vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                     eel_avg, eel_std, eel_std / sqrt_frames,
                                     pol_avg, pol_std, pol_std / sqrt_frames,
                                     sas_avg, sas_std, sas_std / sqrt_frames,
                                     tot_avg, tot_std, tot_std / sqrt_frames])
                    else:
                        text.append(f"{res:14s} | {res2:14s} "
                                    f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                    f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                    f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                    f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                    f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                    f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            if _output_format:
                text.append([])
            else:
                text.append('')

        return text if _output_format else '\n'.join(text) + '\n\n'


def _get_cpptraj_surf(fname):
    """
    This function will parse out the surface areas printed out by cpptraj in a
    standard data file and return it as an EnergyVector instance.
    """
    f = open(fname, 'r')
    vec = EnergyVector()
    for line in f:
        if line.startswith('#'):
            continue
        vec = vec.append(float(line.split()[1]))

    return vec
