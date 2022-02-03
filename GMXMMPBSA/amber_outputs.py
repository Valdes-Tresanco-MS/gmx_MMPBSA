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
from math import sqrt
from GMXMMPBSA.exceptions import (OutputError, LengthError, DecompError)
from GMXMMPBSA.utils import get_std
import h5py
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
        self.com_std = getattr(obj, 'com_stdev', None)

    def stdev(self):
        return self.com_std or self.std()

    def append(self, values):
        return EnergyVector(np.append(self, values))

    def avg(self):
        return np.average(self)

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


class AmberOutput(dict):
    """
    Base Amber output class. It takes a basename as a file name and parses
    through all of the thread-specific output files (assumed to have the suffix
    .# where # spans from 0 to num_files - 1
    """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL']
    chamber_keys = ['CMAP', 'IMP', 'UB']
    # Dictionary that maps each data key to their respective composite keys
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
                      }
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EPOL': 1}

    def __init__(self, mol: str, INPUT, chamber=False, **kwargs):
        super(AmberOutput, self).__init__(**kwargs)
        self.mol = mol
        self.INPUT = INPUT
        self.chamber = chamber
        self.basename = None
        self.num_files = None
        self.is_read = False
        self.apbs = INPUT['sander_apbs']

        if self.chamber:
            for key in self.chamber_keys:
                if key not in self.data_keys:
                    self.data_keys.insert(3, key)

    def parse_from_file(self, basename, num_files=1):
        self.num_files = num_files
        self.basename = basename

        for key in self.data_keys:
            self[key] = EnergyVector()
        for key in self.composite_keys:
            self[key] = EnergyVector()
        AmberOutput._read(self)
        self._fill_composite_terms()

    def parse_from_h5(self, d: dict):
        for key in d:
            if key in ['']:
                continue
            self[key] = EnergyVector(d[key][()])
        self.is_read = True
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
        c = self.INPUT['startframe']
        for i in range(len(self[print_keys[0]])):
            csvwriter.writerow([c] + [round(self[key][i], 4) for key in print_keys])
            c += self.INPUT['interval']

    def set_frame_range(self, start=0, end=None, interval=1):
        for key in self.data_keys:
            self[key] = self[key][start:end:interval]
        self._fill_composite_terms()

    def summary(self, output_format: str = 'ascii'):
        """ Returns a formatted string that can be printed directly to the
            output file
        """
        if not self.is_read:
            raise OutputError('Cannot print summary before reading output files')

        _output_format = 0 if output_format == 'ascii' else 1

        text = []

        if _output_format:
            text.extend([[self.mol.capitalize() + ':'],
                         ['Energy Component', 'Average', 'SD(Prop.)', 'SD', 'SEM', 'SD(Prop.)']])
        else:
            text.extend([self.mol.capitalize() + ':',
                         'Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM',
                         '-------------------------------------------------------------------------------'])

        for key in self.data_keys:
            # Skip terms we don't want to print
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys:
                continue
            # Skip any terms that have zero as every single element (i.e. EDISPER)
            if self.INPUT['sander_apbs'] and key == 'EDISPER':
                continue
            if self[key] == 0:
                continue
            avg = self[key].mean()
            stdev = self[key].stdev()
            semp = stdev / sqrt(len(self[key]))
            std = self[key].std()
            sem = std / sqrt(len(self[key]))
            if _output_format:
                text.append([key, avg, stdev, std, semp, sem])
            else:
                text.append(f'{key:16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}')

        text.append('')
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'TOTAL' and not _output_format:
                text.append('')
            avg = self[key].mean()
            stdev = self[key].stdev()
            semp = stdev / sqrt(len(self[key]))
            std = self[key].std()
            sem = std / sqrt(len(self[key]))
            if _output_format:
                text.append([key, avg, stdev, std, semp, sem])
            else:
                text.append(f'{key:16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}')

        return text if _output_format else '\n'.join(text) + '\n\n'

    def _read(self):
        """
        Internal reading function. This should be called at the end of __init__.
        It loops through all of the output files to populate the arrays
        """
        if self.is_read:
            return None  # don't read through them twice

        for fileno in range(self.num_files):
            with open('%s.%d' % (self.basename, fileno)) as output_file:
                # output_file = open('%s.%d' % (self.basename, fileno), 'r')
                self._get_energies(output_file)
            # output_file.close()
            # If we have to get energies elsewhere (e.g., with GB and ESURF), do
            # that here. This is an empty function when unnecessary
            self._extra_reading(fileno)

        self.is_read = True

    def _extra_reading(self, fileno):
        pass

    def _fill_composite_terms(self):
        """
        Fills in the composite terms WITHOUT adding in terms we're not printing.
        This should be called after the final verbosity level has been set (based
        on whether or not certain terms need to be added in)
        """

        for key in self.composite_keys:
            self[key] = EnergyVector(len(self['EEL']))

        for key in self.data_keys:
            for component in self.data_key_owner[key]:
                self[component] = self[key] + self[component]


class IEout(dict):
    """
    Interaction Entropy output
    """
    def __init__(self, INPUT, **kwargs):
        super(IEout, self).__init__(**kwargs)
        self.INPUT = INPUT

    def parse_from_dict(self, model, d: dict):
        self[model] = {}
        for term, dat in d.items():
            self[model][term] = EnergyVector(dat) if term in ['data', 'iedata'] else dat

    def parse_from_h5(self, d):
        for model in d:
            self[model] = {}
            for term in d[model]:
                if term in ['data', 'iedata']:
                    self[model][term] = EnergyVector(d[model][term][()])
                else:
                    self[model][term] = d[model][term][()]

    def _print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        for model in self:
            csvwriter.writerow([f'Model {model}'])
            csvwriter.writerow(['Frame #', 'Interaction Entropy'])
            f = self.INPUT['startframe']
            for d in self[model]['data']:
                csvwriter.writerow([f] + [d])
                f += self.INPUT['interval']
            csvwriter.writerow([])

    def summary(self, output_format: str = 'ascii'):
        """ Formatted summary of Interaction Entropy results """

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.append(['Model', 'σ(Int. Energy)', 'Average', 'Std. Dev.', 'Std. Err. of Mean'])
        else:
            text.append('Model           σ(Int. Energy)      Average       Std. Dev.   Std. Err. of Mean\n' +
                        '-------------------------------------------------------------------------------')

        for model in self:
            avg = self[model]['iedata'].mean()
            stdev = self[model]['iedata'].stdev()
            if _output_format:
                text.append([model, self[model]['sigma'], avg, stdev, stdev / sqrt(len(self[model]['iedata']))])
            else:
                text.append('%-14s %10.3f %16.3f %15.3f %19.3f\n' % (model, self[model]['sigma'], avg, stdev,
                                                                     stdev / sqrt(len(self[model]['iedata']))))
        return text if _output_format else '\n'.join(text)


class C2out(dict):
    """
    C2 Entropy output
    """

    def __init__(self, **kwargs):
        super(C2out, self).__init__(**kwargs)

    def parse_from_h5(self, d):
        for model in d:
            self[model] = {}
            for term in d[model]:
                self[model][term] = d[model][term][()]

    def summary(self, output_format: str = 'ascii'):
        """ Formatted summary of C2 Entropy results """

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.append(['Model', 'σ(Int. Energy)', 'C2 Value', 'Std. Dev.', 'Conf. Interv. (95%)'])
        else:
            text.append('Model           σ(Int. Energy)    C2 Value         Std. Dev.   Conf. Interv. (95%)\n' +
                        '-------------------------------------------------------------------------------')

        for model in self:
            if _output_format:
                text.append([model, self[model]['sigma'], self[model]['c2data'], self[model]['c2_std'],
                             f"{self[model]['c2_ci'][0]}-{self[model]['c2_ci'][1]}"])
            else:
                text.append('%-14s %10.3f %15.3f %15.3f %13.3f-%5.3f\n' % (model, self[model]['sigma'],
                                                                           self[model]['c2data'],
                                                                           self[model]['c2_std'],
                                                                           self[model]['c2_ci'][0],
                                                                           self[model]['c2_ci'][1]))
        return text if _output_format else '\n'.join(text)


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

    def parse_from_h5(self, d):
        # key: complex, receptor, ligand, delta
        for key in d:
            # key1: Translational, Rotational, Vibrational, Total
            self[key] = {}
            for key1 in d[key]:
                self[key][key1] = d[key][key1][()]

    def summary(self, output_format: str = 'ascii'):
        """ Formatted summary of quasi-harmonic results """

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.append(['', 'Translational', 'Rotational', 'Vibrational', 'Total'])
            text.append(['Complex', self['complex']['Translational'], self['complex']['Rotational'],
                         self['complex']['Vibrational'], self['complex']['Total']])
            if not self.stability:
                text.extend([['Receptor', self['receptor']['Translational'], self['receptor']['Rotational'],
                              self['receptor']['Vibrational'], self['receptor']['Total']],
                             ['Ligand', self['ligand']['Translational'], self['ligand']['Rotational'],
                              self['ligand']['Vibrational'], self['ligand']['Total']]])
                text.append([])
                text.append(['-TΔS', self['delta']['Translational'], self['delta']['Rotational'],
                             self['delta']['Vibrational'], self['delta']['Total']])
            return text
        else:
            text.append('           Translational      Rotational      Vibrational           Total')
            text.append('Complex   %13.4f %15.4f %16.4f %15.4f\n' % (self['complex']['Translational'],
                                                                     self['complex']['Rotational'],
                                                                     self['complex']['Vibrational'],
                                                                     self['complex']['Total']))
            if not self.stability:
                text.extend(['Receptor  %13.4f %15.4f %16.4f %15.4f\n' % (self['receptor']['Translational'],
                                                                          self['receptor']['Rotational'],
                                                                          self['receptor']['Vibrational'],
                                                                          self['receptor']['Total']),
                             'Ligand    %13.4f %15.4f %16.4f %15.4f\n' % (self['ligand']['Translational'],
                                                                          self['ligand']['Rotational'],
                                                                          self['ligand']['Vibrational'],
                                                                          self['ligand']['Total'])])
                text.append('')
                text.append('-TΔS   %13.4f %15.4f %16.4f %15.4f\n' % (self['delta']['Translational'],
                                                                      self['delta']['Rotational'],
                                                                      self['delta']['Vibrational'],
                                                                      self['delta']['Total']))
            return '\n'.join(text) + '\n\n'

    def _read(self):
        """ Parses the output files and fills the data arrays """
        with open(self.filename, 'r') as output:
            rawline = output.readline()
            self['complex'] = {}
            self['receptor'] = {'Total': 0}
            self['ligand'] = {'Total': 0}
            self['delta'] = {}
            comdone = False  # if we've done the complex yet (filled in self.com)
            recdone = False  # if we've done the receptor yet (filled in self.rec)

            # Try to fill in all found entropy values. If we can only find 1 set,
            # we're doing stability calculations
            while rawline:
                if rawline[:6] == " Total":
                    if not comdone:
                        self['complex']['Total'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['Translational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['Rotational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['complex']['Vibrational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        comdone = True
                    elif not recdone:
                        self['receptor']['Total'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['Translational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['Rotational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['receptor']['Vibrational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        recdone = True
                    else:
                        self['ligand']['Total'] = (float(rawline.split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['Translational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['Rotational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        self['ligand']['Vibrational'] = (
                                float(output.readline().split()[3]) * self.temperature / 1000 * -1)
                        break
                rawline = output.readline()
        # end while rawline
        self.stability = not recdone

        # fill the delta if not stability
        if not self.stability:
            self['delta']['Total'] = self['complex']['Total'] - self['receptor']['Total'] - self['ligand']['Total']
            self['delta']['Translational'] = (self['complex']['Translational'] - self['receptor']['Translational'] -
                                              self['ligand']['Translational'])
            self['delta']['Rotational'] = (self['complex']['Rotational'] - self['receptor']['Rotational'] -
                                           self['ligand']['Rotational'])
            self['delta']['Vibrational'] = (self['complex']['Vibrational'] - self['receptor']['Vibrational'] -
                                            self['ligand']['Vibrational'])


class NMODEout(dict):
    """ Normal mode entropy approximation output class """
    # Ordered list of keys in the data dictionary
    data_keys = ['Translational', 'Rotational', 'Vibrational']

    # Other aspects of AmberOutputs, which are just blank arrays
    composite_keys = ['Total']
    data_key_owner = {'Translational': ['Total'], 'Rotational': ['Total'], 'Vibrational': ['Total']}
    print_levels = {'Translational': 1, 'Rotational': 1, 'Vibrational': 1, 'Total': 1}

    def __init__(self, mol, **kwargs):
        super(NMODEout, self).__init__(**kwargs)
        self.mol = mol
        self.is_read = False

    def parse_from_file(self, basename, INPUT, num_files=1, chamber=False):
        self.basename = basename

        for key in self.data_keys:
            self[key] = EnergyVector()

        if chamber:
            logging.warning('nmode is incompatible with chamber topologies!')
        self.temperature = INPUT['temperature']
        self.num_files = num_files
        self.is_read = False
        self._read()
        self._fill_composite_terms()

    def parse_from_h5(self, d):
        # key: complex, receptor, ligand, delta
        for key in d:
            # key1: Translational, Rotational, Vibrational, Total
            self[key] = {}
            for key1 in d[key]:
                self[key][key1] = EnergyVector(d[key][key1][()])

    def print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        # print header
        csvwriter.writerow(['Frame #'] + self.data_keys)

        # print data
        for i in range(len(self[self.data_keys[0]])):
            csvwriter.writerow([i] + [self[key][i] for key in self.data_keys])

    def summary(self, output_format: str = 'ascii'):
        """ Returns the formatted string of output summary """

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[self.mol.capitalize() + ':'],
                         ['Entropy Term', 'Average', 'Std. Dev.', 'Std. Err. of the Mean']])
        else:
            text.append(self.mol.capitalize() + ':\n'
                        'Entropy Term                Average              Std. Dev.   Std. Err. of Mean\n' +
                        '-------------------------------------------------------------------------------')

        for key in self.data_keys:
            avg = self[key].mean()
            stdev = self[key].stdev()
            if _output_format:
                text.append([key, avg, stdev, stdev / sqrt(len(self[key]))])
            else:
                text.append('%-14s %20.4f %21.4f %19.4f\n' % (key, avg, stdev, stdev / sqrt(len(self[key]))))
        return text if _output_format else '\n'.join(text) + '\n\n'

    def _read(self):
        """ Internal reading function to populate the data arrays """

        if self.is_read:
            return   # don't read through again

        # Loop through all filenames
        for fileno in range(self.num_files):
            with open('%s.%d' % (self.basename, fileno), 'r') as output_file:
                self._get_energies(output_file)
        self.is_read = True

    def _get_energies(self, outfile):
        """ Parses the energy terms from the output file. This will parse 1 line
            at a time in order to minimize the memory requirements (we should only
            have to store a single line at a time in addition to the arrays of
            data)
        """
        rawline = outfile.readline()

        while rawline:
            if rawline[:35] == '   |---- Entropy not Calculated---|':
                sys.stderr.write('Not all frames minimized within tolerance')

            if rawline[:6] == 'Total:':
                self['Total'] = self['Total'].append(float(rawline.split()[3]) * self.temperature / 1000 * -1)
                self['Translational'] = self['Translational'].append(
                    float(outfile.readline().split()[3]) * self.temperature / 1000 * -1)
                self['Rotational'] = self['Rotational'].append(
                    float(outfile.readline().split()[3]) * self.temperature / 1000 * -1)
                self['Vibrational'] = self['Vibrational'].append(
                    float(outfile.readline().split()[3]) * self.temperature / 1000 * -1)

            rawline = outfile.readline()

    def _fill_composite_terms(self):
        for key in self.composite_keys:
            self[key] = EnergyVector(len(self['Translational']))

        for key in self.data_keys:
            for component in self.data_key_owner[key]:
                self[component] = self[key] + self[component]


class GBout(AmberOutput):
    """ Amber output class for normal generalized Born simulations """
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF']
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1,
                    'ESURF': 1}

    # Ordered list of keys in the data dictionary

    def __init__(self, mol, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber, **kwargs)

    def _get_energies(self, outfile):
        """ Parses the mdout files for the GB potential terms """
        rawline = outfile.readline()
        while rawline:
            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'] = self['BOND'].append(float(words[2]))
                self['ANGLE'] = self['ANGLE'].append(float(words[5]))
                self['DIHED'] = self['DIHED'].append(float(words[8]))
                words = outfile.readline().split()
                if self.chamber:
                    self['UB'] = self['UB'].append(float(words[2]))
                    self['IMP'] = self['IMP'].append(float(words[5]))
                    self['CMAP'] = self['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()
                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[5]))
                self['EGB'] = self['EGB'].append(float(words[8]))
                words = outfile.readline().split()
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[3]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[7]))
            rawline = outfile.readline()

    def _extra_reading(self, fileno):
        # Load the ESURF data from the cpptraj output
        fname = '%s.%d' % (self.basename, fileno)
        fname = fname.replace('gb.mdout', 'gb_surf.dat')
        surf_data = _get_cpptraj_surf(fname)
        self['ESURF'] = self['ESURF'].append((surf_data * self.INPUT['surften']) + self.INPUT['surfoff'])


class PBout(AmberOutput):
    # Ordered list of keys in the data dictionary
    # FIXME: include Non linear PB
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EPB', 'ENPOLAR']

    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EPB': 1, 'ENPOLAR': 1, 'EDISPER': 1}

    def __init__(self, mol, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber, **kwargs)

        if self.INPUT['inp'] == 2 and 'EDISPER' not in self.data_keys:
            self.data_keys += ['EDISPER']

    def _get_energies(self, outfile):
        """ Parses the energy values from the output files """
        rawline = outfile.readline()
        while rawline:
            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'] = self['BOND'].append(float(words[2]))
                self['ANGLE'] = self['ANGLE'].append(float(words[5]))
                self['DIHED'] = self['DIHED'].append(float(words[8]))
                words = outfile.readline().split()
                if self.chamber:
                    self['UB'] = self['UB'].append(float(words[2]))
                    self['IMP'] = self['IMP'].append(float(words[5]))
                    self['CMAP'] = self['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()
                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[5]))
                self['EPB'] = self['EPB'].append(float(words[8]))
                words = outfile.readline().split()
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[3]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                self['ENPOLAR'] = self['ENPOLAR'].append(float(words[2]))
                if self.INPUT['inp'] == 2 and not self.apbs:
                    self['EDISPER'] = self['EDISPER'].append(float(words[5]))
            rawline = outfile.readline()


class RISMout(AmberOutput):
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'ERISM']
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'ERISM': 1}

    def __init__(self, mol, INPUT, chamber=False, solvtype=0, **kwargs):
        AmberOutput.__init__(self, mol, INPUT, chamber)
        # FIXME: deault value
        self.solvtype = solvtype

    def _get_energies(self, outfile):
        """ Parses the RISM output file for energy terms """

        # Getting the RISM solvation energies requires some decision-making.
        # There are 2 possibilities (right now):
        #
        # 1. Standard free energy (solvtype==0)
        # 2. GF free energy (solvtype==1)

        rawline = outfile.readline()

        while rawline:

            if re.match(r'(solute_epot|solutePotentialEnergy)', rawline):
                words = rawline.split()
                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[3]))
                self['BOND'] = self['BOND'].append(float(words[4]))
                self['ANGLE'] = self['ANGLE'].append(float(words[5]))
                self['DIHED'] = self['DIHED'].append(float(words[6]))
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[7]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[8]))

            elif self.solvtype == 0 and re.match(r'(rism_exchem|rism_excessChemicalPotential)\s', rawline):
                self['ERISM'] = self['ERISM'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(r'(rism_exchGF|rism_excessChemicalPotentialGF)\s', rawline):
                self['ERISM'] = self['ERISM'].append(float(rawline.split()[1]))

            rawline = outfile.readline()


class RISM_std_Out(RISMout):
    """ No polar decomp RISM output file for standard free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 0)


class RISM_gf_Out(RISMout):
    """ No polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 1)


class PolarRISMout(RISMout):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'ERISM', 'POLAR SOLV', 'APOLAR SOLV']
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'POLAR SOLV': 1, 'APOLAR SOLV': 1}

    def _get_energies(self, outfile):
        """ Parses the RISM output file for energy terms """
        # Getting the RISM solvation energies requires some decision-making.
        # There are 2 possibilities (right now):
        #
        # 1. Standard free energy (solvtype==0)
        # 2. GF free energy (solvtype==1)
        rawline = outfile.readline()

        while rawline:

            if re.match(r'(solute_epot|solutePotentialEnergy)', rawline):
                words = rawline.split()
                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[3]))
                self['BOND'] = self['BOND'].append(float(words[4]))
                self['ANGLE'] = self['ANGLE'].append(float(words[5]))
                self['DIHED'] = self['DIHED'].append(float(words[6]))
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[8]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[8]))

            elif self.solvtype == 0 and re.match(
                    r'(rism_polar|rism_polarExcessChemicalPotential)\s', rawline):
                self['POLAR SOLV'] = self['POLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 0 and re.match(
                    r'(rism_apolar|rism_apolarExcessChemicalPotential)\s', rawline):
                self['APOLAR SOLV'] = self['APOLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(
                    r'(rism_polGF|rism_polarExcessChemicalPotentialGF)\s', rawline):
                self['POLAR SOLV'] = self['POLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(
                    r'(rism_apolGF|rism_apolarExcessChemicalPotentialGF)\s', rawline):
                self['APOLAR SOLV'] = self['APOLAR SOLV'].append(float(rawline.split()[1]))

            rawline = outfile.readline()


class PolarRISM_std_Out(PolarRISMout):
    """ Polar decomp RISM output file for standard free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 0)


class PolarRISM_gf_Out(PolarRISMout):
    """ Polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, mol, INPUT, chamber=False):
        RISMout.__init__(self, mol, INPUT, chamber, 1)


class QMMMout(GBout):
    """ Class for QM/MM GBSA output files """
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'ESCF']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1, 'ESURF': 1, 'ESCF': 1}

    def _get_energies(self, outfile):
        """ Parses the energies from a QM/MM output file. NOTE, however, that a
            QMMMout *could* just be a GBout with ESCF==0 if the QM region lies
            entirely outside this system
        """

        rawline = outfile.readline()

        while rawline:

            if rawline[:5] == ' BOND':
                words = rawline.split()
                self['BOND'] = self['BOND'].append(float(words[2]))
                self['ANGLE'] = self['ANGLE'].append(float(words[5]))
                self['DIHED'] = self['DIHED'].append(float(words[8]))
                words = outfile.readline().split()

                if self.chamber:
                    self['UB'] = self['UB'].append(float(words[2]))
                    self['IMP'] = self['IMP'].append(float(words[5]))
                    self['CMAP'] = self['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()

                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[5]))
                self['EGB'] = self['EGB'].append(float(words[8]))
                words = outfile.readline().split()
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[3]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                # This is where ESCF will be. Since ESCF can differ based on which
                # qmtheory was chosen, we just check to see if it's != ESURF:
                if words[0] == 'minimization':
                    self['ESCF'] = self['ESCF'].append(0.0)
                elif words[0].endswith('='):
                    self['ESCF'] = self['ESCF'].append(float(words[1]))
                else:
                    self['ESCF'] = self['ESCF'].append(float(words[2]))
            rawline = outfile.readline()


class BindingStatistics(dict):
    """ Base class for compiling the binding statistics """
    st_null = ['BOND', 'ANGLE', 'DIHED', '1-4 VDW', '1-4 EEL']

    def __init__(self, com, rec, lig, chamber=False, traj_protocol='STP', **kwargs):
        super(BindingStatistics, self).__init__(**kwargs)
        self.com = com
        self.rec = rec
        self.lig = lig
        self.chamber = chamber
        self.traj_protocol = traj_protocol
        self.inconsistent = False
        self.missing_terms = False

        self.data_keys = self.com.data_keys
        self.composite_keys = []


        self.print_levels = self.com.print_levels
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
        if self.traj_protocol == 'STP':
            TINY = 0.005
            for key in self.st_null:
                diff = self.com[key] - self.rec[key] - self.lig[key]
                if diff.abs_gt(TINY):
                    self.inconsistent = True
                    break

        for key in self.com.data_keys:
            self[key] = self.com[key] - self.rec[key] - self.lig[key]
        for key in self.com.composite_keys:
            self[key] = EnergyVector(len(self.com['VDWAALS']))
            self.composite_keys.append(key)
        for key in self.com.data_keys:
            if self.traj_protocol == 'STP' and key in self.st_null:
                continue
            for component in self.com.data_key_owner[key]:
                self[component] = self[key] + self[component]

    def _print_vectors(self, csvwriter):
        """ Output all of the energy terms including the differences if we're
            doing a single trajectory simulation and there are no missing terms
        """
        csvwriter.writerow(['Complex Energy Terms'])
        self.com._print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow(['Receptor Energy Terms'])
        self.rec._print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow(['Ligand Energy Terms'])
        self.lig._print_vectors(csvwriter)
        csvwriter.writerow([])

        csvwriter.writerow(['Delta Energy Terms'])
        print_keys = list(self.data_keys)
        # Add on the composite keys
        print_keys += self.composite_keys

        # write the header
        csvwriter.writerow(['Frame #'] + print_keys)

        # write out each frame
        c = self.com.INPUT['startframe']
        for i in range(len(self[print_keys[0]])):
            csvwriter.writerow([c] + [round(self[key][i], 4) for key in print_keys])
            c += self.com.INPUT['interval']
        csvwriter.writerow([])

    def summary(self, output_format: str = 'ascii'):
        """ Returns a string printing the summary of the binding statistics """

        _output_format = 0 if output_format == 'ascii' else 1
        text = []

        if self.inconsistent:
            if _output_format:
                text.append(['WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL' +
                             'TERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE'])
            else:
                text.append('WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL' +
                            '\nTERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE\n')
        if self.com.INPUT['verbose']:
            if _output_format:
                text.extend(self.com.summary(output_format))
                text.extend(self.rec.summary(output_format))
                text.extend(self.lig.summary(output_format))
            else:
                text.append(self.com.summary())
                text.append(self.rec.summary())
                text.append(self.lig.summary())

        if isinstance(self.com, NMODEout):
            col_name = '%-16s' % 'Entropy Term'
        else:
            col_name = '%-16s' % 'Energy Component'

        if _output_format:
            text.extend([['Delta (Complex - Receptor - Ligand):'],
                         [col_name] + ['Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']])
        else:
            text.append('Delta (Complex - Receptor - Ligand):\n' + f'{col_name:16s}' +
                        '       Average     SD(Prop.)         SD   SEM(Prop.)        SEM\n' +
                        '-------------------------------------------------------------------------------')

        for key in self.data_keys:
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys:
                continue
            # # Skip chamber terms if we aren't using chamber prmtops
            # if not self.chamber and key in ['UB', 'IMP', 'CMAP']:
            #     continue
            # if self.com.INPUT['sander_apbs'] and key == 'EDISPER':
            #     continue
            # Catch special case of NMODEout classes
            if isinstance(self.com, NMODEout) and key == 'Total':
                printkey = '\n-TΔS binding ='
            else:
                printkey = key
            # Now print out the stats
            stdev = self[key].stdev()
            avg = self[key].mean()
            std = self[key].std()
            if not self.missing_terms:
                num_frames = len(self[key])
            else:
                num_frames = min(len(self.com[key]), len(self.rec[key]), len(self.lig[key]))
            semp = stdev / sqrt(num_frames)
            sem = std / sqrt(num_frames)

            if _output_format:
                text.append(['Δ' + printkey, avg, stdev, std, semp, sem])
            else:
                text.append(f"{'Δ' + printkey:16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}")

        if self.composite_keys:
            text.append('')
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'TOTAL':
                text.append('')
            stdev = self[key].stdev()
            avg = self[key].mean()
            std = self[key].std()
            if not self.missing_terms:
                num_frames = len(self[key])
            else:
                num_frames = min(len(self.com[key]), len(self.rec[key]), len(self.lig[key]))
            semp = stdev / sqrt(num_frames)
            sem = std / sqrt(num_frames)
                # num_frames is the same as the one from above
            if _output_format:
                text.append(['Δ' + key, avg, stdev, std, semp, sem])
            else:
                text.append(f"{'Δ' + key:16s} {avg:13.2f} {stdev:13.2f} {std:10.2f} {semp:12.2f} {sem:10.2f}")

        return text if _output_format else '\n'.join(text) + '\n'


class DecompOut(dict):
    """ Class for decomposition output file to collect statistics and output them """
    indicator = "                    PRINT DECOMP - TOTAL ENERGIES"
    descriptions = {'TDC': 'Total Energy Decomposition:',
                    'SDC': 'Sidechain Energy Decomposition:',
                    'BDC': 'Backbone Energy Decomposition:'}

    def __init__(self, mol: str, **kwargs):
        super(DecompOut, self).__init__(**kwargs)

        self.mut = None
        self.mol = mol
        self.resnums = None
        self.decfile = None
        self.num_terms = None
        self.allowed_tokens = tuple(['TDC'])
        self.verbose = None
        self.num_files = None
        self.resl = None
        self.basename = None
        self.csvwriter = None
        self.surften = None
        self.current_file = 0  # File counter
        self.get_next_term = self._get_next_term

    def parse_from_file(self, basename, resl, INPUT, surften, csvwriter, num_files=1, mut=False):
        self.basename = basename  # base name of output files
        self.resl = resl
        self.mut = mut

        self.num_files = num_files  # how many MPI files we created
        self.INPUT = INPUT
        self.verbose = INPUT['dec_verbose']
        self.surften = surften  # explicitly defined since is for GB and PB models
        # Set the term-extractor based on whether we want to dump the values to
        # a CSV file or just get the next term
        if csvwriter:
            self.get_next_term = self._get_next_term_csv

        if self.verbose in [1, 3]:
            self.allowed_tokens = tuple(['TDC', 'SDC', 'BDC'])

        # Create a separate csvwriter for each of the different token types,
        # and store them in a dictionary
        if csvwriter:
            self.csvwriter = {}
            for tok in self.allowed_tokens:
                self.csvwriter[tok] = writer(open(csvwriter + '.' + tok + '.csv', 'w'))
                self.csvwriter[tok].writerow([self.descriptions[tok]])
                self._write_header(self.csvwriter[tok])

        try:
            self.num_terms = int(self._get_num_terms())
        except TypeError:
            raise OutputError('DecompOut: Not a decomp output file')
        for token in self.allowed_tokens:
            self.resnums = [[0 for i in range(self.num_terms)], [0 for i in range(self.num_terms)]]
            self[token] = {}
        self.decfile = open(basename + '.0', 'r')
        self._fill_all_terms()
        self._fill_composite_terms()

    def parse_from_h5(self, d):
        for term in d:
            self[term] = {}
            for res in d[term]:
                self[term][res] = {}
                for res_e in d[term][res]:
                    if isinstance(d[term][res][res_e], np.ndarray):
                        self[term][res][res_e] = EnergyVector(d[term][res][res_e][()])
                    else:
                        self[term][res][res_e] = {}
                        for res2 in d[term][res][res_e]:
                            self[term][res][res_e][res2] = EnergyVector(d[term][res][res_e][res2][()])
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
            raise TypeError("{}.{} have 0 TDC starts".format(self.basename, 0))
        return num_terms

    def _get_next_term(self, expected_type, framenum=1):
        """ Gets the next energy term from the output file(s) """
        line = self.decfile.readline()
        if expected_type and expected_type not in self.allowed_tokens:
            raise OutputError('BUGBUG: expected_type must be in %s' % self.allowed_tokens)
        while line[:3] not in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1:
                    return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename, self.current_file), 'r')
                line = self.decfile.readline()
        # Return [res #, internal, vdw, eel, pol, sas]
        if expected_type and expected_type != line[:3]:
            raise OutputError(('Expecting %s type, but got %s type. Re-run ' +
                               'gmx_MMPBSA with the correct dec_verbose') % (expected_type, line[:3]))
        if self.mut and self.resl[int(line[4:10]) - 1].is_mutant():
            resnum = self.resl[int(line[4:10]) - 1].mutant_string
        else:
            resnum = self.resl[int(line[4:10]) - 1].string
        internal = float(line[11:20])
        vdw = float(line[21:30])
        eel = float(line[31:40])
        pol = float(line[41:50])
        sas = float(line[51:60]) * self.surften

        if resnum not in self[line[:3]]:
            self[line[:3]][resnum] = {
                'int': EnergyVector(),
                'vdw': EnergyVector(),
                'eel': EnergyVector(),
                'pol': EnergyVector(),
                'sas': EnergyVector()
            }

        self[line[:3]][resnum]['int'] = self[line[:3]][resnum]['int'].append(internal)
        self[line[:3]][resnum]['vdw'] = self[line[:3]][resnum]['vdw'].append(vdw)
        self[line[:3]][resnum]['eel'] = self[line[:3]][resnum]['eel'].append(eel)
        self[line[:3]][resnum]['pol'] = self[line[:3]][resnum]['pol'].append(pol)
        self[line[:3]][resnum]['sas'] = self[line[:3]][resnum]['sas'].append(sas)
        return [resnum, internal, vdw, eel, pol, sas]

    def _fill_all_terms(self):
        """
        This is for stability calculations -- just get all of the terms to
        fill up the data arrays.
        """
        token_counter = 0
        searched_type = self.allowed_tokens[0]
        framenum = self.INPUT['startframe']
        frames = 0
        com_token = self.get_next_term(self.allowed_tokens[0], framenum)
        while com_token:
            # Get all of the tokens
            for _ in range(1, self.num_terms):
                com_token = self.get_next_term(searched_type, framenum)
            token_counter += 1
            searched_type = self.allowed_tokens[token_counter % len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0:
                framenum += self.INPUT['interval']
                frames += 1
            com_token = self.get_next_term(searched_type, framenum)

        self.numframes = frames

    def _fill_composite_terms(self):
        for term in self:
            for res in self[term]:
                item = self[term][res][list(self[term][res].keys())[0]]
                # pair decomp scheme
                if isinstance(item, dict):
                    for res2 in self[term][res]:
                        tot = EnergyVector(len(self[term][res][res2]['int']))
                        for e in self[term][res][res2]:
                            tot = tot + self[term][res][res2][e]
                        self[term][res][res2]['tot'] = tot
                else:
                    tot = EnergyVector(len(self[term][res]['int']))
                    for e in self[term][res]:
                        tot = tot + self[term][res][e]
                    self[term][res]['tot'] = tot

    def summary(self, output_format: str = 'ascii'):
        """ Writes the summary in ASCII format to and open output_file """

        _output_format = 0 if output_format == 'ascii' else 1

        text = []
        if _output_format:
            text.append([self.mol.capitalize() + ':'])
        else:
            text.append(self.mol.capitalize() + ':')

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
                sqrt_frames = sqrt(len(self[term][res]['int']))

                if _output_format:
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

    def _get_next_term_csv(self, expected_type, framenum=1):
        """ Gets the next term and prints data to csv file """
        mydat = self._get_next_term(expected_type)
        if mydat:
            self.csvwriter[expected_type].writerow([framenum] + mydat + [sum(mydat[-5:])])
        else:
            self.csvwriter[expected_type].writerow([])
        return mydat

    def _write_header(self, csvwriter):
        """ Writes a table header to a passed CSV file """
        csvwriter.writerow(['Frame #', 'Residue', 'Internal', 'van der Waals',
                            'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])


class PairDecompOut(DecompOut):
    """ Same as DecompOut, but for Pairwise decomposition """
    indicator = "                    PRINT PAIR DECOMP - TOTAL ENERGIES"

    def _get_next_term(self, expected_type=None, framenum=1):
        """ Gets the next energy term from the output file(s) """
        line = self.decfile.readline()
        if expected_type and expected_type not in self.allowed_tokens:
            raise OutputError('BUGBUG: expected_type must be in %s' %
                              self.allowed_tokens)
        while line[:3] not in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1:
                    return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename, self.current_file), 'r')
                line = self.decfile.readline()
        # Return [res #, internal, vdw, eel, pol, sas]
        if expected_type and expected_type != line[:3]:
            raise OutputError(('Expecting %s type, but got %s type. Re-run ' +
                               'gmx_MMPBSA with the correct dec_verbose') % (expected_type, line[:3]))

        if self.mut and self.resl[int(line[4:11]) - 1].is_mutant():
            resnum = self.resl[int(line[4:11]) - 1].mutant_string
        else:
            resnum = self.resl[int(line[4:11]) - 1].string
        if self.mut and self.resl[int(line[13:20]) - 1].is_mutant():
            resnum2 = self.resl[int(line[13:20]) - 1].mutant_string
        else:
            resnum2 = self.resl[int(line[13:20]) - 1].string

        internal = float(line[21:33])
        vdw = float(line[34:46])
        eel = float(line[47:59])
        pol = float(line[60:72])
        sas = float(line[73:85]) * self.surften

        if resnum not in self[line[:3]]:
            self[line[:3]][resnum] = {}

        if resnum2 not in self[line[:3]][resnum]:
            self[line[:3]][resnum][resnum2] = {
                'int': EnergyVector(),
                'vdw': EnergyVector(),
                'eel': EnergyVector(),
                'pol': EnergyVector(),
                'sas': EnergyVector()}

        self[line[:3]][resnum][resnum2]['int'] = self[line[:3]][resnum][resnum2]['int'].append(internal)
        self[line[:3]][resnum][resnum2]['vdw'] = self[line[:3]][resnum][resnum2]['vdw'].append(vdw)
        self[line[:3]][resnum][resnum2]['eel'] = self[line[:3]][resnum][resnum2]['eel'].append(eel)
        self[line[:3]][resnum][resnum2]['pol'] = self[line[:3]][resnum][resnum2]['pol'].append(pol)
        self[line[:3]][resnum][resnum2]['sas'] = self[line[:3]][resnum][resnum2]['sas'].append(sas)
        return [resnum, resnum2, internal, vdw, eel, pol, sas]

    def summary(self, output_format: str = 'ascii'):
        """ Writes the summary in ASCII format to and open output_file """
        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[self.mol.capitalize() + ':']])
        else:
            text.extend([self.mol.capitalize() + ':'])
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

    def _write_header(self, csvwriter):
        """ Writes a table header to the csvwriter """
        csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])


class DecompBinding(dict):
    """ Class for decomposition binding (per-residue) """

    def __init__(self, com, rec, lig, INPUT, csvwriter, desc, **kwargs):
        """
        output should be an open file and csvfile should be a csv.writer class. If
        the output format is specified as csv, then output should be a csv.writer
        class as well.
        """
        super(DecompBinding, self).__init__(**kwargs)
        from csv import writer
        self.com, self.rec, self.lig = com, rec, lig
        self.num_terms = self.com.num_terms
        self.numframes = 0  # frame counter
        self.desc = desc  # Description
        self.INPUT = INPUT
        self.idecomp = INPUT['idecomp']
        self.verbose = INPUT['dec_verbose']
        # Check to see if output is a csv.writer or if it's a file. The invoked
        # method, "parse_all", is set based on whether we're doing a csv output
        # or an ascii output
        # if type(output).__name__ == 'writer':  # yuck... better way?
        #     self.parse_all = self._parse_all_csv
        # else:
        #     self.parse_all = self._parse_all_ascii
        # Set up the data for the DELTAs
        if self.verbose in [1, 3]:
            self.allowed_tokens = tuple(['TDC', 'SDC', 'BDC'])
        else:
            self.allowed_tokens = tuple(['TDC'])
        # Open up a separate CSV writer for all of the allowed tokens
        if csvwriter:
            self.csvwriter = {}
            for tok in self.allowed_tokens:
                self.csvwriter[tok] = writer(open(csvwriter + '.' + tok + '.csv', 'w'))
                self.csvwriter[tok].writerow([DecompOut.descriptions[tok]])
                self._write_header(self.csvwriter[tok])
        else:
            self.csvwriter = None
        for token in self.allowed_tokens:
            self[token] = {}
            self.resnums = {}

    def _write_header(self, csvwriter):
        """ Writes the header to the CSV file (legend at top of chart) """
        csvwriter.writerow(['Frame #', 'Residue', 'Location', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

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

            if self.csvwriter:
                f = self.com.INPUT['startframe']
                for i in range(self.com.numframes):
                    for res in self.com[term]:
                        self.csvwriter[term].writerow([f, res] + [self[term][res][x][i] for x in self[term][res]])
                    f += self.com.INPUT['interval']

    def summary(self, output_format: str = 'ascii'):
        # Parse everything
        self._parse_all_begin()

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[idecompString[self.idecomp]], [self.desc], []])
        else:
            text.append(idecompString[self.idecomp])
            text.append(self.desc)
        if self.verbose > 1:
            if _output_format:
                text.extend(self.com.summary(output_format))
                text.extend(self.com.summary(output_format))
                text.extend(self.com.summary(output_format))
            else:
                text.append(self.com.summary())
                text.append(self.com.summary())
                text.append(self.com.summary())
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

    def _write_header(self, csvwriter):
        """ Writes the header to the CSV file (legend at top of chart) """
        csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

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

            if self.csvwriter:
                f = self.com.INPUT['startframe']
                for i in range(self.com.numframes):
                    for res in self.com[term]:
                        for res2 in self.com[term][res]:
                            self.csvwriter[term].writerow([f, res, res2] +
                                                          [self[term][res][res2][x][i] for x in self[term][res][res2]])
                    f += self.com.INPUT['interval']

    def summary(self, output_format: str = 'ascii'):
        # Parse everything
        self._parse_all_begin()

        _output_format = 0 if output_format == 'ascii' else 1
        text = []
        if _output_format:
            text.extend([[idecompString[self.idecomp]], [self.desc], []])
        else:
            text.append(idecompString[self.idecomp])
            text.append(self.desc)
            text.append('')
        if self.verbose > 1:
            if _output_format:
                text.extend(self.com.summary(output_format))
                text.extend(self.rec.summary(output_format))
                text.extend(self.lig.summary(output_format))
            else:
                text.append(self.com.summary())
                text.append(self.rec.summary())
                text.append(self.lig.summary())
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


class H5Output:
    def __init__(self, fname):
        self.h5f = h5py.File(fname, 'r')
        self.app_namespace = SimpleNamespace(INPUT={}, FILES=SimpleNamespace(), INFO={})
        self.calc_types = SimpleNamespace(normal={}, mutant={}, decomp_normal={}, decomp_mutant={})

        for key in self.h5f:
            if key in ['normal', 'mutant']:
                self._h52e(key)
            elif key in ['decomp_normal', 'decomp_mutant']:
                self._h52decomp(key)
            elif key in ['INFO', 'INPUT', 'FILES']:
                self._h52app_namespace(key)
        self.h5f.close()

    def _h52app_namespace(self, key):
        for x in self.h5f[key]:
            tvar = self.h5f[key][x][()]
            if isinstance(tvar, bytes):
                cvar = tvar.decode()
            elif isinstance(tvar, np.float):
                cvar = None if np.isnan(tvar) else tvar
            elif isinstance(tvar, np.ndarray):
                cvar = [x.decode() if isinstance(x, bytes) else x for x in tvar if isinstance(x, bytes)]
            else:
                cvar = tvar
            if key == 'INPUT':
                self.app_namespace.INPUT[x] = cvar
            elif key == 'FILES':
                setattr(self.app_namespace.FILES, x, cvar)
            else:
                self.app_namespace.INFO[x] = cvar

    def _h52e(self, key):

        GBClass = QMMMout if self.app_namespace.INPUT['ifqnt'] else GBout
        # Determine which kind of RISM output class we are based on std/gf and
        # polardecomp
        if self.app_namespace.INPUT['polardecomp']:
            RISM_GF = PolarRISM_gf_Out
            RISM_Std = PolarRISM_std_Out
        else:
            RISM_GF = RISM_gf_Out
            RISM_Std = RISM_std_Out
        energy_outkeys = {'nmode': NMODEout, 'gb': GBClass, 'pb': PBout, 'rism std': RISM_Std, 'rism gf': RISM_GF}
        ent_outkeys = {'ie': IEout, 'c2': C2out}
        # key: normal or mutant
        calc_types = getattr(self.calc_types, key)
        # key  Energy: [gb, pb, rism std, rism gf], Decomp: [gb, pb], Entropy: [nmode, qh, ie, c2]
        for key1 in self.h5f[key]:
            # if key in ['gb', 'pb', 'rism std', 'rism gf', 'nmode', 'qh', 'ie', 'c2']:
            if key1 == 'qh':
                calc_types[key1] = QHout()
                calc_types[key1].parse_from_h5(self.h5f[key][key1])
                continue
            elif key1 in ent_outkeys:
                calc_types[key1] = ent_outkeys[key1]()
                calc_types[key1].parse_from_h5(self.h5f[key][key1])
                continue
            else:
                calc_types[key1] = {}
            # key2 is complex, receptor, ligand, delta or model for ie and c2
            for key2 in self.h5f[key][key1]:
                if key1 not in energy_outkeys:
                    continue
                calc_types[key1][key2] = energy_outkeys[key1](key2, self.app_namespace.INPUT,
                                                              self.app_namespace.INFO['using_chamber'])
                calc_types[key1][key2].parse_from_h5(self.h5f[key][key1][key2])

    def _h52decomp(self, key):
        DecompClass = DecompOut if self.app_namespace.INPUT['idecomp'] in [1, 2] else PairDecompOut
        calc_types = getattr(self.calc_types, key)
        for key1 in self.h5f[key]:
            # model
            calc_types[key1] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in self.h5f[key][key1]:
                calc_types[key1][key2] = DecompClass(key2)
                calc_types[key1][key2].parse_from_h5(self.h5f[key][key1][key2])


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
