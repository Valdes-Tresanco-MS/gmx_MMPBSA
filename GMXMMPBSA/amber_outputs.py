"""
This module contains all of the classes and code to collect data and calculate
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

from math import sqrt
from GMXMMPBSA.exceptions import (OutputError, LengthError, DecompError, InternalError)
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
        if obj is None: return
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
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EPOL', 'ENPOL']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'UB': ['GGAS', 'TOTAL'],
                      'IMP': ['GGAS', 'TOTAL'], 'CMAP': ['GGAS', 'TOTAL'],
                      'VDWAALS': ['GGAS', 'TOTAL'], 'EEL': ['GGAS', 'TOTAL'],
                      '1-4 VDW': ['GGAS', 'TOTAL'], '1-4 EEL': ['GGAS', 'TOTAL'],
                      'EPOL': ['GSOLV', 'TOTAL'], 'ENPOL': ['GSOLV', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'UB': 2, 'IMP': 2, 'CMAP': 2,
                    'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EPOL': 1,
                    'ENPOL': 1}

    def __init__(self, INPUT, chamber=False, **kwargs):
        super(AmberOutput, self).__init__(**kwargs)
        self.INPUT = INPUT
        self.chamber = chamber
        self.basename = None
        self.num_files = None
        self.is_read = False

    def parse_from_file(self, basename, num_files=1):
        self.num_files = num_files
        self.basename = basename

        for key in self.data_keys:
            self[key] = EnergyVector()
        for key in self.composite_keys:
            self[key] = EnergyVector()
        AmberOutput._read(self)
        self.fill_composite_terms()

    def parse_from_dict(self, d: dict):
        for key in d:
            if key in ['']:
                continue
            self[key] = EnergyVector(d[key])
        self.is_read = True
        self.fill_composite_terms()

    def print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        print_keys = [key for key in self.data_keys]
        # Add on the composite keys
        print_keys += self.composite_keys

        # write the header
        csvwriter.writerow(['Frame #'] + print_keys)

        # write out each frame
        c = self.INPUT['startframe']
        for i in range(len(self[print_keys[0]])):
            csvwriter.writerow([c] + [round(self[key][i], 4) for key in print_keys])
            c += self.INPUT['interval']

    def print_summary_csv(self, csvwriter):
        """ Prints the summary in CSV format """
        # print the header
        csvwriter.writerow(['Energy Component', 'Average', 'Std. Dev.',
                            'Std. Err. of Mean'])

        for key in self.data_keys:
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Skip any terms that have zero as every single element (i.e. EDISPER)
            if self[key] == 0: continue
            stdev = self[key].stdev()
            avg = self[key].mean()
            csvwriter.writerow([key, avg, stdev, stdev / sqrt(len(self[key]))])

        for key in self.composite_keys:
            # Now print out the composite terms
            stdev = self[key].stdev()
            avg = self[key].mean()
            csvwriter.writerow([key, avg, stdev, stdev / sqrt(len(self[key]))])

    def print_summary(self, mol: str = None):
        """ Returns a formatted string that can be printed directly to the
            output file
        """
        if not self.is_read:
            raise OutputError('Cannot print summary before reading output files')

        ret_str = ''
        if mol:
            ret_str = mol.capitalize() + '\n'
        ret_str += 'Energy Component            Average              Std. Dev.   Std. Err. of Mean\n'
        ret_str += '-------------------------------------------------------------------------------\n'

        for key in self.data_keys:
            # Skip terms we don't want to print
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Skip any terms that have zero as every single element (i.e. EDISPER)
            if self[key] == 0: continue
            stdev = self[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key, self[key].mean(), stdev, stdev / sqrt(len(
                self[key])))

        ret_str += '\n'
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'TOTAL': ret_str += '\n'
            stdev = self[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key, self[key].mean(), stdev, stdev / sqrt(len(
                self[key])))

        return ret_str + '\n\n'

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

    def fill_composite_terms(self):
        """
        Fills in the composite terms WITHOUT adding in terms we're not printing.
        This should be called after the final verbosity level has been set (based
        on whether or not certain terms need to be added in)
        """

        for key in self.composite_keys:
            self[key] = EnergyVector(len(self['EEL']))

        for key in self.data_keys:
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']:
                continue
            for component in self.data_key_owner[key]:
                self[component] = self[key] + self[component]


class IEout(dict):
    """
    Interaction Entropy output
    """

    def __init__(self, **kwargs):
        super(IEout, self).__init__(**kwargs)

    def print_summary_csv(self, csvwriter):
        """ Output summary of quasi-harmonic results in CSV format """
        csvwriter.writerow([f"Iteration Entropy calculation from last {self['ieframes']} frames..."])
        csvwriter.writerow(['Iteration Entropy:', '{:.2f}'.format(self.avg())])

    def sum(self, other, model, key1, key2):
        """
        Takes the sum between 2 keys of 2 different BindingStatistics
        classes and returns the average and standard deviation of that diff.
        """
        if len(self[model][key1]) != len(other[key2]):
            return (self[model][key1].mean() + other[key2].mean(),
                    sqrt(self[model][key1].stdev() ** 2 + other[key2].stdev() ** 2))
        mydiff = self[model][key1] + other[key2]
        return mydiff.mean(), mydiff.stdev()

    # def print_vectors(self, csvwriter):
    #     """ Prints the energy vectors to a CSV file for easy viewing
    #         in spreadsheets
    #     """
    #     csvwriter.writerow(['Frame #', 'Interaction Entropy'])
    #     for f, d in zip(self['frames'], self['data']):
    #         csvwriter.writerow([f] + [d])
    #
    def print_summary(self):
        """ Formatted summary of Interaction Entropy results """

        ret_str = 'Model           σ(Int. Energy)      Average       Std. Dev.   Std. Err. of Mean\n'
        ret_str += '-------------------------------------------------------------------------------\n'
        for model in self:
            self[model]['iedata'] = EnergyVector(self[model]['iedata'])
            stdev = self[model]['iedata'].stdev()
            avg = self[model]['iedata'].mean()
            ret_str += '%-14s %10.3f %16.3f %15.3f %19.3f\n' % (model, self[model]['sigma'], avg, stdev,
                                                                stdev / sqrt(len(self[model]['iedata'])))
        return ret_str


class C2out(dict):
    """
    Interaction Entropy output
    """

    def __init__(self, **kwargs):
        super(C2out, self).__init__(**kwargs)

    def print_summary_csv(self, csvwriter):
        """ Output summary of C2 results in CSV format """
        csvwriter.writerow([f"C2 Entropy calculation from last {self['ieframes']} frames..."])
        csvwriter.writerow(['C2 Entropy:', '{:.2f}'.format(self.avg())])

    def sum(self, other, model, key1, key2):
        """
        Takes the sum between 2 keys of 2 different BindingStatistics
        classes and returns the average and standard deviation of that diff.
        """
        return self[model][key1] + other[key2].mean(), other[key2].stdev()

    # def print_vectors(self, csvwriter):
    #     """ Prints the energy vectors to a CSV file for easy viewing
    #         in spreadsheets
    #     """
    #     csvwriter.writerow(['Frame #', 'Interaction Entropy'])
    #     for f, d in zip(self['frames'], self['data']):
    #         csvwriter.writerow([f] + [d])
    #
    def print_summary(self):
        """ Formatted summary of C2 Entropy results """

        ret_str = 'Model           σ(Int. Energy)      Value         Std. Dev.   Conf. Interv. (95%)\n'
        ret_str += '-------------------------------------------------------------------------------\n'
        for model in self:
            ret_str += '%-14s %10.3f %15.3f %15.3f %13.3f-%5.3f\n' % (model, self[model]['sigma'],
                                                                      self[model]['c2data'],
                                                                      self[model]['c2_std'],
                                                                      self[model]['c2_ci'][0],
                                                                      self[model]['c2_ci'][1])
        return ret_str


class QHout(dict):
    """ Quasi-harmonic output file class. QH output files are strange so we won't
        derive from AmberOutput
    """

    def __init__(self, filename, temp=298.15, **kwargs):
        super(QHout, self).__init__(**kwargs)
        self.filename = filename
        self.temperature = temp
        self.stability = False
        self._read()

    def print_summary_csv(self, csvwriter):
        """ Output summary of quasi-harmonic results in CSV format """
        csvwriter.writerow(['System', 'Translational', 'Rotational', 'Vibrational',
                            'Total'])
        csvwriter.writerow(['Complex:', self.com[1], self.com[2], self.com[3], self.com[0]])
        if not self.stability:
            csvwriter.writerow(['Receptor:', self.rec[1], self.rec[2], self.rec[3], self.rec[0]])
            csvwriter.writerow(['Ligand:', self.lig[1], self.lig[2], self.lig[3], self.lig[0]])
            csvwriter.writerow(['Delta S:',
                                self.com[1] - self.rec[1] - self.lig[1],
                                self.com[2] - self.rec[2] - self.lig[2],
                                self.com[3] - self.rec[3] - self.lig[3],
                                self.com[0] - self.rec[0] - self.lig[0]])

    def print_summary(self):
        """ Formatted summary of quasi-harmonic results """
        ret_str = '           Translational      Rotational      Vibrational           Total\n'

        ret_str += 'Complex:   %13.4f %15.4f %16.4f %15.4f\n' % (self.com[1], self.com[2], self.com[3], self.com[0])

        if not self.stability:
            ret_str += 'Receptor:  %13.4f %15.4f %16.4f %15.4f\n' % (self.rec[1], self.rec[2], self.rec[3], self.rec[0])
            ret_str += 'Ligand:    %13.4f %15.4f %16.4f %15.4f\n' % (self.lig[1], self.lig[2], self.lig[3], self.lig[0])
            ret_str += '\nTΔS:   %13.4f %15.4f %16.4f %15.4f\n' % (
                self.com[1] - self.rec[1] - self.lig[1],
                self.com[2] - self.rec[2] - self.lig[2],
                self.com[3] - self.rec[3] - self.lig[3],
                self.com[0] - self.rec[0] - self.lig[0])

        return ret_str

    def total_avg(self):
        """ Returns the average of the total """
        return self.com[0] - self.rec[0] - self.lig[0]

    def _read(self):
        """ Parses the output files and fills the data arrays """
        output = open(self.filename, 'r')
        rawline = output.readline()
        self.com = EnergyVector(4)
        self.rec = EnergyVector(4)
        self.lig = EnergyVector(4)
        comdone = False  # if we've done the complex yet (filled in self.com)
        recdone = False  # if we've done the receptor yet (filled in self.rec)

        # Try to fill in all found entropy values. If we can only find 1 set,
        # we're doing stability calculations
        while rawline:
            if rawline[0:6] == " Total":
                if not comdone:
                    self.com[0] = (float(rawline.split()[3]) * self.temperature / 1000)
                    self.com[1] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.com[2] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.com[3] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    comdone = True
                elif not recdone:
                    self.rec[0] = (float(rawline.split()[3]) * self.temperature / 1000)
                    self.rec[1] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.rec[2] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.rec[3] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    recdone = True
                else:
                    self.lig[0] = (float(rawline.split()[3]) * self.temperature / 1000)
                    self.lig[1] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.lig[2] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    self.lig[3] = (float(output.readline().split()[3]) *
                                   self.temperature / 1000)
                    break
            rawline = output.readline()
        # end while rawline

        self.stability = not recdone

        output.close()


class NMODEout(dict):
    """ Normal mode entropy approximation output class """
    # Ordered list of keys in the data dictionary
    data_keys = ['Translational', 'Rotational', 'Vibrational', 'Total']

    # Other aspects of AmberOutputs, which are just blank arrays
    composite_keys = []
    data_key_owner = {}
    print_levels = {'Translational': 1, 'Rotational': 1, 'Vibrational': 1, 'Total': 1}

    def __init__(self, basename, INPUT, num_files=1, chamber=False, **kwargs):
        super(NMODEout, self).__init__(**kwargs)
        import warnings
        self.basename = basename

        for key in self.data_keys:
            self[key] = EnergyVector()

        if chamber:
            warnings.warn('nmode is incompatible with chamber topologies!')
        self.temp = INPUT['temperature']
        self.num_files = num_files
        self.is_read = False

        self._read()

    def print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        # print header
        csvwriter.writerow(['Frame #'] + self.data_keys)

        # print data
        for i in range(len(self[self.data_keys[0]])):
            csvwriter.writerow([i] + [self[key][i] for key in self.data_keys])

    def print_summary_csv(self, csvwriter):
        """ Writes summary in CSV format """
        csvwriter.writerow(['Entropy Term', 'Average', 'Std. Dev.',
                            'Std. Err. of the Mean'])

        for key in self.data_keys:
            stdev = self[key].stdev()
            avg = self[key].mean()
            csvwriter.writerow([key, avg, stdev, stdev / sqrt(len(self[key]))])

    def print_summary(self):
        """ Returns the formatted string of output summary """

        ret_str = ('Entropy Term                Average              ' +
                   'Std. Dev.   Std. Err. of Mean\n')
        ret_str += ('------------------------------------------------' +
                    '-------------------------------\n')
        for key in self.data_keys:
            stdev = self[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key, self[key].mean(), stdev, stdev / sqrt(len(self[key])))

        return ret_str + '\n'

    def _read(self):
        """ Internal reading function to populate the data arrays """

        if self.is_read: return None  # don't read through again

        # Loop through all filenames
        for fileno in range(self.num_files):
            output_file = open('%s.%d' % (self.basename, fileno), 'r')
            self._get_energies(output_file)
            output_file.close()

        self.is_read = True

    def _get_energies(self, outfile):
        """ Parses the energy terms from the output file. This will parse 1 line
            at a time in order to minimize the memory requirements (we should only
            have to store a single line at a time in addition to the arrays of
            data)
        """

        rawline = outfile.readline()

        while rawline:
            if rawline[0:35] == '   |---- Entropy not Calculated---|':
                sys.stderr.write('Not all frames minimized within tolerance')

            if rawline[0:6] == 'Total:':
                self['Total'] = self['Total'].append(float(rawline.split()[3]) * self.temp / 1000)
                self['Translational'] = self['Translational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)
                self['Rotational'] = self['Rotational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)
                self['Vibrational'] = self['Vibrational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)

            rawline = outfile.readline()

    def fill_composite_terms(self):
        """ No-op for this class """
        pass


class GBout(AmberOutput):
    """ Amber output class for normal generalized Born simulations """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'UB': ['GGAS', 'TOTAL'],
                      'IMP': ['GGAS', 'TOTAL'], 'CMAP': ['GGAS', 'TOTAL'],
                      'VDWAALS': ['GGAS', 'TOTAL'], 'EEL': ['GGAS', 'TOTAL'],
                      '1-4 VDW': ['GGAS', 'TOTAL'], '1-4 EEL': ['GGAS', 'TOTAL'],
                      'EGB': ['GSOLV', 'TOTAL'], 'ESURF': ['GSOLV', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'UB': 2, 'IMP': 2, 'CMAP': 2,
                    'VDWAALS': 1, 'EEL': 1, '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1,
                    'ESURF': 1}

    # Ordered list of keys in the data dictionary

    def __init__(self, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, INPUT, chamber, **kwargs)
        # self.surften = INPUT['surften']
        # self.surfoff = INPUT['surfoff']
        # if read:
        #     AmberOutput._read(self)

    def _get_energies(self, outfile):
        """ Parses the mdout files for the GB potential terms """
        rawline = outfile.readline()
        while rawline:
            if rawline[0:5] == ' BOND':
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
                else:
                    self['UB'] = self['UB'].append(0.0)
                    self['IMP'] = self['IMP'].append(0.0)
                    self['CMAP'] = self['CMAP'].append(0.0)
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
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EPB', 'ENPOLAR', 'EDISPER']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'UB': ['GGAS', 'TOTAL'],
                      'IMP': ['GGAS', 'TOTAL'], 'CMAP': ['GGAS', 'TOTAL'],
                      'VDWAALS': ['GGAS', 'TOTAL'], 'EEL': ['GGAS', 'TOTAL'],
                      '1-4 VDW': ['GGAS', 'TOTAL'], '1-4 EEL': ['GGAS', 'TOTAL'],
                      'EPB': ['GSOLV', 'TOTAL'], 'ENPOLAR': ['GSOLV', 'TOTAL'],
                      'EDISPER': ['GSOLV', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EPB': 1, 'ENPOLAR': 1, 'UB': 2,
                    'IMP': 2, 'CMAP': 2, 'EDISPER': 1}

    def __init__(self, INPUT, chamber=False, **kwargs):
        AmberOutput.__init__(self, INPUT, chamber, **kwargs)
        self.apbs = INPUT['sander_apbs']
        AmberOutput._read(self)
        # FIXME:
        if self.apbs: self.print_levels['EDISPER'] = 3  # never print this for APBS

    def _get_energies(self, outfile):
        """ Parses the energy values from the output files """
        rawline = outfile.readline()
        while rawline:
            if rawline[0:5] == ' BOND':
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
                else:
                    self['UB'] = self['UB'].append(0.0)
                    self['IMP'] = self['IMP'].append(0.0)
                    self['CMAP'] = self['CMAP'].append(0.0)
                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[5]))
                self['EPB'] = self['EPB'].append(float(words[8]))
                words = outfile.readline().split()
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[3]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                self['ENPOLAR'] = self['ENPOLAR'].append(float(words[2]))
                if not self.apbs:
                    self['EDISPER'] = self['EDISPER'].append(float(words[5]))
                else:
                    self['EDISPER'] = self['EDISPER'].append(0.0)
            rawline = outfile.readline()


class RISMout(AmberOutput):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW',
                 '1-4 EEL', 'ERISM']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'VDWAALS': ['GGAS', 'TOTAL'],
                      'EEL': ['GGAS', 'TOTAL'], '1-4 VDW': ['GGAS', 'TOTAL'],
                      '1-4 EEL': ['GGAS', 'TOTAL'], 'ERISM': ['GSOLV', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'ERISM': 1}

    def __init__(self, INPUT, chamber=False, solvtype=0, **kwargs):
        AmberOutput.__init__(self, INPUT, chamber)
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

    def __init__(self, INPUT, chamber=False):
        RISMout.__init__(self, INPUT, chamber, 0)


class RISM_gf_Out(RISMout):
    """ No polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, INPUT, chamber=False):
        RISMout.__init__(self, INPUT, chamber, 1)


class PolarRISMout(RISMout):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW',
                 '1-4 EEL', 'POLAR SOLV', 'APOLAR SOLV']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'VDWAALS': ['GGAS', 'TOTAL'],
                      'EEL': ['GGAS', 'TOTAL'], '1-4 VDW': ['GGAS', 'TOTAL'],
                      '1-4 EEL': ['GGAS', 'TOTAL'],
                      'POLAR SOLV': ['GSOLV', 'TOTAL'],
                      'APOLAR SOLV': ['GSOLV', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
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

            if re.match(r'(solute_epot|solutePotentialEnergy)',
                        rawline):
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

    def __init__(self, INPUT, chamber=False):
        RISMout.__init__(self, INPUT, chamber, 0)


class PolarRISM_gf_Out(PolarRISMout):
    """ Polar decomp RISM output file for Gaussian Fluctuation free energy """

    def __init__(self, INPUT, chamber=False):
        RISMout.__init__(self, INPUT, chamber, 1)


class QMMMout(GBout):
    """ Class for QM/MM GBSA output files """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'ESCF']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND': ['GGAS', 'TOTAL'], 'ANGLE': ['GGAS', 'TOTAL'],
                      'DIHED': ['GGAS', 'TOTAL'], 'UB': ['GGAS', 'TOTAL'],
                      'IMP': ['GGAS', 'TOTAL'], 'CMAP': ['GGAS', 'TOTAL'],
                      'VDWAALS': ['GGAS', 'TOTAL'], 'EEL': ['GGAS', 'TOTAL'],
                      '1-4 VDW': ['GGAS', 'TOTAL'], '1-4 EEL': ['GGAS', 'TOTAL'],
                      'EGB': ['GSOLV', 'TOTAL'], 'ESURF': ['GSOLV', 'TOTAL'],
                      'ESCF': ['TOTAL']}
    # Which of those keys are composite
    composite_keys = ['GGAS', 'GSOLV', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND': 2, 'ANGLE': 2, 'DIHED': 2, 'VDWAALS': 1, 'EEL': 1,
                    '1-4 VDW': 2, '1-4 EEL': 2, 'EGB': 1, 'ESURF': 1, 'ESCF': 1,
                    'UB': 2, 'IMP': 2, 'CMAP': 2}

    def _get_energies(self, outfile):
        """ Parses the energies from a QM/MM output file. NOTE, however, that a
            QMMMout *could* just be a GBout with ESCF==0 if the QM region lies
            entirely outside this system
        """

        rawline = outfile.readline()

        while rawline:

            if rawline[0:5] == ' BOND':
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
                else:
                    self['UB'] = self['UB'].append(0.0)
                    self['IMP'] = self['IMP'].append(0.0)
                    self['CMAP'] = self['CMAP'].append(0.0)

                self['VDWAALS'] = self['VDWAALS'].append(float(words[2]))
                self['EEL'] = self['EEL'].append(float(words[5]))
                self['EGB'] = self['EGB'].append(float(words[8]))
                words = outfile.readline().split()
                self['1-4 VDW'] = self['1-4 VDW'].append(float(words[3]))
                self['1-4 EEL'] = self['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                # This is where ESCF will be. Since ESCF can differ based on which
                # qmtheory was chosen, we just check to see if it's != ESURF:
                if words[0] != 'minimization':
                    # It's possible that there is no space between ***ESCF and the =.
                    # If not, the ESCF variable will be the second, not the 3rd word
                    if words[0].endswith('='):
                        self['ESCF'] = self['ESCF'].append(float(words[1]))
                    else:
                        self['ESCF'] = self['ESCF'].append(float(words[2]))
                else:
                    self['ESCF'] = self['ESCF'].append(0.0)

            rawline = outfile.readline()

            # end if rawline[0:5] == ' BOND':

        # end while rawline


class BindingStatistics(dict):
    """ Base class for compiling the binding statistics """

    def __init__(self, com, rec, lig, chamber=False, **kwargs):
        super(BindingStatistics, self).__init__(**kwargs)
        self.com = com
        self.rec = rec
        self.lig = lig
        self.chamber = chamber
        self.inconsistent = False
        self.missing_terms = False

        if type(self.com) != type(self.rec) or type(self.com) != type(self.lig):
            raise TypeError('Binding statistics requires identical types')

        self.data_keys = self.com.data_keys
        self.composite_keys = []

        for key in self.com.composite_keys:
            self.composite_keys.append('DELTA ' + key)
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
        TINY = 0.005
        # First thing we do is check to make sure that all of the terms that
        # should *not* be printed actually cancel out (i.e. bonded terms)
        # for key in self.com.print_levels:
        #     if self.com.print_levels[key] > 1:
        #         diff = self.com[key] - self.rec[key] - self.lig[key]
        #         if diff.abs_gt(SingleTrajBinding.TINY):
        #             self.inconsistent = True
        #             # Now we have to print out everything
        #             self.com.verbose = 2
        #             self.rec.verbose = 2
        #             self.lig.verbose = 2
        #             break

        # FIXME: ya lo hice?
        # self.com.fill_composite_terms()
        # self.rec.fill_composite_terms()
        # self.lig.fill_composite_terms()

        for key in self.com.data_keys:
            self[key] = self.com[key] - self.rec[key] - self.lig[key]
        for key in self.com.composite_keys:
            self['DELTA ' + key] = self.com[key] - self.rec[key] - self.lig[key]

    def print_vectors(self, csvwriter):
        """ Output all of the energy terms including the differences if we're
            doing a single trajectory simulation and there are no missing terms
        """
        csvwriter.writerow(['Complex Energy Terms'])
        self.com.print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow(['Receptor Energy Terms'])
        self.rec.print_vectors(csvwriter)
        csvwriter.writerow([])
        csvwriter.writerow(['Ligand Energy Terms'])
        self.lig.print_vectors(csvwriter)
        csvwriter.writerow([])

        csvwriter.writerow(['Delta Energy Terms'])
        print_keys = [key for key in self.data_keys]
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

    def print_summary_csv(self, csvwriter):
        """ Prints formatted output in CSV format """
        if self.inconsistent:
            csvwriter.writerow(['WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL ' +
                                'POTENTIAL TERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY ' +
                                'QUESTIONABLE'])

        self.csvwriter.writerow(['Complex:'])
        self.com.print_summary_csv(csvwriter)
        self.csvwriter.writerow(['Receptor:'])
        self.rec.print_summary_csv(csvwriter)
        self.csvwriter.writerow(['Ligand:'])
        self.lig.print_summary_csv(csvwriter)

        csvwriter.writerow(['Differences (Complex - Receptor - Ligand):'])
        csvwriter.writerow(['Energy Component', 'Average', 'Std. Dev.', 'Std. Err. of Mean'])
        # Set verbose level. If verbose==0, that means we don't print com/rec/lig
        # but we print the differences as though verbose==1

        for key in self.data_keys:
            # Skip terms we don't want to print
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Catch special case of NMODEout classes
            if isinstance(self.com, NMODEout) and key == 'Total':
                printkey = '\nTΔS binding ='
            else:
                printkey = key
            # Now print out the stats
            if not self.missing_terms:
                stdev = self[key].stdev()
                avg = self[key].mean()
                num_frames = len(self[key])
            else:
                stdev = self[key][1]
                avg = self[key][0]
                num_frames = min(len(self.com[key]), len(self.rec[key]),
                                 len(self.lig[key]))
            csvwriter.writerow([printkey, avg, stdev, stdev / sqrt(num_frames)])

        for key in self.composite_keys:
            # Now print out the composite terms
            if not self.missing_terms:
                stdev = self[key].stdev()
                avg = self[key].mean()
                num_frames = len(self[key])
            else:
                stdev = self[key][1]
                avg = self[key][0]
                # num_frames is the same as the one from above
            csvwriter.writerow([key, avg, stdev, stdev / sqrt(len(self[key]))])

    def print_summary(self):
        """ Returns a string printing the summary of the binding statistics """

        if self.inconsistent:
            ret_str = ('WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL' +
                       '\nTERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE\n')
        else:
            ret_str = ''

        ret_str += 'Complex:\n' + self.com.print_summary()
        ret_str += 'Receptor:\n' + self.rec.print_summary()
        ret_str += 'Ligand:\n' + self.lig.print_summary()

        if isinstance(self.com, NMODEout):
            col_name = '%-16s' % 'Entropy Term'
        else:
            col_name = '%-16s' % 'Energy Component'
        ret_str += 'Differences (Complex - Receptor - Ligand):\n'
        ret_str += col_name + '            Average              Std. Dev.   Std. Err. of Mean\n'
        ret_str += '-------------------------------------------------------------------------------\n'

        for key in self.data_keys:
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys:
                continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']:
                continue
            # Catch special case of NMODEout classes
            if isinstance(self.com, NMODEout) and key == 'Total':
                printkey = '\nTΔS binding ='
            else:
                printkey = key
            # Now print out the stats
            stdev = self[key].stdev()
            avg = self[key].mean()
            if not self.missing_terms:
                num_frames = len(self[key])
            else:
                num_frames = min(len(self.com[key]), len(self.rec[key]), len(self.lig[key]))
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (printkey, avg, stdev, stdev / sqrt(num_frames))

        if self.composite_keys:
            ret_str += '\n'
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'DELTA TOTAL':
                ret_str += '\n'
            stdev = self[key].stdev()
            avg = self[key].mean()
            if not self.missing_terms:
                num_frames = len(self[key])
            else:
                num_frames = min(len(self.com[key]), len(self.rec[key]), len(self.lig[key]))
                # num_frames is the same as the one from above
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key, avg, stdev, stdev / sqrt(num_frames))

        return ret_str + '\n\n'


class DecompOut(dict):
    " Class for decomposition output file to collect statistics and output them "
    indicator = "                    PRINT DECOMP - TOTAL ENERGIES"
    descriptions = {'TDC': 'Total Energy Decomposition:',
                    'SDC': 'Sidechain Energy Decomposition:',
                    'BDC': 'Backbone Energy Decomposition:'}

    def __init__(self, **kwargs):
        super(DecompOut, self).__init__(**kwargs)

        self.resnums = None
        self.decfile = None
        self.num_terms = None
        self.allowed_tokens = tuple(['TDC'])
        self.verbose = None
        self.num_files = None
        self.prmtop = None
        self.basename = None
        self.csvwriter = None
        self.current_file = 0  # File counter
        self.get_next_term = self._get_next_term

    def parse_from_file(self, basename, prmtop, INPUT, csvwriter, num_files=1):
        self.basename = basename  # base name of output files
        self.prmtop = prmtop  # AmberParm prmtop object
        self.num_files = num_files  # how many MPI files we created
        self.INPUT = INPUT
        self.verbose = INPUT['dec_verbose']
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
        self.fill_composite_terms()

    def parse_from_dict(self, d: dict):
        for term, v in d.items():
            self[term] = {}
            for res, v1 in v.items():
                self[term][res] = {}
                for res_e, v2 in v1.items():
                    if isinstance(v2, dict):
                        self[term][res][res_e] = {_k: EnergyVector(_v) for _k, _v in v2.items()}
                    else:
                        self[term][res][res_e] = EnergyVector(v2)
        self.fill_composite_terms()

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
        while line[0:3] not in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1: return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename, self.current_file), 'r')
                line = self.decfile.readline()
        # Return [res #, internal, vdw, eel, pol, sas]
        if expected_type and expected_type != line[0:3]:
            raise OutputError(('Expecting %s type, but got %s type. Re-run ' +
                               'gmx_MMPBSA with the correct dec_verbose') % (expected_type, line[0:3]))
        resnum = int(line[4:10])
        internal = float(line[11:20])
        vdw = float(line[21:30])
        eel = float(line[31:40])
        pol = float(line[41:50])
        sas = float(line[51:60]) * self.INPUT['surften']

        if resnum not in self[line[0:3]]:
            self[line[0:3]][resnum] = {'int': EnergyVector(), 'vdw': EnergyVector(), 'eel': EnergyVector(),
                                       'pol': EnergyVector(), 'sas': EnergyVector()}

        self[line[0:3]][resnum]['int'] = self[line[0:3]][resnum]['int'].append(internal)
        self[line[0:3]][resnum]['vdw'] = self[line[0:3]][resnum]['vdw'].append(vdw)
        self[line[0:3]][resnum]['eel'] = self[line[0:3]][resnum]['eel'].append(eel)
        self[line[0:3]][resnum]['pol'] = self[line[0:3]][resnum]['pol'].append(pol)
        self[line[0:3]][resnum]['sas'] = self[line[0:3]][resnum]['sas'].append(sas)
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
            for i in range(1, self.num_terms):
                com_token = self.get_next_term(searched_type, framenum)
            token_counter += 1
            searched_type = self.allowed_tokens[token_counter % len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0:
                framenum += self.INPUT['interval']
                frames += 1
            com_token = self.get_next_term(searched_type, framenum)

        self.numframes = frames

    def fill_composite_terms(self):
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

    def print_summary(self):
        """ Writes the summary in ASCII format to and open output_file """
        ref_str = ''
        for term in self:
            ref_str += self.descriptions[term] + '\n'
            ref_str += ('Residue |       Internal      |    ' +
                        'van der Waals    |    Electrostatic    |   Polar Solvation  ' +
                        ' |   Non-Polar Solv.   |       TOTAL\n')
            ref_str += ('-------------------------------------------' +
                        '--------------------------------------------------------------' +
                        '----------------------------------\n')
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
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res]
                res_str = '%3s%4d' % (resnm, res)
                ref_str += (f"{str(res_str):7s} "
                            f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                            f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                            f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                            f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                            f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                            f"|{tot_avg:9.3f} +/- {tot_std:6.3f}\n")

        print(ref_str)

    def _get_next_term_csv(self, expected_type, framenum=1):
        """ Gets the next term and prints data to csv file """
        mydat = self._get_next_term(expected_type)
        if mydat:
            resnm = self.prmtop.parm_data['RESIDUE_LABEL'][mydat[0] - 1]
            res_str = '%3s%4d' % (resnm, mydat[0])
            if len(mydat) == 7:
                resnm2 = self.prmtop.parm_data['RESIDUE_LABEL'][mydat[1] - 1]
                res_str2 = '%3s%4d' % (resnm2, mydat[1])
                self.csvwriter[expected_type].writerow([framenum] + [res_str, res_str2] + mydat[2:] + [sum(mydat[-5:])])
            else:
                self.csvwriter[expected_type].writerow([framenum] + [res_str] + mydat[1:] + [sum(mydat[-5:])])
        else:
            self.csvwriter[expected_type].writerow([])
        return mydat

    def _write_header(self, csvwriter):
        """ Writes a table header to a passed CSV file """
        csvwriter.writerow(['Frame #', 'Residue', 'Internal', 'van der Waals',
                            'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    def write_summary(self, output_file):
        """ Writes the summary in ASCII format to and open output_file """
        for term in self:
            output_file.writeline(self.descriptions[term])
            output_file.writeline('Residue |       Internal      |    ' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation  ' +
                                  ' |   Non-Polar Solv.   |       TOTAL')
            output_file.writeline('-------------------------------------------' +
                                  '--------------------------------------------------------------' +
                                  '----------------------------------')
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

                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res]
                res_str = '%3s%4d' % (resnm, res)
                output_file.writeline(f"{res_str} "
                                      f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                      f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                      f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                      f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                      f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                      f"|{tot_avg:9.3f} +/- {tot_std:6.3f}\n")
            output_file.writeline('')

    def write_summary_csv(self, csvwriter):
        """ Writes the summary to a CSV file """
        for term in self.allowed_tokens:
            csvwriter.writerow([self.descriptions[term], ])
            csvwriter.writerow(['Residue', 'Internal', '', '', 'van der Waals', '',
                                '', 'Electrostatic', '', '', 'Polar Solvation', '',
                                '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''])
            csvwriter.writerow([''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6)
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

                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res - 1]
                res_str = '%3s%4d' % (resnm, res)
                sqrt_frames = sqrt(len(self[term][res]['int']))

                csvwriter.writerow([res_str, int_avg, int_std, int_std / sqrt_frames,
                                    vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                    eel_avg, eel_std, eel_std / sqrt_frames,
                                    pol_avg, pol_std, pol_std / sqrt_frames,
                                    sas_avg, sas_std, sas_std / sqrt_frames,
                                    tot_avg, tot_std, tot_std / sqrt_frames])
        csvwriter.writerow([])


class PairDecompOut(DecompOut):
    """ Same as DecompOut, but for Pairwise decomposition """
    indicator = "                    PRINT PAIR DECOMP - TOTAL ENERGIES"

    def _get_next_term(self, expected_type=None, framenum=1):
        """ Gets the next energy term from the output file(s) """
        line = self.decfile.readline()
        if expected_type and not expected_type in self.allowed_tokens:
            raise OutputError('BUGBUG: expected_type must be in %s' %
                              self.allowed_tokens)
        while not line[0:3] in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1: return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename, self.current_file), 'r')
                line = self.decfile.readline()
        # Return [res #, internal, vdw, eel, pol, sas]
        if expected_type and expected_type != line[0:3]:
            raise OutputError(('Expecting %s type, but got %s type. Re-run ' +
                               'gmx_MMPBSA with the correct dec_verbose') % (expected_type, line[0:3]))
        resnum = int(line[4:11])
        resnum2 = int(line[13:20])
        internal = float(line[21:33])
        vdw = float(line[34:46])
        eel = float(line[47:59])
        pol = float(line[60:72])
        sas = float(line[73:85]) * self.INPUT['surften']

        if resnum not in self[line[0:3]]:
            self[line[0:3]][resnum] = {}

        if resnum2 not in self[line[:3]][resnum]:
            self[line[:3]][resnum][resnum2] = {'int': EnergyVector(), 'vdw': EnergyVector(), 'eel': EnergyVector(),
                                               'pol': EnergyVector(), 'sas': EnergyVector()}

        self[line[0:3]][resnum][resnum2]['int'] = self[line[0:3]][resnum][resnum2]['int'].append(internal)
        self[line[0:3]][resnum][resnum2]['vdw'] = self[line[0:3]][resnum][resnum2]['vdw'].append(vdw)
        self[line[0:3]][resnum][resnum2]['eel'] = self[line[0:3]][resnum][resnum2]['eel'].append(eel)
        self[line[0:3]][resnum][resnum2]['pol'] = self[line[0:3]][resnum][resnum2]['pol'].append(pol)
        self[line[0:3]][resnum][resnum2]['sas'] = self[line[0:3]][resnum][resnum2]['sas'].append(sas)
        return [resnum, resnum2, internal, vdw, eel, pol, sas]

    def print_summary(self):
        """ Writes the summary in ASCII format to and open output_file """
        ref_str = ''
        for term in self:
            ref_str += self.descriptions[term] + '\n'
            ref_str += ('Resid 1 | Resid 2 |       Internal      |    ' +
                        'van der Waals    |    Electrostatic    |   Polar Solvation   |  ' +
                        'Non-Polar Solv.   |       TOTAL\n')
            ref_str += ('------------------------------------------------------------------------------------------'
                        '-----------------------------------------------------------\n')
            for res in self[term]:
                for res2 in self[term][res]:
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
                    resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res]
                    res_str0 = '%3s%4d' % (resnm, res)
                    resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res2]
                    res_str1 = '%3s%4d' % (resnm, res2)
                    ref_str += (f"{res_str0} | {res_str1}"
                                f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                f"|{tot_avg:9.3f} +/- {tot_std:6.3f}\n")
        print(ref_str)

    def write_summary(self, output_file):
        """ Writes the summary in ASCII format to and open output_file """
        for term in self.allowed_tokens:
            output_file.writeline(self.descriptions[term])
            output_file.writeline('Resid 1 | Resid 2 |       Internal      |    ' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation   |  ' +
                                  'Non-Polar Solv.   |       TOTAL')
            output_file.writeline('---------------------------------------------' +
                                  '----------------------------------------------------------------' +
                                  '----------------------------------------')
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
                    resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res - 1]
                    res_str0 = '%3s%4d' % (resnm, res)
                    resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res2 - 1]
                    res_str1 = '%3s%4d' % (resnm, res2)
                    output_file.writeline(f"{res_str0} | {res_str1} "
                                          f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                          f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                          f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                          f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                          f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                          f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            output_file.writeline('')

    def _write_header(self, csvwriter):
        """ Writes a table header to the csvwriter """
        csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    def write_summary_csv(self, csvwriter):
        """ Writes the summary to a CSV file """
        for term in self:
            csvwriter.writerow([self.descriptions[term]])
            csvwriter.writerow(['Resid 1', 'Resid 2', 'Internal', '', '',
                                'van der Waals', '', '', 'Electrostatic', '', '',
                                'Polar Solvation', '', '', 'Non-Polar Solv.', '',
                                '', 'TOTAL', '', ''])
            csvwriter.writerow([''] * 2 + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6)
            for res in self[term]:
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res]
                res_str0 = '%3s%4d' % (resnm, res)
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

                    resnm = self.prmtop.parm_data['RESIDUE_LABEL'][res2]
                    res_str1 = '%3s%4d' % (resnm, res2)
                    sqrt_frames = sqrt(len(self[term][res][res2]['int']))
                    csvwriter.writerow([res_str0, res_str1,
                                        int_avg, int_std, int_std / sqrt_frames,
                                        vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                        eel_avg, eel_std, eel_std / sqrt_frames,
                                        pol_avg, pol_std, pol_std / sqrt_frames,
                                        sas_avg, sas_std, sas_std / sqrt_frames,
                                        tot_avg, tot_std, tot_std / sqrt_frames])
            csvwriter.writerow([])


class DecompBinding(dict):
    """ Class for decomposition binding (per-residue) """

    def __init__(self, com, rec, lig, prmtop_system, idecomp, verbose, output,
                 csvwriter, desc, **kwargs):
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
        self.output = output
        self.desc = desc  # Description
        self.prmtop_system = prmtop_system
        self.idecomp = idecomp
        self.verbose = verbose
        # Check to see if output is a csv.writer or if it's a file. The invoked
        # method, "parse_all", is set based on whether we're doing a csv output
        # or an ascii output
        if type(output).__name__ == 'writer':  # yuck... better way?
            self.parse_all = self._parse_all_csv
        else:
            self.parse_all = self._parse_all_ascii
        # Set up the data for the DELTAs
        if verbose in [1, 3]:
            self.allowed_tokens = tuple(['TDC', 'SDC', 'BDC'])
        else:
            self.allowed_tokens = tuple(['TDC'])
        # Open up a separate CSV writer for all of the allowed tokens
        if csvwriter:
            self.csvwriter = {}
            for tok in self.allowed_tokens:
                self.csvwriter[tok] = writer(open(csvwriter + '.' + tok + '.csv', 'w'))
                self.csvwriter[tok].writerow(['DELTA', DecompOut.descriptions[tok]])
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
                res_label = self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][res - 1]

                if self.prmtop_system.res_list[res].receptor_number:
                    rec_num = self.prmtop_system.res_list[res].receptor_number - 1
                    self.resnums[res] = ['%3s%4d' % (res_label, res), f"R {res_label:3s}{rec_num:4d}"]
                    other_token = self.rec[term][rec_num]
                else:
                    lig_num = self.prmtop_system.res_list[res].ligand_number - 1
                    self.resnums[res] = ['%3s%4d' % (res_label, res), f"L {res_label:3s}{lig_num:4d}"]
                    other_token = self.lig[term][lig_num]

                for e in self.com[term][res]:
                    self[term][res][e] = self.com[term][res][e] - other_token[e]

            if self.csvwriter:
                f = self.com.INPUT['startframe']
                for i in range(self.com.numframes):
                    for res in self.com[term]:
                        self.csvwriter[term].writerow([f] + self.resnums[res] +
                                                      [self[term][res][x][i] for x in self[term][res]])
                    f += self.com.INPUT['interval']
                # if self.csvwriter:
                #     self.csvwriter[term].writerow([framenum, res, self.resnums[res]] +
                #                                   [self[term][res][x].mean() for x in self[term][res]])

    def _parse_all_csv(self):
        """ Specifically parses everything and dumps avg results to a CSV file """
        # Parse everything
        self._parse_all_begin()
        # Now we have all of our DELTAs, time to print them to the CSV
        # Print the header
        self.output.writerow([idecompString[self.idecomp]])
        self.output.writerow([self.desc])
        if self.verbose > 1:
            self.output.writerow(['Complex:'])
            self.com.write_summary_csv(self.output)
            self.output.writerow(['Receptor:'])
            self.rec.write_summary_csv(self.output)
            self.output.writerow(['Ligand:'])
            self.lig.write_summary_csv(self.output)
        # Now write the DELTAs

        self.output.writerow(['DELTAS:'])
        for term in self:
            self.output.writerow([DecompOut.descriptions[term]])
            self.output.writerow(['Residue', 'Location', 'Internal', '', '',
                                  'van der Waals', '', '', 'Electrostatic', '', '',
                                  'Polar Solvation', '', '', 'Non-Polar Solv.',
                                  '', '', 'TOTAL', '', ''])
            self.output.writerow(['', ''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6)
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

                resnm = self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][res - 1]
                res_str = '%3s%4d' % (resnm, res)
                sqrt_frames = sqrt(len(self[term][res]['int']))

                self.output.writerow(self.resnums[res] +
                                     [int_avg, int_std, int_std / sqrt_frames,
                                      vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                      eel_avg, eel_std, eel_std / sqrt_frames,
                                      pol_avg, pol_std, pol_std / sqrt_frames,
                                      sas_avg, sas_std, sas_std / sqrt_frames,
                                      tot_avg, tot_std, tot_std / sqrt_frames])
            self.output.writerow([])

    def _parse_all_ascii(self):
        """ Parses all output files and prints to ASCII output format """
        # Parse everything
        self._parse_all_begin()
        # Now we have all of our DELTAs, time to print them to the CSV
        # Print the header
        self.output.writeline(idecompString[self.idecomp])
        self.output.writeline(self.desc)
        self.output.writeline('')
        if self.verbose > 1:
            self.output.writeline('')
            self.output.writeline('Complex:')
            self.com.write_summary(self.output)
            self.output.writeline('Receptor:')
            self.rec.write_summary(self.output)
            self.output.writeline('Ligand:')
            self.lig.write_summary(self.output)
        # Now write the DELTAs

        self.output.writeline('DELTAS:')

        for term in self:
            self.output.writeline(DecompOut.descriptions[term])
            self.output.writeline('Residue |  Location |       Internal      |    ' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation   |' +
                                  '    Non-Polar Solv.  |       TOTAL')
            self.output.writeline('---------------------------------------------' +
                                  '--------------------------------------------------------------' +
                                  '--------------------------------------------')
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
                self.output.writeline(f"{self.resnums[res][0]} | {self.resnums[res][1]} "
                                      f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                      f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                      f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                      f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                      f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                      f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            self.output.writeline('')


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
                res_label = self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][res - 1]
                for res2 in self.com[term][res]:
                    self[term][res][res2] = {}
                    res2_label = self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][res2 - 1]
                    if self.prmtop_system.res_list[res].receptor_number:
                        res_num = self.prmtop_system.res_list[res].receptor_number - 1
                        res_name = f"R {res_label:3s}{res_num:4d}"
                        if self.prmtop_system.res_list[res2].receptor_number:
                            rec_num = self.prmtop_system.res_list[res2].receptor_number - 1
                            res_name2 = f"R {res2_label:3s}{rec_num:4d}"
                            # Both residues are in the receptor -- pull the next one
                            other_token = self.rec[term][res_num][rec_num]
                        else:
                            lig_num = self.prmtop_system.res_list[res2].ligand_number - 1
                            res_name2 = f"L {res2_label:3s}{lig_num:4d}"
                            other_token = {}
                    else:
                        res_num = self.prmtop_system.res_list[res].ligand_number - 1
                        res_name = f"L {res_label:3s}{res_num:4d}"
                        if self.prmtop_system.res_list[res2].ligand_number:
                            lig_num = self.prmtop_system.res_list[res2].ligand_number - 1
                            res_name2 = f"L {res2_label:3s}{lig_num:4d}"
                            # Both residues are in the ligand -- pull the next one
                            other_token = self.lig[term][res_num][lig_num]
                        else:
                            rec_num = self.prmtop_system.res_list[res2].receptor_number - 1
                            res_name2 = f"R {res2_label:3s}{rec_num:4d}"
                            other_token = {}
                    self.resnums[(res, res2)] = (res_name, res_name2)
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
                            self.csvwriter[term].writerow([f, self.resnums[(res, res2)][0],
                                                           self.resnums[(res, res2)][1]] +
                                                          [self[term][res][res2][x][i] for x in self[term][res][res2]])
                    f += self.com.INPUT['interval']

    def _parse_all_csv(self):
        """ Specifically parses everything and dumps avg results to a CSV file """
        # Parse everything
        self._parse_all_begin()
        # Now we have all of our DELTAs, time to print them to the CSV
        # Print the header
        self.output.writerow([idecompString[self.idecomp]])
        self.output.writerow([self.desc])
        self.output.writerow([])
        if self.verbose > 1:
            self.output.writerow(['Complex:'])
            self.com.write_summary_csv(self.output)
            self.output.writerow(['Receptor:'])
            self.rec.write_summary_csv(self.output)
            self.output.writerow(['Ligand:'])
            self.lig.write_summary_csv(self.output)
        # Now write the DELTAs

        self.output.writerow(['DELTAS:'])
        for term in self:
            self.output.writerow([DecompOut.descriptions[term]])
            self.output.writerow(['Resid 1', 'Resid 2', 'Internal', '', '',
                                  'van der Waals', '', '', 'Electrostatic', '', '',
                                  'Polar Solvation', '', '', 'Non-Polar Solv.',
                                  '', '', 'TOTAL', '', ''])
            self.output.writerow(['', ''] + ['Avg.', 'Std. Dev.', 'Std. Err. of Mean'] * 6)
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
                    self.output.writerow([self.resnums[(res, res2)][0], self.resnums[(res, res2)][1],
                                          int_avg, int_std, int_std / sqrt_frames,
                                          vdw_avg, vdw_std, vdw_std / sqrt_frames,
                                          eel_avg, eel_std, eel_std / sqrt_frames,
                                          pol_avg, pol_std, pol_std / sqrt_frames,
                                          sas_avg, sas_std, sas_std / sqrt_frames,
                                          tot_avg, tot_std, tot_std / sqrt_frames])
            self.output.writerow([])

    def _parse_all_ascii(self):
        """ Parses all output files and prints to ASCII output format """
        # Parse everything
        self._parse_all_begin()
        # Now we have all of our DELTAs, time to print them to the CSV
        # Print the header
        self.output.writeline(idecompString[self.idecomp])
        self.output.writeline(self.desc)
        self.output.writeline('')
        if self.verbose > 1:
            self.output.writeline('')
            self.output.writeline('Complex:')
            self.com.write_summary(self.output)
            self.output.writeline('Receptor:')
            self.rec.write_summary(self.output)
            self.output.writeline('Ligand:')
            self.lig.write_summary(self.output)
        # Now write the DELTAs

        self.output.writeline('DELTAS:')
        for term in self:
            self.output.writeline(DecompOut.descriptions[term])
            self.output.writeline('Resid 1   | Resid 2   |       Internal      |    van der Waals    |    '
                                  'Electrostatic    |   Polar Solvation   |    Non-Polar Solv.  |       TOTAL')
            self.output.writeline('-------------------------------------------------------------------------------'
                                  '--------------------------------------------------------------------------')
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
                    print(res, res2, self.resnums[(res, res2)])
                    self.output.writeline(f"{self.resnums[(res, res2)][0]} | {self.resnums[(res, res2)][1]} "
                                          f"|{int_avg:9.3f} +/- {int_std:6.3f} "
                                          f"|{vdw_avg:9.3f} +/- {vdw_std:6.3f} "
                                          f"|{eel_avg:9.3f} +/- {eel_std:6.3f} "
                                          f"|{pol_avg:9.3f} +/- {pol_std:6.3f} "
                                          f"|{sas_avg:9.3f} +/- {sas_std:6.3f} "
                                          f"|{tot_avg:9.3f} +/- {tot_std:6.3f}")
            self.output.writeline('')


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
        # key: normal or mutant
        calc_types = getattr(self.calc_types, key)
        # key  Energy: [gb, pb, rism std, rism gf], Decomp: [gb, pb], Entropy: [nmode, qh, ie, c2]
        for key1 in self.h5f[key]:
            # if key in ['gb', 'pb', 'rism std', 'rism gf', 'nmode', 'qh', 'ie', 'c2']:
            calc_types[key1] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in self.h5f[key][key1]:
                calc_types[key1][key2] = {}
                # Energetic terms
                for key3 in self.h5f[key][key1][key2]:
                    calc_types[key1][key2][key3] = self.h5f[key][key1][key2][key3][()]

    def _h52decomp(self, key):
        calc_types = getattr(self.calc_types, key)
        for key1 in self.h5f[key]:
            # model
            calc_types[key1] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in self.h5f[key][key1]:
                calc_types[key1][key2] = {}
                # TDC, SDC, BDC
                for key3 in self.h5f[key][key1][key2]:
                    # residue first level
                    for key4 in self.h5f[key][key1][key2][key3]:
                        for key5 in self.h5f[key][key1][key2][key3][key4]:
                            if isinstance(self.h5f[key][key1][key2][key3][key4], h5py.Group):
                                # residue sec level
                                for key6 in self.h5f[key][key1][key2][key3][key4][key5]:
                                    calc_types[key][key2][(key3, key4, key5, key6)] = \
                                        self.h5f[key][key1][key2][key3][key4][key5][key6][()]
                            else:
                                # energy terms
                                for key5 in self.h5f[key][key1][key2][key3][key4]:
                                    calc_types[key][key2][(key3, key4, key5)] = \
                                        self.h5f[key][key1][key2][key3][key4][key5][()]


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
