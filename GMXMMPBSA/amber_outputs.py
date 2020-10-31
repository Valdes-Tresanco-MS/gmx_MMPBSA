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
from GMXMMPBSA.exceptions import (OutputError, LengthError, DecompError,
                                  InternalError)
import sys

idecompString = ['idecomp = 0: No decomposition analysis',
                 'idecomp = 1: Per-residue decomp adding 1-4 interactions to Internal.',
                 'idecomp = 2: Per-residue decomp adding 1-4 interactions to EEL and VDW.',
                 'idecomp = 3: Pairwise decomp adding 1-4 interactions to Internal.',
                 'idecomp = 4: Pairwise decomp adding 1-4 interactions to EEL and VDW.']

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

def _std_dev(sum_squares, running_sum, num):
    """ Returns a propagated result for the std. dev. using sum of squares """
    return sqrt(abs(sum_squares/num - (running_sum/num) * (running_sum/num)))

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class EnergyVector(list):
    """ list-derived class for holding energies that defines add and subtract
        methods to add and subtract different vectors term by term
    """

    #==================================================

    def __add__(self, other):
        # If other is a scalar, do scalar addition
        if isinstance(other, float) or isinstance(other, int):
            return EnergyVector([x + other for x in self])
        # Otherwise, vector addition
        if len(self) != len(other):
            raise LengthError('length mismatch in energy vectors')
        return EnergyVector([self[i] + other[i] for i in range(len(self))])

    #==================================================

    def __sub__(self, other):
        # If other is a scalar, do scalar subtraction
        if isinstance(other, float) or isinstance(other, int):
            return EnergyVector([x - other for x in self])
        # Otherwise, vector subtraction
        if len(self) != len(other):
            raise LengthError('length mismatch in energy vectors')
        return EnergyVector([self[i] - other[i] for i in range(len(self))])

    #==================================================

    def __mul__(self, scalar):
        """
        Return a copy of the vector whose every element is multiplied by the
        scalar
        """
        if not isinstance(scalar, float) and not isinstance(scalar, int):
            raise TypeError('Cannot multiply a vector by a non-numeric scalar')
        return EnergyVector([x * scalar for x in self])

    #==================================================

    def __imul__(self, scalar):
        """ In-place multiplication -- no allocation of another EnergyVector """
        if not isinstance(scalar, float) and not isinstance(scalar, int):
            raise TypeError('Cannot multiply a vector by a non-numeric scalar')
        for i, x in enumerate(self):
            self[i] = x * scalar

    #==================================================

    def __iadd__(self, other):
        """ In-place addition -- no allocation of another EnergyVector """
        # Scalar addition?
        if isinstance(other, int) or isinstance(other, float):
            for i in range(len(self)):
                self[i] += other
            return
        # Otherwise, vector addition
        if len(self) != len(other):
            raise LengthError('length mismatch in energy vectors')
        for i in range(len(self)):
            self[i] += other[i]

    #==================================================

    def __isub__(self, other):
        """ In-place subtraction -- no allocation of another EnergyVector """
        # Scalar subtraction?
        if isinstance(other, int) or isinstance(other, float):
            for i in range(len(self)):
                self[i] -= other
            return
        # Otherwise, vector subtraction
        if len(self) != len(other):
            raise LengthError('length mismatch in energy vectors')
        for i in range(len(self)):
            self[i] -= other[i]

    #==================================================

    def __lt__(self, other):
        """
        This allows us to compare a list term-by-term, or we can compare every
        term to a single integer or floating point number
        """
        try:
            for val in self:
                if val > float(other): return False
            return True
        except TypeError:
            if len(self) != len(other):
                return self.avg() < other.avg()
            for i, s in enumerate(self):
                if s > other[i]: return False
            return True

        raise InternalError('Should not be here!')

    #==================================================

    def __eq__(self, other):
        """
        This allows a term-by-term comparison or compare every term to a single
        integer/float.
        """
        try:
            for val in self:
                if val != float(other): return False
            return True
        except TypeError:
            if len(self) != len(other):
                return self.avg() == other.avg()
            else:
                for i, s in enumerate(self):
                    if s != other[i]: return False
                return True

        raise InternalError('Should not be here!')

    #==================================================

    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)

    #==================================================

    def __gt__(self, other):
        """
        Opposite of __lt__
        """
        return not self.__lt__(other) and not self.__eq__(other)

    #==================================================

    def __ge__(self, other):
        return not self.__lt__(other)

    #==================================================

    def avg(self):
        return (sum(self) / len(self))

    #==================================================

    def stdev(self):
        rsum2 = 0
        for num in self:
            rsum2 += num * num
        rsum2 /= len(self)
        avg = self.avg()
        return sqrt(abs(rsum2 - avg * avg))

    #==================================================

    def abs_gt(self, val):
        """ If any element's absolute value is greater than a # """
        for myval in self:
            if abs(myval) > val: return True
        return False

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class AmberOutput(object):
    """
    Base Amber output class. It takes a basename as a file name and parses
    through all of the thread-specific output files (assumed to have the suffix
    .# where # spans from 0 to num_files - 1
    """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EPOL', 'ENPOL']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'UB':['G gas', 'TOTAL'],
                      'IMP':['G gas', 'TOTAL'], 'CMAP':['G gas', 'TOTAL'],
                      'VDWAALS':['G gas', 'TOTAL'], 'EEL':['G gas', 'TOTAL'],
                      '1-4 VDW':['G gas', 'TOTAL'], '1-4 EEL':['G gas', 'TOTAL'],
                      'EPOL':['G solv', 'TOTAL'], 'ENPOL':['G solv', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'UB':2, 'IMP':2, 'CMAP':2,
                    'VDWAALS':1, 'EEL':1, '1-4 VDW':2, '1-4 EEL':2, 'EPOL':1,
                    'ENPOL':1}

    #==================================================

    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        self.verbose = INPUT['verbose']
        self.num_files = num_files
        self.basename = basename
        self.chamber = chamber
        self.data = {}
        for key in self.data_keys:
            self.data[key] = EnergyVector()
        for key in self.composite_keys:
            self.data[key] = EnergyVector()

        self.is_read = False

    #==================================================

    def print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        # Determine which keys we want to print
        print_keys = []
        for key in self.data_keys:
            if self.print_levels[key] > self.verbose: continue
            print_keys.append(key)
        # Add on the composite keys
        print_keys += self.composite_keys

        # write the header
        csvwriter.writerow(['Frame #'] + print_keys)

        # write out each frame
        for i in range(len(self.data[print_keys[0]])):
            csvwriter.writerow([i] + [self.data[key][i] for key in print_keys])

    #==================================================

    def print_summary_csv(self, csvwriter):
        """ Prints the summary in CSV format """
        # print the header
        csvwriter.writerow(['Energy Component','Average','Std. Dev.',
                            'Std. Err. of Mean'])

        for key in self.data_keys:
            # Skip terms we don't want to print
            if self.verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Skip any terms that have zero as every single element (i.e. EDISPER)
            if self.data[key] == 0: continue
            stdev = self.data[key].stdev()
            avg = self.data[key].avg()
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(len(self.data[key]))])

        for key in self.composite_keys:
            # Now print out the composite terms
            stdev = self.data[key].stdev()
            avg = self.data[key].avg()
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(len(self.data[key]))])

    #==================================================

    def print_summary(self):
        """ Returns a formatted string that can be printed directly to the
            output file
        """
        if not self.is_read:
            raise OutputError('Cannot print summary before reading output files')

        ret_str = ('Energy Component            Average              ' +
                   'Std. Dev.   Std. Err. of Mean\n')
        ret_str += ('------------------------------------------------' +
                    '-------------------------------\n')

        for key in self.data_keys:
            # Skip terms we don't want to print
            if self.verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Skip any terms that have zero as every single element (i.e. EDISPER)
            if self.data[key] == 0: continue
            stdev = self.data[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key,
                                                         self.data[key].avg(), stdev, stdev / sqrt(len(self.data[key])))

        ret_str += '\n'
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'TOTAL': ret_str += '\n'
            stdev = self.data[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key,
                                                         self.data[key].avg(), stdev, stdev / sqrt(len(self.data[key])))

        return ret_str + '\n\n'

    #==================================================

    def _read(self):
        """
        Internal reading function. This should be called at the end of __init__.
        It loops through all of the output files to populate the arrays
        """
        if self.is_read: return None # don't read through them twice

        for fileno in range(self.num_files):
            output_file = open('%s.%d' % (self.basename, fileno), 'r')
            self._get_energies(output_file)
            output_file.close()
            # If we have to get energies elsewhere (e.g., with GB and ESURF), do
            # that here. This is an empty function when unnecessary
            self._extra_reading(fileno)

        self.is_read = True

    #==================================================

    def _extra_reading(self, fileno):
        pass

    #==================================================

    def fill_composite_terms(self):
        """
        Fills in the composite terms WITHOUT adding in terms we're not printing.
        This should be called after the final verbosity level has been set (based
        on whether or not certain terms need to be added in)
        """

        for key in self.composite_keys:
            self.data[key]=EnergyVector([0 for i in range(len(self.data['EEL']))])

        for key in self.data_keys:
            if self.print_levels[key] > self.verbose: continue
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            for component in self.data_key_owner[key]:
                self.data[component] = self.data[key] + self.data[component]

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class QHout(object):
    """ Quasi-harmonic output file class. QH output files are strange so we won't
        derive from AmberOutput
    """

    #==================================================

    def __init__(self, filename, temp=298.15):
        self.filename = filename
        self.temperature = temp
        self.stability = False
        self._read()

    #==================================================

    def print_summary_csv(self, csvwriter):
        """ Output summary of quasi-harmonic results in CSV format """
        csvwriter.writerow(['System','Translational','Rotational','Vibrational',
                            'Total'])
        csvwriter.writerow(['Complex:',self.com[1],self.com[2],self.com[3],
                            self.com[0]])
        if not self.stability:
            csvwriter.writerow(['Receptor:',self.rec[1],self.rec[2],self.rec[3],
                                self.rec[0]])
            csvwriter.writerow(['Ligand:',self.lig[1],self.lig[2],self.lig[3],
                                self.lig[0]])
            csvwriter.writerow(['Delta S:',
                                self.com[1] - self.rec[1] - self.lig[1],
                                self.com[2] - self.rec[2] - self.lig[2],
                                self.com[3] - self.rec[3] - self.lig[3],
                                self.com[0] - self.rec[0] - self.lig[0] ])

    #==================================================

    def print_summary(self):
        """ Formatted summary of quasi-harmonic results """
        ret_str =  ('           Translational      Rotational      ' +
                    'Vibrational           Total\n')

        ret_str += 'Complex:   %13.4f %15.4f %16.4f %15.4f\n' % (self.com[1],
                                                                 self.com[2], self.com[3], self.com[0])

        if not self.stability:
            ret_str += 'Receptor:  %13.4f %15.4f %16.4f %15.4f\n' % (self.rec[1],
                                                                     self.rec[2], self.rec[3], self.rec[0])
            ret_str += 'Ligand:    %13.4f %15.4f %16.4f %15.4f\n' % (self.lig[1],
                                                                     self.lig[2], self.lig[3], self.lig[0])
            ret_str += '\nDELTA S:   %13.4f %15.4f %16.4f %15.4f\n' % (
                self.com[1] - self.rec[1] - self.lig[1],
                self.com[2] - self.rec[2] - self.lig[2],
                self.com[3] - self.rec[3] - self.lig[3],
                self.com[0] - self.rec[0] - self.lig[0] )

        return ret_str

    #==================================================

    def total_avg(self):
        """ Returns the average of the total """
        return self.com[0] - self.rec[0] - self.lig[0]

    #==================================================

    def _read(self):
        """ Parses the output files and fills the data arrays """
        output = open(self.filename, 'r')
        rawline = output.readline()
        self.com = EnergyVector([0,0,0,0])
        self.rec = EnergyVector([0,0,0,0])
        self.lig = EnergyVector([0,0,0,0])
        comdone = False # if we've done the complex yet (filled in self.com)
        recdone = False # if we've done the receptor yet (filled in self.rec)

        # Try to fill in all found entropy values. If we can only find 1 set,
        # we're doing stability calculations
        while rawline:
            if rawline[0:6] == " Total":
                if not comdone:
                    self.com[0] = (float(rawline.split()[3]) * self.temperature/1000)
                    self.com[1] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.com[2] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.com[3] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    comdone = True
                elif not recdone:
                    self.rec[0] = (float(rawline.split()[3]) * self.temperature/1000)
                    self.rec[1] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.rec[2] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.rec[3] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    recdone = True
                else:
                    self.lig[0] = (float(rawline.split()[3]) * self.temperature/1000)
                    self.lig[1] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.lig[2] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    self.lig[3] = (float(output.readline().split()[3]) *
                                   self.temperature/1000)
                    break
            rawline = output.readline()
        # end while rawline

        self.stability = not recdone

        output.close()

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class NMODEout(object):
    """ Normal mode entropy approximation output class """
    # Ordered list of keys in the data dictionary
    data_keys = ['Translational', 'Rotational', 'Vibrational', 'Total']

    # Other aspects of AmberOutputs, which are just blank arrays
    composite_keys = []
    data_key_owner = {}
    print_levels = {'Translational':1,'Rotational':1,'Vibrational':1,'Total':1}

    #==================================================

    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        import warnings
        self.basename = basename

        self.data = {}
        for key in self.data_keys:
            self.data[key] = EnergyVector()

        if chamber:
            warnings.warn('nmode is incompatible with chamber topologies!')
        self.temp = INPUT['temp']
        self.num_files = num_files
        self.is_read = False

        self._read()

    #==================================================

    def print_vectors(self, csvwriter):
        """ Prints the energy vectors to a CSV file for easy viewing
            in spreadsheets
        """
        # print header
        csvwriter.writerow(['Frame #'] + self.data_keys)

        # print data
        for i in range(len(self.data[self.data_keys[0]])):
            csvwriter.writerow([i] + [self.data[key][i] for key in self.data_keys])

    #==================================================

    def print_summary_csv(self, csvwriter):
        """ Writes summary in CSV format """
        csvwriter.writerow(['Entropy Term','Average','Std. Dev.',
                            'Std. Err. of the Mean'])

        for key in self.data_keys:
            stdev = self.data[key].stdev()
            avg = self.data[key].avg()
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(len(self.data[key]))])

    #==================================================

    def print_summary(self):
        """ Returns the formatted string of output summary """

        ret_str = ('Entropy Term                Average              ' +
                   'Std. Dev.   Std. Err. of Mean\n')
        ret_str += ('------------------------------------------------' +
                    '-------------------------------\n')
        for key in self.data_keys:
            stdev = self.data[key].stdev()
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key,
                                                         self.data[key].avg(), stdev, stdev/sqrt(len(self.data[key])))

        return ret_str + '\n'

    #==================================================

    def _read(self):
        """ Internal reading function to populate the data arrays """

        if self.is_read: return None # don't read through again

        # Loop through all filenames
        for fileno in range(self.num_files):
            output_file = open('%s.%d' % (self.basename, fileno), 'r')
            self._get_energies(output_file)
            output_file.close()

        self.is_read = True

    #==================================================

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
                self.data['Total'].append(float(rawline.split()[3]) *
                                          self.temp / 1000)
                self.data['Translational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)
                self.data['Rotational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)
                self.data['Vibrational'].append(
                    float(outfile.readline().split()[3]) * self.temp / 1000)

            rawline = outfile.readline()

    #==================================================

    def fill_composite_terms(self):
        """ No-op for this class """
        pass

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class GBout(AmberOutput):
    """ Amber output class for normal generalized Born simulations """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'UB':['G gas', 'TOTAL'],
                      'IMP':['G gas', 'TOTAL'], 'CMAP':['G gas', 'TOTAL'],
                      'VDWAALS':['G gas', 'TOTAL'], 'EEL':['G gas', 'TOTAL'],
                      '1-4 VDW':['G gas', 'TOTAL'], '1-4 EEL':['G gas', 'TOTAL'],
                      'EGB':['G solv', 'TOTAL'], 'ESURF':['G solv', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'UB':2, 'IMP':2, 'CMAP':2,
                    'VDWAALS':1, 'EEL':1, '1-4 VDW':2, '1-4 EEL':2, 'EGB':1,
                    'ESURF':1}
    # Ordered list of keys in the data dictionary

    #==================================================

    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        AmberOutput.__init__(self, basename, INPUT, num_files, chamber)
        self.surften = INPUT['surften']
        self.surfoff = INPUT['surfoff']
        AmberOutput._read(self)

    #==================================================

    def _get_energies(self, outfile):
        """ Parses the mdout files for the GB potential terms """
        rawline = outfile.readline()

        while rawline:

            if rawline[0:5] == ' BOND':
                words = rawline.split()
                self.data['BOND'].append(float(words[2]))
                self.data['ANGLE'].append(float(words[5]))
                self.data['DIHED'].append(float(words[8]))
                words = outfile.readline().split()

                if self.chamber:
                    self.data['UB'].append(float(words[2]))
                    self.data['IMP'].append(float(words[5]))
                    self.data['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()
                else:
                    self.data['UB'].append(0.0)
                    self.data['IMP'].append(0.0)
                    self.data['CMAP'].append(0.0)

                self.data['VDWAALS'].append(float(words[2]))
                self.data['EEL'].append(float(words[5]))
                self.data['EGB'].append(float(words[8]))
                words = outfile.readline().split()
                self.data['1-4 VDW'].append(float(words[3]))
                self.data['1-4 EEL'].append(float(words[7]))
            # end if rawline[0:5] == 'BOND'

            rawline = outfile.readline()

        # end while rawline

    #==================================================

    def _extra_reading(self, fileno):
        # Load the ESURF data from the cpptraj output
        fname = '%s.%d' % (self.basename, fileno)
        fname = fname.replace('gb.mdout','gb_surf.dat')
        surf_data = _get_cpptraj_surf(fname)
        self.data['ESURF'].extend((surf_data * self.surften) + self.surfoff)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class PBout(AmberOutput):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EPB', 'ENPOLAR', 'EDISPER']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'UB':['G gas', 'TOTAL'],
                      'IMP':['G gas', 'TOTAL'], 'CMAP':['G gas', 'TOTAL'],
                      'VDWAALS':['G gas', 'TOTAL'], 'EEL':['G gas', 'TOTAL'],
                      '1-4 VDW':['G gas', 'TOTAL'], '1-4 EEL':['G gas', 'TOTAL'],
                      'EPB':['G solv', 'TOTAL'], 'ENPOLAR':['G solv', 'TOTAL'],
                      'EDISPER':['G solv', 'TOTAL'] }
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'VDWAALS':1, 'EEL':1,
                    '1-4 VDW':2, '1-4 EEL':2, 'EPB':1, 'ENPOLAR':1, 'UB':2,
                    'IMP':2,'CMAP':2, 'EDISPER':1}

    #==================================================

    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        AmberOutput.__init__(self, basename, INPUT, num_files, chamber)
        self.apbs = INPUT['sander_apbs']
        AmberOutput._read(self)
        if self.apbs: self.print_levels['EDISPER'] = 3 # never print this for APBS

    #==================================================

    def _get_energies(self, outfile):
        """ Parses the energy values from the output files """

        rawline = outfile.readline()

        while rawline:

            if rawline[0:5] == ' BOND':
                words = rawline.split()
                self.data['BOND'].append(float(words[2]))
                self.data['ANGLE'].append(float(words[5]))
                self.data['DIHED'].append(float(words[8]))
                words = outfile.readline().split()

                if self.chamber:
                    self.data['UB'].append(float(words[2]))
                    self.data['IMP'].append(float(words[5]))
                    self.data['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()
                else:
                    self.data['UB'].append(0.0)
                    self.data['IMP'].append(0.0)
                    self.data['CMAP'].append(0.0)

                self.data['VDWAALS'].append(float(words[2]))
                self.data['EEL'].append(float(words[5]))
                self.data['EPB'].append(float(words[8]))
                words = outfile.readline().split()
                self.data['1-4 VDW'].append(float(words[3]))
                self.data['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                self.data['ENPOLAR'].append(float(words[2]))
                if not self.apbs: self.data['EDISPER'].append(float(words[5]))
                else: self.data['EDISPER'].append(0.0)
            # end if rawline == ' BOND'

            rawline = outfile.readline()

        # end while rawline

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
import re
class RISMout(AmberOutput):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW',
                 '1-4 EEL', 'ERISM']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'VDWAALS':['G gas', 'TOTAL'],
                      'EEL':['G gas', 'TOTAL'], '1-4 VDW':['G gas', 'TOTAL'],
                      '1-4 EEL':['G gas', 'TOTAL'], 'ERISM':['G solv', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'VDWAALS':1, 'EEL':1,
                    '1-4 VDW':2, '1-4 EEL':2, 'ERISM':1}

    #==================================================

    def __init__(self, basename, INPUT, num_files=1, chamber=False, solvtype=0):
        AmberOutput.__init__(self, basename, INPUT, num_files, chamber)
        self.solvtype = solvtype
        AmberOutput._read(self)

    #==================================================

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
                self.data['VDWAALS'].append(float(words[2]))
                self.data['EEL'].append(float(words[3]))
                self.data['BOND'].append(float(words[4]))
                self.data['ANGLE'].append(float(words[5]))
                self.data['DIHED'].append(float(words[6]))
                self.data['1-4 VDW'].append(float(words[7]))
                self.data['1-4 EEL'].append(float(words[8]))

            elif self.solvtype == 0 and re.match(
                    r'(rism_exchem|rism_excessChemicalPotential)\s',rawline):
                self.data['ERISM'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(
                    r'(rism_exchGF|rism_excessChemicalPotentialGF)\s',rawline):
                self.data['ERISM'].append(float(rawline.split()[1]))

            rawline = outfile.readline()

    #==================================================

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class RISM_std_Out(RISMout):
    """ No polar decomp RISM output file for standard free energy """
    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        RISMout.__init__(self, basename, INPUT, num_files, chamber, 0)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class RISM_gf_Out(RISMout):
    """ No polar decomp RISM output file for Gaussian Fluctuation free energy """
    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        RISMout.__init__(self, basename, INPUT, num_files, chamber, 1)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class PolarRISMout(RISMout):
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW',
                 '1-4 EEL', 'POLAR SOLV', 'APOLAR SOLV']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'VDWAALS':['G gas', 'TOTAL'],
                      'EEL':['G gas', 'TOTAL'], '1-4 VDW':['G gas', 'TOTAL'],
                      '1-4 EEL':['G gas', 'TOTAL'],
                      'POLAR SOLV':['G solv', 'TOTAL'],
                      'APOLAR SOLV':['G solv', 'TOTAL']}
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # Which of those keys belong to the gas phase energy contributions
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'VDWAALS':1, 'EEL':1,
                    '1-4 VDW':2, '1-4 EEL':2, 'POLAR SOLV':1, 'APOLAR SOLV':1}

    #==================================================

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
                self.data['VDWAALS'].append(float(words[2]))
                self.data['EEL'].append(float(words[3]))
                self.data['BOND'].append(float(words[4]))
                self.data['ANGLE'].append(float(words[5]))
                self.data['DIHED'].append(float(words[6]))
                self.data['1-4 VDW'].append(float(words[8]))
                self.data['1-4 EEL'].append(float(words[8]))

            elif self.solvtype == 0 and re.match(
                    r'(rism_polar|rism_polarExcessChemicalPotential)\s',rawline):
                self.data['POLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 0 and re.match(
                    r'(rism_apolar|rism_apolarExcessChemicalPotential)\s',rawline):
                self.data['APOLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(
                    r'(rism_polGF|rism_polarExcessChemicalPotentialGF)\s',rawline):
                self.data['POLAR SOLV'].append(float(rawline.split()[1]))
            elif self.solvtype == 1 and re.match(
                    r'(rism_apolGF|rism_apolarExcessChemicalPotentialGF)\s',rawline):
                self.data['APOLAR SOLV'].append(float(rawline.split()[1]))

            rawline = outfile.readline()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PolarRISM_std_Out(PolarRISMout):
    """ Polar decomp RISM output file for standard free energy """
    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        RISMout.__init__(self, basename, INPUT, num_files, chamber, 0)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class PolarRISM_gf_Out(PolarRISMout):
    """ Polar decomp RISM output file for Gaussian Fluctuation free energy """
    def __init__(self, basename, INPUT, num_files=1, chamber=False):
        RISMout.__init__(self, basename, INPUT, num_files, chamber, 1)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

class QMMMout(GBout):
    """ Class for QM/MM GBSA output files """
    # Ordered list of keys in the data dictionary
    data_keys = ['BOND', 'ANGLE', 'DIHED', 'UB', 'IMP', 'CMAP', 'VDWAALS', 'EEL',
                 '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'ESCF']
    # Dictionary that maps each data key to their respective composite keys
    data_key_owner = {'BOND':['G gas', 'TOTAL'], 'ANGLE':['G gas', 'TOTAL'],
                      'DIHED':['G gas', 'TOTAL'], 'UB':['G gas', 'TOTAL'],
                      'IMP':['G gas', 'TOTAL'], 'CMAP':['G gas', 'TOTAL'],
                      'VDWAALS':['G gas', 'TOTAL'], 'EEL':['G gas', 'TOTAL'],
                      '1-4 VDW':['G gas', 'TOTAL'], '1-4 EEL':['G gas', 'TOTAL'],
                      'EGB':['G solv', 'TOTAL'], 'ESURF':['G solv', 'TOTAL'],
                      'ESCF':['TOTAL']}
    # Which of those keys are composite
    composite_keys = ['G gas', 'G solv', 'TOTAL']
    # What the value of verbosity must be to print out this data
    print_levels = {'BOND':2, 'ANGLE':2, 'DIHED':2, 'VDWAALS':1, 'EEL':1,
                    '1-4 VDW':2, '1-4 EEL':2, 'EGB':1, 'ESURF':1, 'ESCF':1,
                    'UB':2, 'IMP':2, 'CMAP':2}

    #==================================================

    def _get_energies(self, outfile):
        """ Parses the energies from a QM/MM output file. NOTE, however, that a
            QMMMout *could* just be a GBout with ESCF==0 if the QM region lies
            entirely outside this system
        """

        rawline = outfile.readline()

        while rawline:

            if rawline[0:5] == ' BOND':
                words = rawline.split()
                self.data['BOND'].append(float(words[2]))
                self.data['ANGLE'].append(float(words[5]))
                self.data['DIHED'].append(float(words[8]))
                words = outfile.readline().split()

                if self.chamber:
                    self.data['UB'].append(float(words[2]))
                    self.data['IMP'].append(float(words[5]))
                    self.data['CMAP'].append(float(words[8]))
                    words = outfile.readline().split()
                else:
                    self.data['UB'].append(0.0)
                    self.data['IMP'].append(0.0)
                    self.data['CMAP'].append(0.0)

                self.data['VDWAALS'].append(float(words[2]))
                self.data['EEL'].append(float(words[5]))
                self.data['EGB'].append(float(words[8]))
                words = outfile.readline().split()
                self.data['1-4 VDW'].append(float(words[3]))
                self.data['1-4 EEL'].append(float(words[7]))
                words = outfile.readline().split()
                # This is where ESCF will be. Since ESCF can differ based on which
                # qmtheory was chosen, we just check to see if it's != ESURF:
                if words[0] != 'minimization':
                    # It's possible that there is no space between ***ESCF and the =.
                    # If not, the ESCF variable will be the second, not the 3rd word
                    if words[0].endswith('='):
                        self.data['ESCF'].append(float(words[1]))
                    else:
                        self.data['ESCF'].append(float(words[2]))
                else:
                    self.data['ESCF'].append(0.0)

            rawline = outfile.readline()

            # end if rawline[0:5] == ' BOND':

        # end while rawline

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class BindingStatistics(object):
    """ Base class for compiling the binding statistics """

    #==================================================

    def __init__(self, com, rec, lig, verbose=1, chamber=False):
        self.com = com
        self.rec = rec
        self.lig = lig
        self.verbose = verbose
        self.chamber = chamber
        self.inconsistent = False
        self.missing_terms = False

        if type(self.com) != type(self.rec) or type(self.com) != type(self.lig):
            raise TypeError('Binding statistics requires identical types')

        self.data_keys = self.com.data_keys
        self.composite_keys = []
        self.data = {}
        for key in self.com.composite_keys:
            self.composite_keys.append('DELTA ' + key)
        self.print_levels = self.com.print_levels
        try:
            self.delta()
            self.missing_terms = False
        except LengthError:
            self.delta2()
            self.missing_terms = True

    #==================================================

    def delta(self):
        """
        Calculates the delta statistics. Should check for any consistencies that
        would cause verbosity levels to change, and it should change them
        accordingly in the child classes
        """
        pass

    #==================================================

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

    #==================================================

    def print_summary_csv(self, csvwriter):
        """ Prints formatted output in CSV format """
        if self.inconsistent:
            csvwriter.writerow(['WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL ' +
                                'POTENTIAL TERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY ' +
                                'QUESTIONABLE'])

        if self.verbose:
            self.csvwriter.writerow(['Complex:'])
            self.com.print_summary_csv(csvwriter)
            self.csvwriter.writerow(['Receptor:'])
            self.rec.print_summary_csv(csvwriter)
            self.csvwriter.writerow(['Ligand:'])
            self.lig.print_summary_csv(csvwriter)

        csvwriter.writerow(['Differences (Complex - Receptor - Ligand):'])
        csvwriter.writerow(['Energy Component','Average','Std. Dev.',
                            'Std. Err. of Mean'])
        # Set verbose level. If verbose==0, that means we don't print com/rec/lig
        # but we print the differences as though verbose==1
        verbose = self.verbose
        if self.verbose == 0: verbose = 1

        for key in self.data_keys:
            # Skip terms we don't want to print
            if verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Catch special case of NMODEout classes
            if isinstance(self.com, NMODEout) and key == 'Total':
                printkey = '\nDELTA S total='
            else:
                printkey = key
            # Now print out the stats
            if not self.missing_terms:
                stdev = self.data[key].stdev()
                avg = self.data[key].avg()
                num_frames = len(self.data[key])
            else:
                stdev = self.data[key][1]
                avg = self.data[key][0]
                num_frames = min(len(self.com.data[key]),len(self.rec.data[key]),
                                 len(self.lig.data[key]))
            csvwriter.writerow([printkey, avg, stdev, stdev/sqrt(num_frames)])

        for key in self.composite_keys:
            # Now print out the composite terms
            if not self.missing_terms:
                stdev = self.data[key].stdev()
                avg = self.data[key].avg()
                num_frames = len(self.data[key])
            else:
                stdev = self.data[key][1]
                avg = self.data[key][0]
                # num_frames is the same as the one from above
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(len(self.data[key]))])

    #==================================================

    def print_summary(self):
        """ Returns a string printing the summary of the binding statistics """

        if self.inconsistent:
            ret_str = ('WARNING: INCONSISTENCIES EXIST WITHIN INTERNAL POTENTIAL' +
                       '\nTERMS. THE VALIDITY OF THESE RESULTS ARE HIGHLY QUESTIONABLE\n')
        else:
            ret_str = ''

        if self.verbose:
            ret_str += 'Complex:\n' + self.com.print_summary()
            ret_str += 'Receptor:\n' + self.rec.print_summary()
            ret_str += 'Ligand:\n' + self.lig.print_summary()

        if isinstance(self.com, NMODEout):
            col_name = '%-16s' % 'Entropy Term'
        else:
            col_name = '%-16s' % 'Energy Component'
        ret_str += 'Differences (Complex - Receptor - Ligand):\n'
        ret_str += (col_name+'            Average              ' +
                    'Std. Dev.   Std. Err. of Mean\n')
        ret_str += ('------------------------------------------------' +
                    '-------------------------------\n')

        # verbose == 0 only suppresses Complex, receptor, and ligand printout
        verbose = self.verbose
        if self.verbose == 0: verbose = 1
        if self.inconsistent: verbose = 2

        for key in self.data_keys:
            # Skip terms we don't want to print
            if verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            # Catch special case of NMODEout classes
            if isinstance(self.com, NMODEout) and key == 'Total':
                printkey = '\nDELTA S total='
            else:
                printkey = key
            # Now print out the stats
            if not self.missing_terms:
                stdev = self.data[key].stdev()
                avg = self.data[key].avg()
                num_frames = len(self.data[key])
            else:
                stdev = self.data[key][1]
                avg = self.data[key][0]
                num_frames = min(len(self.com.data[key]),len(self.rec.data[key]),
                                 len(self.lig.data[key]))
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (printkey, avg, stdev,
                                                         stdev/sqrt(num_frames))

        if self.composite_keys: ret_str += '\n'
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'DELTA TOTAL': ret_str += '\n'
            if not self.missing_terms:
                stdev = self.data[key].stdev()
                avg = self.data[key].avg()
                num_frames = len(self.data[key])
            else:
                stdev = self.data[key][1]
                avg = self.data[key][0]
                # num_frames is the same as the one from above
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key, avg, stdev,
                                                         stdev/sqrt(num_frames))

        return ret_str + '\n\n'

    #==================================================

    def diff(self, other, key1, key2):
        """
        Takes the difference between 2 keys of 2 different BindingStatistics
        classes and returns the average and standard deviation of that diff.
        """
        if len(self.data[key1]) != len(other.data[key2]):
            return (self.data[key1].avg() - other.data[key2].avg(),
                    sqrt(self.data[key1].stdev()**2 + other.data[key2].stdev()**2))

        mydiff = self.data[key1] - other.data[key2]

        return mydiff.avg(), mydiff.stdev()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class SingleTrajBinding(BindingStatistics):
    """ Statistics calculated from a single trajectory binding calculation """

    # QM/MM runs (at least the test) have ~0.001 error in the 1-4s. This should
    # be a safe value.
    TINY = 0.005

    #==================================================

    def delta(self):
        """ Calculates the delta statistics """
        # First thing we do is check to make sure that all of the terms that
        # should *not* be printed actually cancel out (i.e. bonded terms)
        for key in self.com.print_levels:
            if self.com.print_levels[key] > 1:
                diff = self.com.data[key] - self.rec.data[key] - self.lig.data[key]
                if diff.abs_gt(SingleTrajBinding.TINY):
                    self.inconsistent = True
                    # Now we have to print out everything
                    self.com.verbose = 2
                    self.rec.verbose = 2
                    self.lig.verbose = 2
                    break

        self.com.fill_composite_terms()
        self.rec.fill_composite_terms()
        self.lig.fill_composite_terms()

        for key in self.com.data_keys:
            self.data[key] = self.com.data[key] - self.rec.data[key] - self.lig.data[key]
        for key in self.com.composite_keys:
            self.data['DELTA ' + key] = self.com.data[key] - self.rec.data[key] - self.lig.data[key]

    #==================================================

    def delta2(self):
        """
        In the off-chance that not every frame was calculated (normal mode calc
        in which not every frame was minimized, for instance), then treat this
        the same way we treat multi-traj binding
        """
        self.com.verbose = 2
        self.rec.verbose = 2
        self.lig.verbose = 2

        self.com.fill_composite_terms()
        self.rec.fill_composite_terms()
        self.lig.fill_composite_terms()

        for key in self.com.data_keys:
            self.data[key] = [self.com.data[key].avg() - self.rec.data[key].avg() -
                              self.lig.data[key].avg(),
                              sqrt(self.com.data[key].stdev() ** 2 +
                                   self.rec.data[key].stdev() ** 2 +
                                   self.lig.data[key].stdev() ** 2) ]

        for key in self.com.composite_keys:
            self.data['DELTA ' + key] = \
                [self.com.data[key].avg() - self.rec.data[key].avg() -
                 self.lig.data[key].avg(),
                 sqrt(self.com.data[key].stdev() ** 2 +
                      self.rec.data[key].stdev() ** 2 +
                      self.lig.data[key].stdev() ** 2) ]

    #==================================================

    def print_vectors(self, csvwriter):
        """ Single trajectory binding may be able to print DELTA terms also """
        BindingStatistics.print_vectors(self, csvwriter)
        if not self.missing_terms:
            csvwriter.writerow(['DELTA Energy Terms'])

            # Determine our print_keys
            print_keys = []
            for key in self.data_keys:
                if self.print_levels[key] > self.verbose: continue
                print_keys.append(key)
            print_keys += self.composite_keys

            # write the header
            csvwriter.writerow(['Frame #'] + print_keys)

            # write out each frame
            for i in range(len(self.data[print_keys[0]])):
                csvwriter.writerow([i]+[self.data[key][i] for key in print_keys])
            csvwriter.writerow([])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MultiTrajBinding(BindingStatistics):
    """ Statistics calculated from multiple trajectories binding calculation """

    #==================================================

    def delta(self):
        """ Calculates the delta statistics """
        # We have to print out *ALL* terms, since different trajectories mean that
        # internal potential terms don't cancel out
        self.com.verbose = 2
        self.rec.verbose = 2
        self.lig.verbose = 2

        self.com.fill_composite_terms()
        self.rec.fill_composite_terms()
        self.lig.fill_composite_terms()

        for key in self.com.data_keys:
            self.data[key] = self.com.data[key] - self.rec.data[key] - self.lig.data[key]
            # self.data[key] = [self.com.data[key].avg() - self.rec.data[key].avg() -
            #                   self.lig.data[key].avg(),
            #                   sqrt(self.com.data[key].stdev() ** 2 +
            #                        self.rec.data[key].stdev() ** 2 +
            #                        self.lig.data[key].stdev() ** 2) ]
        for key in self.com.composite_keys:
            self.data['DELTA ' + key] = self.com.data[key] - self.rec.data[key]- self.lig.data[key]
            # self.data['DELTA ' + key] = \
            #     [self.com.data[key].avg() - self.rec.data[key].avg() -
            #      self.lig.data[key].avg(),
            #      sqrt(self.com.data[key].stdev() ** 2 +
            #           self.rec.data[key].stdev() ** 2 +
            #           self.lig.data[key].stdev() ** 2) ]
    #==================================================

    def print_summary_csv(self, csvwriter):
        """ Prints the summary of the binding statistics in CSV format """
        if self.verbose:
            csvwriter.writerow(['Complex:'])
            self.com.print_summary_csv(csvwriter)
            csvwriter.writerow(['Receptor:'])
            self.rec.print_summary_csv(csvwriter)
            csvwriter.writerow(['Ligand:'])
            self.lig.print_summary_csv(csvwriter)

        csvwriter.writerow(['Differences (Complex - Receptor - Ligand):'])

        # verbose == 0 means don't print com/rec/lig, but print diffs as though
        # verbose == 2 for multi traj
        if self.verbose < 2: verbose = 2
        else: verbose = self.verbose

        for key in self.data_keys:
            # Skip terms we don't want to print
            if verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            stdev = self.data[key][1]
            avg = self.data[key][0]
            num_frames = min(len(self.com.data[key]),len(self.rec.data[key]),
                             len(self.lig.data[key]))
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(num_frames)])

        for key in self.composite_keys:
            # Now print out the composite terms
            stdev = self.data[key][1]
            avg = self.data[key][0]
            csvwriter.writerow([key, avg, stdev, stdev/sqrt(num_frames)])

    #==================================================

    def print_summary(self):
        """ Returns a string printing the summary of the binding statistics """

        if self.verbose:
            ret_str = 'Complex:\n' + self.com.print_summary()
            ret_str += 'Receptor:\n' + self.rec.print_summary()
            ret_str += 'Ligand:\n' + self.lig.print_summary()

        if isinstance(self.com, NMODEout):
            col_name = '%-16s' % 'Entropy Term'
        else:
            col_name = '%-16s' % 'Energy Component'
        ret_str += 'Differences (Complex - Receptor - Ligand):\n'
        ret_str += (col_name+'            Average              ' +
                    'Std. Dev.   Std. Err. of Mean\n')
        ret_str += ('------------------------------------------------' +
                    '-------------------------------\n')

        # verbose == 0 means don't print com/rec/lig, but print diffs as though
        # verbose == 2 for multi traj
        if self.verbose < 2: verbose = 2
        else: verbose = self.verbose

        for key in self.data_keys:
            # Skip terms we don't want to print
            if verbose < self.print_levels[key]: continue
            # Skip the composite terms, since we print those at the end
            if key in self.composite_keys: continue
            # Skip chamber terms if we aren't using chamber prmtops
            if not self.chamber and key in ['UB', 'IMP', 'CMAP']: continue
            stdev = self.data[key][1]
            num_frames = min(len(self.com.data[key]),len(self.rec.data[key]),
                             len(self.lig.data[key]))
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key,
                                                         self.data[key][0], stdev, stdev / sqrt(num_frames))

        ret_str += '\n'
        for key in self.composite_keys:
            # Now print out the composite terms
            if key == 'DELTA TOTAL': ret_str += '\n'
            stdev = self.data[key][1]
            ret_str += '%-14s %20.4f %21.4f %19.4f\n' % (key,
                                                         self.data[key][0], stdev, stdev / sqrt(num_frames))

        return ret_str + '\n\n'

    #==================================================

    def diff(self, other, key1, key2):
        """
        Takes the difference between 2 keys of 2 different BindingStatistics
        classes and returns the average and standard deviation of that diff.
        """
        return (self.data[key1][0] - other.data[key2][0],
                sqrt(self.data[key1][1]**2 + other.data[key2][1]**2))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class DecompOut(object):
    " Class for decomposition output file to collect statistics and output them "

    indicator = "                    PRINT DECOMP - TOTAL ENERGIES"
    descriptions = { 'TDC' : 'Total Energy Decomposition:',
                     'SDC' : 'Sidechain Energy Decomposition:',
                     'BDC' : 'Backbone Energy Decomposition:' }

    #==================================================

    def __init__(self, basename, prmtop, surften, csvwriter, num_files=1,
                 verbose=1):
        from csv import writer
        self.basename = basename # base name of output files
        self.prmtop = prmtop # AmberParm prmtop object
        self.num_files = num_files # how many MPI files we created
        self.termnum = 0 # which term number we are on
        self.surften = surften # surface tension to multiply SAS by
        self.verbose = verbose
        # Set the term-extractor based on whether we want to dump the values to
        # a CSV file or just get the next term
        if csvwriter:
            self.get_next_term = self._get_next_term_csv
        else:
            self.get_next_term = self._get_next_term
        if verbose in [1,3]:
            self.allowed_tokens = tuple(['TDC','SDC','BDC'])
        else:
            self.allowed_tokens = tuple(['TDC'])

        # Create a separate csvwriter for each of the different token types,
        # and store them in a dictionary
        if csvwriter:
            self.csvwriter = {}
            for tok in self.allowed_tokens:
                self.csvwriter[tok] = writer(open(csvwriter+'.'+tok+'.csv', 'w'))
                self.csvwriter[tok].writerow([self.descriptions[tok]])
                self._write_header(self.csvwriter[tok])
        else:
            self.csvwriter = None

        self.current_file = 0 # File counter
        try:
            self.num_terms = int(self.get_num_terms())
        except TypeError:
            raise OutputError('DecompOut: Not a decomp output file')
        self.data = {}
        for token in self.allowed_tokens:
            self.data[token] = {'int': [EnergyVector(), EnergyVector()],
                                'vdw': [EnergyVector(), EnergyVector()],
                                'eel': [EnergyVector(), EnergyVector()],
                                'pol': [EnergyVector(), EnergyVector()],
                                'sas': [EnergyVector(), EnergyVector()],
                                'tot': [EnergyVector(), EnergyVector()] }
            self.data[token]['int'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['int'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['vdw'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['vdw'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['eel'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['eel'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['pol'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['pol'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['sas'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['sas'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['tot'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['tot'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.resnums               = [[0 for i in range(self.num_terms)],
                                          [0 for i in range(self.num_terms)]]
        self.decfile = open(basename + '.0', 'r')

    #==================================================

    def get_num_terms(self):
        """ Gets the number of terms in the output file """
        decfile = open('%s.%d' % (self.basename, 0), 'r')
        num_terms = 0
        line = decfile.readline()
        while line:
            if not line.startswith(self.indicator):
                line = decfile.readline()
                continue
            while line[0:3] != 'TDC':
                line = decfile.readline()
            while line[0:3] == 'TDC':
                line = decfile.readline()
                num_terms += 1
            # We've now gotten to the end of the Total Decomp Contribution,
            # so we know how many terms we have
            break
        decfile.close()
        return num_terms

    #==================================================

    def _get_next_term(self, expected_type, framenum=1):
        """ Gets the next energy term from the output file(s) """
        line = self.decfile.readline()
        if expected_type and not expected_type in self.allowed_tokens:
            raise OutputError('BUGBUG: expected_type must be in %s' %
                              self.allowed_tokens)
        while not line[0:3] in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            self.termnum = 0
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1: return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename,self.current_file),'r')
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
        sas = float(line[51:60]) * self.surften
        tot = internal + vdw + eel + pol + sas
        self.resnums[0][self.termnum] = resnum
        self.data[line[0:3]]['int'][0][self.termnum] += internal
        self.data[line[0:3]]['int'][1][self.termnum] += internal * internal
        self.data[line[0:3]]['vdw'][0][self.termnum] += vdw
        self.data[line[0:3]]['vdw'][1][self.termnum] += vdw * vdw
        self.data[line[0:3]]['eel'][0][self.termnum] += eel
        self.data[line[0:3]]['eel'][1][self.termnum] += eel * eel
        self.data[line[0:3]]['pol'][0][self.termnum] += pol
        self.data[line[0:3]]['pol'][1][self.termnum] += pol * pol
        self.data[line[0:3]]['sas'][0][self.termnum] += sas
        self.data[line[0:3]]['sas'][1][self.termnum] += sas * sas
        self.data[line[0:3]]['tot'][0][self.termnum] += tot
        self.data[line[0:3]]['tot'][1][self.termnum] += tot * tot
        self.termnum += 1
        return [resnum, internal, vdw, eel, pol, sas, tot]

    #==================================================

    def fill_all_terms(self):
        """
        This is for stability calculations -- just get all of the terms to
        fill up the data arrays.
        """
        token_counter = 0
        searched_type = self.allowed_tokens[0]
        framenum = 1
        com_token = self.get_next_term(self.allowed_tokens[0], framenum)
        while com_token:
            # Get all of the tokens
            for i in range(1, self.num_terms):
                com_token = self.get_next_term(searched_type, framenum)

            token_counter += 1
            searched_type = self.allowed_tokens[token_counter %
                                                len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0: framenum += 1
            com_token = self.get_next_term(searched_type, framenum)

        self.numframes = framenum - 1

    #==================================================

    def _get_next_term_csv(self, expected_type, framenum=1):
        """ Gets the next term and prints data to csv file """
        mydat = self._get_next_term(expected_type)
        if mydat: self.csvwriter[expected_type].writerow([framenum]+mydat)
        else: self.csvwriter[expected_type].writerow([])
        return mydat

    #==================================================

    def _write_header(self, csvwriter):
        """ Writes a table header to a passed CSV file """
        csvwriter.writerow(['Frame #', 'Residue', 'Internal', 'van der Waals',
                            'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    #==================================================

    def write_summary(self, numframes, output_file):
        """ Writes the summary in ASCII format to and open output_file """
        for term in self.allowed_tokens:
            output_file.writeline(self.descriptions[term])
            output_file.writeline('Residue |       Internal      |    ' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation  ' +
                                  ' |   Non-Polar Solv.   |       TOTAL')
            output_file.writeline('-------------------------------------------' +
                                  '--------------------------------------------------------------' +
                                  '----------------------------------')
            # Now loop over all of the terms.
            for i in range(self.num_terms):
                int_avg = self.data[term]['int'][0][i] / numframes
                int_std = sqrt(abs(self.data[term]['int'][1][i]/numframes -
                                   int_avg**2))
                vdw_avg = self.data[term]['vdw'][0][i] / numframes
                vdw_std = sqrt(abs(self.data[term]['vdw'][1][i]/numframes -
                                   vdw_avg**2))
                eel_avg = self.data[term]['eel'][0][i] / numframes
                eel_std = sqrt(abs(self.data[term]['eel'][1][i]/numframes -
                                   eel_avg**2))
                pol_avg = self.data[term]['pol'][0][i] / numframes
                pol_std = sqrt(abs(self.data[term]['pol'][1][i]/numframes -
                                   pol_avg**2))
                sas_avg = self.data[term]['sas'][0][i] / numframes
                sas_std = sqrt(abs(self.data[term]['sas'][1][i]/numframes -
                                   sas_avg**2))
                tot_avg = self.data[term]['tot'][0][i] / numframes
                tot_std = sqrt(abs(self.data[term]['tot'][1][i]/numframes -
                                   tot_avg**2))
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[0][i]-1]
                res_str = '%3s%4d' % (resnm, self.resnums[0][i])
                output_file.writeline(('%s |%9.3f +/- %6.3f |%9.3f +/- %6.3f ' +
                                       '|%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f' +
                                       ' +/- %6.3f') % (res_str, int_avg, int_std, vdw_avg, vdw_std,
                                                        eel_avg, eel_std, pol_avg, pol_std, sas_avg, sas_std,
                                                        tot_avg, tot_std))
            output_file.writeline('')

    #==================================================

    def write_summary_csv(self, numframes, csvwriter):
        """ Writes the summary to a CSV file """
        for term in self.allowed_tokens:
            csvwriter.writerow([self.descriptions[term],])
            csvwriter.writerow(['Residue', 'Internal', '', '', 'van der Waals', '',
                                '', 'Electrostatic', '', '', 'Polar Solvation', '',
                                '', 'Non-Polar Solv.', '', '', 'TOTAL', '', ''])
            csvwriter.writerow([''] + ['Avg.','Std. Dev.', 'Std. Err. of Mean']*6)
            for i in range(self.num_terms):
                sqrt_frames = sqrt(numframes)
                int_avg = self.data[term]['int'][0][i] / numframes
                int_std = sqrt(abs(self.data[term]['int'][1][i]/numframes -
                                   int_avg**2))
                vdw_avg = self.data[term]['vdw'][0][i] / numframes
                vdw_std = sqrt(abs(self.data[term]['vdw'][1][i]/numframes -
                                   vdw_avg**2))
                eel_avg = self.data[term]['eel'][0][i] / numframes
                eel_std = sqrt(abs(self.data[term]['eel'][1][i]/numframes -
                                   eel_avg**2))
                pol_avg = self.data[term]['pol'][0][i] / numframes
                pol_std = sqrt(abs(self.data[term]['pol'][1][i]/numframes -
                                   pol_avg**2))
                sas_avg = self.data[term]['sas'][0][i] / numframes
                sas_std = sqrt(abs(self.data[term]['sas'][1][i]/numframes -
                                   sas_avg**2))
                tot_avg = self.data[term]['tot'][0][i] / numframes
                tot_std = sqrt(abs(self.data[term]['tot'][1][i]/numframes -
                                   tot_avg**2))
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[0][i]-1]
                res_str = '%3s%4d' % (resnm, self.resnums[0][i])
                csvwriter.writerow([res_str, int_avg, int_std, int_std/sqrt_frames,
                                    vdw_avg, vdw_std, vdw_std/sqrt_frames,
                                    eel_avg, eel_std, eel_std/sqrt_frames,
                                    pol_avg, pol_std, pol_std/sqrt_frames,
                                    sas_avg, sas_std, sas_std/sqrt_frames,
                                    tot_avg, tot_std, tot_std/sqrt_frames])
        csvwriter.writerow([])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PairDecompOut(DecompOut):
    """ Same as DecompOut, but for Pairwise decomposition """

    indicator = "                    PRINT PAIR DECOMP - TOTAL ENERGIES"

    #==================================================

    def _get_next_term(self, expected_type=None, framenum=1):
        """ Gets the next energy term from the output file(s) """
        line = self.decfile.readline()
        if expected_type and not expected_type in self.allowed_tokens:
            raise OutputError('BUGBUG: expected_type must be in %s' %
                              self.allowed_tokens)
        while not line[0:3] in self.allowed_tokens:
            # We only get in here if we've gone off the end of a block, so our
            # current term number is 0 now.
            self.termnum = 0
            line = self.decfile.readline()
            if not line:
                self.decfile.close()
                if self.current_file == self.num_files - 1: return []
                self.current_file += 1
                self.decfile = open('%s.%d' % (self.basename,self.current_file),'r')
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
        sas = float(line[73:85]) * self.surften
        tot = internal + vdw + eel + pol + sas
        self.data[line[0:3]]['int'][0][self.termnum] += internal
        self.data[line[0:3]]['int'][1][self.termnum] += internal * internal
        self.data[line[0:3]]['vdw'][0][self.termnum] += vdw
        self.data[line[0:3]]['vdw'][1][self.termnum] += vdw * vdw
        self.data[line[0:3]]['eel'][0][self.termnum] += eel
        self.data[line[0:3]]['eel'][1][self.termnum] += eel * eel
        self.data[line[0:3]]['pol'][0][self.termnum] += pol
        self.data[line[0:3]]['pol'][1][self.termnum] += pol * pol
        self.data[line[0:3]]['sas'][0][self.termnum] += sas
        self.data[line[0:3]]['sas'][1][self.termnum] += sas * sas
        self.data[line[0:3]]['tot'][0][self.termnum] += tot
        self.data[line[0:3]]['tot'][1][self.termnum] += tot * tot
        self.resnums[0][self.termnum] = resnum
        self.resnums[1][self.termnum] = resnum2
        self.termnum += 1
        return [resnum, resnum2, internal, vdw, eel, pol, sas, tot]

    #==================================================

    def write_summary(self, numframes, output_file):
        """ Writes the summary in ASCII format to and open output_file """
        for term in self.allowed_tokens:
            output_file.writeline(self.descriptions[term])
            output_file.writeline('Resid 1 | Resid 2 |       Internal      |    ' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation   |  ' +
                                  'Non-Polar Solv.   |       TOTAL')
            output_file.writeline('---------------------------------------------' +
                                  '----------------------------------------------------------------' +
                                  '----------------------------------------')
            # Now loop over all of the terms.
            for i in range(self.num_terms):
                int_avg = self.data[term]['int'][0][i] / numframes
                int_std = sqrt(abs(self.data[term]['int'][1][i]/numframes -
                                   int_avg**2))
                vdw_avg = self.data[term]['vdw'][0][i] / numframes
                vdw_std = sqrt(abs(self.data[term]['vdw'][1][i]/numframes -
                                   vdw_avg**2))
                eel_avg = self.data[term]['eel'][0][i] / numframes
                eel_std = sqrt(abs(self.data[term]['eel'][1][i]/numframes -
                                   eel_avg**2))
                pol_avg = self.data[term]['pol'][0][i] / numframes
                pol_std = sqrt(abs(self.data[term]['pol'][1][i]/numframes -
                                   pol_avg**2))
                sas_avg = self.data[term]['sas'][0][i] / numframes
                sas_std = sqrt(abs(self.data[term]['sas'][1][i]/numframes -
                                   sas_avg**2))
                tot_avg = self.data[term]['tot'][0][i] / numframes
                tot_std = sqrt(abs(self.data[term]['tot'][1][i]/numframes -
                                   tot_avg**2))
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[0][i]-1]
                res_str0 = '%3s%4d' % (resnm, self.resnums[0][i])
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[1][i]-1]
                res_str1 = '%3s%4d' % (resnm, self.resnums[1][i])
                output_file.writeline(('%s | %s |%9.3f +/- %6.3f |%9.3f +/- %6.3f' +
                                       ' |%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f +/- %6.3f ' +
                                       '|%9.3f +/- %6.3f') % (res_str0, res_str1,
                                                              int_avg, int_std, vdw_avg, vdw_std, eel_avg, eel_std,
                                                              pol_avg, pol_std, sas_avg, sas_std, tot_avg, tot_std))
            output_file.writeline('')

    #==================================================

    def _write_header(self, csvwriter):
        """ Writes a table header to the csvwriter """
        csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    #==================================================

    def write_summary_csv(self, numframes, csvwriter):
        """ Writes the summary to a CSV file """
        for term in self.allowed_tokens:
            csvwriter.writerow([self.descriptions[term]])
            csvwriter.writerow(['Resid 1', 'Resid 2', 'Internal', '', '',
                                'van der Waals', '', '', 'Electrostatic', '', '',
                                'Polar Solvation', '', '', 'Non-Polar Solv.', '',
                                '', 'TOTAL', '', ''])
            csvwriter.writerow(['']*2+['Avg.','Std. Dev.', 'Std. Err. of Mean']*6)
            for i in range(self.num_terms):
                sqrt_frames = sqrt(numframes)
                int_avg = self.data[term]['int'][0][i] / numframes
                int_std = sqrt(abs(self.data[term]['int'][1][i]/numframes -
                                   int_avg**2))
                vdw_avg = self.data[term]['vdw'][0][i] / numframes
                vdw_std = sqrt(abs(self.data[term]['vdw'][1][i]/numframes -
                                   vdw_avg**2))
                eel_avg = self.data[term]['eel'][0][i] / numframes
                eel_std = sqrt(abs(self.data[term]['eel'][1][i]/numframes -
                                   eel_avg**2))
                pol_avg = self.data[term]['pol'][0][i] / numframes
                pol_std = sqrt(abs(self.data[term]['pol'][1][i]/numframes -
                                   pol_avg**2))
                sas_avg = self.data[term]['sas'][0][i] / numframes
                sas_std = sqrt(abs(self.data[term]['sas'][1][i]/numframes -
                                   sas_avg**2))
                tot_avg = self.data[term]['tot'][0][i] / numframes
                tot_std = sqrt(abs(self.data[term]['tot'][1][i]/numframes -
                                   tot_avg**2))
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[0][i]-1]
                res_str0 = '%3s%4d' % (resnm, self.resnums[0][i])
                resnm = self.prmtop.parm_data['RESIDUE_LABEL'][self.resnums[1][i]-1]
                res_str1 = '%3s%4d' % (resnm, self.resnums[1][i])
                csvwriter.writerow([res_str0, res_str1,
                                    int_avg, int_std, int_std/sqrt_frames,
                                    vdw_avg, vdw_std, vdw_std/sqrt_frames,
                                    eel_avg, eel_std, eel_std/sqrt_frames,
                                    pol_avg, pol_std, pol_std/sqrt_frames,
                                    sas_avg, sas_std, sas_std/sqrt_frames,
                                    tot_avg, tot_std, tot_std/sqrt_frames])
            csvwriter.writerow([])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class DecompBinding(object):
    """ Class for decomposition binding (per-residue) """

    #==================================================

    def __init__(self, com, rec, lig, prmtop_system, idecomp, verbose, output,
                 csvwriter, desc):
        """
        output should be an open file and csvfile should be a csv.writer class. If
        the output format is specified as csv, then output should be a csv.writer
        class as well.
        """
        from csv import writer
        (self.com, self.rec, self.lig) = com, rec, lig
        self.num_terms = self.com.num_terms
        self.numframes = 0 # frame counter
        self.output = output
        self.desc = desc # Description
        self.prmtop_system = prmtop_system
        self.idecomp = idecomp
        self.verbose = verbose
        # Check to see if output is a csv.writer or if it's a file. The invoked
        # method, "parse_all", is set based on whether we're doing a csv output
        # or an ascii output
        if type(output).__name__ == 'writer': # yuck... better way?
            self.parse_all = self._parse_all_csv
        else:
            self.parse_all = self._parse_all_ascii
        # Set up the data for the DELTAs
        if verbose in [1,3]:
            self.allowed_tokens = tuple(['TDC','SDC','BDC'])
        else:
            self.allowed_tokens = tuple(['TDC'])
        # Open up a separate CSV writer for all of the allowed tokens
        if csvwriter:
            self.csvwriter = {}
            for tok in self.allowed_tokens:
                self.csvwriter[tok] = writer(open(csvwriter+'.'+tok+'.csv', 'w'))
                self.csvwriter[tok].writerow(['DELTA', DecompOut.descriptions[tok]])
                self._write_header(self.csvwriter[tok])
        else:
            self.csvwriter = None
        self.data = {}
        self.data_stats = {}
        for token in self.allowed_tokens:
            self.data[token] = {'int': [EnergyVector(), EnergyVector()],
                                'vdw': [EnergyVector(), EnergyVector()],
                                'eel': [EnergyVector(), EnergyVector()],
                                'pol': [EnergyVector(), EnergyVector()],
                                'sas': [EnergyVector(), EnergyVector()],
                                'tot': [EnergyVector(), EnergyVector()] }
            self.data_stats[token] = {'int': [EnergyVector(), EnergyVector()],
                                      'vdw': [EnergyVector(), EnergyVector()],
                                      'eel': [EnergyVector(), EnergyVector()],
                                      'pol': [EnergyVector(), EnergyVector()],
                                      'sas': [EnergyVector(), EnergyVector()],
                                      'tot': [EnergyVector(), EnergyVector()] }
            self.data[token]['int'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['int'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['vdw'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['vdw'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['eel'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['eel'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['pol'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['pol'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['sas'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['sas'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['tot'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data[token]['tot'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['int'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['int'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['vdw'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['vdw'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['eel'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['eel'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['pol'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['pol'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['sas'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['sas'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['tot'][0] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.data_stats[token]['tot'][1] = EnergyVector(
                [0 for i in range(self.num_terms)])
            self.resnums               = [[0 for i in range(self.num_terms)],
                                          [0 for i in range(self.num_terms)]]

    #==================================================

    def _write_header(self, csvwriter):
        """ Writes the header to the CSV file (legend at top of chart) """
        csvwriter.writerow(['Frame #', 'Residue', 'Location', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    #==================================================

    def _parse_all_begin(self):
        """ Parses through all of the terms in all of the frames, but doesn't
            do any printing
        """
        # For per-residue decomp, we need terms to match up
        if self.com.num_terms != (self.rec.num_terms + self.lig.num_terms):
            raise DecompError('Mismatch in number of decomp terms!')
        # Get the first complex term, then parse through the rest of them
        token_counter = 0
        searched_token = self.allowed_tokens[0]
        framenum = 1
        com_token = self.com.get_next_term(searched_token, framenum)
        while com_token:
            # Figure out our resnum and location
            self.resnums[0][0] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                             parm_data['RESIDUE_LABEL'][com_token[0]-1], com_token[0])
            if self.prmtop_system.res_list[com_token[0]-1].receptor_number:
                self.resnums[1][0] = 'R %3s%4d' % (self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                                                        com_token[0]-1],
                                                   self.prmtop_system.res_list[com_token[0]-1].receptor_number)
                other_token = self.rec.get_next_term(searched_token, framenum)
            else:
                self.resnums[1][0] = 'L %3s%4d' % (self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                                                        com_token[0]-1],
                                                   self.prmtop_system.res_list[com_token[0]-1].ligand_number)
                other_token = self.lig.get_next_term(searched_token, framenum)

            # Fill the data array
            self.data[searched_token]['int'][0][0] += (com_token[1] - other_token[1])
            self.data[searched_token]['int'][1][0] += (com_token[1] - other_token[1]) * (com_token[1] - other_token[1])
            self.data[searched_token]['vdw'][0][0] += (com_token[2] - other_token[2])
            self.data[searched_token]['vdw'][1][0] += (com_token[2] - other_token[2]) * (com_token[2] - other_token[2])
            self.data[searched_token]['eel'][0][0] += (com_token[3] - other_token[3])
            self.data[searched_token]['eel'][1][0] += (com_token[3] - other_token[3]) * (com_token[3] - other_token[3])
            self.data[searched_token]['pol'][0][0] += (com_token[4] - other_token[4])
            self.data[searched_token]['pol'][1][0] += (com_token[4] - other_token[4]) * (com_token[4] - other_token[4])
            self.data[searched_token]['sas'][0][0] += (com_token[5] - other_token[5])
            self.data[searched_token]['sas'][1][0] += (com_token[5] - other_token[5]) * (com_token[5] - other_token[5])
            self.data[searched_token]['tot'][0][0] += (com_token[6] - other_token[6])
            self.data[searched_token]['tot'][1][0] += (com_token[6] - other_token[6]) * (com_token[6] - other_token[6])
            if self.csvwriter:
                self.csvwriter[searched_token].writerow([framenum,
                                                         self.resnums[0][0],
                                                         self.resnums[1][0],
                                                         com_token[1]-other_token[1],
                                                         com_token[2]-other_token[2],
                                                         com_token[3]-other_token[3],
                                                         com_token[4]-other_token[4],
                                                         com_token[5]-other_token[5],
                                                         com_token[6]-other_token[6] ])
            for i in range(1, self.num_terms):
                com_token = self.com.get_next_term(searched_token, framenum)
                # First see if it's in the receptor
                if self.prmtop_system.res_list[com_token[0]-1].receptor_number:
                    self.resnums[1][i] = 'R %3s%4d' % (
                        self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                            com_token[0]-1],
                        self.prmtop_system.res_list[com_token[0]-1].receptor_number)
                    other_token = self.rec.get_next_term(searched_token, framenum)
                else:
                    self.resnums[1][i] = 'L %3s%4d' % (
                        self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                            com_token[0]-1],
                        self.prmtop_system.res_list[com_token[0]-1].ligand_number)
                    other_token = self.lig.get_next_term(searched_token, framenum)
                # Figure out which residue numbers we are and where we are mapped
                self.resnums[0][i] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                                 parm_data['RESIDUE_LABEL'][com_token[0]-1], com_token[0])
                if self.prmtop_system.res_list[com_token[0]-1].receptor_number:
                    self.resnums[1][i] = 'R %3s%4d' % (
                        self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                            com_token[0]-1],
                        self.prmtop_system.res_list[com_token[0]-1].receptor_number)
                else:
                    self.resnums[1][i] = 'L %3s%4d' % (
                        self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                            com_token[0]-1],
                        self.prmtop_system.res_list[com_token[0]-1].ligand_number)
                # Now fill the arrays
                self.data[searched_token]['int'][0][i] += (com_token[1] - other_token[1])
                self.data[searched_token]['int'][1][i] += (com_token[1] - other_token[1])\
                                                          * (com_token[1] - other_token[1])
                self.data[searched_token]['vdw'][0][i] += (com_token[2] - other_token[2])
                self.data[searched_token]['vdw'][1][i] += (com_token[2] - other_token[2])\
                                                          * (com_token[2] - other_token[2])
                self.data[searched_token]['eel'][0][i] += (com_token[3] - other_token[3])
                self.data[searched_token]['eel'][1][i] += (com_token[3] - other_token[3])\
                                                          * (com_token[3] - other_token[3])
                self.data[searched_token]['pol'][0][i] += (com_token[4] - other_token[4])
                self.data[searched_token]['pol'][1][i] += (com_token[4] - other_token[4])\
                                                          * (com_token[4] - other_token[4])
                self.data[searched_token]['sas'][0][i] += (com_token[5] - other_token[5])
                self.data[searched_token]['sas'][1][i] += (com_token[5] - other_token[5])\
                                                          * (com_token[5] - other_token[5])
                self.data[searched_token]['tot'][0][i] += (com_token[6] - other_token[6])
                self.data[searched_token]['tot'][1][i] += (com_token[6] - other_token[6])\
                                                          * (com_token[6] - other_token[6])
                if self.csvwriter:
                    self.csvwriter[searched_token].writerow([framenum,
                                                             self.resnums[0][i],
                                                             self.resnums[1][i],
                                                             com_token[1]-other_token[1],
                                                             com_token[2]-other_token[2],
                                                             com_token[3]-other_token[3],
                                                             com_token[4]-other_token[4],
                                                             com_token[5]-other_token[5],
                                                             com_token[6]-other_token[6] ])
            # end for i in range(self.num_terms)
            token_counter += 1
            searched_token = self.allowed_tokens[token_counter %
                                                 len(self.allowed_tokens)]
            # We are going on to the next frame
            if token_counter % len(self.allowed_tokens) == 0:
                framenum += 1
            # Get the first com_token of the next searched_token
            com_token = self.com.get_next_term(searched_token, framenum)

        # Now figure out how many frames we calculated -- this is just how many
        # total tokens we counted // number of distinct tokens we have
        self.numframes = framenum - 1
        self.num_com_frames = self.numframes
        self.num_rec_frames = self.numframes
        self.num_lig_frames = self.numframes
        self._calc_avg_stdev()

    #==================================================

    def _calc_avg_stdev(self):
        """ Calculates the averages and standard deviations of all of the data """
        # Do population averages and such
        for key1 in list(self.data.keys()):
            for key2 in list(self.data[key1].keys()):
                for i in range(self.num_terms):
                    myavg = self.data[key1][key2][0][i] / self.numframes
                    myavg2 = self.data[key1][key2][1][i] / self.numframes
                    self.data_stats[key1][key2][0][i] = myavg
                    self.data_stats[key1][key2][1][i] = sqrt(abs(myavg2-myavg*myavg))

    #==================================================

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
            self.com.write_summary_csv(self.num_com_frames, self.output)
            self.output.writerow(['Receptor:'])
            self.rec.write_summary_csv(self.num_rec_frames, self.output)
            self.output.writerow(['Ligand:'])
            self.lig.write_summary_csv(self.num_lig_frames, self.output)
        # Now write the DELTAs

        self.output.writerow(['DELTAS:'])
        for term in self.allowed_tokens:
            self.output.writerow([DecompOut.descriptions[term]])
            self.output.writerow(['Residue', 'Location', 'Internal', '', '',
                                  'van der Waals', '', '', 'Electrostatic', '', '',
                                  'Polar Solvation', '', '', 'Non-Polar Solv.',
                                  '', '', 'TOTAL', '', ''])
            self.output.writerow(['', ''] +
                                 ['Avg.','Std. Dev.', 'Std. Err. of Mean']*5)
            for i in range(self.num_terms):
                sqrt_frames = sqrt(self.num_com_frames)
                int_avg = self.data_stats[term]['int'][0][i]
                int_std = self.data_stats[term]['int'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                eel_avg = self.data_stats[term]['eel'][0][i]
                eel_std = self.data_stats[term]['eel'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                pol_avg = self.data_stats[term]['pol'][0][i]
                pol_std = self.data_stats[term]['pol'][1][i]
                sas_avg = self.data_stats[term]['sas'][0][i]
                sas_std = self.data_stats[term]['sas'][1][i]
                tot_avg = self.data_stats[term]['tot'][0][i]
                tot_std = self.data_stats[term]['tot'][1][i]
                self.output.writerow([self.resnums[0][i], self.resnums[1][i],
                                      int_avg, int_std, int_std/sqrt_frames,
                                      vdw_avg, vdw_std, vdw_std/sqrt_frames,
                                      eel_avg, eel_std, eel_std/sqrt_frames,
                                      pol_avg, pol_std, pol_std/sqrt_frames,
                                      sas_avg, sas_std, sas_std/sqrt_frames,
                                      tot_avg, tot_std, tot_std/sqrt_frames])
            self.output.writerow([])

    #==================================================

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
            self.com.write_summary(self.num_com_frames, self.output)
            self.output.writeline('Receptor:')
            self.rec.write_summary(self.num_rec_frames, self.output)
            self.output.writeline('Ligand:')
            self.lig.write_summary(self.num_lig_frames, self.output)
        # Now write the DELTAs

        self.output.writeline('DELTAS:')
        for term in self.allowed_tokens:
            self.output.writeline(DecompOut.descriptions[term])
            self.output.writeline(
                'Residue |  Location |       Internal      |    ' +
                'van der Waals    |    Electrostatic    |   Polar Solvation   |' +
                '    Non-Polar Solv.  |       TOTAL')
            self.output.writeline('---------------------------------------------' +
                                  '--------------------------------------------------------------' +
                                  '--------------------------------------------')
            for i in range(self.num_terms):
                int_avg = self.data_stats[term]['int'][0][i]
                int_std = self.data_stats[term]['int'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                eel_avg = self.data_stats[term]['eel'][0][i]
                eel_std = self.data_stats[term]['eel'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                pol_avg = self.data_stats[term]['pol'][0][i]
                pol_std = self.data_stats[term]['pol'][1][i]
                sas_avg = self.data_stats[term]['sas'][0][i]
                sas_std = self.data_stats[term]['sas'][1][i]
                tot_avg = self.data_stats[term]['tot'][0][i]
                tot_std = self.data_stats[term]['tot'][1][i]
                self.output.writeline(('%s | %s |%9.3f +/- %6.3f |%9.3f +/- %6.3f' +
                                       ' |%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f ' +
                                       '+/- %6.3f') % (self.resnums[0][i], self.resnums[1][i],
                                                       int_avg, int_std, vdw_avg, vdw_std, eel_avg, eel_std,
                                                       pol_avg, pol_std, sas_avg, sas_std, tot_avg, tot_std))
            self.output.writeline('')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PairDecompBinding(DecompBinding):
    """ Class for decomposition binding (pairwise) """

    #==================================================

    def _write_header(self, csvwriter):
        """ Writes the header to the CSV file (legend at top of chart) """
        csvwriter.writerow(['Frame #', 'Resid 1', 'Resid 2', 'Internal',
                            'van der Waals', 'Electrostatic', 'Polar Solvation',
                            'Non-Polar Solv.', 'TOTAL'])

    #==================================================

    def _parse_all_begin(self):
        """ Parses through all of the terms in all of the frames, but doesn't
            do any printing
        """
        # Get the first complex term, then parse through the rest of them
        token_counter = 0
        searched_token = self.allowed_tokens[0]
        framenum = 1
        com_token = self.com.get_next_term(searched_token, framenum)
        while com_token:
            # Figure out our resnum and location
            self.resnums[0][0] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                             parm_data['RESIDUE_LABEL'][com_token[0]-1], com_token[0])
            self.resnums[1][0] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                             parm_data['RESIDUE_LABEL'][com_token[1]-1], com_token[1])
            if self.prmtop_system.res_list[com_token[0]-1].receptor_number:
                if self.prmtop_system.res_list[com_token[1]-1].receptor_number:
                    # Both residues are in the receptor -- pull the next one
                    other_token = self.rec.get_next_term(searched_token, framenum)
                else:
                    other_token = [0,0,0,0,0,0,0,0]
            else:
                if self.prmtop_system.res_list[com_token[1]-1].ligand_number:
                    # Both residues are in the ligand -- pull the next one
                    other_token = self.lig.get_next_term(searched_token, framenum)
                else:
                    other_token = [0,0,0,0,0,0,0,0]

            # Fill the data array
            self.data[searched_token]['int'][0][0] += (com_token[2] - other_token[2])
            self.data[searched_token]['int'][1][0] += (com_token[2] - other_token[2]) * (com_token[2] - other_token[2])
            self.data[searched_token]['vdw'][0][0] += (com_token[3] - other_token[3])
            self.data[searched_token]['vdw'][1][0] += (com_token[3] - other_token[3]) * (com_token[3] - other_token[3])
            self.data[searched_token]['eel'][0][0] += (com_token[4] - other_token[4])
            self.data[searched_token]['eel'][1][0] += (com_token[4] - other_token[4]) * (com_token[4] - other_token[4])
            self.data[searched_token]['pol'][0][0] += (com_token[5] - other_token[5])
            self.data[searched_token]['pol'][1][0] += (com_token[5] - other_token[5]) * (com_token[5] - other_token[5])
            self.data[searched_token]['sas'][0][0] += (com_token[6] - other_token[6])
            self.data[searched_token]['sas'][1][0] += (com_token[6] - other_token[6]) * (com_token[6] - other_token[6])
            self.data[searched_token]['tot'][0][0] += (com_token[7] - other_token[7])
            self.data[searched_token]['tot'][1][0] += (com_token[7] - other_token[7]) * (com_token[7] - other_token[7])
            if self.csvwriter:
                self.csvwriter[searched_token].writerow([framenum,
                                                         self.resnums[0][0],
                                                         self.resnums[1][0],
                                                         com_token[2]-other_token[2],
                                                         com_token[3]-other_token[3],
                                                         com_token[4]-other_token[4],
                                                         com_token[5]-other_token[5],
                                                         com_token[6]-other_token[6],
                                                         com_token[7]-other_token[7]  ])
            for i in range(1, self.num_terms):
                com_token = self.com.get_next_term(searched_token, framenum)

                # Figure out our resnum and location
                self.resnums[0][i] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                                 parm_data['RESIDUE_LABEL'][com_token[0]-1], com_token[0])
                self.resnums[1][i] = '%3s%4d' % (self.prmtop_system.complex_prmtop.
                                                 parm_data['RESIDUE_LABEL'][com_token[1]-1], com_token[1])
                if self.prmtop_system.res_list[com_token[0]-1].receptor_number:
                    if self.prmtop_system.res_list[com_token[1]-1].receptor_number:
                        # Both residues are in the receptor -- pull the next one
                        other_token = self.rec.get_next_term(searched_token, framenum)
                    else:
                        other_token = [0,0,0,0,0,0,0,0]
                else:
                    if self.prmtop_system.res_list[com_token[1]-1].ligand_number:
                        # Both residues are in the ligand -- pull the next one
                        other_token = self.lig.get_next_term(searched_token, framenum)
                    else:
                        other_token = [0,0,0,0,0,0,0,0]

                # Now fill the arrays
                self.data[searched_token]['int'][0][i] += (com_token[2] - other_token[2])
                self.data[searched_token]['int'][1][i] += (com_token[2] - other_token[2])\
                                                            * (com_token[2] - other_token[2])
                self.data[searched_token]['vdw'][0][i] += (com_token[3] - other_token[3])
                self.data[searched_token]['vdw'][1][i] += (com_token[3] - other_token[3])\
                                                            * (com_token[3] - other_token[3])
                self.data[searched_token]['eel'][0][i] += (com_token[4] - other_token[4])
                self.data[searched_token]['eel'][1][i] += (com_token[4] - other_token[4])\
                                                            * (com_token[4] - other_token[4])
                self.data[searched_token]['pol'][0][i] += (com_token[5] - other_token[5])
                self.data[searched_token]['pol'][1][i] += (com_token[5] - other_token[5])\
                                                            * (com_token[5] - other_token[5])
                self.data[searched_token]['sas'][0][i] += (com_token[6] - other_token[6])
                self.data[searched_token]['sas'][1][i] += (com_token[6] - other_token[6])\
                                                            * (com_token[6] - other_token[6])
                self.data[searched_token]['tot'][0][i] += (com_token[7] - other_token[7])
                self.data[searched_token]['tot'][1][i] += (com_token[7] - other_token[7])\
                                                            * (com_token[7] - other_token[7])
                if self.csvwriter:
                    self.csvwriter[searched_token].writerow([framenum,
                                                             self.resnums[0][i],
                                                             self.resnums[1][i],
                                                             com_token[2]-other_token[2],
                                                             com_token[3]-other_token[3],
                                                             com_token[4]-other_token[4],
                                                             com_token[5]-other_token[5],
                                                             com_token[6]-other_token[6],
                                                             com_token[7]-other_token[7] ])
            # end for i in range(self.num_terms)
            token_counter += 1
            searched_token = self.allowed_tokens[token_counter %
                                                 len(self.allowed_tokens)]
            # Get the first com_token of the next searched_token
            if token_counter % len(self.allowed_tokens) == 0: framenum += 1
            com_token = self.com.get_next_term(searched_token, framenum)

        # Now figure out how many frames we calculated -- this is just how many
        # total tokens we counted // number of distinct tokens we have
        self.numframes = framenum - 1
        self.num_com_frames = self.numframes
        self.num_rec_frames = self.numframes
        self.num_lig_frames = self.numframes
        self._calc_avg_stdev()

    #==================================================

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
            self.com.write_summary_csv(self.num_com_frames, self.output)
            self.output.writerow(['Receptor:'])
            self.rec.write_summary_csv(self.num_rec_frames, self.output)
            self.output.writerow(['Ligand:'])
            self.lig.write_summary_csv(self.num_lig_frames, self.output)
        # Now write the DELTAs

        self.output.writerow('DELTAS:')
        for term in self.allowed_tokens:
            self.output.writerow(DecompOut.descriptions[term])
            self.output.writerow(['Resid 1', 'Resid 2', 'Internal', '', '',
                                  'van der Waals', '', '', 'Electrostatic', '', '',
                                  'Polar Solvation', '', '', 'Non-Polar Solv.',
                                  '', '', 'TOTAL', '', ''])
            self.output.writerow(['',''] +
                                 ['Avg.','Std. Dev.', 'Std. Err. of Mean']*5)
            for i in range(self.num_terms):
                sqrt_frames = sqrt(self.num_com_frames)
                int_avg = self.data_stats[term]['int'][0][i]
                int_std = self.data_stats[term]['int'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                eel_avg = self.data_stats[term]['eel'][0][i]
                eel_std = self.data_stats[term]['eel'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                pol_avg = self.data_stats[term]['pol'][0][i]
                pol_std = self.data_stats[term]['pol'][1][i]
                sas_avg = self.data_stats[term]['sas'][0][i]
                sas_std = self.data_stats[term]['sas'][1][i]
                tot_avg = self.data_stats[term]['tot'][0][i]
                tot_std = self.data_stats[term]['tot'][1][i]
                self.output.writerow([self.resnums[0][i], self.resnums[1][i],
                                      int_avg, int_std, int_std/sqrt_frames,
                                      vdw_avg, vdw_std, vdw_std/sqrt_frames,
                                      eel_avg, eel_std, eel_std/sqrt_frames,
                                      pol_avg, pol_std, pol_std/sqrt_frames,
                                      sas_avg, sas_std, sas_std/sqrt_frames,
                                      tot_avg, tot_std, tot_std/sqrt_frames])
            self.output.writerow([])

    #==================================================

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
            self.com.write_summary(self.num_com_frames, self.output)
            self.output.writeline('Receptor:')
            self.rec.write_summary(self.num_rec_frames, self.output)
            self.output.writeline('Ligand:')
            self.lig.write_summary(self.num_lig_frames, self.output)
        # Now write the DELTAs

        self.output.writeline('DELTAS:')
        for term in self.allowed_tokens:
            self.output.writeline(DecompOut.descriptions[term])
            self.output.writeline('Resid 1 | Resid 2 |       Internal      |' +
                                  'van der Waals    |    Electrostatic    |   Polar Solvation   |' +
                                  'Non-Polar Solv.   |       TOTAL')
            self.output.writeline('---------------------------------------------' +
                                  '--------------------------------------------------------------' +
                                  '------------------------------------------')
            for i in range(self.num_terms):
                int_avg = self.data_stats[term]['int'][0][i]
                int_std = self.data_stats[term]['int'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                eel_avg = self.data_stats[term]['eel'][0][i]
                eel_std = self.data_stats[term]['eel'][1][i]
                vdw_avg = self.data_stats[term]['vdw'][0][i]
                vdw_std = self.data_stats[term]['vdw'][1][i]
                pol_avg = self.data_stats[term]['pol'][0][i]
                pol_std = self.data_stats[term]['pol'][1][i]
                sas_avg = self.data_stats[term]['sas'][0][i]
                sas_std = self.data_stats[term]['sas'][1][i]
                tot_avg = self.data_stats[term]['tot'][0][i]
                tot_std = self.data_stats[term]['tot'][1][i]
                self.output.writeline(('%s | %s |%9.3f +/- %6.3f |%9.3f +/- ' +
                                       '%6.3f |%9.3f +/- %6.3f |%9.3f +/- %6.3f |%9.3f +/- %6.3f  ' +
                                       '|%9.3f +/- %6.3f') % (self.resnums[0][i], self.resnums[1][i],
                                                              int_avg, int_std, vdw_avg, vdw_std, eel_avg, eel_std,
                                                              pol_avg, pol_std, sas_avg, sas_std, tot_avg, tot_std))
            self.output.writeline('')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MultiTrajDecompBinding(DecompBinding):
    """ Same as DecompBinding class, except for multiple trajectories """

    #==================================================

    def _parse_all_begin(self):
        """ Parses all of the files """
        self.num_com_frames = 0
        self.num_rec_frames = 0
        self.num_lig_frames = 0
        token_counter = 0
        framenum = 1
        searched_token = self.allowed_tokens[0]
        my_term = self.com.get_next_term(searched_token, framenum)
        while my_term:
            for i in range(1, self.com.num_terms):
                my_term = self.com.get_next_term(searched_token, framenum)

            token_counter += 1
            searched_token = self.allowed_tokens[token_counter %
                                                 len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0: framenum += 1
            my_term = self.com.get_next_term(searched_token, framenum)
        self.num_com_frames = framenum - 1

        token_counter = 0
        searched_token = self.allowed_tokens[0]
        framenum = 1
        my_term = self.rec.get_next_term(searched_token, framenum)
        while my_term:
            for i in range(1, self.rec.num_terms):
                my_term = self.rec.get_next_term(searched_token, framenum)

            token_counter += 1
            searched_token = self.allowed_tokens[token_counter %
                                                 len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0: framenum += 1
            my_term = self.rec.get_next_term(searched_token, framenum)
        self.num_rec_frames = framenum - 1

        token_counter = 0
        searched_token = self.allowed_tokens[0]
        framenum = 1
        my_term = self.lig.get_next_term(searched_token, framenum)
        while my_term:
            for i in range(1, self.lig.num_terms):
                my_term = self.lig.get_next_term(searched_token, framenum)

            token_counter += 1
            searched_token = self.allowed_tokens[token_counter %
                                                 len(self.allowed_tokens)]
            if token_counter % len(self.allowed_tokens) == 0: framenum += 1
            my_term = self.lig.get_next_term(searched_token, framenum)
        self.num_lig_frames = framenum - 1

        # Fill the self.resnums
        for i in range(self.com.num_terms):
            resnm = self.com.resnums[0][i]-1
            resnam = self.prmtop_system.complex_prmtop.parm_data['RESIDUE_LABEL'][
                resnm]
            self.resnums[0][i] = '%3s%4d' % (resnam, resnm+1)
            if self.prmtop_system.res_list[resnm].receptor_number:
                self.resnums[1][i] = 'R %3s%4d' % (resnam,
                                                   self.prmtop_system.res_list[resnm].receptor_number)
            else:
                self.resnums[1][i] = 'L %3s%4d' % (resnam,
                                                   self.prmtop_system.res_list[resnm].ligand_number)
        self._calc_avg_stdev()

    #==================================================

    def _calc_avg_stdev(self):
        """ Calculates standard deviation and averages """
        # Use error propagation and deltas of averages
        for key1 in list(self.data.keys()):
            for key2 in list(self.data[key1].keys()):
                num_rec_terms, num_lig_terms = 0, 0
                for i in range(self.com.num_terms):
                    # Get the value of the other term -- either in ligand
                    # or receptor
                    if self.prmtop_system.res_list[self.com.resnums[0][i]
                                                   -1].receptor_number:
                        other = self.rec
                        numframes = self.num_rec_frames
                        oi = num_rec_terms
                        num_rec_terms += 1
                    else:
                        other = self.lig
                        numframes = self.num_lig_frames
                        oi = num_lig_terms
                        num_lig_terms += 1
                    cavg = self.com.data[key1][key2][0][i] / self.num_com_frames
                    oavg = other.data[key1][key2][0][oi] / numframes
                    cvar = abs(self.com.data[key1][key2][1][i]/self.num_com_frames -
                               cavg * cavg )
                    ovar = abs(other.data[key1][key2][1][oi] / numframes -
                               oavg * oavg )
                    self.data_stats[key1][key2][0][i] = cavg - oavg
                    self.data_stats[key1][key2][1][i] = sqrt(cvar + ovar)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class MultiTrajPairDecompBinding(MultiTrajDecompBinding, PairDecompBinding):
    """ Same as PairDecompBinding, but for multiple trajectories """

    #==================================================

    def _calc_avg_stdev(self):
        """ Calculates standard deviation and averages """
        # Use error propagation and deltas of averages
        for key1 in list(self.data.keys()):
            for key2 in list(self.data[key1].keys()):
                num_rec_terms, num_lig_terms = 0, 0
                for i in range(self.num_terms):
                    # Get the value of the other term -- either in ligand
                    # or receptor
                    if self.prmtop_system.res_list[self.com.resnums[0][i]-1]. \
                            receptor_number:
                        if self.prmtop_system.res_list[self.com.resnums[1][i]-1]. \
                                receptor_number:
                            alone = False
                            other = self.rec
                            numframes = self.num_rec_frames
                            oi = num_rec_terms
                            num_rec_terms += 1
                        else: alone = True # Only complex has this interaction
                    else:
                        if self.prmtop_system.res_list[self.com.resnums[1][i]-1]. \
                                ligand_number:
                            alone = False
                            other = self.lig
                            numframes = self.num_lig_frames
                            oi = num_lig_terms
                            num_lig_terms += 1
                        else: alone = True # Only complex has this interaction

                    if alone:
                        self.data_stats[key1][key2][0][i] = \
                            self.com.data[key1][key2][0][i]
                        self.data_stats[key1][key2][1][i] = \
                            self.com.data[key1][key2][1][i]
                    else:
                        cavg = self.com.data[key1][key2][0][i] / self.num_com_frames
                        oavg = other.data[key1][key2][0][oi] / numframes
                        cvar = abs(self.com.data[key1][key2][1][i]/self.num_com_frames
                                   - cavg * cavg )
                        ovar = abs(other.data[key1][key2][1][oi] / numframes
                                   - oavg * oavg )
                        self.data_stats[key1][key2][0][i] = cavg - oavg
                        self.data_stats[key1][key2][1][i] = sqrt(cvar + ovar)

    #==================================================

    # Some methods should inherit from PairDecompBinding instead of
    # MultiTrajDecompBinding

    _parse_all_csv = PairDecompBinding._parse_all_csv

    _parse_all_ascii = PairDecompBinding._parse_all_ascii

    # The rest are inherited from MultiTrajDecompBinding

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

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
        vec.append(float(line.split()[1]))

    return vec
