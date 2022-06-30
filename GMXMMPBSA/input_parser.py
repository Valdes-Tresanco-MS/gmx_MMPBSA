"""
This is a module that contains functions responsible for parsing the
input file for gmx_MMPBSA. It must be included with gmx_MMPBSA to
ensure proper functioning.
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

from GMXMMPBSA.exceptions import InputError, InternalError
from GMXMMPBSA import __version__
import re


class Variable(object):
    """
    Base variable class. It has a name and a single value
    """

    def __init__(self, varname, dat_type=int, default=None, description='', int_dat_type=str):
        """ Initializes the variable type. Sets the default value as well as
          specifying how many characters are required by the parser to trigger
          recognition
        """
        # Catch illegalities
        if dat_type not in (int, str, float, list, tuple):
            raise InputError('Variable has unknown data type %s' % dat_type.__name__)
        if int_dat_type not in (int, str, float):
            raise InputError('Variable has unknown internal data type %s' % int_dat_type)

        self.name = varname
        self.datatype = dat_type
        self.int_datatype = int_dat_type
        if default is None:
            self.value = None
        elif self.datatype is str:
            self.value = default.replace("'", '').replace('"', '')
        elif self.datatype in [list, tuple]:
            if isinstance(default, str):
                self.value = [self.int_datatype(x.strip()) for x in re.split(';\s*|,\s*', default.replace('"',''))]
            else:
                self.value = default
        else:
            self.value = self.datatype(default)
        self.description = description

    def __str__(self):
        """ Prints statistics for the variable """
        string = 'Variable name:  %s\n' % self.name
        string += 'Variable type:  %s\n' % self.datatype
        string += 'Variable value: %s\n' % self.value
        string += 'Description:    %s\n' % self.description
        return string

    def help_str(self):
        """ returns the string [<name> = <value>.... # description] """
        if self.datatype is str:
            valstring = f'{self.name:20s} = "{self.value:s}"'
        elif self.datatype in [list, tuple]:
            if self.int_datatype == str:
                v = ','.join(self.value)
                valstring = f'{self.name:20s} = "{v}"'
            else:
                v = ','.join(map(str, self.value))
                valstring = f'{self.name:20s} = {v}'
        else:
            valstring = f'{self.name:20s} = {self.value}'
        length = 70
        valstring += ' ' + ' ' * (length - len(valstring) - 2) + ' '
        return valstring + '# %s' % self.description

    def __eq__(self, teststring):
        """ Determines if a variable string matches this variable """
        return self.name == teststring

    def __ne__(self, teststring):
        """ Not equal """
        return not self.__eq__(teststring)

    def SetValue(self, value):
        """ Sets the value of the variable """
        if self.datatype is str:
            self.value = value.replace('"', '').replace("'", '')
        elif self.datatype in [list, tuple]:
            data = value.replace('"', '').replace("'", '')
            self.value = [self.int_datatype(x.strip()) for x in re.split(';\s*|,\s*', data)]
        else:
            self.value = self.datatype(value)


class Namelist(object):
    """ Sets up a namelist. This holds many different Variables, and these
       variables will only be recognized when parsing this particular namelist.
       Set up to mimic the behavior of a Fortran namelist (so that the input is
       similar to the rest of Amber). Some of the known deficiencies:

         o the body of the namelist cannot start on the same line as the start
           or end of the namelist

         o the end of the namelist must be &end or / and must be on its own line

         o It will not (yet) recognize array lengths -- those are simply parsed
           as strings
   """

    def __init__(self, trigger, full_name, to_match=3):
        """ Sets up the list of variables, which is just the trigger for now. The
          trigger is a logical variable that gets set to true if this namelist
          is parsed. Open() trips the trigger if it exists. It can be passed in
          as anything that evaluates to False if it doesn't exist.
      """
        self.trigger = trigger
        self.variables = {}
        if self.trigger is not None:
            self.variables = {self.trigger: False}
        self.open = False
        self.full_name = full_name
        self.to_match = to_match

    def __eq__(self, nml):
        """ A namelist is equal if the name matches properly """
        return nml == self.full_name[:len(nml)] and len(nml) >= min(self.to_match, len(self.full_name))

    def __ne__(self, nml):
        """ Not equal """
        return not self.__eq__(nml)

    def addVariable(self, varname, datatype, default=None, description=None, int_dat_type=str):
        """ Adds a variable to this namelist. It checks to make sure that it's
          going to create a conflict with an existing variable.
        """
        if varname in list(self.variables.keys()):
            raise InternalError('Duplicated variable %s in Namelist' % varname)
        self.variables[varname] = Variable(varname, datatype, default, description, int_dat_type)

    def Open(self):
        """ Signifies that the namelist is open """
        if self.open:
            raise InputError('Namelist already open. Cannot open before closing')

        if self.trigger: self.variables[self.trigger] = True
        self.open = True

    def __str__(self):
        """
        Prints out the full contents of this namelist in the Fortran namelist
        format
        """
        retstr = '&%s\n' % self.full_name
        for variable in self.variables:
            if variable is self.trigger: continue
            retstr += '  %s\n' % self.variables[variable].help_str()
        return retstr + '/'


class InputFile(object):
    """ Defines the Input File and parses it. You have to add stuff to the parser
       via addNamelist. Use it as follows:

       input = InputFile()

       input.addNamelist('gb', 'gb', [['saltcon', float, 0], ...],
                         trigger='gbrun')
       input.addNamelist('ala', 'alanine_scanning', [['mutant', int, 0]],
                         trigger='alarun')

       INPUT = input.Parse('mmpbsa.in')
   """

    def __init__(self):
        """ Initializes the input file, sets up empty arrays/dictionaries """
        self.ordered_namelist_keys = []
        self.namelists = {}
        self.text = ''  # text of the given input file

    def __str__(self):
        """ Prints out the input file """
        if not self.text:
            return

        ret_text = self.text.replace('\n', '\n|')  # Add | to start of each line

        return ('|Input file:\n|---------------------------------------' +
                '-----------------------\n|' + ret_text +
                '-----------------------------------------------------' +
                '---------\n')

    def print_contents(self, destination, calc_list=None):
        """ Prints the contents of the input file """
        # Open a file to write to if need be
        # section description
        sd = {'general': '# General namelist variables',
              'gb': '# (AMBER) Generalized-Born namelist variables',
              'pb': '# (AMBER) Possion-Boltzmann namelist variables',
              'rism': '# 3D-RISM namelist variables',
              'decomp': '# Decomposition namelist variables',
              'ala': '# Alanine scanning namelist variables',
              'nmode': '# Normal Modes Entropy namelist variables'}

        dest = destination if hasattr(destination, 'write') else open(destination, 'w')
        if calc_list:
            dest.write(f'Input file generated by gmx_MMPBSA ({__version__})\n'
                       f'Be careful with the variables you modify, some can have severe consequences on the results '
                       f'you obtain.\n\n')
        for namelist in self.ordered_namelist_keys:
            if calc_list and namelist in calc_list or not calc_list:
                dest.write(f'{sd[namelist]}\n')
                dest.write('%s\n\n' % self.namelists[namelist])

        # Close the file if we opened it.
        if dest is not destination:
            dest.close()

    def addNamelist(self, name, full_name, variable_list, trigger=None):
        """ Adds a namelist to the input file that will be parsed. Variable list
          must be an array of arrays. Each array element must be an array that
          has the information [varname, datatype, default, chars to match]. The
          'trigger' is the logical variable that gets set to true if this
          namelist is specified in the input file.
        """

        if name in self.ordered_namelist_keys:
            raise InputError('Namelist %s defined multiple times' % name)

        self.ordered_namelist_keys.append(name)
        self.namelists[name] = Namelist(trigger, full_name)

        for var in variable_list:

            if not isinstance(var, (list, tuple)) or len(var) not in [4, 5]:
                raise InputError('variables in variable_list must be lists of ' +
                                 'length 4 or 5. [varname, datatype, default, description, internal_datatype ('
                                 'Optional)]')
            if len(var) == 4:
                self.namelists[name].addVariable(var[0], var[1], var[2], var[3])
            else:
                self.namelists[name].addVariable(var[0], var[1], var[2], var[3], var[4])

    def _full_namelist_name(self, nml):
        """ Determines what the full namelist name is. We try to make as many
          allowances as possible. We will match the first 3 characters and
          replace all _'s with
      """
        nml = nml.replace(' ', '_')  # replaces spaces with _'s
        for key in self.ordered_namelist_keys:
            if self.namelists[key] == nml: return key

        raise InputError('Unrecognized namelist %s' % nml)

    def Parse(self, filename):
        """
        This subroutine parses the input file. Only data in namelists are
          parsed, and all namelists must be set prior to calling this routine.

          It will create a dictionary of Input variables for all variables in
          all namelists. They all flood the same namespace. If there are any
          conflicts between variables in namelists, an error will be raised.
          Make sure all input variables are unique!
        """
        from pathlib import Path

        # Make sure our file exists

        if filename is None:
            raise InputError("No input file was provided!")
        if not Path(filename).exists():
            raise InputError("Can't find input file (%s)" % filename)

        # Load the whole thing into memory. This should be plenty short enough.
        lines = open(filename, 'r').readlines()

        # Save the text of the input file so we can echo it back later
        self.text = ''.join(lines)

        # We will loop through the input file three times:
        #
        # 1st: Load all of the data into an array (namelist_fields)
        # 2nd: Combine multi-element values (to allow commas in input variables)
        # 3rd: Loop through the fields to change the values of the variables.

        declared_namelists = []  # keep track of the namelists we've found so far
        namelist_fields = []  # entries in a given namelist
        innml = False  # are we in a namelist now? Don't enter multiple

        # split up the input file into separate fields by comma

        for line in lines:
            # Skip title lines (we are flexible here) and comments
            if not innml and not line.strip().startswith('&'):
                continue
            if line.strip().startswith('#') or line.strip().startswith('!'):
                continue

            # Catch some errors
            if innml and line.strip().startswith('&'):
                raise InputError('Invalid input. Terminate each namelist prior to starting another one.')

            # End of a namelist
            elif innml and line.strip() in ['/', '&end']:
                innml = False

            # Now if we finally find a namelist
            elif not innml and line.strip().startswith('&'):
                innml = True
                namelist = line.strip()[1:].lower()
                namelist = self._full_namelist_name(namelist)

                if namelist in declared_namelists:
                    raise InputError('Namelist %s specified multiple times' %
                                     namelist)

                self.namelists[namelist].Open()
                declared_namelists.append(namelist)
                namelist_fields.append([])

            # We are in a namelist here, now fill in the fields
            elif innml:
                line = line[:line.strip().index('#')] if '#' in line else line
                items = line.strip().split(',')
                # Screen any blank fields
                j = 0
                while j < len(items):
                    items[j] = items[j].strip()
                    if len(items[j]) == 0:
                        items.pop(j)
                    else:
                        j += 1
                namelist_fields[-1].extend(items)
            # end if [elif innml]
        # # end for line in lines

        # # Combine any multi-element fields into the last field that has a = in it
        begin_field = -1
        for i in range(len(namelist_fields)):
            for j in range(len(namelist_fields[i])):
                if '=' in namelist_fields[i][j]:
                    begin_field = j
                elif begin_field == -1:
                    raise ('Invalid input file! Error reading namelist %s' % declared_namelists[i])
                else:
                    namelist_fields[i][begin_field] += ',%s' % namelist_fields[i][j]
        # Now parse through the items to add them to the master dictionary. Note
        # that thanks to the last step, all data in namelist_fields will be
        # contained within fields that have a '='. All others can be ignored
        for i in range(len(namelist_fields)):
            for j in range(len(namelist_fields[i])):
                if '=' not in namelist_fields[i][j]:
                    continue
                var = namelist_fields[i][j].split('=')
                var[0] = var[0].strip()
                var[1] = var[1].strip()

                # Now we have to loop through all variables in that namelist to
                # see if this is the variable we want.
                found = False
                for key in list(self.namelists[declared_namelists[i]].variables.keys()):
                    if self.namelists[declared_namelists[i]].variables[key] == var[0]:
                        self.namelists[declared_namelists[i]].variables[key].SetValue(var[1])
                        found = True
                        break

                if not found:
                    raise InputError('Unknown variable %s in &%s' % (var[0], declared_namelists[i]))
        # Now it's time to fill the INPUT dictionary
        INPUT = {}
        for nml in self.ordered_namelist_keys:
            for var in list(self.namelists[nml].variables.keys()):
                # Here, the triggers are just bool types, so protect from accessing
                # an attribute that doesn't exist! We only allow Variable types and
                # bool types
                var_object = self.namelists[nml].variables[var]
                try:
                    INPUT[var] = self.namelists[nml].variables[var].value
                except AttributeError:
                    if isinstance(var_object, bool):
                        INPUT[var] = var_object
                    else:
                        raise InputError('Disallowed namelist variable type')
        return INPUT


# Define the MM/PBSA input file here
input_file = InputFile()

# Add namelists with a list of variables. The variables are added to the
# namelists in lists. The entries are:
# [<variable name> <variable type> <default value> <# of characters to match>]

input_file.addNamelist('general', 'general',
                       [
                            # Basic options
                           ['sys_name', str, '', 'System name'],
                           ['startframe', int, 1, 'First frame to analyze'],
                           ['endframe', int, 9999999, 'Last frame to analyze'],
                           ['interval', int, 1, 'Number of frames between adjacent frames analyzed'],

                            # Parameters options
                           ['forcefields', list, 'oldff/leaprc.ff99SB, leaprc.gaff', 'Define the force field to build '
                                                                                     'the Amber topology'],
                           ['ions_parameters', int, 1, 'Define ions parameters to build the Amber topology'],
                           ['PBRadii', int, 3, 'Define PBRadii to build amber topology from GROMACS files'],
                           ['temperature', float, 298.15, 'Temperature'],

                            # Entropy options
                           ['qh_entropy', int, 0, 'Do quasi-harmonic calculation'],
                           ['interaction_entropy', int, 0, 'Do Interaction Entropy calculation'],
                           ['ie_segment', int, 25, 'Trajectory segment to calculate interaction entropy'],
                           ['c2_entropy', int, 0, 'Do C2 Entropy calculation'],

                            # Miscellaneous options
                           ['assign_chainID', int, 0, 'Assign chains ID'],
                           ['exp_ki', list, [0.0], 'Experimental Ki in nM', float],
                           ['full_traj', int, 0, 'Print a full traj. AND the thread trajectories'],
                           ['gmx_path', str, '', 'Force to use this path to get GROMACS executable'],
                           ['keep_files', int, 2, 'How many files to keep after successful completion'],

                           ['netcdf', int, 0, 'Use NetCDF intermediate trajectories'],
                           ['solvated_trajectory', int, 1, 'Define if it is necessary to cleanup the trajectories'],
                           ['verbose', int, 1, 'How many energy terms to print in the final output']
                       ], trigger=None)

input_file.addNamelist('gb', 'gb',
                       [
                           ['igb', int, 5, 'GB model to use'],
                           ['intdiel', float, 1.0, 'Internal dielectric constant for sander'],
                           ['extdiel', float, 78.5, 'External dielectric constant for sander'],

                           ['saltcon', float, 0, 'Salt concentration (M)'],
                           ['surften', float, 0.0072, 'Surface tension'],
                           ['surfoff', float, 0.0, 'Surface tension offset'],
                           ['molsurf', int, 0, 'Use Connelly surface (\'molsurf\' program)'],
                           ['msoffset', float, 0.0, 'Offset for molsurf calculation'],
                           ['probe', float, 1.4, 'Solvent probe radius for surface area calc'],

                            # Options for QM
                           ['ifqnt', int, 0, 'Use QM on part of the system'],
                           ['qm_theory', str, '', 'Semi-empirical QM theory to use'],
                           ['qm_residues', str, '', 'Residues to treat with QM'],

                           # TODO: deprecated since 1.5.0. Automatic charge assignment
                           ['qmcharge_com', int, 0, 'Charge of QM region in complex'],
                           ['qmcharge_lig', int, 0, 'Charge of QM region in ligand'],
                           ['qmcharge_rec', int, 0, 'Charge of QM region in receptor'],
                           ['qmcut', float, 9999, 'Cutoff in the QM region'],
                           ['scfconv', float, 1.0e-8, 'Convergence criteria for the SCF calculation, in kcal/mol'],
                           ['peptide_corr', int, 0, 'Apply MM correction to peptide linkages'],
                           ['writepdb', int, 1, 'Write a PDB file of the selected QM region'],
                           ['verbosity', int, 0, 'Controls the verbosity of QM/MM related output'],

                           # Options for alpb
                           ['alpb', int, 0, 'Use Analytical Linearized Poisson-Boltzmann (ALPB)'],
                           ['arad_method', int, 1, 'Selected method to estimate the effective electrostatic size']
                       ], trigger='gbrun')

input_file.addNamelist('pb', 'pb',
                       [
                            # Basic input options
                           ['ipb', int, 2, 'Dielectric model for PB'],
                           ['inp', int, 2, 'Nonpolar solvation method'],
                           ['sander_apbs', int, 0, 'Use sander.APBS?'],

                           # Options to define the physical constants
                           ['indi', float, 1, 'Internal dielectric constant'],
                           ['exdi', float, 80, 'External dielectric constant'],
                           ['emem', float, 4.0, 'Membrane dielectric constant'],
                           ['smoothopt', int, 1, 'Set up dielectric values for finite-difference grid edges that are '
                                                 'located across the solute/solvent dielectric boundary'],
                           ['istrng', float, 0.0, 'Ionic strength (M)'],
                           ['radiopt', int, 1, 'Use optimized radii?'],
                           ['prbrad', float, 1.4, 'Probe radius'],
                           ['iprob', float, 2.0, 'Mobile ion probe radius (Angstroms) for ion accessible surface used '
                                                 'to define the Stern layer'],
                           ['sasopt', int, 0, 'Molecular surface in PB implict model'],
                           ['arcres', float, 0.25, 'The resolution (Å) to compute solvent accessible arcs'],

                           # Options for Implicit Membranes
                           ['memopt', int, 0, 'Use PB optimization for membrane'],
                           ['mprob', float, 2.70, 'Membrane probe radius in Å'],
                           ['mthick', float, 40.0, 'Membrane thickness'],
                           ['mctrdz', float, 0.0, 'Distance to offset membrane in Z direction'],
                           ['poretype', int, 1, 'Use exclusion region for channel proteins'],

                           # Options to select numerical procedures
                           ['npbopt', int, 0, 'Use NonLinear PB solver?'],
                           ['solvopt', int, 1, 'Select iterative solver'],
                           ['accept', float, 0.001, 'Sets the iteration convergence criterion (relative to the initial '
                                                    'residue)'],
                           ['linit', int, 1000, 'Number of SCF iterations'],
                           ['fillratio', float, 4, 'Ratio between the longest dimension of the rectangular '
                                                   'finite-difference grid and that of the solute'],
                           ['scale', float, 2.0, '1/scale = grid spacing for the finite difference solver (default = '
                                                 '1/2 Å)'],
                           ['nbuffer', float, 0, 'Sets how far away (in grid units) the boundary of the finite '
                                                 'difference grid is away from the solute surface'],
                           ['nfocus', int, 2, 'Electrostatic focusing calculation'],
                           ['fscale', int, 8, 'Set the ratio between the coarse and fine grid spacings in an '
                                              'electrostatic focussing calculation'],
                           ['npbgrid', int, 1, 'Sets how often the finite-difference grid is regenerated'],

                            # Options to compute energy and forces
                           ['bcopt', int, 5, 'Boundary condition option'],
                           ['eneopt', int, 2, 'Compute electrostatic energy and forces'],
                           ['frcopt', int, 0, 'Output for computing electrostatic forces'],
                           ['scalec', int, 0, 'Option to compute reaction field energy and forces'],
                           ['cutfd', float, 5.0, 'Cutoff for finite-difference interactions'],
                           ['cutnb', float, 0.0, 'Cutoff for nonbonded interations'],
                           ['nsnba', int, 1, 'Sets how often atom-based pairlist is generated'],

                            # Options to select a non-polar solvation treatment
                           ['decompopt', int, 2, 'Option to select different decomposition schemes when INP = 2'],
                           ['use_rmin', int, 1, 'The option to set up van der Waals radii'],
                           ['sprob', float, 0.557, 'Solvent probe radius for SASA used to compute the dispersion term'],
                           ['vprob', float, 1.300, 'Solvent probe radius for molecular volume (the volume enclosed by '
                                                   'SASA)'],
                           ['rhow_effect', float, 1.129, 'Effective water density used in the non-polar dispersion '
                                                         'term calculation'],
                           ['use_sav', int, 1, 'Use molecular volume (the volume enclosed by SASA) for cavity term '
                                               'calculation'],
                           ['cavity_surften', float, 0.0378, 'Surface tension'],
                           ['cavity_offset', float, -0.5692, 'Offset for nonpolar solvation calc'],
                           ['maxsph', int, 400, 'Approximate number of dots to represent the maximum atomic solvent '
                                                'accessible surface'],
                           ['maxarcdot', int, 1500, 'Number of dots used to store arc dots per atom'],

                           # Options for output
                           ['npbverb', int, 0, 'Option to turn on verbose mode']
                       ], trigger='pbrun')

input_file.addNamelist('rism', 'rism',
                       [
                           ['closure', list, ['kh'], 'Closure equation to use'],
                           ['gfcorrection', int, 0, 'Compute the Gaussian fluctuation excess chemical potential '
                                                    'functional'],
                           ['pcpluscorrection', int, 0, 'Compute the PC+/3D-RISM excess chemical potential functional'],
                           ['noasympcorr', int, 1, 'Turn off long range asymptotic corrections for thermodynamic '
                                                   'output only'],
                           ['buffer', float, 14, 'Distance between solute and edge of grid'],
                           ['solvcut', float, -1, 'Cutoff of the box'],
                           ['grdspc', list, [0.5, 0.5, 0.5], 'Grid spacing', float],
                           ['ng', list, [-1, -1, -1], 'Number of grid points', int],
                           ['solvbox', list, [-1, -1, -1], 'Box limits', int],
                           ['tolerance', list, [1.0e-5], 'Convergence tolerance', float],
                           ['ljTolerance', float, -1.0, 'Determines the Lennard-Jones cutoff distance based on the '
                                                        'desired accuracy of the calculation'],
                           ['asympKSpaceTolerance', float, -1.0, 'Determines the reciprocal space long range '
                                                                 'asymptotics cutoff distance based on the desired '
                                                                 'accuracy of the calculation'],
                           ['treeDCF', int, 1, 'Use the treecode approximation to calculate the direct '
                                               'correlation function (DCF) long-range asymptotic correction'],
                           ['treeTCF', int, 1, 'Use the treecode approximation to calculate the total '
                                               'correlation function (TCF) long-range asymptotic correction'],
                           ['treeCoulomb', int, 0, 'Use direct sum or the treecode approximation to calculate the '
                                                   'Coulomb potential energy'],
                           ['treeDCFMAC', float, 0.1, 'Treecode multipole acceptance criterion for the DCF long-range '
                                                      'asymptotic correction'],
                           ['treeTCFMAC', float, 0.1, 'Treecode multipole acceptance criterion for the TCF long-range '
                                                      'asymptotic correction'],
                           ['treeCoulombMAC', float, 0.1, 'Treecode multipole acceptance criterion for the Coulomb '
                                                          'potential energy'],
                           ['treeDCFOrder', int, 2, 'Treecode Taylor series order for the DCF long-range asymptotic '
                                                    'correction'],
                           ['treeTCFOrder', int, 2, 'Treecode Taylor series order for the TCF long-range asymptotic '
                                                    'correction'],
                           ['treeCoulombOrder', int, 2, 'Treecode Taylor series order for the Coulomb potential '
                                                        'energy'],
                           ['treeDCFN0', int, 500, 'Maximum number of grid points contained within the treecode leaf '
                                                   'clusters for the DCF'],
                           ['treeTCFN0', int, 500, 'Maximum number of grid points contained within the treecode leaf '
                                                   'clusters for the  TCF'],
                           ['treeCoulombN0', int, 500, 'Maximum number of grid points contained within the treecode '
                                                       'leaf clusters for the Coulomb potential energy'],
                           ['mdiis_del', float, 0.7, 'MDIIS step size'],
                           ['mdiis_nvec', int, 5, 'Number of previous iterations MDIIS uses to predict a new solution'],
                           ['mdiis_restart', float, 10.0, 'Use lowest residual solution in memory if '
                                                          'current residual is mdiis_restart times larger than '
                                                          'the smallest residual in memory'],
                           ['maxstep', int, 10000, 'Maximum number of iterative steps per solution'],
                           ['npropagate', int, 5, 'Number of previous solutions to use in predicting a new solution'],
                           ['polardecomp', int, 0, 'Break solv. energy into polar and nonpolar terms'],
                           # TODO: work with entropicDecomp? need more tests...
                           ['entropicdecomp', int, 0, 'Decomposes solvation free energy into energy and entropy '
                                                      'components'],
                           # ['centering', int, 1, 'Select how solute is centered in the solvent box'],
                           ['rism_verbose', int, 0, 'Control how much 3D-RISM info to print']
                       ], trigger='rismrun')

input_file.addNamelist('ala', 'alanine_scanning',
                       [
                           ['mutant_res', str, '', 'Which residue will be mutated'],
                           ['mutant', str, 'ALA', 'Defines if Alanine or Glycine scanning will be performed'],
                           ['mutant_only', int, 0, 'Only compute mutant energies'],
                           ['cas_intdiel', int, 0, 'Change the intdiel value based on which aa is mutated'],
                           ['intdiel_nonpolar', int, 1, 'intdiel for nonpolar residues'],
                           ['intdiel_polar', int, 3, 'intdiel for polar residues'],
                           ['intdiel_positive', int, 5, 'intdiel for positive charged residues'],
                           ['intdiel_negative', int, 5, 'intdiel for negative charged residues']
                       ], trigger='alarun')

input_file.addNamelist('decomp', 'decomposition',
                       [
                           ['idecomp', int, 0, 'Which type of decomposition analysis to do'],
                           ['dec_verbose', int, 0, 'Control energy terms are printed to the output'],
                           ['print_res', str, 'within 6', 'Which residues to print decomposition data for'],
                           ['csv_format', int, 1, 'Write decomposition data in CSV format']
                       ], trigger='decomprun')

input_file.addNamelist('nmode', 'nmode',
                       [
                            # Basic Options
                           ['nmstartframe', int, 1, 'First frame to analyze for normal modes'],
                           ['nmendframe', int, 1000000, 'Last frame to analyze for normal modes'],
                           ['nminterval', int, 1, 'Interval to take snapshots for normal mode analysis'],
                            # Parameters options
                           ['nmode_igb', int, 1, 'GB model to use for normal mode calculation'],
                           ['nmode_istrng', float, 0, 'Ionic strength for GB model (M)'],
                           ['dielc', float, 1, 'Dielectric constant'],
                            # Minimization options
                           ['drms', float, 0.001, 'Minimization gradient cutoff'],
                           ['maxcyc', int, 10000, 'Maximum number of minimization cycles'],
                       ], trigger='nmoderun')
