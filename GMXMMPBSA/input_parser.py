"""
This is a module that contains functions responsible for parsing the
input file for gmx_MMPBSA. It must be included with gmx_MMPBSA to
ensure proper functioning.
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

from GMXMMPBSA.exceptions import InputError, InternalError
from GMXMMPBSA import __version__
import re


DEFAULT_QM_THEORY = 'PM6-DH+'
SUPPORTED_QM_THEORIES = (
    'PM3', 'AM1', 'RM1', 'MNDO', 'PM3-PDDG', 'PM3-PDDG_08', 'MNDO-PDDG',
    'PM3-CARB1', 'PM3-ZNB', 'PM3-MAIS', 'AM1-D*', 'AM1-DH+', 'MNDO/D',
    'AM1/D', 'PM6', 'PM6-D', 'PM6-DH+', 'DFTB', 'DFTB2', 'DFTB3',
)


PB_MEMBRANE_TEMPLATE = {
    'memopt': 1,  #Use a heterogeneous membrane dielectric constant in a slab-like implicit membrane
    'emem': 7.0,
    'indi': 4.0,
    'mctrdz': 'automatic',
    'mthick': 'automatic',
    'poretype': 1,
    'radiopt': 0,
    'istrng': 0.150,
    'fillratio': 1.25,
    'inp': 2,
    'sasopt': 0,
    'solvopt': 2,
    'ipb': 1,
    'bcopt': 10,
    'nfocus': 1,
    'linit': 1000,
    'eneopt': 1,
    'cutfd': 7.0,
    'cutnb': 99.0,
    'maxarcdot': 15000,
    'npbverb': 1,
    # Automatic membrane detection uses phosphorus unless the user changes it.
    'membrane_atoms': 'P',
}


_USE_CURRENT_VALUE = object()


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
        self.allow_none = default is None and self.datatype is not str
        if default is None:
            self.value = None
        elif self.datatype is str:
            self.value = default.replace("'", '').replace('"', '')
        elif self.datatype in [list, tuple]:
            if isinstance(default, str):
                self.value = [self.int_datatype(x.strip()) for x in re.split(r';\s*|,\s*', default.replace('"',''))]
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

    def help_str(self, value=_USE_CURRENT_VALUE):
        """ returns the string [<name> = <value>.... # description] """
        name_width = 30
        if value is _USE_CURRENT_VALUE:
            value = self.value
        # ``None`` represents an optional value that should be left to the
        # downstream program.  It is not a valid Fortran namelist literal for
        # numeric variables, so keep it visible in generated templates but
        # comment it out instead of emitting e.g. ``ndiis_attempts = None``.
        if value is None:
            valstring = f'# {self.name:{name_width}s} = None'
        elif self.datatype is str:
            valstring = f'{self.name:{name_width}s} = "{value:s}"'
        elif self.datatype in [list, tuple]:
            v = ','.join(map(str, value))
            if self.int_datatype == str:
                valstring = f'{self.name:{name_width}s} = "{v}"'
            else:
                valstring = f'{self.name:{name_width}s} = {v}'
        else:
            valstring = f'{self.name:{name_width}s} = {value}'
        length = 72
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
        value = value.strip()
        if self.allow_none and value.lower() == 'none':
            self.value = None
            return

        if self.datatype is str:
            self.value = value.replace('"', '').replace("'", '')
        elif self.datatype in [list, tuple]:
            data = value.replace('"', '').replace("'", '')
            self.value = [self.int_datatype(x.strip()) for x in re.split(r';\s*|,\s*', data)]
        else:
            try:
                self.value = self.datatype(value)
            except (TypeError, ValueError) as exc:
                raise InputError(
                    f'Invalid value {value!r} for {self.name}; expected {self.datatype.__name__}'
                ) from exc


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
        if varname in self.variables:
            raise InternalError(f'Duplicated variable {varname} in Namelist')
        self.variables[varname] = Variable(varname, datatype, default, description, int_dat_type)

    def Open(self):
        """ Signifies that the namelist is open """
        if self.open:
            raise InputError('Namelist already open. Cannot open before closing')

        if self.trigger: self.variables[self.trigger] = True
        self.open = True

    def contents(self, overrides=None):
        """
        Prints out the full contents of this namelist in the Fortran namelist
        format
        """
        overrides = overrides or {}
        retstr = '&%s\n' % self.full_name
        for variable in self.variables:
            if variable is self.trigger: continue
            value = overrides.get(variable, _USE_CURRENT_VALUE)
            retstr += '  %s\n' % self.variables[variable].help_str(value)
        return f'{retstr}/'

    def __str__(self):
        return self.contents()


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
              'gbnsr6': '# GBNSR6 namelist variables',
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
        membrane_template = calc_list and 'pb_mem' in calc_list
        for namelist in self.ordered_namelist_keys:
            selected = (not calc_list or namelist in calc_list or
                        (namelist == 'pb' and membrane_template))
            if selected:
                dest.write(f'{sd[namelist]}\n')
                if namelist == 'pb' and membrane_template:
                    contents = self.namelists[namelist].contents(PB_MEMBRANE_TEMPLATE)
                else:
                    contents = str(self.namelists[namelist])
                dest.write('%s\n\n' % contents)

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

        # Parse marks namelists open while reading; always release that state so
        # shared InputFile instances remain reusable after success or failure.
        try:
            # split up the input file into separate fields by comma

            for line in lines:
                # Skip title lines (we are flexible here) and comments
                if not innml and not line.strip().startswith('&'):
                    continue
                if line.strip().startswith('#') or line.strip().startswith('!'):
                    continue

                # Catch some errors
                if innml and line.strip().startswith('&') and line.strip() != '&end':
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
                        raise InputError('Namelist %s specified multiple times' % namelist)

                    self.namelists[namelist].Open()
                    declared_namelists.append(namelist)
                    namelist_fields.append([])

                # We are in a namelist here, now fill in the fields
                elif innml:
                    line = line[:line.strip().index('#')] if '#' in line else line.strip('\n')
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
            # # Combine any multi-element fields into the last field that has a = in it
            begin_field = -1
            for i, _ in enumerate(namelist_fields):
                for j, _ in enumerate(namelist_fields[i]):
                    if '=' in namelist_fields[i][j]:
                        begin_field = j
                    elif begin_field == -1:
                        raise InputError(f'Invalid input file! Error reading namelist {declared_namelists[i]}')
                    else:
                        namelist_fields[i][begin_field] += f',{namelist_fields[i][j]}'
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
                    for key in self.namelists[declared_namelists[i]].variables:
                        if self.namelists[declared_namelists[i]].variables[key] == var[0]:
                            self.namelists[declared_namelists[i]].variables[key].SetValue(var[1])
                            found = True
                            break

                    if not found:
                        raise InputError(f'Unknown variable {var[0]} in &{declared_namelists[i]}')
            # Now it's time to fill the INPUT dictionary
            INPUT = {}
            for nml in self.ordered_namelist_keys:
                INPUT[nml] = {}
                for var in self.namelists[nml].variables:
                    # Here, the triggers are just bool types, so protect from accessing
                    # an attribute that doesn't exist! We only allow Variable types and
                    # bool types
                    var_object = self.namelists[nml].variables[var]
                    try:
                        INPUT[nml][var] = self.namelists[nml].variables[var].value
                    except AttributeError:
                        if isinstance(var_object, bool):
                            INPUT[nml][var] = var_object
                        else:
                            raise InputError('Disallowed namelist variable type')
            return INPUT
        finally:
            for namelist in self.namelists.values():
                namelist.open = False


# Define the MM/PBSA input file here
input_file = InputFile()

# Add namelists with a list of variables. The variables are added to the
# namelists in lists. The entries are:
# [<variable name> <variable type> <default value> <# of characters to match>]

input_file.addNamelist('general', 'general',
                       [
                           # Basic options
                           ['sys_name', str, '', 'System name; e.g. "complex"'],
                           ['startframe', int, 1, 'First frame; e.g. 1'],
                           ['endframe', int, 9999999, 'Last frame; e.g. 100'],
                           ['interval', int, 1, 'Frame interval; e.g. 1'],

                           # Parameters options
                           ['forcefields', list, 'oldff/leaprc.ff99SB, leaprc.gaff',
                            'Force fields; e.g. "leaprc.protein.ff14SB"'],
                           ['ions_parameters', int, 1, 'Ion params; e.g. 1'],
                           ['PBRadii', int, 4, 'PB radii set; 1-7'],
                           ['temperature', float, 298.15, 'Temperature (K); e.g. 298.15'],

                           # Entropy options
                           ['qh_entropy', int, 0, 'Legacy QH output reader; new calculations reject 1'],
                           ['interaction_entropy', int, 0, 'Run IE entropy; 0/1'],
                           ['ie_segment', int, 25, 'IE segment length (%); e.g. 25'],
                           ['c2_entropy', int, 0, 'Run C2 entropy; 0/1'],

                           # Miscellaneous options
                           ['assign_chainID', int, 0, 'Assign chain IDs; 0/1'],
                           ['exp_ki', list, [0.0], 'Experimental Ki (nM); e.g. 0.0', float],
                           ['full_traj', int, 0, 'Write full trajectory; 0/1'],
                           ['gmx_path', str, '', 'GROMACS path; e.g. "/usr/bin"'],
                           ['keep_files', int, 2, 'Files to keep; 0-2'],

                           ['netcdf', int, 0, 'Use NetCDF; 0/1'],
                           ['solvated_trajectory', int, 1, 'Clean solvated traj.; 0/1'],
                           ['explicit_waters', int, 0, 'Explicit waters; e.g. 10'],
                           ['explicit_waters_mask', str, '',
                            'Water reference; e.g. ":1-10", "within 4", "dASA"'],
                           ['explicit_waters_group', str, '',
                            'Solvent group; e.g. "TIP3"'],
                           ['explicit_waters_dasa_cutoff', float, 0.5,
                            'dASA cutoff; e.g. 0.5'],
                           ['explicit_waters_as', str, 'receptor', 'Water owner; e.g. "receptor"'],
                           ['explicit_waters_extra_points', str, 'error',
                            'Virtual sites; "error" or "strip"'],
                           ['verbose', int, 1, 'Output verbosity; 0-2']
                       ], trigger=None)

input_file.addNamelist('gb', 'gb',
                       [
                           ['igb', int, 8, 'GB model, e.g. 2 or 8'],
                           ['intdiel', float, 1.0, 'Internal dielectric; e.g. 1.0'],
                           ['extdiel', float, 78.5, 'External dielectric; e.g. 78.5'],

                           ['saltcon', float, 0, 'Salt conc. (M); e.g. 0.150'],
                           ['surften', float, 0.0072, 'Surface tension; e.g. 0.0072'],
                           ['surfoff', float, 0.0, 'Surface offset; e.g. 0.0'],
                           ['molsurf', int, 0, 'Use molsurf; 0/1'],
                           ['msoffset', float, 0.0, 'Molsurf offset; e.g. 0.0'],
                           ['probe', float, 1.4, 'Probe radius (A); e.g. 1.4'],

                            # Options for QM
                           ['ifqnt', int, 0, 'Enable QM/MM; 0/1'],
                           ['qm_theory', str, DEFAULT_QM_THEORY, 'QM theory; e.g. "PM6-DH+"'],
                           ['qm_residues', str, '', 'QM residues; e.g. ":1-5"'],

                           ['com_qmmask', str, '', 'Complex QM mask; e.g. ":1-5"'],
                           ['rec_qmmask', str, '', 'Receptor QM mask; e.g. ":1-5"'],
                           ['lig_qmmask', str, '', 'Ligand QM mask; e.g. ":1"'],

                           # deprecated since 1.5.0. Automatic charge assignment except when using user defined masks
                           ['qmcharge_com', int, 0, 'Complex QM charge; e.g. 0'],
                           ['qmcharge_lig', int, 0, 'Ligand QM charge; e.g. 0'],
                           ['qmcharge_rec', int, 0, 'Receptor QM charge; e.g. 0'],

                           ['qmcut', float, 9999, 'QM cutoff (A); e.g. 9999'],
                           ['scfconv', float, 1.0e-8, 'SCF convergence; e.g. 1.0e-8'],
                           ['itrmax', int, 1000, 'Maximum SCF iterations; e.g. 5000'],
                           ['ndiis_attempts', int, None,
                            'Maximum DIIS attempts per SCF cycle; e.g. 700'],
                           ['peptide_corr', int, 0, 'Peptide correction; 0/1'],
                           ['writepdb', int, 1, 'Write QM PDB; 0/1'],
                           ['verbosity', int, 0, 'QM/MM verbosity; 0-5'],

                           # Options for alpb
                           ['alpb', int, 0, 'Use ALPB; 0/1'],
                           ['arad_method', int, 1, 'ALPB size method; e.g. 1']
                       ], trigger='gbrun')

input_file.addNamelist('gbnsr6', 'gbnsr6',
                       [
                           ['b', float, 0.028, 'GBNSR6 offset; e.g. 0.028'],
                           ['alpb', int, 1, 'Use ALPB; 0/1'],
                           ['epsin', float, 1.0, 'Solute dielectric; e.g. 1.0'],
                           ['epsout', float, 78.5, 'Solvent dielectric; e.g. 78.5'],
                           # FIXME: convert to M
                           ['istrng', float, 0.0, 'Ionic strength (M); e.g. 0.150'],
                           ['rs', float, 0.52, 'Boundary shift; e.g. 0.52'],
                           ['dprob', float, 1.4, 'Probe radius (A); e.g. 1.4'],
                           ['space', float, 0.5, 'Grid spacing (A); e.g. 0.5'],
                           ['arcres', float, 0.2, 'Arc resolution; e.g. 0.2'],
                           ['radiopt', int, 0, 'Radii option; e.g. 0'],
                           ['chagb', int, 0, 'Use CHAGB; 0/1'],
                           ['roh', int, 1, 'RzOH value; e.g. 1'],
                           ['tau', float, 1.47, 'CHAGB tau; e.g. 1.47'],
                           ['cavity_surften', float, 0.005, 'Cavity surften; e.g. 0.005'],
                       ], trigger='gbnsr6run')

input_file.addNamelist('pb', 'pb',
                       [
                           # Basic input options
                           ['ipb', int, 2, 'PB model; e.g. 2'],
                           ['inp', int, 1, 'Nonpolar method; 1 or 2'],
                           ['sander_apbs', int, 0, 'Use sander.APBS; 0/1'],

                           # Options to define the physical constants
                           ['indi', float, 1, 'Internal dielectric; e.g. 1.0'],
                           ['exdi', float, 78.5, 'External dielectric; e.g. 78.5'],
                           ['emem', float, 4.0, 'Membrane dielectric; e.g. 4.0'],
                           ['smoothopt', int, 1, 'Dielectric smoothing; 0-2'],
                           ['istrng', float, 0.0, 'Ionic strength (M); e.g. 0.150'],
                           ['radiopt', int, 1, 'Use optimized radii; 0/1'],
                           ['prbrad', float, 1.4, 'Probe radius (A); e.g. 1.4'],
                           ['iprob', float, 2.0, 'Ion probe (A); e.g. 2.0'],
                           ['sasopt', int, 0, 'PB surface option; 0/1'],
                           ['arcres', float, 0.25, 'Arc resolution (A); e.g. 0.25'],

                           # Options for Implicit Membranes
                           ['memopt', int, 0, 'Use membrane PB; 0/1'],
                           ['mprob', float, 2.70, 'Membrane probe (A); e.g. 2.7'],
                           ['mthick', str, 'automatic', 'Membrane thickness (A), or automatic'],
                           ['mctrdz', str, 'automatic', 'Membrane Z offset (A), or automatic'],
                           ['membrane_atoms', str, 'P',
                            'Atom names for automatic membrane parameters; semicolon-separated'],
                           ['poretype', int, 1, 'Pore type; 1 or 2'],

                           # Options to select numerical procedures
                           ['npbopt', int, 0, 'Use nonlinear PB; 0/1'],
                           ['solvopt', int, 1, 'PB solver; e.g. 1'],
                           ['accept', float, 0.001, 'Convergence; e.g. 0.001'],
                           ['linit', int, 1000, 'SCF iterations; e.g. 1000'],
                           ['fillratio', float, 4, 'Grid fill ratio; e.g. 4'],
                           ['scale', float, 2.0, 'Grid scale; e.g. 2'],
                           ['nbuffer', float, 0, 'Grid buffer; e.g. 0'],
                           ['nfocus', int, 2, 'Focus levels; e.g. 2'],
                           ['fscale', int, 8, 'Focus scale; e.g. 8'],
                           ['npbgrid', int, 1, 'Grid update freq.; e.g. 1'],

                           # Options to compute energy and forces
                           ['bcopt', int, 5, 'Boundary condition; e.g. 5'],
                           ['eneopt', int, 2, 'Energy option; e.g. 2'],
                           ['frcopt', int, 0, 'Force output; e.g. 0'],
                           ['scalec', int, 0, 'Reaction field option; e.g. 0'],
                           ['cutfd', float, 5.0, 'FD cutoff (A); e.g. 5'],
                           ['cutnb', float, 0.0, 'Nonbonded cutoff (A); e.g. 0'],
                           ['nsnba', int, 1, 'Pairlist frequency; e.g. 1'],

                           # Options to select a non-polar solvation treatment
                           ['decompopt', int, 2, 'Decomp scheme; 1 or 2'],
                           ['use_rmin', int, 1, 'Use Rmin radii; 0/1'],
                           ['sprob', float, 0.557, 'SASA probe (A); e.g. 0.557'],
                           ['vprob', float, 1.300, 'Volume probe (A); e.g. 1.3'],
                           ['rhow_effect', float, 1.129, 'Water density; e.g. 1.129'],
                           ['use_sav', int, 1, 'Use SAV cavity; 0/1'],
                           ['cavity_surften', float, 0.0378, 'Cavity surften; e.g. 0.0378'],
                           ['cavity_offset', float, -0.5692, 'Cavity offset; e.g. -0.5692'],
                           ['maxsph', int, 400, 'Max surface dots; e.g. 400'],
                           ['maxarcdot', int, 1500, 'Max arc dots; e.g. 1500'],

                           # Options for output
                           ['npbverb', int, 0, 'PB verbosity; 0/1']
                       ], trigger='pbrun')

input_file.addNamelist('rism', 'rism',
                       [
                           ['closure', list, ['kh'], 'Closure equation; e.g. "kh"'],
                           ['gfcorrection', int, 0, 'GF correction; 0/1'],
                           ['pcpluscorrection', int, 0, 'PC+ correction; 0/1'],
                           ['noasympcorr', int, 1, 'Disable asymptotic corr.; 0/1'],
                           ['buffer', float, 14, 'Grid buffer (A); e.g. 14'],
                           ['solvcut', float, -1, 'Solvent cutoff (A); e.g. -1'],
                           ['grdspc', list, [0.5, 0.5, 0.5], 'Grid spacing; e.g. 0.5,0.5,0.5', float],
                           ['ng', list, [-1, -1, -1], 'Grid points; e.g. -1,-1,-1', int],
                           ['solvbox', list, [-1, -1, -1], 'Solvent box; e.g. -1,-1,-1', int],
                           ['tolerance', list, [1.0e-5], 'Convergence tol.; e.g. 1.0e-5', float],
                           ['ljTolerance', float, -1.0, 'LJ tolerance; e.g. -1.0'],
                           ['asympKSpaceTolerance', float, -1.0, 'K-space tolerance; e.g. -1.0'],
                           ['treeDCF', int, 1, 'Use DCF treecode; 0/1'],
                           ['treeTCF', int, 1, 'Use TCF treecode; 0/1'],
                           ['treeCoulomb', int, 0, 'Use Coulomb treecode; 0/1'],
                           ['treeDCFMAC', float, 0.1, 'DCF MAC; e.g. 0.1'],
                           ['treeTCFMAC', float, 0.1, 'TCF MAC; e.g. 0.1'],
                           ['treeCoulombMAC', float, 0.1, 'Coulomb MAC; e.g. 0.1'],
                           ['treeDCFOrder', int, 2, 'DCF tree order; e.g. 2'],
                           ['treeTCFOrder', int, 2, 'TCF tree order; e.g. 2'],
                           ['treeCoulombOrder', int, 2, 'Coulomb tree order; e.g. 2'],
                           ['treeDCFN0', int, 500, 'DCF leaf size; e.g. 500'],
                           ['treeTCFN0', int, 500, 'TCF leaf size; e.g. 500'],
                           ['treeCoulombN0', int, 500, 'Coulomb leaf size; e.g. 500'],
                           ['mdiis_del', float, 0.7, 'MDIIS step size; e.g. 0.7'],
                           ['mdiis_nvec', int, 5, 'MDIIS vectors; e.g. 5'],
                           ['mdiis_restart', float, 10.0, 'MDIIS restart; e.g. 10.0'],
                           ['maxstep', int, 10000, 'Max iterations; e.g. 10000'],
                           ['npropagate', int, 5, 'Propagation history; e.g. 5'],
                           ['polardecomp', int, 0, 'Polar decomposition; 0/1'],
                           # TODO: work with entropicDecomp? need more tests...
                           ['entropicdecomp', int, 0, 'Entropic decomposition; 0/1'],
                           # ['centering', int, 1, 'Select how solute is centered in the solvent box'],
                           ['rism_verbose', int, 0, 'RISM verbosity; 0-2']
                       ], trigger='rismrun')

input_file.addNamelist('ala', 'alanine_scanning',
                       [
                           ['mutant_res', str, '', 'Residue to mutate; e.g. "A/23"'],
                           ['mutant', str, 'ALA', 'Mutation target; "ALA" or "GLY"'],
                           ['mutant_only', int, 0, 'Mutant energies only; 0/1'],
                           ['cas_intdiel', int, 0, 'Set intdiel by residue; 0/1'],
                           ['intdiel_nonpolar', int, 1, 'Nonpolar intdiel; e.g. 1'],
                           ['intdiel_polar', int, 3, 'Polar intdiel; e.g. 3'],
                           ['intdiel_positive', int, 5, 'Positive intdiel; e.g. 5'],
                           ['intdiel_negative', int, 5, 'Negative intdiel; e.g. 5']
                       ], trigger='alarun')

input_file.addNamelist('decomp', 'decomposition',
                       [
                           ['idecomp', int, 2, 'Decomp mode; 0-4'],
                           ['dec_verbose', int, 1, 'Decomp verbosity; 0-3'],
                           ['print_res', str, 'within 6',
                            'Residues to print; e.g. "all", "within 6", "A/2-10"'],
                           ['csv_format', int, 1, 'Write CSV output; 0/1']
                       ], trigger='decomprun')

input_file.addNamelist('nmode', 'nmode',
                       [
                           # Basic Options
                           ['nmstartframe', int, 1, 'First NM frame; e.g. 1'],
                           ['nmendframe', int, 1000000, 'Last NM frame; e.g. 100'],
                           ['nminterval', int, 1, 'NM frame stride; e.g. 1'],
                           # Parameters options
                           ['nmode_igb', int, 1, 'GB model, e.g. 1'],
                           ['nmode_istrng', float, 0, 'NM ionic strength (M); e.g. 0.0'],
                           ['dielc', float, 1, 'NM dielectric; e.g. 1.0'],
                           # Minimization options
                           ['drms', float, 0.001, 'Min. gradient cutoff; e.g. 0.001'],
                           ['maxcyc', int, 10000, 'Max minimization cycles; e.g. 10000'],
                       ], trigger='nmoderun')
