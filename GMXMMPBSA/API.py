"""
This module provides a means for users to take advantage of gmx_MMPBSA's parsing
ability. It exposes the free energy data (optionally to numpy arrays) so that
users can write a simple script to carry out custom data analyses, leveraging
the full power of Python's extensions, if they want (e.g., numpy, scipy, etc.)
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

from copy import deepcopy
from GMXMMPBSA import infofile, main, amber_outputs
from GMXMMPBSA.exceptions import SetupError, NoFileExists
from GMXMMPBSA.fake_mpi import MPI
import os
import warnings
# Try importing numpy to see if we have it available
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    from array import array
    HAS_NUMPY = False

__all__ = ['load_mmpbsa_info']

# Abstract the array-making routines
if HAS_NUMPY:
    make_array = lambda x: np.fromiter(x, float)
    make_array_len = lambda x: np.zeros(x, float)
else:
    make_array = lambda x: array('d', x)
    make_array_len = lambda x: array('d', [0 for i in range(x)])

class mmpbsa_data(dict):
    """ Main class that holds all of the Free Energy data """
    def __init__(self, app):
        """ Load data from an info object """
        global HAS_NUMPY
        if not isinstance(app, main.MMPBSA_App):
            raise TypeError('mmpbsa_data can only take an MMPBSA_App!')
        # Loop through all of the data
        if not hasattr(app, 'calc_types'):
            raise SetupError('Output files have not yet been parsed!')
        # Now load the data into the dict
        has_mutant = False
        # See if we are doing stability
        self.stability = app.stability
        # Now load the data
        for key in app.calc_types:
            if key == 'mutant' or key =='qh':
                has_mutant = True
                continue
            self[key] = {}
            tmpdict = {}
            for dkey in app.calc_types[key]['complex'].data:
                tmpdict[dkey] = make_array(app.calc_types[key]['complex'].data[dkey])
            self[key]['complex'] = tmpdict
            if not self.stability:
                tmpdict = {}
                for dkey in app.calc_types[key]['receptor'].data:
                    tmpdict[dkey] = make_array(app.calc_types[key]['receptor'].data[dkey])
                self[key]['receptor'] = tmpdict
                tmpdict = {}
                for dkey in app.calc_types[key]['ligand'].data:
                    tmpdict[dkey] = make_array(app.calc_types[key]['ligand'].data[dkey])
                self[key]['ligand'] = tmpdict
                tmpdict = {}
                for dkey in app.calc_types[key]['delta'].data:
                    tmpdict[dkey] = make_array(app.calc_types[key]['delta'].data[dkey])
                self[key]['delta'] = tmpdict

        # Are we doing a mutant?
        if has_mutant:
            self.mutant = {}
            for key in app.calc_types['mutant']:
                if key == 'qh': continue
                self.mutant[key] = {}
                tmpdict = {}
                for dkey in app.calc_types['mutant'][key]['complex'].data:
                    tmpdict[dkey] = make_array(
                        app.calc_types['mutant'][key]['complex'].data[dkey])
                self.mutant[key]['complex'] = tmpdict
                if not self.stability:
                    self.mutant[key]['receptor'] = {}
                    tmpdict = {}
                    for dkey in app.calc_types['mutant'][key]['receptor'].data:
                        tmpdict[dkey] = make_array(
                            app.calc_types['mutant'][key]['receptor'].data[dkey])
                    self.mutant[key]['receptor'] = tmpdict
                    tmpdict = {}
                    for dkey in app.calc_types['mutant'][key]['ligand'].data:
                        tmpdict[dkey] = make_array(
                            app.calc_types['mutant'][key]['ligand'].data[dkey])
                    self.mutant[key]['ligand'] = tmpdict
                    tmpdict = {}
                    for dkey in app.calc_types['mutant'][key]['delta'].data:
                        tmpdict[dkey] = make_array(
                            app.calc_types['mutant'][key]['delta'].data[dkey])
                    self.mutant[key]['delta'] = tmpdict


    def __iadd__(self, other):
        """
        Adding one to another extends every array. The way we do this depends on
        whether we're using numpy arrays or
        """
        global HAS_NUMPY
        if HAS_NUMPY:
            return self._add_numpy(other)
        else:
            return self._add_nonumpy(other)

    def _add_numpy(self, other):
        """
        If we have numpy available, we need to extend every array in a numpy-valid
        way
        """
        used_keys = []
        for key in self:
            used_keys.append(key)
            try:
                for dkey in self[key]['complex']:
                    _combine_np_arrays(self[key]['complex'][dkey],
                                       other[key]['complex'][dkey])
                for dkey in self[key]['receptor']:
                    _combine_np_arrays(self[key]['receptor'][dkey],
                                       other[key]['receptor'][dkey])
                for dkey in self[key]['ligand']:
                    _combine_np_arrays(self[key]['ligand'][dkey],
                                       other[key]['ligand'][dkey])
            except KeyError:
                pass
        for key in other:
            if key in used_keys:
                continue
            # If we didn't have a particular calc type, copy that array in here
            self[key] = deepcopy(other[key])
        # Check mutant statuses. If the other has mutant and I don't, copy other
        # If we both have mutant, combine.  If only I do, already done
        if self.mutant and not other.mutant:
            self.mutant = deepcopy(other.mutant)
        elif self.mutant and other.mutant:
            used_keys_mutant = []
            for key in self.mutant:
                used_keys_mutant.append(key)
                try:
                    for dkey in self[key]['complex']:
                        _combine_np_arrays(self[key]['complex'][dkey],
                                           other[key]['complex'][dkey])
                    for dkey in self[key]['receptor']:
                        _combine_np_arrays(self[key]['receptor'][dkey],
                                           other[key]['receptor'][dkey])
                    for dkey in self[key]['ligand']:
                        _combine_np_arrays(self[key]['ligand'][dkey],
                                           other[key]['ligand'][dkey])
                except KeyError:
                    pass
            for key in other.mutant:
                if key in used_keys_mutant:
                    continue
                self.mutant[key] = deepcopy(other.mutant[key])

    def _add_nonumpy(self, other):
        """ Adds up 2 array objects (just use 'extend' method) """
        used_keys = []
        for key in self:
            used_keys.append(key)
            for dkey in self[key]:
                try:
                    self[key]['complex'][dkey].extend(other[key]['complex'][dkey])
                    self[key]['receptor'][dkey].extend(other[key]['receptor'][dkey])
                    self[key]['ligand'][dkey].extend(other[key]['ligand'][dkey])
                except KeyError:
                    pass
        for key in other:
            if key in used_keys:
                continue
            # If we didn't have a particular calc type, copy that array in here
            self[key] = deepcopy(other[key])
        # Check mutant statuses. If the other has mutant and I don't, copy other
        # If we both have mutant, combine.  If only I do, already done
        if self.mutant and not other.mutant:
            self.mutant = deepcopy(other.mutant)
        elif self.mutant and other.mutant:
            used_keys_mutant = []
            for key in self.mutant:
                used_keys_mutant.append(key)
                for dkey in self.mutant[key]:
                    try:
                        self.mutant[key]['complex'].extend(
                            other.mutant[key]['complex'])
                        self.mutant[key]['receptor'].extend(
                            other.mutant[key]['receptor'])
                        self.mutant[key]['ligand'].extend(
                            other.mutant[key]['ligand'])
                    except KeyError:
                        pass
            for key in other.mutant:
                if key in used_keys_mutant:
                    continue
                self.mutant[key] = deepcopy(other.mutant[key])


def _combine_np_arrays(nparray1, nparray2):
    origsize = nparray1.shape[0]
    nparray1.resize(origsize + nparray2.shape[0])
    for i in range(nparray2.shape[0]):
        nparray1[origsize + i] = nparray2[i]

    def load_topologies(self):
        """
        Loads the topology files so we have residue information for decomp
        analyses
        """
        self.app.loadcheck_prmtops()

class APIDecompOut(amber_outputs.DecompOut):

    def __init__(self, basename, surften, num_files, verbose, nframes):
        amber_outputs.DecompOut.__init__(self, basename, None, surften, False,
                                         num_files, verbose)
        self.array_data = {}
        # Make a new dict for all printed tokens (TDC,SDC,BDC)
        for key in self.allowed_tokens:
            self.array_data[key] = {}
            for i in range(nframes):
                for j in range(self.num_terms):
                    rnum, internal, vdw, eel, pol, sas, tot = self.get_next_term(key)
                    try:
                        self.array_data[key][rnum]['int'][i] = internal
                        self.array_data[key][rnum]['vdw'][i] = vdw
                        self.array_data[key][rnum]['eel'][i] = eel
                        self.array_data[key][rnum]['pol'][i] = pol
                        self.array_data[key][rnum]['sas'][i] = sas
                        self.array_data[key][rnum]['tot'][i] = tot
                    except KeyError:
                        # This is the first frame, we don't have the rnum dict yet, so
                        # make that dict here and create the arrays, then fill the
                        # first term
                        self.array_data[key][rnum] = {}
                        for k in ('int', 'vdw', 'eel', 'pol', 'sas', 'tot'):
                            self.array_data[key][rnum][k] = make_array_len(nframes)
                        self.array_data[key][rnum]['int'][i] = internal
                        self.array_data[key][rnum]['vdw'][i] = vdw
                        self.array_data[key][rnum]['eel'][i] = eel
                        self.array_data[key][rnum]['pol'][i] = pol
                        self.array_data[key][rnum]['sas'][i] = sas
                        self.array_data[key][rnum]['tot'][i] = tot

class APIPairDecompOut(amber_outputs.PairDecompOut):

    def __init__(self, basename, surften, num_files, verbose, nframes):
        amber_outputs.DecompOut.__init__(self, basename, None, surften, False,
                                         num_files, verbose)
        self.array_data = {}
        # Make a new dict for all printed tokens (TDC,SDC,BDC)
        for key in self.allowed_tokens:
            self.array_data[key] = {}

        for i in range(nframes):
            for key in self.allowed_tokens:
                for j in range(self.num_terms):
                    rnum, rnum2, internal, vdw, eel, pol, sas, tot = \
                        self.get_next_term(key)
                    dkey = '%d-%d' % (rnum, rnum2)
                    try:
                        self.array_data[key][dkey]['int'][i] = internal
                        self.array_data[key][dkey]['vdw'][i] = vdw
                        self.array_data[key][dkey]['eel'][i] = eel
                        self.array_data[key][dkey]['pol'][i] = pol
                        self.array_data[key][dkey]['sas'][i] = sas
                        self.array_data[key][dkey]['tot'][i] = tot
                    except KeyError:
                        # This is the first frame, we don't have the rnum dict yet, so
                        # make that dict here and create the arrays, then fill the
                        # first term
                        self.array_data[key][dkey] = {}
                        for k in ('int', 'vdw', 'eel', 'pol', 'sas', 'tot'):
                            self.array_data[key][dkey][k] = make_array_len(nframes)
                        self.array_data[key][dkey]['int'][i] = internal
                        self.array_data[key][dkey]['vdw'][i] = vdw
                        self.array_data[key][dkey]['eel'][i] = eel
                        self.array_data[key][dkey]['pol'][i] = pol
                        self.array_data[key][dkey]['sas'][i] = sas
                        self.array_data[key][dkey]['tot'][i] = tot

class APIDecompOutGMX(amber_outputs.DecompOut):

    def __init__(self, basename, surften, num_files, verbose, nframes, lig_num=None):
        amber_outputs.DecompOut.__init__(self, basename, None, surften, False,
                                         num_files, verbose)
        self.array_data = {}
        # Make a new dict for all printed tokens (TDC,SDC,BDC)
        for key in self.allowed_tokens:
            self.array_data[key] = {}
        for i in range(nframes):
            for key in self.allowed_tokens:
                n = 0
                for j in range(self.num_terms):
                    rnum, internal, vdw, eel, pol, sas, tot = self.get_next_term(key)
                    if lig_num:
                        if n == len(lig_num):
                            n = 0
                        rnum = lig_num[n]
                        n += 1
                    if rnum not in self.array_data[key]:
                        self.array_data[key][rnum] = {}
                        for k in ('int', 'vdw', 'eel', 'pol', 'sas', 'tot'):
                            self.array_data[key][rnum][k] = make_array_len(nframes)
                    self.array_data[key][rnum]['int'][i] = internal
                    self.array_data[key][rnum]['vdw'][i] = vdw
                    self.array_data[key][rnum]['eel'][i] = eel
                    self.array_data[key][rnum]['pol'][i] = pol
                    self.array_data[key][rnum]['sas'][i] = sas
                    self.array_data[key][rnum]['tot'][i] = tot

class APIPairDecompOutGMX(amber_outputs.PairDecompOut):

    def __init__(self, basename, surften, num_files, verbose, nframes, lig_num=None):
        amber_outputs.DecompOut.__init__(self, basename, None, surften, False,
                                         num_files, verbose)
        self.array_data = {}
        # Make a new dict for all printed tokens (TDC,SDC,BDC)
        for key in self.allowed_tokens:
            self.array_data[key] = {}

        for i in range(nframes):
            for key in self.allowed_tokens:
                m = 0
                n = 0
                for j in range(self.num_terms):
                    rnum, rnum2, internal, vdw, eel, pol, sas, tot = self.get_next_term(key)
                    if lig_num: # rename lig residues like complex lig names
                        if n == len(lig_num):
                            m += 1
                            n = 0
                        rnum = lig_num[m]
                        rnum2 = lig_num[n]
                        n += 1
                    if rnum not in self.array_data[key]:
                        self.array_data[key][rnum] = {}
                    if rnum2 not in self.array_data[key][rnum]:
                        self.array_data[key][rnum][rnum2] =  {}
                        for k in ('int', 'vdw', 'eel', 'pol', 'sas', 'tot'):
                            self.array_data[key][rnum][rnum2][k] = make_array_len(nframes)
                    self.array_data[key][rnum][rnum2]['int'][i] = internal
                    self.array_data[key][rnum][rnum2]['vdw'][i] = vdw
                    self.array_data[key][rnum][rnum2]['eel'][i] = eel
                    self.array_data[key][rnum][rnum2]['pol'][i] = pol
                    self.array_data[key][rnum][rnum2]['sas'][i] = sas
                    self.array_data[key][rnum][rnum2]['tot'][i] = tot

def get_delta_decomp(app, decomp_calc_type, data):
    """

    :param calc_type: gb or pb
    :param data:
    :return:
    """
    data_out = {decomp_calc_type: {'delta': {}}}
    # data_out[decomp_calc_type] = {'delta': {}}
    tempdict = {}
    # complex, if stability: receptor, ligand
    com = data[decomp_calc_type]['complex']
    rec = data[decomp_calc_type]['receptor']
    lig = data[decomp_calc_type]['ligand']
    for p in com:
        tempdict[p] = {}
        for res in com[p]:
            if res not in tempdict[p]:
                tempdict[p][res] = {}
            if app.INPUT['idecomp'] in [1, 2]:
                if rec and res in rec[p]:
                    for para in com[p][res]:
                        d = com[p][res][para] - rec[p][res][para]
                        tempdict[p][res][para] = d
                elif lig and res in lig[p]:
                    for para in com[p][res]:
                        d = com[p][res][para] - lig[p][res][para]
                        tempdict[p][res][para] = d
            else:
                if rec and res in rec[p]:
                    for resp, resp1 in zip(com[p][res], rec[p][res]):
                        tempdict[p][res][resp] = {}
                        for para in com[p][res][resp]:
                            d = com[p][res][resp][para] - rec[p][res][resp1][para]
                            tempdict[p][res][resp][para] = d
                    for resp in com[p][res]:
                        tempdict[p][res][resp] = com[p][res][resp]
                elif lig and res in lig[p]:
                    for resp, resp1 in zip(list(com[p][res])[(len(com[p][res]) - len(lig[p][res])):],
                                           lig[p][res]):
                        tempdict[p][res][resp] = {}
                        for para in com[p][res][resp]:
                            d = com[p][res][resp][para] - lig[p][res][resp1][para]
                            tempdict[p][res][resp][para] = d
                    for resp in com[p][res]:
                        tempdict[p][res][resp] = com[p][res][resp]
    return tempdict

def load_mmpbsa_info(fname):
    """
    Loads up an MMPBSA.py info file and returns a mmpbsa_data instance with all
    of the data available in numpy arrays if numpy is available. The returned
    object is a mmpbsa_data instance.

    mmpbsa_data attributes:
    -----------------------
       o  Derived from "dict"
       o  Each solvent model is a dictionary key for a numpy array (if numpy is
          available) or array.array (if numpy is unavailable) for each of the
          species (complex, receptor, ligand) present in the calculation.
       o  The alanine scanning mutant data is under another dict denoted by the
          'mutant' key.

    Data Layout:
    ------------
       Solvent Model     |  Dictionary Key    |  Data Keys Available
       -------------------------------------------------------------------
       Generalized Born  |  'gb'              |  EGB, ESURF, *
       Poisson-Boltzmann |  'pb'              |  EPB, EDISPER, ECAVITY, *
       3D-RISM (GF)      |  'rism gf'         |
       3D-RISM (Standard)|  'rism std'        |
       Normal Mode       |  'nmode'           |
       Quasi-harmonic    |  'qh'              |

    * == TOTAL, VDW, EEL, 1-4 EEL, 1-4 VDW, BOND, ANGLE, DIHED

    The keys above are entries for the main dict as well as the sub-dict whose
    key is 'mutant' in the main dict.  Each entry in the main (and mutant sub-)
    dict is, itself, a dict with 1 or 3 keys; 'complex', 'receptor', 'ligand';
    where 'receptor' and 'ligand' are missing for stability calculations.
    If numpy is available, all data will be numpy.ndarray instances.  Otherwise,
    all data will be array.array instances.

    All of the objects referenced by the listed 'Dictionary Key's are dicts in
    which the listed 'Data Keys Available' are keys to the data arrays themselves

    Examples:
    ---------
       # Load numpy for our analyses (optional)
       import numpy as np

       # Load the _MMPBSA_info file:
       mydata = load_mmpbsa_info('_MMPBSA_info')

       # Access the complex GB data structure and calculate the autocorr. fcn.
       autocorr = np.correlate(mydata['gb']['complex']['TOTAL'],
                               mydata['gb']['complex']['TOTAL'])

       # Calculate the standard deviation of the alanine mutant receptor in PB
       print mydata.mutant['pb']['receptor']['TOTAL'].std()
    """
    if not HAS_NUMPY:
        warnings.warn('numpy was not found. Data will be packed in normal Python '
                      'arrays. Install numpy for more efficient array handling.',
                      NotImplemented)

    if not os.path.exists(fname):
        raise NoFileExists("cannot find %s!" % fname)
    app = main.MMPBSA_App(MPI)
    info = infofile.InfoFile(app)
    info.read_info(fname)
    app.normal_system = app.mutant_system = None
    app.parse_output_files()
    return_data = mmpbsa_data(app)
    # Since Decomp data is parsed in a memory-efficient manner (by not storing
    # all of the data in arrays, but rather by printing each data point as it's
    # parsed), we need to handle the decomp data separately here
    if app.INPUT['decomprun']:
        # Simplify the decomp class instance creation
        if app.INPUT['idecomp'] in (1, 2):
            DecompClass = lambda x, y: APIDecompOut(x, y, app.mpi_size,
                                                    app.INPUT['dec_verbose'], app.numframes)
        else:
            DecompClass = lambda x, y: APIPairDecompOut(x, y, app.mpi_size,
                                                        app.INPUT['dec_verbose'], app.numframes)

        if not app.INPUT['mutant_only']:
            # Do normal GB
            if app.INPUT['gbrun']:
                return_data['decomp'] = {'gb' : {}}
                return_data['decomp']['gb']['complex'] = DecompClass(
                    app.FILES.prefix + 'complex_gb.mdout',
                    app.INPUT['surften']).array_data
                if not app.stability:
                    return_data['decomp']['gb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'receptor_gb.mdout',
                        app.INPUT['surften']).array_data
                    return_data['decomp']['gb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'ligand_gb.mdout',
                        app.INPUT['surften']).array_data
            # Do normal PB
            if app.INPUT['pbrun']:
                return_data['decomp'] = {'pb' : {}}
                return_data['decomp']['pb']['complex'] = DecompClass(
                    app.FILES.prefix + 'complex_pb.mdout',
                    app.INPUT['surften']).array_data
                if not app.stability:
                    return_data['decomp']['pb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'receptor_pb.mdout',
                        app.INPUT['surften']).array_data
                    return_data['decomp']['pb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'ligand_pb.mdout',
                        app.INPUT['surften']).array_data
        if app.INPUT['alarun']:
            # Do mutant GB
            if app.INPUT['gbrun']:
                return_data.mutant['decomp'] = {'gb' : {}}
                return_data.mutant['decomp']['gb']['complex'] = DecompClass(
                    app.FILES.prefix + 'mutant_complex_gb.mdout',
                    app.INPUT['surften']).array_data
                if not app.stability:
                    return_data.mutant['decomp']['gb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'mutant_receptor_gb.mdout',
                        app.INPUT['surften']).array_data
                    return_data.mutant['decomp']['gb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'mutant_ligand_gb.mdout',
                        app.INPUT['surften']).array_data
            # Do mutant PB
            if app.INPUT['pbrun']:
                return_data.mutant['decomp'] = {'pb' : {}}
                return_data.mutant['decomp']['pb']['complex'] = DecompClass(
                    app.FILES.prefix + 'mutant_complex_pb.mdout',
                    app.INPUT['surften']).array_data
                if not app.stability:
                    return_data.mutant['decomp']['pb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'mutant_receptor_pb.mdout',
                        app.INPUT['surften']).array_data
                    return_data.mutant['decomp']['pb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'mutant_ligand_pb.mdout',
                        app.INPUT['surften']).array_data
        else:
            return_data.mutant = None

    return return_data

def load_gmxmmpbsa_info(fname):
    """
    Loads up an gmx_MMPBSA info file and returns a mmpbsa_data instance with all
    of the data available in numpy arrays if numpy is available. The returned
    object is a mmpbsa_data instance.

    change the structure to get more easy way to graph per residue

    mmpbsa_data attributes:
    -----------------------
       o  Derived from "dict"
       o  Each solvent model is a dictionary key for a numpy array (if numpy is
          available) or array.array (if numpy is unavailable) for each of the
          species (complex, receptor, ligand) present in the calculation.
       o  The alanine scanning mutant data is under another dict denoted by the
          'mutant' key.

    Data Layout:
    ------------
       Solvent Model     |  Dictionary Key    |  Data Keys Available
       -------------------------------------------------------------------
       Generalized Born  |  'gb'              |  EGB, ESURF, *
       Poisson-Boltzmann |  'pb'              |  EPB, EDISPER, ECAVITY, *
       3D-RISM (GF)      |  'rism gf'         |
       3D-RISM (Standard)|  'rism std'        |
       Normal Mode       |  'nmode'           |
       Quasi-harmonic    |  'qh'              |

    * == TOTAL, VDW, EEL, 1-4 EEL, 1-4 VDW, BOND, ANGLE, DIHED

    The keys above are entries for the main dict as well as the sub-dict whose
    key is 'mutant' in the main dict.  Each entry in the main (and mutant sub-)
    dict is, itself, a dict with 1 or 3 keys; 'complex', 'receptor', 'ligand';
    where 'receptor' and 'ligand' are missing for stability calculations.
    If numpy is available, all data will be numpy.ndarray instances.  Otherwise,
    all data will be array.array instances.

    All of the objects referenced by the listed 'Dictionary Key's are dicts in
    which the listed 'Data Keys Available' are keys to the data arrays themselves

    Examples:
    ---------
       # Load numpy for our analyses (optional)
       import numpy as np

       # Load the _MMPBSA_info file:
       mydata = load_mmpbsa_info('_MMPBSA_info')

       # Access the complex GB data structure and calculate the autocorr. fcn.
       autocorr = np.correlate(mydata['gb']['complex']['TOTAL'],
                               mydata['gb']['complex']['TOTAL'])

       # Calculate the standard deviation of the alanine mutant receptor in PB
       print mydata.mutant['pb']['receptor']['TOTAL'].std()
    """
    if not HAS_NUMPY:
        warnings.warn('numpy was not found. Data will be packed in normal Python '
                      'arrays. Install numpy for more efficient array handling.',
                      NotImplemented)

    if not os.path.exists(fname):
        raise NoFileExists("cannot find %s!" % fname)
    app = main.MMPBSA_App(MPI)
    info = infofile.InfoFile(app)
    info.read_info(fname)
    app.normal_system = app.mutant_system = None
    app.parse_output_files()
    return_data = mmpbsa_data(app)
    # Since Decomp data is parsed in a memory-efficient manner (by not storing
    # all of the data in arrays, but rather by printing each data point as it's
    # parsed), we need to handle the decomp data separately here
    if not app.INPUT['alarun']:
        return_data.mutant = None
    lig_res = None
    if app.INPUT['decomprun']:
        # Simplify the decomp class instance creation
        if app.INPUT['idecomp'] in (1, 2):
            DecompClass = lambda x, y, z=None: APIDecompOutGMX(x, y, app.mpi_size,
                                                    app.INPUT['dec_verbose'], app.numframes, lig_num=z)
        else:
            DecompClass = lambda x, y, z=None: APIPairDecompOutGMX(x, y, app.mpi_size,
                                                        app.INPUT['dec_verbose'], app.numframes, lig_num=z)

        if not app.INPUT['mutant_only']:
            # Do normal GB
            if app.INPUT['gbrun']:
                return_data['decomp'] = {'gb' : {}}
                return_data['decomp']['gb']['complex'] = DecompClass(app.FILES.prefix + 'complex_gb.mdout',
                                                                     app.INPUT['surften']).array_data
                com_res = []
                for k in return_data['decomp']['gb']['complex']['TDC']:
                    if k not in com_res:
                        com_res.append(k)

                if not app.stability:
                    return_data['decomp']['gb']['receptor'] = DecompClass(app.FILES.prefix + 'receptor_gb.mdout',
                                                                          app.INPUT['surften']).array_data
                    rec_res = []
                    for k in return_data['decomp']['gb']['receptor']['TDC']:
                        if k not in rec_res:
                            rec_res.append(k)
                    lig_res = com_res[len(rec_res):]
                    return_data['decomp']['gb']['ligand'] = DecompClass(app.FILES.prefix + 'ligand_gb.mdout',
                                                                        app.INPUT['surften'], lig_res).array_data
                    return_data['decomp']['gb']['delta'] = get_delta_decomp(app, 'gb', return_data['decomp'])
            # Do normal PB
            if app.INPUT['pbrun']:
                return_data['decomp'] = {'pb' : {}}
                return_data['decomp']['pb']['complex'] = DecompClass(app.FILES.prefix + 'complex_pb.mdout',
                                                                     app.INPUT['surften']).array_data
                if not app.stability:
                    return_data['decomp']['pb']['receptor'] = DecompClass(app.FILES.prefix + 'receptor_pb.mdout',
                                                                          app.INPUT['surften']).array_data
                    return_data['decomp']['pb']['ligand'] = DecompClass(app.FILES.prefix + 'ligand_pb.mdout',
                                                                        app.INPUT['surften'], lig_res).array_data
                    return_data['decomp']['pb']['delta'] = get_delta_decomp(app, 'pb', return_data['decomp'])
        if app.INPUT['alarun']:
            # Do mutant GB
            if app.INPUT['gbrun']:
                return_data.mutant['decomp'] = {'gb' : {}}
                return_data.mutant['decomp']['gb']['complex'] = DecompClass(app.FILES.prefix +
                                                                             'mutant_complex_gb.mdout',
                                                                            app.INPUT['surften']).array_data
                if not app.stability:
                    return_data.mutant['decomp']['gb']['receptor'] = DecompClass(app.FILES.prefix +
                                                                                  'mutant_receptor_gb.mdout',
                                                                                 app.INPUT['surften']).array_data
                    return_data.mutant['decomp']['gb']['ligand'] = DecompClass(app.FILES.prefix +
                                                                                'mutant_ligand_gb.mdout',
                                                                               app.INPUT['surften'], lig_res).array_data
                    return_data.mutant['decomp']['gb']['delta'] = get_delta_decomp(app, 'gb',
                                                                                   return_data.mutant['decomp'])
            # Do mutant PB
            if app.INPUT['pbrun']:
                return_data.mutant['decomp'] = {'pb' : {}}
                return_data.mutant['decomp']['pb']['complex'] = DecompClass(app.FILES.prefix +
                                                                             'mutant_complex_pb.mdout',
                                                                            app.INPUT['surften']).array_data
                if not app.stability:
                    return_data.mutant['decomp']['pb']['receptor'] = DecompClass(app.FILES.prefix +
                                                                                  'mutant_receptor_pb.mdout',
                                                                                 app.INPUT['surften']).array_data
                    return_data.mutant['decomp']['pb']['ligand'] = DecompClass(app.FILES.prefix +
                                                                                'mutant_ligand_pb.mdout',
                                                                               app.INPUT['surften'], lig_res).array_data
                    return_data.mutant['decomp']['pb']['delta'] = get_delta_decomp(app, 'pb',
                                                                                   return_data.mutant['decomp'])
        else:
            return_data.mutant = None

    return return_data, app

