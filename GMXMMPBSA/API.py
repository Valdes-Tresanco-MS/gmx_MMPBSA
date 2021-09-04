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

from pathlib import Path
import parmed
import os
import numpy as np
from types import SimpleNamespace

__all__ = ['load_gmxmmpbsa_info']


class DataStore(dict):
    def __init__(self):
        super(DataStore, self).__init__()
        self.mutant = {}
        self.decomp = type('calc_types', (dict,), {'mutant': {}})()


class H52Data:
    def __init__(self, fname):
        self.h5f = h5py.File(fname, 'r')
        self.app_namespace = SimpleNamespace(INPUT={}, FILES=SimpleNamespace(), INFO={})
        self.calc_types = DataStore()

        for key in self.h5f:
            if key in ['INFO', 'INPUT', 'FILES']:
                self._h52app_namespace(key)
            elif key == 'decomp':
                self._h52decomp(self.h5f[key])
            elif key == 'mutant':
                for mkey in self.h5f[key]:
                    self._h52e(self.h5f[key], mkey, True)
                    if mkey == 'decomp':
                        self._h52decomp(self.h5f[key][mkey], True)
            else:
                self._h52e(self.h5f, key)
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

    def _h52e(self, d, key, mut=False):

        calc_types = self.calc_types.mutant if mut else self.calc_types
        # key  Energy: [gb, pb, rism std, rism gf], Decomp: [gb, pb], Entropy: [nmode, qh, ie, c2]
        if key in ['gb', 'pb', 'rism std', 'rism gf']:
            calc_types[key] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in d[key]:
                calc_types[key][key2] = {}
                # complex, receptor, etc., is a class and the data is contained in the attribute data
                for key3 in d[key][key2]:
                    calc_types[key][key2][key3] = d[key][key2][key3][()]
        elif key in ['nmode', 'qh']:
            calc_types[key] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in d[key]:
                calc_types[key][key2] = {}
                # vibrational, translational, rotational, total
                for key3 in d[key][key2]:
                    calc_types[key][key2][key3] = d[key][key2][key3][()]

        elif key in ['ie', 'c2']:
            calc_types[key] = {}
            # key2 is PB, GB or RISM?
            for key2 in d[key]:
                calc_types[key][key2] = {}
                for key3 in d[key][key2]:
                    calc_types[key][key2][key3] = d[key][key2][key3][()]

    def _h52decomp(self, d, mut=False):
        calc_types = self.calc_types.decomp.mutant if mut else self.calc_types.decomp
        for key in d:
            # model
            calc_types[key] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in d[key]:
                calc_types[key][key2] = {}
                # TDC, SDC, BDC
                for key3 in d[key][key2]:
                    # residue first level
                    for key4 in d[key][key2][key3]:
                        if isinstance(d[key][key2][key3][key4], h5py.Group):
                            # residue sec level
                            for key5 in d[key][key2][key3][key4]:
                                # calc_types terms
                                for key6 in d[key][key2][key3][key4][key5]:
                                    calc_types[key][key2][(key3, key4, key5, key6)] = d[key][key2][key3][key4][
                                        key5][key6][()]
                        else:
                            # energy terms
                            for key5 in d[key][key2][key3][key4]:
                                calc_types[key][key2][(key3, key4, key5)] = d[key][key2][key3][key4][key5][()]

    mut_com_res_info = []
    mut_com = mut_rec.copy()
    mut_com.update(mut_lig)
    for key, value in sorted(mut_com.items()):
        mut_com_res_info.append(value)

    mut_rec_res_info = []
    for key, value in mut_rec.items():
        mut_rec_res_info.append(value)

    mut_lig_res_info = []
    for key, value in mut_lig.items():
        mut_lig_res_info.append(value)

    if not app.INPUT['alarun']:
        return_data.mutant = {}
    if app.INPUT['decomprun']:
        # Simplify the decomp class instance creation
        if app.INPUT['idecomp'] in (1, 2):
            DecompClass = lambda x, part: APIDecompOut(x, part, app)
        else:
            DecompClass = lambda x, part: APIPairDecompOut(x, part, app)

        if not app.INPUT['mutant_only']:
            # Do normal GB
            if app.INPUT['gbrun']:
                return_data['decomp'] = {'gb' : {}}
                return_data['decomp']['gb']['complex'] = DecompClass(app.FILES.prefix + 'complex_gb.mdout',
                                                                     com_res_info).array_data

                if not app.stability:
                    return_data['decomp']['gb']['receptor'] = DecompClass(app.FILES.prefix + 'receptor_gb.mdout',
                                                                          rec_res_info).array_data
                    return_data['decomp']['gb']['ligand'] = DecompClass(app.FILES.prefix + 'ligand_gb.mdout',
                                                                        lig_res_info).array_data
                    return_data['decomp']['gb']['delta'] = get_delta_decomp(app, 'gb', return_data['decomp'])
            # Do normal PB
            if app.INPUT['pbrun']:
                return_data['decomp'] = {'pb' : {}}
                return_data['decomp']['pb']['complex'] = DecompClass(app.FILES.prefix + 'complex_pb.mdout',
                                                                     com_res_info).array_data
                if not app.stability:
                    return_data['decomp']['pb']['receptor'] = DecompClass(app.FILES.prefix + 'receptor_pb.mdout',
                                                                          rec_res_info).array_data
                    return_data['decomp']['pb']['ligand'] = DecompClass(app.FILES.prefix + 'ligand_pb.mdout',
                                                                        lig_res_info).array_data
                    return_data['decomp']['pb']['delta'] = get_delta_decomp(app, 'pb', return_data['decomp'])
        if app.INPUT['alarun']:
            # Do mutant GB
            if app.INPUT['gbrun']:
                return_data.mutant['decomp'] = {'gb' : {}}
                return_data.mutant['decomp']['gb']['complex'] = DecompClass(
                    app.FILES.prefix + 'mutant_complex_gb.mdout', mut_com_res_info).array_data
                if not app.stability:
                    return_data.mutant['decomp']['gb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'mutant_receptor_gb.mdout', mut_rec_res_info).array_data
                    return_data.mutant['decomp']['gb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'mutant_ligand_gb.mdout', mut_lig_res_info).array_data
                    return_data.mutant['decomp']['gb']['delta'] = get_delta_decomp(app, 'gb',
                                                                                   return_data.mutant['decomp'])
            # Do mutant PB
            if app.INPUT['pbrun']:
                return_data.mutant['decomp'] = {'pb' : {}}
                return_data.mutant['decomp']['pb']['complex'] = DecompClass(
                    app.FILES.prefix + 'mutant_complex_pb.mdout', mut_com_res_info).array_data
                if not app.stability:
                    return_data.mutant['decomp']['pb']['receptor'] = DecompClass(
                        app.FILES.prefix + 'mutant_receptor_pb.mdout', mut_rec_res_info).array_data
                    return_data.mutant['decomp']['pb']['ligand'] = DecompClass(
                        app.FILES.prefix + 'mutant_ligand_pb.mdout', mut_lig_res_info).array_data
                    return_data.mutant['decomp']['pb']['delta'] = get_delta_decomp(app, 'pb',
                                                                                   return_data.mutant['decomp'])
        else:
            return_data.mutant = None

    app_namespace = SimpleNamespace(FILES=app.FILES, INPUT=app.INPUT, numframes=app.numframes,
                                    numframes_nmode=app.numframes_nmode)

    return return_data, app_namespace

