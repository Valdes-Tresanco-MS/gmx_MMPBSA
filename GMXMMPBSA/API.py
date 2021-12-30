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

from typing import Union
from GMXMMPBSA import infofile, main
from GMXMMPBSA.exceptions import NoFileExists
from GMXMMPBSA.fake_mpi import MPI
import pandas as pd
from pathlib import Path
import os
import numpy as np
import h5py
from types import SimpleNamespace

__all__ = ['load_gmxmmpbsa_info']


class DataStore(dict):
    def __init__(self, *args):
        super(DataStore, self).__init__(*args)


class H52Data:
    def __init__(self, fname):
        self.h5f = h5py.File(fname, 'r')
        self.app_namespace = SimpleNamespace(INPUT={}, FILES=SimpleNamespace(), INFO={})
        self.calc_types = DataStore()
        self.calc_types.mutant = DataStore()
        self.calc_types.decomp = DataStore()
        self.calc_types.decomp.mutant = DataStore()


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
        if key in ['gb', 'pb', 'rism std', 'rism gf', 'nmode', 'qh', 'ie', 'c2']:
            calc_types[key] = {}
            # key2 is complex, receptor, ligand, delta
            for key2 in d[key]:
                calc_types[key][key2] = {}
                # complex, receptor, etc., is a class and the data is contained in the attribute data
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
                        for key5 in d[key][key2][key3][key4]:
                            if isinstance(d[key][key2][key3][key4], h5py.Group):
                                # residue sec level
                                for key6 in d[key][key2][key3][key4][key5]:
                                    calc_types[key][key2][(key3, key4, key5, key6)] = d[key][key2][key3][key4][
                                        key5][key6][()]
                            else:
                                # energy terms
                                for key5 in d[key][key2][key3][key4]:
                                    calc_types[key][key2][(key3, key4, key5)] = d[key][key2][key3][key4][key5][()]


class DataMMPBSA:
    """ Main class that holds all of the Free Energy data """

    def __init__(self):
        self.frames = []
        self.app_namespace = SimpleNamespace()
        self.data = DataStore()
        self.data.mutant = DataStore()

    def get_fromH5(self, h5file):

        h5file = H52Data(h5file)
        self.app_namespace = h5file.app_namespace
        self._get_data(h5file)

    def get_fromApp(self, ifile):

        app = main.MMPBSA_App(MPI)
        info = infofile.InfoFile(app)
        info.read_info(ifile)
        app.normal_system = app.mutant_system = None
        app.parse_output_files()
        self.app_namespace = self._get_namespace(app)
        self._get_data(app)

    @staticmethod
    def _get_namespace(app):

        input_file_text = ('|Input file:\n|--------------------------------------------------------------\n|'
                           + ''.join(open(app.FILES.input_file).readlines()).replace('\n', '\n|') +
                           '--------------------------------------------------------------\n')
        INFO = {'COM_PDB': ''.join(open(app.FILES.complex_fixed).readlines()),
                'input_file': input_file_text,
                'mut_str': app.mut_str,
                'numframes': app.numframes,
                'numframes_nmode': app.numframes_nmode,
                # FIXME: arreglar compatibilidad con la versiones anteriores, algunas variables no existen
                'output_file': ''.join(open(app.FILES.output_file).readlines()),
                'size': app.mpi_size,
                'using_chamber': app.using_chamber}
        if app.INPUT['decomprun']:
            INFO['decomp_output_file'] = ''.join(open(app.FILES.decompout).readlines())

        return SimpleNamespace(FILES=app.FILES, INPUT=app.INPUT, INFO=INFO)

    def _get_data(self, data_object: Union[H52Data, main.MMPBSA_App]):
        # check the data_object type
        h5 = isinstance(data_object, H52Data)
        INPUT = self.app_namespace.INPUT
        numframes = self.app_namespace.INFO['numframes']
        numframes_nmode = self.app_namespace.INFO['numframes_nmode']
        # See if we are doing stability
        self.stability = self.app_namespace.FILES.stability

        # Now load the data into the dict
        self.frames = list(
            range(
                INPUT['startframe'],
                INPUT['startframe'] + numframes * INPUT['interval'],
                INPUT['interval'],
            )
        )

        INPUT['endframe'] = self.frames[-1]

        self.nmode_frames = list(
            range(
                INPUT['nmstartframe'],
                INPUT['nmstartframe'] + numframes_nmode * INPUT['interval'],
                INPUT['interval'],
            )
        )

        if numframes_nmode:
            INPUT['nmendframe'] = self.nmode_frames[-1]
        # Now load the data
        if not INPUT['mutant_only']:
            self._get_edata(data_object.calc_types, h5)
        # Are we doing a mutant?
        if data_object.calc_types.mutant:
            self._get_edata(data_object.calc_types.mutant, h5, True)

        # Now we load the decomp data. Avoid to get the decomp data in first place in gmx_MMPBSA_ana
        if INPUT['decomprun']:
            if not INPUT['mutant_only']:
                self.data['decomp'] = self._get_ddata(data_object.calc_types.decomp, h5)
            if data_object.calc_types.mutant:
                self.data.mutant['decomp'] = self._get_ddata(data_object.calc_types.decomp.mutant, h5)

    def _get_edata(self, calc_types, h5=False, mut=False):

        data = self.data.mutant if mut else self.data
        for key in calc_types:
            if key == 'decomp':
                continue
            if key in ['ie', 'c2']:
                # since the model data object in MMPBSA_App contain the data in the attribute data and H5 not,
                # we need to define a conditional object
                calc_type_data = calc_types[key] if h5 else calc_types[key].data
                data[key] = {}
                if key == 'ie':
                    for iekey in calc_type_data:
                        data[key][iekey] = {
                            'data': pd.DataFrame(calc_type_data[iekey]['data'], columns=['data'],
                                                 index=self.frames),
                            'iedata': pd.DataFrame(calc_type_data[iekey]['iedata'], columns=['iedata'],
                                                   index=self.frames[-calc_type_data[iekey]['ieframes']:]),
                            'ieframes': calc_type_data[iekey]['ieframes'],
                            'sigma': calc_type_data[iekey]['sigma']}
                else:
                    for c2key in calc_type_data:
                        data[key][c2key] = {'c2data': calc_type_data[c2key]['c2data'],
                                            'sigma': calc_type_data[c2key]['sigma'],
                                            'c2_std': calc_type_data[c2key]['c2_std'],
                                            'c2_ci': calc_type_data[c2key]['c2_ci']}
            else:
                # Since the energy models have the same structure as nmode and qh, we only need to worry about
                # correctly defining the Dataframe index, that is, the frames
                if key == 'nmode':
                    cframes = self.nmode_frames
                elif key == 'qh':
                    cframes = [0]
                else:
                    cframes = self.frames
                # since the model data object in MMPBSA_App contain the data in the attribute data and H5 not,
                # we need to define a conditional object
                com_calc_type_data = (calc_types[key]['complex'] if h5
                                      else calc_types[key]['complex'].data)
                df = complex = pd.DataFrame({dkey: com_calc_type_data[dkey] for dkey in com_calc_type_data},
                                            index=cframes)
                if not self.stability:
                    rec_calc_type_data = (calc_types[key]['receptor'] if h5
                                          else calc_types[key]['receptor'].data)
                    receptor = pd.DataFrame({dkey: rec_calc_type_data[dkey] for dkey in rec_calc_type_data},
                                            index=cframes)
                    lig_calc_type_data = (calc_types[key]['ligand'] if h5
                                          else calc_types[key]['ligand'].data)
                    ligand = pd.DataFrame({dkey: lig_calc_type_data[dkey] for dkey in lig_calc_type_data},
                                          index=cframes)
                    delta = complex - receptor - ligand
                    df = pd.concat([complex, receptor, ligand, delta], axis=1,
                                   keys=['complex', 'receptor', 'ligand', 'delta'])
                data[key] = df

    def _get_ddata(self, calc_types, h5=False):
        data = {}
        # Take the decomp data
        for key in calc_types:
            # since the model data object in MMPBSA_App contain the data in the attribute data and H5 not,
            # we need to define a conditional object. Also, the decomp data must be re-structured for multiindex
            # Dataframe
            com_calc_type_data = (calc_types[key]['complex'] if h5
                                  else self._transform_from_lvl_decomp(calc_types[key]['complex']))
            df = complex = pd.DataFrame(com_calc_type_data, index=self.frames)
            if not self.stability:
                rec_calc_type_data = (calc_types[key]['receptor'] if h5
                                      else self._transform_from_lvl_decomp(calc_types[key]['receptor']))
                receptor = pd.DataFrame(rec_calc_type_data, index=self.frames)

                lig_calc_type_data = (calc_types[key]['ligand'] if h5
                                      else self._transform_from_lvl_decomp(calc_types[key]['ligand']))
                ligand = pd.DataFrame(lig_calc_type_data, index=self.frames)

                delta = complex.subtract(pd.concat([receptor, ligand], axis=1)).combine_first(
                    complex).reindex_like(df)
                df = pd.concat([complex, receptor, ligand, delta], axis=1,
                               keys=['complex', 'receptor', 'ligand', 'delta'])
            data[key] = df
        return data

    @staticmethod
    def _transform_from_lvl_decomp(nd):
        data = {}
        for k2, v2 in nd.items():  # TDC, SDC, BDC
            for k3, v3 in v2.items():  # residue
                for k4, v4 in v3.items():  # residue in per-wise or terms in per-res
                    if isinstance(v4, dict):  # per-wise
                        for k5, v5 in v4.items():
                            data[(k2, k3, k4, k5)] = v5
                    else:
                        data[(k2, k3, k4)] = v4
        return data


def load_gmxmmpbsa_info(fname: Union[Path, str]):
    """
    Loads up an gmx_MMPBSA info or h5 file and returns a dictionary with all
    energy terms and  a namespace with some variables needed in gmx_MMPBSA_ana

    Depending on the data type store the variable can be a Dataframe, string or
    scalar

    Data attributes:
    -----------------------
       o  All attributes from dict
       o  Each solvent model is a dictionary key for a numpy array (if numpy is
          available) or array.array (if numpy is unavailable) for each of the
          species (complex, receptor, ligand) present in the calculation.
       o  The alanine scanning mutant data is under another dict denoted by the
          'mutant' key.

    Data Layout:
    ------------
               Model       | Dict Key    |  Data Keys Available
       --------------------------------------------------------------
       Generalized Born    | 'gb'        |  EGB, ESURF, *
       Poisson-Boltzmann   | 'pb'        |  EPB, EDISPER, ECAVITY, *
       3D-RISM (GF)        | 'rism gf'   |  POLAR SOLV, APOLAR SOLV, *
       3D-RISM (Standard)  | 'rism std'  |
       Normal Mode         | 'nmode'     |
       Quasi-harmonic      | 'qh'        |
       Interaction entropy | 'ie'        |
       C2 Entropy          | 'c2'        |


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
    if not isinstance(fname, Path):
        fname = Path(fname)

    if not fname.exists():
        raise NoFileExists("cannot find %s!" % fname)
    os.chdir(fname.parent)
    d_mmpbsa = DataMMPBSA()

    if fname.suffix == '.h5':
        d_mmpbsa.get_fromH5(fname)
    else:
        d_mmpbsa.get_fromApp(fname)

    return d_mmpbsa.data, d_mmpbsa.app_namespace
