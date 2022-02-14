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
from copy import copy, deepcopy
from typing import Union

from GMXMMPBSA.calculation import InteractionEntropyCalc, C2EntropyCalc

from GMXMMPBSA import infofile, main
from GMXMMPBSA.exceptions import NoFileExists
from GMXMMPBSA.fake_mpi import MPI
from GMXMMPBSA.amber_outputs import (H5Output, BindingStatistics, IEout, C2out, DeltaBindingStatistics,
                                     DecompBinding, PairDecompBinding)
import pandas as pd
from pathlib import Path
import os
import numpy as np
import h5py
from types import SimpleNamespace


class MMPBSA_API():
    """ Main class that holds all the Free Energy data """

    def __init__(self):
        self.app_namespace = SimpleNamespace()
        self.raw_energy = None
        self.data = {}

        self.set_config()

    def set_config(self, starttime=0, timestep=0, timeunit='ps'):
        self.starttime = starttime
        self.timestep = timestep
        self.timeunit = timeunit

    def load_file(self, fname: Union[Path, str]):
        """
        Load the info or h5 file and extract the info

        Args:
            fname: String- or Path-like file reference

        Returns:

        """
        if not isinstance(fname, Path):
            fname = Path(fname)
        if not fname.exists():
            raise NoFileExists("cannot find %s!" % fname)
        os.chdir(fname.parent)

        if fname.suffix == '.h5':
            self._get_fromH5(fname)
        else:
            self._get_fromApp(fname)

    def get_info(self):
        """
        Get all variables in the INFO dictionary
        Returns:

        """
        return self.app_namespace.INFO

    def get_input(self):
        """
        Get all variables defined in the INPUT dictionary
        Returns:

        """
        return self.app_namespace.INPUT

    def get_files(self):
        """
        Get all variables in the FILES dictionary
        Returns:

        """
        return self.app_namespace.FILES

    def get_raw_energy(self):
        """
        Get the energy data for normal and mutant in structured dict
        Returns:

        """
        return self.data

    def _get_df(self, raw_energy, key, delta=False):
        df_models = {}
        terms = {'energy': [], 'entropy': []}
        if self.timestep:
            index = pd.Series(list(self.frames.values()), name=f'Time ({self.timeunit})')
        else:
            index = pd.Series(list(self.frames.keys()), name='Frames')
        for model, data in raw_energy[key].items():
            if model in ['gb', 'pb', 'rism std', 'rism gf', 'nmode', 'qh']:
                if model in ['nmode', 'qh']:
                    terms['entropy'].append(model)
                elif model in ['gb', 'pb', 'rism std', 'rism gf']:
                    terms['energy'].append(model)
                df_models[model] = pd.DataFrame(self._energy2flatdict(data), index=index)
            elif model == 'ie':
                terms['entropy'].append(model)
                temp = {}
                for m, iedata in data.items():
                    df = pd.DataFrame({'data': iedata['data']}, index=index)
                    df1 = pd.DataFrame({'iedata': iedata['iedata'], 'sigma': iedata['sigma']},
                                       index=index[-iedata['ieframes']:])
                    temp[m] = pd.concat([df, df1], axis=1)
                df_models[model] = pd.concat(temp.values(), axis=1, keys=temp.keys())
            elif model == 'c2':
                terms['entropy'].append(model)
                temp = {m: pd.DataFrame({'c2data': c2data['c2data'], 'sigma': c2data['sigma']}, index=[0])
                        for m, c2data in data.items()}
                df_models[model] = pd.concat(temp.values(), axis=1, keys=temp.keys())

        # Calculate binding
        if terms['energy'] and terms['entropy']:
            df_models['binding'] = {}
            for m in terms['energy']:
                for e in terms['entropy']:
                    total = df_models[m]['delta']['TOTAL']
                    total.name = 'ΔH'
                    if e in ['nmode', 'qh']:
                        df_models['binding'][f"{m}+{e}"] = pd.concat(df_models[m], axis=1)
                        ent = df_models[e]['delta']['Total']
                    else:
                        k = 'iedata' if e == 'ie' else 'c2data'
                        ent = df_models[e][m][k]
                    ent.name = '-TΔS'
                    temp = pd.concat([total, ent], axis=1)
                    dg = pd.Series(total.mean() + ent.mean(), index=[total.index.tolist()[-1]], name='ΔG')
                    df_models['binding'][f"{m}+{e}"] = pd.concat([temp, dg], axis=1)
        return df_models

    def get_energy(self):

        raw_energy = self.get_raw_energy()

        energy = {}
        if raw_energy['normal']:
            energy['normal'] = self._get_df(raw_energy, 'normal')
        if raw_energy['mutant']:
            energy['mutant'] = self._get_df(raw_energy, 'mutant')
        if raw_energy['mutant-normal']:
            energy['mutant-normal'] = self._get_df(raw_energy, 'mutant-normal', True)

        if raw_energy['decomp_normal']:
            energy['decomp_normal'] = self._get_df(raw_energy, 'decomp_normal')
        if raw_energy['decomp_mutant']:
            energy['decomp_mutant'] = self._get_df(raw_energy, 'decomp_mutant')
            # print(energy['decomp_normal'])

        return energy

    def update_energy_frames(self, startframe, endframe, interval=1):

        self._get_frames()
        start = list(self.frames.keys()).index(startframe)
        end = list(self.frames.keys()).index(endframe) + 1

        # 6- get the number of frames
        self._update_frames(startframe, endframe, interval)

        # TODO:
        # 1- Iterar sobre cada clase que contiene la data y crear una copia (sino no se puede recobrar la data original)
        self.data = deepcopy(self._oringin)
        for key, v in self.data.items():
            # key: normal, mutant, decomp_normal, decomp_mutant
            if key in ['normal', 'mutant']:
                for model, v1 in v.items():
                    if model in ['gb', 'pb', 'rism std', 'rism gf']:
        # 2- Actualizar los frames y recalcular los composites
                        v1['complex'].set_frame_range(start, end, interval)
                        if not self.stability:
                            v1['receptor'].set_frame_range(start, end, interval)
                            v1['ligand'].set_frame_range(start, end, interval)
        # 3- Recalcular los deltas
                            v1['delta'] = BindingStatistics(v1['complex'], v1['receptor'], v1['ligand'],
                                                               self.app_namespace.INFO['using_chamber'],
                                                               self.traj_protocol)
        # 4- Recalcular las entropías basadas en GGAS
                            if self.app_namespace.INPUT['interaction_entropy']:
                                edata = v1['delta']['GGAS']
                                ie = InteractionEntropyCalc(edata, self.app_namespace.INPUT)
                                ieout = IEout(self.app_namespace.INPUT)
                                v['ie'][model] = ieout.parse_from_dict(key, {'data': ie.data, 'iedata': ie.iedata,
                                                                                  'ieframes': ie.ieframes,
                                                                                  'sigma': ie.ie_std})
                            if self.app_namespace.INPUT['c2_entropy']:
                                edata = v1['delta']['GGAS']
                                c2 = C2EntropyCalc(edata, self.app_namespace.INPUT)
                                v['c2'][model] = {'c2data': c2.c2data, 'sigma': c2.ie_std, 'c2_std': c2.c2_std,
                                                  'c2_ci': c2.c2_ci}
            elif key in ['decomp_normal', 'decomp_mutant']:
                for model, v1 in v.items():
                    v1['complex'].set_frame_range(start, end, interval)
                    if not self.stability:
                        v1['receptor'].set_frame_range(start, end, interval)
                        v1['ligand'].set_frame_range(start, end, interval)
                        # 3- Recalcular los deltas
                        if self.app_namespace.INPUT['idecomp'] in [1, 2]:
                            Decomp_Delta = DecompBinding
                        else:
                            Decomp_Delta = PairDecompBinding
                        v1['delta'] = Decomp_Delta(v1['complex'], v1['receptor'], v1['ligand'],
                                                   self.app_namespace.INPUT)


        # 5- Recalcular los delta delta sin CAS
        if 'normal' in self.data and self.data['normal'] and 'mutant' in self.data and self.data['mutant']:
            for model in self.data['normal']:
                if model in ['gb', 'pb', 'rism std', 'rism gf']:
                    self.data['mutant-normal'][model] = DeltaBindingStatistics(self.data['mutant'][model]['delta'],
                                                                     self.data['normal'][model]['delta'])


        # 6- Recalcular los summary
        self.get_summary()
        # 7- actualizar la data de salida
        # return self.data

    def get_summary(self):
        summary = {}

        for m, v in self.data.items():
            if m not in ['mutant-normal', 'mutant', 'normal']:
                continue
            for model, v1 in v.items():
                if m == 'mutant-normal':
                    # summary[m][model] = [v1.mean(), v1.std(), v1.std() /len(self.frames)]
                    continue
                else:
                    summary[(m, model)] = {}
                for mol, v2 in v1.items():
                    try:
                        if v2:
                            summ = v2.summary('csv')
                            summary[(m, model, (mol,))] = pd.DataFrame(summ[2:], columns=summ[1])
                            # summary[m][model][(mol,)] = v2.summary('csv')[1:]

                        else:
                            print(mol)
            #             summary[m][model][mol][(mol,)] = [v2.mean(), v2.stdev(), v2.std(), v2.stdev() / len(self.frames),
            #                                               v2.std() /len(self.frames)]
                    except:
                        pass
        # df = pd.DataFrame(summary['gb']['complex'], columns=['Category', 'Name', 'Marks'])

        return summary

    def get_gb_energy(self):
        energy = self.get_energy()
        gb_energy = {}
        if 'normal' in energy and 'gb' in energy['normal']:
            gb_energy['normal'] = energy['normal']['gb']
        if 'mutant' in energy and 'gb' in energy['mutant']:
            gb_energy['mutant'] = energy['mutant']['gb']
        return gb_energy

    def get_normal_gb_energy(self):
        energy = self.get_gb_energy()
        if 'normal' in energy:
            return energy['normal']

    def get_mutant_gb_energy(self):
        energy = self.get_gb_energy()
        if 'mutant' in energy:
            return energy['mutant']

    def _get_fromH5(self, h5file):
        h5file = H5Output(h5file)
        self.app_namespace = h5file.app_namespace
        self._oringin = {'normal': h5file.calc_types.normal, 'mutant': h5file.calc_types.mutant,
                         'decomp_normal': h5file.calc_types.decomp_normal, 'decomp_mutant': h5file.calc_types.decomp_mutant,
                         'mutant-normal': h5file.calc_types.mut_norm,
                         'decomp_mutant-normal': h5file.calc_types.decomp_mut_norm}
        self.data = copy(self._oringin)
        self._get_frames()

    def _get_fromApp(self, ifile):

        app = main.MMPBSA_App(MPI)
        info = infofile.InfoFile(app)
        info.read_info(ifile)
        app.normal_system = app.mutant_system = None
        app.parse_output_files()
        self.app_namespace = self._get_namespace(app)
        self._oringin = {'normal': app.calc_types.normal, 'mutant': app.calc_types.mutant,
                         'decomp_normal': app.calc_types.decomp_normal, 'decomp_mutant': app.calc_types.decomp_mutant,
                         'mutant-normal': app.calc_types.mut_norm,
                         'decomp_mutant-normal': app.calc_types.decomp_mut_norm}
        self.data = copy(self._oringin)
        self._get_frames()
        self._get_data(None)

    @staticmethod
    def _get_namespace(app):

        # input_file_text = ('|Input file:\n|--------------------------------------------------------------\n|'
        #                    + ''.join(open(app.FILES.input_file).readlines()).replace('\n', '\n|') +
        #                    '--------------------------------------------------------------\n')
        INFO = {'COM_PDB': ''.join(open(app.FILES.complex_fixed).readlines()),
                'input_file': app.input_file_text,
                'mutant_index': app.mutant_index,
                'mut_str': app.resl[app.mutant_index].mutant_label if app.mutant_index else '',
                'numframes': app.numframes,
                'numframes_nmode': app.numframes_nmode,
                'output_file': ''.join(open(app.FILES.output_file).readlines()),
                'size': app.mpi_size,
                'using_chamber': app.using_chamber}
        if app.INPUT['decomprun']:
            INFO['decomp_output_file'] = ''.join(open(app.FILES.decompout).readlines())

        return SimpleNamespace(FILES=app.FILES, INPUT=app.INPUT, INFO=INFO)

    def _get_frames(self):

        INPUT = self.app_namespace.INPUT
        numframes = self.app_namespace.INFO['numframes']

        frames_list = list(range(INPUT['startframe'], INPUT['startframe'] + numframes * INPUT['interval'],
                                 INPUT['interval']))
        INPUT['endframe'] = frames_list[-1]

        ts = 1 if not self.timestep else self.timestep
        time_step_list = list(range(self.starttime, self.starttime + len(frames_list) * ts, ts * INPUT['interval']))
        self.frames = dict(zip(frames_list, time_step_list))

    def _update_frames(self, startframe, endframe, interval):

        # get the original frames list
        self._get_frames()


        curr_frames = list(range(startframe, endframe + interval, interval))
        frames = list(self.frames.keys())
        for x in frames:
            if x not in curr_frames:
                self.frames.pop(x)

        print('updated frames', self.frames)

    def get_com(self):
        return copy(self.data['normal']['gb']['complex'])

    def _get_data(self, calc_types):
        """

        Args:
            calc_types:

        Returns:

        """
        # check the calc_types type
        INPUT = self.app_namespace.INPUT
        numframes = self.app_namespace.INFO['numframes']
        numframes_nmode = self.app_namespace.INFO['numframes_nmode']
        # See if we are doing stability
        self.stability = self.app_namespace.FILES.stability
        self.mutant_only = self.app_namespace.INPUT['mutant_only']
        self.traj_protocol = ('MTP' if self.app_namespace.FILES.receptor_trajs or
                                       self.app_namespace.FILES.ligand_trajs else 'STP')


        # Now load the data
        # if not INPUT['mutant_only']:
        #     self._get_edata(calc_types.normal)

    # def _get_edata(self, calc_type_data):
    #     from GMXMMPBSA.amber_outputs import GBout
    #     for key in calc_type_data:
    #
    #         if key == 'gb':
    #             com = GBout(None, self.app_namespace.INPUT, read=False)
    #             com.parse_from_h5(calc_type_data[key]['complex'])
    #             # com.get_energies_fromdict(calc_type_data[key]['complex'])
    #         # if key in ['ie', 'c2']:
    #         #     data[key] = {}
    #         #     if key == 'ie':
    #         #         for iekey in calc_type_data:
    #         #             data[key][iekey] = {
    #         #                 'data': pd.DataFrame(calc_type_data[iekey]['data'], columns=['data'],
    #         #                                      index=pd.Index(self.frames, names='frames')),
    #         #                 'iedata': pd.DataFrame(calc_type_data[iekey]['iedata'], columns=['iedata'],
    #         #                                        index=self.frames[-calc_type_data[iekey]['ieframes']:]),
    #         #                 'ieframes': calc_type_data[iekey]['ieframes'],
    #         #                 'sigma': calc_type_data[iekey]['sigma']}
    #         #     else:
    #         #         for c2key in calc_type_data:
    #         #             data[key][c2key] = {'c2data': calc_type_data[c2key]['c2data'],
    #         #                                 'sigma': calc_type_data[c2key]['sigma'],
    #         #                                 'c2_std': calc_type_data[c2key]['c2_std'],
    #         #                                 'c2_ci': calc_type_data[c2key]['c2_ci']}
    #         # else:
    #     #         # Since the energy models have the same structure as nmode and qh, we only need to worry about
    #     #         # correctly defining the Dataframe index, that is, the frames
    #     #         if key == 'nmode':
    #     #             cframes = self.nmode_frames
    #     #         elif key == 'qh':
    #     #             cframes = [0]
    #     #         else:
    #     #             cframes = self.frames
    #     #         # since the model data object in MMPBSA_App contain the data in the attribute data and H5 not,
    #     #         # we need to define a conditional object
    #     #         com_calc_type_data = (calc_types[key]['complex'] if h5
    #     #                               else calc_types[key]['complex'].data)
    #     #         df = complex = pd.DataFrame({dkey: com_calc_type_data[dkey] for dkey in com_calc_type_data},
    #     #                                     index=pd.Series(cframes, name='frames'))
    #     #         if not self.stability:
    #     #             rec_calc_type_data = (calc_types[key]['receptor'] if h5
    #     #                                   else calc_types[key]['receptor'].data)
    #     #             receptor = pd.DataFrame({dkey: rec_calc_type_data[dkey] for dkey in rec_calc_type_data},
    #     #                                     index=pd.Series(cframes, name='frames'))
    #     #             lig_calc_type_data = (calc_types[key]['ligand'] if h5
    #     #                                   else calc_types[key]['ligand'].data)
    #     #             ligand = pd.DataFrame({dkey: lig_calc_type_data[dkey] for dkey in lig_calc_type_data},
    #     #                                   index=pd.Series(cframes, name='frames'))
    #     #             delta = complex - receptor - ligand
    #     #             df = pd.concat([complex, receptor, ligand, delta], axis=1,
    #     #                            keys=['complex', 'receptor', 'ligand', 'delta'])
    #     #         data[key] = df
    #
    # def _get_ddata(self, calc_types, h5=False):
    #     data = {}
    #     # Take the decomp data
    #     for key in calc_types:
    #         # since the model data object in MMPBSA_App contain the data in the attribute data and H5 not,
    #         # we need to define a conditional object. Also, the decomp data must be re-structured for multiindex
    #         # Dataframe
    #         com_calc_type_data = (calc_types[key]['complex'] if h5
    #                               else self._transform_from_lvl_decomp(calc_types[key]['complex']))
    #         df = complex = pd.DataFrame(com_calc_type_data, index=self.frames)
    #         if not self.stability:
    #             rec_calc_type_data = (calc_types[key]['receptor'] if h5
    #                                   else self._transform_from_lvl_decomp(calc_types[key]['receptor']))
    #             receptor = pd.DataFrame(rec_calc_type_data, index=self.frames)
    #
    #             lig_calc_type_data = (calc_types[key]['ligand'] if h5
    #                                   else self._transform_from_lvl_decomp(calc_types[key]['ligand']))
    #             ligand = pd.DataFrame(lig_calc_type_data, index=self.frames)
    #
    #             delta = complex.subtract(pd.concat([receptor, ligand], axis=1)).combine_first(
    #                 complex).reindex_like(df)
    #             df = pd.concat([complex, receptor, ligand, delta], axis=1,
    #                            keys=['complex', 'receptor', 'ligand', 'delta'])
    #         data[key] = df
    #     return data

    @staticmethod
    def _energy2flatdict(nd, omit=None):
        if not omit:
            omit = []
        data = {}
        for k1, v1 in nd.items():  # Complex, receptor, ligand and delta
            if k1 in omit:
                continue
            for k2, v2 in v1.items():  # TDC, SDC, BDC or terms
                if isinstance(v2, dict):  # per-wise
                    for k3, v3 in v2.items():  # residue
                        for k4, v4 in v3.items():  # residue in per-wise or terms in per-res
                            if isinstance(v4, dict):  # per-wise
                                for k5, v5 in v4.items():
                                    data[(k1, k2, k3, k4, k5)] = v5
                            else:
                                data[(k1, k2, k3, k4)] = v4
                else:
                    data[(k1, k2)] = v2
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
       print mydata.mutant['pb']['receptor']['TOTAL'].stdev()
    """
    if not isinstance(fname, Path):
        fname = Path(fname)

    if not fname.exists():
        raise NoFileExists("cannot find %s!" % fname)
    os.chdir(fname.parent)
    d_mmpbsa = MMPBSA_API()
    d_mmpbsa.load_file(fname)
    return d_mmpbsa.get_energy(), d_mmpbsa.app_namespace
