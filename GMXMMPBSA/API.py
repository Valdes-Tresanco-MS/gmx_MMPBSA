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
import logging
import math
import pickle
import shutil
from multiprocessing.pool import ThreadPool
from copy import copy
from typing import Union

from GMXMMPBSA.calculation import InteractionEntropyCalc, C2EntropyCalc

from GMXMMPBSA import infofile, main, utils
from GMXMMPBSA.exceptions import NoFileExists
from GMXMMPBSA.fake_mpi import MPI
from GMXMMPBSA.amber_outputs import IEout, C2out
import pandas as pd
from pathlib import Path
import os
from types import SimpleNamespace

from GMXMMPBSA.utils import emapping, flatten, mask2list


def _remove_empty_charts(data):
    if isinstance(data, pd.Series):
        if (data.abs() < 0.01).all():
            return True
    elif (data.abs() < 0.01).all().all():
        return True


def _itemdata_properties(data, decomp=False):
    """
    Pre-processing the items data.
    Get the following properties:
    - separable: if contains subcategories (DH [energetics components], Per-residue[receptor and ligand])

    Also, remove empty terms and charts according to selected options
    @param data:
    @return:
    """
    groups = {}
    if not decomp:
        _extracted_from__itemdata_properties_5(data, groups)
    else:
        groups['Receptor'] = []
        groups['Ligand'] = []
        for k in data.columns:
            if k[0].startswith('R:') and k[0] not in groups['Receptor']:
                groups['Receptor'].append(k[0])
            elif k[0].startswith('L:') and k[0] not in groups['Ligand']:
                groups['Ligand'].append(k[0])
    return groups


# TODO Rename this here and in `_itemdata_properties`
def _extracted_from__itemdata_properties_5(data, groups):
    sep_ggas_keys = []
    sep_gsolv_keys = []
    # remove empty charts? (BOND, ANGLE and DIHEDRAL for STP)
    # FIXME: NLPBsolver ?
    ggas_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'UB', 'IMP', 'CMAP', 'ESCF']
    gsolv_keys = ['EGB', 'ESURF', 'EPB', 'ENPOLAR', 'EDISPER', 'POLAR SOLV', 'APOLAR SOLV', 'ERISM']
    for k in data.columns:
        if k in ggas_keys:
            sep_ggas_keys.append(k)
        elif k in gsolv_keys:
            sep_gsolv_keys.append(k)
    if sep_ggas_keys:
        groups['GGAS'] = sep_ggas_keys
        groups['GSOLV'] = sep_gsolv_keys
        groups['TOTAL'] = ['GGAS', 'GSOLV', 'TOTAL']


def _setup_data(data: pd.DataFrame, level=0, iec2=False, name=None, index=None,
                memory: dict = None, id=None):
    # this variable show if the data changed or not. At first time, it is true, then when plotting become in false
    change = True

    inmemory: bool = memory.get('inmemory', False)
    temp_path: Path = memory.get('temp_path')
    parquet_file = temp_path.joinpath('_'.join(id) + '_%s').as_posix()

    cont = {'line_plot_data': None, 'bar_plot_data': None, 'heatmap_plot_data': None}
    if level == 0:
        options = {'iec2': iec2}
        data.name = data.name[-1]
        line_plot_data = data[:-3].to_frame()
        if inmemory:
            cont['line_plot_data'] = [line_plot_data, options, change]
        else:
            line_plot_data.to_parquet(parquet_file % 'lp', compression=None)
            cont['line_plot_data'] = [parquet_file % 'lp', options, change]
    elif level == 1:
        options = ({'iec2': True} if iec2 else {}) | dict(groups=_itemdata_properties(data))
        bar_plot_data = data[-3:].reindex(columns=index)
        if inmemory:
            cont['bar_plot_data'] = [bar_plot_data, options, change]
        else:
            bar_plot_data.to_parquet(parquet_file % 'bp', compression=None)
            cont['bar_plot_data'] = [parquet_file % 'bp', options, change]
    elif level == 1.5:
        options = {'iec2': True}
        line_plot_data = data[['AccIntEnergy', 'ie']][:-3]
        bar_plot_data = data[['ie', 'sigma']][-3:]
        if inmemory:
            cont['line_plot_data'] = [line_plot_data, options, change]
            cont['bar_plot_data'] = [bar_plot_data, options, change]
        else:
            line_plot_data.to_parquet(parquet_file % 'lp', compression=None)
            bar_plot_data.to_parquet(parquet_file % 'bp', compression=None)
            cont['line_plot_data'] = [parquet_file % 'lp', options, change]
            cont['bar_plot_data'] = [parquet_file % 'bp', options, change]

    elif level == 2:
        tempdf = data.loc[:, data.columns.get_level_values(1) == 'tot'].droplevel(level=1, axis=1)
        temp_data = tempdf.reindex(columns=index)
        line_plot_data = temp_data[:-3].sum(axis=1).rename(name).to_frame()
        bar_plot_data = tempdf[-3:]
        heatmap_plot_data = tempdf[:-3].T
        if memory:
            cont['line_plot_data'] = [line_plot_data, {}, change]
            cont['bar_plot_data'] = [bar_plot_data, dict(groups=_itemdata_properties(bar_plot_data)), change]
            cont['heatmap_plot_data'] = [heatmap_plot_data, {}, change]
        else:
            line_plot_data.to_parquet(parquet_file % 'lp', compression=None)
            bar_plot_data.to_parquet(parquet_file % 'bp', compression=None)
            heatmap_plot_data.to_parquet(parquet_file % 'hp', compression=None)
            cont['line_plot_data'] = [parquet_file % 'lp', {}, change]
            cont['bar_plot_data'] = [parquet_file % 'bp', dict(groups=_itemdata_properties(bar_plot_data)), change]
            cont['heatmap_plot_data'] = [parquet_file % 'hp', {}, change]
    elif level == 3:
        # Select only the "tot" column, remove the level, change first level of columns to rows and remove the mean
        # index
        tempdf = data.loc[:, data.columns.get_level_values(2) == 'tot']
        line_plot_data = tempdf[:-3].groupby(axis=1, level=0, sort=False).sum().reindex(
            columns=index).sum(axis=1).rename(name).to_frame()
        bar_plot_data = tempdf[:-3].groupby(axis=1, level=0, sort=False).sum().agg(
            [lambda x: x.mean(), lambda x: x.std(ddof=0), lambda x: x.std(ddof=0)/ math.sqrt(len(x))]
            ).reindex(columns=index)
        bar_plot_data.index = ['Average', 'SD', 'SEM']
        heatmap_plot_data = tempdf.loc[["Average"]].droplevel(level=2, axis=1).stack().droplevel(level=0).reindex(
            columns=index, index=index)
        if memory:
            cont['line_plot_data'] = [line_plot_data, {}, change]
            cont['bar_plot_data'] = [bar_plot_data, dict(groups=_itemdata_properties(bar_plot_data)), change]
            cont['heatmap_plot_data'] = [heatmap_plot_data, {}, change]
        else:
            line_plot_data.to_parquet(parquet_file % 'lp', compression=None)
            bar_plot_data.to_parquet(parquet_file % 'bp', compression=None)
            heatmap_plot_data.to_parquet(parquet_file % 'hp', compression=None)
            cont['line_plot_data'] = [parquet_file % 'lp', {}, change]
            cont['bar_plot_data'] = [parquet_file % 'bp', dict(groups=_itemdata_properties(bar_plot_data)), change]
            cont['heatmap_plot_data'] = [parquet_file % 'hp', {}, change]

    return cont

def calculatestar(arg):
    func, args, id, t = arg
    data = args.get('data')
    level = args.get('level', 0)
    iec2 = args.get('iec2', False)
    name = args.get('name')
    index = args.get('index')
    memory = args.get('memory')
    return t, id, func(data, level, iec2, name, index, memory, id)


class MMPBSA_API():
    """ Main class that holds all the Free Energy data """

    def __init__(self, settings: dict = None, t=None):
        self.print_keys = None
        self.app_namespace = None
        self.raw_energy = None
        self.data = {}
        self.fname = None
        self.settings = settings
        self.temp_folder = None

        self._settings = {'data_on_disk': False, 'create_temporal_data': True, 'use_temporal_data': True, 
                          'overwrite_temp': True}
        if settings:
            not_keys = []
            for k, v in settings.items():
                if k not in self._settings:
                    not_keys.append(k)
                    continue
                self._settings[k] = v
            if not_keys:
                logging.warning(f'Not keys {not_keys}. Will be ignored')

    def setting_time(self, timestart=0, timestep=0, timeunit='ps'):
        self.starttime = timestart
        self.timestep = timestep
        self.timeunit = timeunit

    def setup_file(self, fname: Union[Path, str]):
        self.fname = fname if isinstance(fname, Path) else Path(fname)
        if not self.fname.exists():
            raise NoFileExists(f"cannot find {self.fname}!")
        os.chdir(self.fname.parent)

    def load_file(self, fname: Union[Path, str] = None):
        """

        Args:
            fname: String- or Path-like file reference

        Returns:

        """
        if not fname and not self.fname:
            raise

        if not self.fname:
            self.fname = fname if isinstance(fname, Path) else Path(fname)
            if not self.fname.exists():
                raise NoFileExists(f"cannot find {self.fname}!")
            os.chdir(self.fname.parent)

        if self.fname.suffix == '.mmxsa':
            self._get_fromBinary(self.fname)
        else:
            self._get_fromApp(self.fname)

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

    def _get_frames_index(self, framestype, startframe, endframe, interval):

        if framestype == 'energy':
            start = list(self.frames.keys()).index(startframe) if startframe else startframe
            end = list(self.frames.keys()).index(endframe) + 1 if endframe else endframe
            index_frames = {f: self.frames[f] for f in list(self.frames.keys())[start:end:interval]}
        else:
            start = list(self.nmframes.keys()).index(startframe) if startframe else startframe
            end = list(self.nmframes.keys()).index(endframe) + 1 if endframe else endframe
            index_frames = {f: self.nmframes[f] for f in list(self.nmframes.keys())[start:end:interval]}
        return (start, end, pd.Series(index_frames.values(), name=f'Time ({self.timeunit})') if self.timestep
                else pd.Series(index_frames.keys(), name='Frames'))

    @staticmethod
    def arg2tuple(arg):
        return arg if isinstance(arg, tuple) else tuple(arg)

    def get_energy(self, etype: tuple = None, model: tuple = None, mol: tuple = None, term: tuple = None,
                    remove_empty_terms=True, threshold=0.01, startframe=None, endframe=None, interval=1):
        """
        Get energy
        Args:
            keys: Energy levels

        Returns: Energy pd.Dataframe

        """
        # get start and end for frames range
        if etype:
            etype = self.arg2tuple(etype)
        if model:
            model = self.arg2tuple(model)
        if mol:
            mol = self.arg2tuple(mol)
        if term:
            term = self.arg2tuple(term)
        s, e, index = self._get_frames_index('energy', startframe, endframe, interval)
        e_map = {}
        energy = {}
        summ_df = {}
        comp_summary = {}
        temp_print_keys = etype or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data)

        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not etype {et} in data')
            elif not self.data[et]:
                continue
            else:
                e_map[et] = {}
                energy[et] = {}
                summ_df[et] = {}
                comp_summary[et] = {}
                temp_model_keys = model or tuple(x for x in self.data[et].keys() if x not in ['nmode', 'qh', 'ie','c2'])

                for m in temp_model_keys:
                    if m not in self.data[et]:
                        print(f'Not model {m} in etype {et}')
                    else:
                        e_map[et][m] = {}
                        model_energy = {}
                        temp_mol_keys = mol or tuple(self.data[et][m].keys())

                        for m1 in temp_mol_keys:
                            if m1 not in self.data[et][m]:
                                print(f'Not mol {m1} in etype {et} > model {m}')
                            else:
                                e_map[et][m][m1] = []
                                model_energy[m1] = {}
                                temp_terms_keys = ([x for x in self.data[et][m][m1].keys() if x in term] if term
                                                   else tuple(self.data[et][m][m1].keys()))
                                valid_terms = []

                                for t in temp_terms_keys:
                                    if t not in self.data[et][m][m1]:
                                        print(f'Not term {t} in etype {et} > model {m} > mol {m1}')
                                    elif (
                                        not remove_empty_terms
                                        or abs(self.data[et][m][m1][t].mean())
                                        >= threshold or t in ['GSOLV', 'GGAS', 'TOTAL']
                                    ):
                                        e_map[et][m][m1].append(t)
                                        model_energy[m1][t] = self.data[et][m][m1][t][s:e:interval]
                                        valid_terms.append(t)
                        energy[et][m], summ_df[et][m] = self._model2df(model_energy, index)

        corr = {}
        for et, v in summ_df.items():
            corr[et] = {}
            for m, v2 in v.items():
                if 'delta' in v2:
                    corr[et][m] = pd.DataFrame(
                        {('ΔGeff', t): v for t, v in v2['delta']['TOTAL'].to_dict().items()}, index=[0])
        return {'map': e_map, 'data': energy, 'summary': summ_df, 'correlation': corr}

    def _model2df(self, energy, index):
        energy_df = pd.DataFrame(flatten(energy), index=index)
        s = pd.concat([energy_df.mean(), energy_df.std(ddof=0), energy_df.std(ddof=0) / math.sqrt(len(index))], axis=1)
        s.columns = ['Average', 'SD', 'SEM']
        summary_df = s.T
        df = pd.concat([energy_df, summary_df])
        df.index.name = index.name
        return df, summary_df

    def get_nmode_entropy(self, nmtype: tuple = None, mol: tuple = None, term: tuple = None,
                          startframe=None, endframe=None, interval=None):

        s, e, index = self._get_frames_index('entropy', startframe, endframe, interval)

        energy = {}
        summ_df = {}
        temp_print_keys = nmtype or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data and
                                          self.data[x])
        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not nmtype {et} in data')
            else:
                energy[et] = {'nmode': {}}
                summ_df[et] = {'nmode': {}}
                model_energy = {}
                temp_mol_keys = mol or tuple(self.data[et]['nmode'].keys())
                for m1 in temp_mol_keys:
                    if m1 not in self.data[et]['nmode']:
                        print(f'Not mol {m1} in etype {et}')
                    else:
                        model_energy[m1] = {}
                        temp_terms_keys = ([x for x in self.data[et]['nmode'][m1].keys() if x in term] if term
                                           else tuple(self.data[et]['nmode'][m1].keys()))
                        valid_terms = []
                        for t in temp_terms_keys:
                            if t not in self.data[et]['nmode'][m1]:
                                print(f'Not term {t} in etype {et} > mol {m1}')
                            else:
                                model_energy[m1][t] = self.data[et]['nmode'][m1][t][s:e:interval]
                                valid_terms.append(t)
                energy[et]['nmode'], summ_df[et]['nmode'] = self._model2df(model_energy, index)

        return {'map': emapping(energy), 'data': energy, 'summary': summ_df}

    def get_qh_entropy(self, qhtype: tuple = None, mol: tuple = None, term: tuple = None):
        temp_print_keys = qhtype or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data)
        entropy = {}
        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not qhtype {et} in data')
            else:
                entropy[et] = {'qh': {}}
                temp_mol_keys = mol or tuple(self.data[et]['qh'].keys())
                for m1 in temp_mol_keys:
                    if m1 not in self.data[et]['qh']:
                        print(f'Not mol {m1} in qhtype {et}')
                    else:
                        entropy[et]['qh'][m1] = {}
                        temp_terms_keys = ([x for x in self.data[et]['qh'][m1].keys() if x in term] if term
                                           else tuple(self.data[et]['qh'][m1].keys()))
                        for t in temp_terms_keys:
                            if t not in self.data[et]['qh'][m1]:
                                print(f'Not term {t} in etype {et} > mol {m1}')
                            else:
                                entropy[et]['qh'][m1][t] = self.data[et]['qh'][m1][t]

        df = pd.DataFrame(flatten(entropy))
        return {'map': emapping(entropy), 'data': df, 'summary': df.xs(('delta', 'TOTAL'), level=[1, 2], axis=1)}

    def get_c2_entropy(self, c2type: tuple = None, startframe=None, endframe=None, interval=None):
        temp_print_keys = c2type or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data and
                                          self.data[x])
        entropy = {}
        entropy_df = {}
        recalc = bool((startframe and startframe != self.app_namespace.INPUT['startframe'] or
                       endframe and endframe != self.app_namespace.INPUT['endframe'] or
                       interval and interval != self.app_namespace.INPUT['interval']))

        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not c2type {et} in data')
            if recalc:
                d = self._recalc_iec2('c2', et, startframe, endframe, interval)
            else:
                d = self.data[et]['c2']
            entropy[et] = {'c2': {x: None for x in ['c2', 'sigma']}}
            entropy_df[et] = {'c2': pd.DataFrame({'c2': [d['c2data'], d['c2_std'], d['c2_std']], 'sigma': [d['sigma'], 0, 0]},
                                                 index=['Average', 'SD', 'SEM'])}

        return {'map': emapping(entropy), 'data': entropy_df, 'summary': entropy_df}

    def get_ie_entropy(self, ietype: tuple = None, startframe=None, endframe=None, interval=None,
                       ie_segment = 25):
        temp_print_keys = ietype or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data and
                                          self.data[x])
        entropy = {}
        summ_df = {}
        entropy_df = {}
        recalc = bool((startframe and startframe != self.app_namespace.INPUT['startframe'] or
                       endframe and endframe != self.app_namespace.INPUT['endframe'] or
                       interval and interval != self.app_namespace.INPUT['interval']))

        s, e, index = self._get_frames_index('energy', startframe, endframe, interval)
        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not ietype {et} in data')
                continue
            if recalc:
                d = self._recalc_iec2('ie', et, startframe, endframe, interval, ie_segment)
            else:
                d = self.data[et]['ie']
            ieframes = math.ceil(len(d['data']) * ie_segment / 100)
            entropy[et] = {'ie': {x: None for x in ['AccIntEnergy', 'ie', 'sigma']}}
            df = pd.DataFrame({'AccIntEnergy': d['data']}, index=index)
            df1 = pd.DataFrame({'ie': d['data'][-ieframes:]}, index=index[-ieframes:])
            df2 = pd.concat([df, df1], axis=1)
            df3 = pd.DataFrame({'ie': [float(d['iedata'].mean()),
                                       float(d['iedata'].std()),
                                       float(d['iedata'].std() / math.sqrt(ieframes))],
                                'sigma': [d['sigma'], 0, 0]}, index=['Average', 'SD', 'SEM'])
            summ_df[et] = {'ie': df3}
            df4 = pd.concat([df2, df3])
            df4.index.name = df.index.name
            entropy_df[et] = {'ie': df4}
        return {'map': emapping(entropy), 'data': entropy_df, 'summary': summ_df}

    @staticmethod
    def _merge_ent(res_dict: dict, d: dict):
        for e1, v1 in d.items():
            if e1 not in res_dict:
                res_dict[e1] = v1
            else:
                res_dict[e1].update(v1)

    def get_entropy(self, etype: tuple = None, model: tuple = None, mol: tuple = None, term: tuple = None,
                    nmstartframe=None, nmendframe=None, nminterval=1, startframe=None, endframe=None, interval=None,
                    ie_segment=25):
        temp_print_keys = etype or tuple(x for x in ['normal', 'mutant', 'mutant-normal'] if x in self.data)
        entropy_keys = []

        for et in temp_print_keys:
            if et not in self.data:
                print(f'Not etype {et} in data')
            else:
                entropy_keys.extend(x for x in self.data[et].keys() if x in ['nmode', 'qh', 'ie','c2'])
        temp_model_keys = tuple(x for x in entropy_keys if x in model) if model else tuple(entropy_keys)

        ent_summ = {}
        ent_map = {}
        ent_data = {}
        for m in temp_model_keys:
            if m == 'nmode':
                _m, _d, _s = self.get_nmode_entropy(etype, mol, term, nmstartframe, nmendframe, nminterval).values()
            elif m == 'qh':
                _m, _d, _s = self.get_qh_entropy(etype, mol, term).values()
            elif m == 'ie':
                _m, _d, _s = self.get_ie_entropy(etype, startframe, endframe, interval, ie_segment).values()
            else:
                _m, _d, _s = self.get_c2_entropy(etype, startframe, endframe, interval).values()
            self._merge_ent(ent_map, _m)
            self._merge_ent(ent_data, _d)
            self._merge_ent(ent_summ, _s)

        return {'map': ent_map, 'data': ent_data, 'summary': ent_summ}

    def _recalc_iec2(self, method, etype, startframe=None, endframe=None, interval=None, ie_segment=25):
        allowed_met = ['gb', 'pb', 'rism std', 'rism gf', 'rism pcplus', 'gbnsr6']
        result = None
        start = list(self.frames.keys()).index(startframe) if startframe else startframe
        end = list(self.frames.keys()).index(endframe) + 1 if endframe else endframe
        for key in allowed_met:
            if key in self.data[etype]:
                edata = self.data[etype][key]['delta']['GGAS'][start:end:interval]
                if method == 'ie':
                    ie = InteractionEntropyCalc(edata,
                                                dict(temperature=self.app_namespace.INPUT['temperature'],
                                                     startframe=startframe, endframe=endframe, interval=interval),
                                                iesegment=ie_segment)
                    result = IEout({})
                    result.parse_from_dict(dict(data=ie.data, sigma=ie.ie_std, iedata=ie.iedata))
                else:
                    c2 = C2EntropyCalc(edata, dict(temperature=self.app_namespace.INPUT['temperature']))
                    result = C2out()
                    result.parse_from_dict(dict(c2data=c2.c2data, c2_std=c2.c2_std, sigma=c2.ie_std, c2_ci=c2.c2_ci))
                break
        return result

    def get_binding(self, energy_summary=None, entropy_summary=None):
        binding  = {}
        b_map = {}
        corr = {}
        if energy_summary and entropy_summary:
            dg_string = {'ie': 'ΔGie', 'c2': 'ΔGc2', 'nmode': 'ΔGnm', 'qh': 'ΔGqh'}
            for et, ev in energy_summary.items():
                binding[et] = {}
                b_map[et] = {}
                corr[et] = {}
                for em, emv in ev.items():
                    binding[et][em] = {}
                    b_map[et][em] = []
                    mol = 'complex' if self.app_namespace.FILES.stability else 'delta'
                    edata = emv[mol]['TOTAL']
                    edata.name = 'ΔH'
                    if et in entropy_summary:
                        if mol == 'delta':
                            corr[et][em] = {}
                        c = {}
                        for ent, etv in entropy_summary[et].items():
                            b_map[et][em].append(ent)
                            if ent == 'ie' and not self.app_namespace.FILES.stability:
                                entdata = etv['ie']
                            elif ent == 'c2' and not self.app_namespace.FILES.stability:
                                entdata = etv['c2']
                            else:
                                entdata = etv[mol]['TOTAL']
                            entdata.name = '-TΔS'
                            dg = edata.loc['Average'] + entdata.loc['Average']
                            std = utils.get_std(edata.loc['SD'], entdata.loc['SD'])
                            dgdata = pd.Series([dg, std, std], index=['Average', 'SD', 'SEM'], name='ΔG')
                            binding[et][em][ent] = pd.concat([edata, entdata, dgdata], axis=1)
                            if mol == 'delta':
                                c[(dg_string[ent], 'Average')] = dg
                                c[(dg_string[ent], 'SD')] = std
                                c[(dg_string[ent], 'SEM')] = std
                        if mol == 'delta':
                            corr[et][em] = pd.DataFrame(c, index=[0])
        return {'map': b_map, 'data': binding, 'correlation': corr}

    def get_decomp_energy(self, etype: tuple = None, model: tuple=None, mol: tuple = None, contribution: tuple = None,
                          res1: tuple = None, res2: tuple = None, term: tuple = None, res_threshold=0.5,
                          startframe=None, endframe=None, interval=None):

        s, e, index = self._get_frames_index('energy', startframe, endframe, interval)
        name = index.name
        index = pd.concat([index, pd.Series(['Average', 'SD', 'SEM'])])
        index.name = name

        temp_print_keys = etype or tuple(x for x in ['decomp_normal', 'decomp_mutant'] if self.data.get(x))
        self.print_keys = []

        decomp_energy = {}
        e_map = {}
        for et in temp_print_keys:
            if not self.data.get(et):
                print(f'Not etype {et} in data')
            else:
                etkey = et.split('_')[1]
                decomp_energy[etkey] = {}
                e_map[etkey] = {}
                temp_model_keys = model or tuple(self.data[et].keys())

                for m in temp_model_keys:
                    if m not in self.data[et]:
                        print(f'Not model {m} in etype {et}')
                    else:
                        model_decomp_energy = {}

                        e_map[etkey][m] = {}
                        temp_mol_keys = mol or tuple(self.data[et][m].keys())

                        for m1 in temp_mol_keys:
                            if m1 not in self.data[et][m]:
                                print(f'Not mol {m1} in etype {et} > model {m}')
                            else:
                                model_decomp_energy[m1] = {}
                                e_map[etkey][m][m1] = {}
                                temp_comp_keys = contribution or tuple(self.data[et][m][m1].keys())
                                for c in temp_comp_keys:
                                    if c not in self.data[et][m][m1]:
                                        print(f'Not component {c} in etype {et} > model {m} > mol {m1}')
                                    else:
                                        model_decomp_energy[m1][c] = {}
                                        e_map[etkey][m][m1][c] = {}
                                        temp_res1_keys = res1 or tuple(self.data[et][m][m1][c].keys())

                                        for r1 in temp_res1_keys:
                                            remove = False
                                            if r1 not in self.data[et][m][m1][c]:
                                                print(f'Not res {r1} in etype {et} > model {m} > mol {m1} > comp {c} ')
                                            else:
                                                temp_energy = {}

                                                if self.app_namespace.INPUT['idecomp'] in [1, 2]:
                                                    temp_emap = []
                                                    temp_terms_keys = term or tuple(self.data[et][m][m1][c][r1].keys())

                                                    for t in temp_terms_keys:
                                                        if t not in self.data[et][m][m1][c][r1]:
                                                            print(f'Not term {t} in etype {et} > model {m} > mol {m1} '
                                                                  f'> comp {c} > res {r1}')
                                                        else:
                                                            temp_energy[t] = self.data[et][m][m1][c][r1][t][s:e:interval]
                                                            temp_emap.append(t)
                                                            mean = temp_energy[t].mean()
                                                            std = temp_energy[t].std()
                                                            sem = std / math.sqrt(len(temp_energy[t]))
                                                            temp_energy[t] = temp_energy[t].append([mean, std, sem])
                                                            if (t == 'tot' and res_threshold > 0 and
                                                                    abs(mean) < res_threshold):
                                                                remove = True
                                                else:
                                                    temp_emap = {}
                                                    temp_res2_keys = res2 or tuple(self.data[et][m][m1][c][r1].keys())
                                                    # for per-wise only since for per-residue we get the tot value
                                                    res1_contrib = 0
                                                    for r2 in temp_res2_keys:
                                                        if r2 not in self.data[et][m][m1][c][r1]:
                                                            print(f'Not res {r2} in etype {et} > model {m} > mol {m1} '
                                                                  f'> comp {c} > res {r1}')
                                                        else:
                                                            temp_energy_r2 = {}
                                                            temp_emap[r2] = []
                                                            # energy[et][m][m1][c][r1][r2] = {}
                                                            temp_terms_keys = term or tuple(self.data[et][m][m1][c][r1][r2].keys())
                                                            for t in temp_terms_keys:
                                                                if t not in self.data[et][m][m1][c][r1][r2]:
                                                                    print(
                                                                        f'Not term {t} in etype {et} > model {m} '
                                                                        f'> mol {m1} > comp {c} > res {r1} > res {r2}')
                                                                else:
                                                                    temp_emap[r2].append(t)
                                                                    temp_energy_r2[t] = self.data[et][m][m1][c][r1][
                                                                                            r2][t][s:e:interval]
                                                                    mean = temp_energy_r2[t].mean()
                                                                    std = temp_energy_r2[t].std()
                                                                    sem = std / math.sqrt(len(temp_energy_r2[t]))
                                                                    temp_energy_r2[t] = temp_energy_r2[t].append([
                                                                        mean, std, sem])
                                                                    if t == 'tot':
                                                                        res1_contrib += mean
                                                            if temp_energy_r2:
                                                                temp_energy[r2] = temp_energy_r2
                                                    if res_threshold > 0 and float(abs(res1_contrib)) < res_threshold:
                                                        remove = True
                                                if not remove:
                                                    model_decomp_energy[m1][c][r1] = temp_energy
                                                    e_map[etkey][m][m1][c][r1] = temp_emap
                        decomp_energy[etkey][m] = pd.DataFrame(flatten(model_decomp_energy), index=index)
        return {'map': e_map, 'data': decomp_energy}

    def get_ana_data(self, energy_options=None, entropy_options=None, decomp_options=None, performance_options=None,
                     correlation=False):
        TASKs = []
        d = {}
        corr = {}
        ecorr = None
        memory_args = dict(temp_path=self.temp_folder)
        if energy_options:
            energy_map, energy, energy_summary, ecorr = self.get_energy(**energy_options).values()
            if energy_map:
                for et, v in ecorr.items():
                    corr[et] = v
                d['enthalpy'] = {'map': energy_map, 'keys': {}, 'summary': energy_summary}
                memory_args['inmemory'] = performance_options.get('energy_memory')
                for level, value in energy_map.items():
                    for level1, value1 in value.items():
                        for level2, value2 in value1.items():
                            TASKs.append([_setup_data, dict(data=energy[level][level1][level2][value2], level=1,
                                                            memory=memory_args),
                                          (level, level1, level2), 'enthalpy'])
                            TASKs.extend([_setup_data, dict(data=energy[level][level1][(level2, level3)], level=0,
                                                            memory=memory_args),
                                          (level, level1, level2,level3), 'enthalpy'] for level3 in value2)

        if entropy_options:
            entropy_map, entropy, entropy_summary = self.get_entropy(**entropy_options).values()
            if entropy_map:
                d['entropy'] = {'map': entropy_map, 'keys': {}, 'summary': entropy_summary}
                memory_args['inmemory'] = performance_options.get('energy_memory')
                for level, value in entropy_map.items():
                    for level1, value1 in value.items():
                        if level1 in ['nmode', 'qh']:
                            for level2, value2 in value1.items():
                                TASKs.append([_setup_data, dict(data=entropy[level][level1][level2], level=1,
                                                                memory=memory_args),
                                              (level, level1, level2), 'entropy'])
                                TASKs.extend([_setup_data, dict(data=entropy[level][level1][(level2, level3)],
                                                                level=0, memory=memory_args),
                                              (level, level1, level2, level3), 'entropy'] for level3 in value2)
                        elif level1 == 'c2':
                            TASKs.append([_setup_data, dict(data=entropy[level][level1], level=1, iec2=True,
                                                            memory=memory_args),
                                          (level, level1), 'entropy'])
                        elif level1 == 'ie':
                            TASKs.append([_setup_data, dict(data=entropy[level][level1], level=1.5,
                                                            memory=memory_args),
                                          (level, level1), 'entropy'])

        if energy_options and entropy_options:
            bind_map, binding, bcorr = self.get_binding(energy_summary, entropy_summary).values()
            if bind_map:
                for et, v in bcorr.items():
                    for m, v1 in v.items():
                        if ecorr:
                            corr[et][m] = pd.concat([ecorr[et][m], v1], axis=1)
                    
                d['binding'] = {'map': bind_map, 'keys': {}, 'summary': binding}
                memory_args['inmemory'] = performance_options.get('energy_memory')
                for level, value in bind_map.items():
                    for level1, value1 in value.items():
                        TASKs.extend([_setup_data, dict(data=binding[level][level1][level2], level=1,
                                                        memory=memory_args),
                                      (level, level1, level2), 'binding'] for level2 in value1)
        if decomp_options:
            decomp_map, decomp = self.get_decomp_energy(**decomp_options).values()
            if decomp_map:
                d['decomposition'] = {'map': decomp_map, 'keys': {}}
                memory_args['inmemory'] = performance_options.get('decomp_memory')
                for level, value in decomp_map.items():
                    for level1, value1 in value.items():
                        for level2, value2 in value1.items():
                            for level3, value3 in value2.items():
                                item_lvl = 2 if self.app_namespace.INPUT['idecomp'] in [1, 2] else 3
                                index = list(value3.keys())
                                for level4, value4 in value3.items():
                                    if item_lvl == 3 and len(index) == 1:
                                        index.append(list(value4.keys()))
                                    item_lvl2 = 1 if self.app_namespace.INPUT['idecomp'] in [1, 2] else 2
                                    if self.app_namespace.INPUT['idecomp'] in [1, 2]:
                                        TASKs.append(
                                            [_setup_data, dict(data=decomp[level][level1][level2][level3][level4],
                                                               level=item_lvl2, name=level4, index=value4,
                                                               memory=memory_args),
                                             (level, level1, level2, level3, level4), 'decomposition'])
                                        TASKs.extend([_setup_data,
                                                      dict(data=decomp[level][level1][level2][level3][level4][level5],
                                                           memory=memory_args),
                                                    (level, level1, level2, level3, level4, level5), 'decomposition']
                                                     for level5 in value4)
                                    else:
                                        TASKs.append(
                                            [_setup_data, dict(data=decomp[level][level1][level2][level3][level4],
                                                               level=item_lvl2, name=level4, index=list(value4.keys()),
                                                               memory=memory_args),
                                             (level, level1, level2, level3, level4), 'decomposition'])
                                        TASKs.extend([_setup_data,
                                                      dict(data=decomp[level][level1][level2][level3][level4][level5],
                                                           level=1, index=value5, memory=memory_args),
                                                      (level, level1, level2, level3, level4, level5), 'decomposition']
                                                     for level5, value5 in value4.items())
                                TASKs.append([_setup_data,
                                              dict(data=decomp[level][level1][level2][level3], level=item_lvl,
                                                   name=level3, index=index, memory=memory_args),
                                          (level, level1, level2, level3), 'decomposition'])

        with ThreadPool(performance_options.get('jobs')) as pool:
            imap_unordered_it = pool.imap_unordered(calculatestar, TASKs)
            for t, id, result in imap_unordered_it:
                d[t]['keys'][id] = result

        if correlation:
            d['correlation'] = corr

        return d

    def _get_fromApp(self, ifile):

        app = main.MMPBSA_App(MPI)
        info = infofile.InfoFile(app)
        info.read_info(ifile)
        app.normal_system = app.mutant_system = None
        app.parse_output_files(from_calc=False)
        self.app_namespace = self._get_namespace(app, 'App')
        self._oringin = {'normal': app.calc_types.normal, 'mutant': app.calc_types.mutant,
                         'decomp_normal': app.calc_types.decomp_normal, 'decomp_mutant': app.calc_types.decomp_mutant,
                         'mutant-normal': app.calc_types.mut_norm,
                         }
        self._finalize_reading(ifile)

    def _get_fromBinary(self, ifile):
        with open(ifile, 'rb') as bf:
            info = pickle.load(bf)
            self.app_namespace = self._get_namespace(info, 'Binary')
            bdata = pickle.load(bf)
            self._oringin = {'normal': bdata.normal, 'mutant': bdata.mutant, 'decomp_normal': bdata.decomp_normal,
                             'decomp_mutant': bdata.decomp_mutant, 'mutant-normal': bdata.mut_norm}
            self._finalize_reading(ifile)

    def _finalize_reading(self, ifile):
        self.temp_folder = ifile.parent.joinpath('.gmx_mmpbsa_temp')
        if self.temp_folder.exists():
            shutil.rmtree(self.temp_folder)
        self.temp_folder.mkdir()
        self.data = copy(self._oringin)
        self._get_frames()

    @staticmethod
    def _get_namespace(app, tfile):

        if tfile == 'App':
            com_pdb = ''.join(open(app.FILES.complex_fixed).readlines())
            input_file = app.input_file_text
            output_file = ''.join(open(app.FILES.output_file).readlines())
            decomp_output_file = ''.join(open(app.FILES.decompout).readlines()) if app.INPUT['decomprun'] else None
            size = app.mpi_size
        else:
            com_pdb = app.COM_PDB
            input_file = app.input_file
            output_file = app.output_file
            decomp_output_file = app.decomp_output_file if app.INPUT['decomprun'] else None
            size = app.size
            app.resl = mask2list(app.FILES.complex_fixed, app.INPUT['receptor_mask'], app.INPUT['ligand_mask'])

        INFO = {'COM_PDB': com_pdb,
                'input_file': input_file,
                'mutant_index': app.mutant_index,
                'mut_str': app.resl[app.mutant_index].mutant_label if app.mutant_index else '',
                'numframes': app.numframes,
                'numframes_nmode': app.numframes_nmode,
                'output_file': output_file,
                'size': size,
                'using_chamber': app.using_chamber,
                'decomp_output_file': decomp_output_file}

        return SimpleNamespace(FILES=app.FILES, INPUT=app.INPUT, INFO=INFO)

    def _get_frames(self):

        ts = self.timestep or 1

        INPUT = self.app_namespace.INPUT
        numframes = self.app_namespace.INFO['numframes']
        nmnumframes = self.app_namespace.INFO['numframes_nmode']

        frames_list = list(range(INPUT['startframe'],
                                 INPUT['startframe'] + numframes * INPUT['interval'],
                                 INPUT['interval']))
        self.frames_list = frames_list
        INPUT['endframe'] = frames_list[-1]
        time_step_list = list(range(self.starttime,
                                    self.starttime + len(frames_list) * ts * INPUT['interval'],
                                    ts * INPUT['interval']))
        self.frames = dict(zip(frames_list, time_step_list))

        if INPUT['nmoderun']:
            nmframes_list = list(range(INPUT['nmstartframe'],
                                       INPUT['nmstartframe'] + nmnumframes * INPUT['nminterval'],
                                       INPUT['interval']))
            INPUT['nmendframe'] = nmframes_list[-1]

            nm_start = (nmframes_list[0] - frames_list[0]) * INPUT['interval']
            nmtime_step_list = list(range(self.starttime + nm_start,
                                          self.starttime + nm_start + len(nmframes_list) * ts * INPUT['nminterval'],
                                          ts * INPUT['nminterval']))
            self.nmframes_list = nmframes_list
            self.nmframes = dict(zip(nmframes_list, nmtime_step_list))

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
       3D-RISM (PCPLUS)    |'rism pcplus'|
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
