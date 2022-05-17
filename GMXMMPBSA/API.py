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

from GMXMMPBSA.utils import emapping, flatten


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


def _setup_data(data, level=0, iec2=False, name=None, index=None):
    # this variable show if the data changed or not. At first time, it is true, then when plotting become in false
    change = True

    cont = {'ie_plot_data': None, 'line_plot_data': None, 'bar_plot_data': None, 'heatmap_plot_data': None}
    if level == 0:
        options = {'iec2': iec2}
        # if isinstance(data, pd.Series):
        data.name = data.name[-1]
        cont['line_plot_data'] = [data[:-2], options, change]
    elif level == 1:
        options = ({'iec2': True} if iec2 else {}) | dict(groups=_itemdata_properties(data))
        cont['bar_plot_data'] = [data[-2:].reindex(columns=index), options, change]
    elif level == 1.5:
        options = {'iec2': True}
        cont['line_plot_data'] = [data[['AccIntEnergy', 'ie']][:-2], options, change]
        cont['bar_plot_data'] = [data[['ie', 'sigma']][-2:], options, change]
    elif level == 2:
        tempdf = data.loc[:, data.columns.get_level_values(1) == 'tot']
        bar_plot_data = tempdf.droplevel(level=1, axis=1).reindex(columns=index)

        cont['line_plot_data'] = [bar_plot_data[:-2].sum(axis=1).rename(name), {}, change]
        cont['bar_plot_data'] = [bar_plot_data[-2:], dict(groups=_itemdata_properties(bar_plot_data)), change]
        cont['heatmap_plot_data'] = [bar_plot_data[:-2].T, {}, change]
    elif level == 3:
        # Select only the "tot" column, remove the level, change first level of columns to rows and remove the mean
        # index
        from icecream import ic

        tempdf = data.loc[:, data.columns.get_level_values(2) == 'tot']
        cont['heatmap_plot_data'] = [
            tempdf.loc[["Average"]].droplevel(level=2, axis=1).stack().droplevel(level=0).reindex(columns=index[0],
                                                                                                  index=index[1]),
            {}, change]
        bar_plot_data = tempdf.groupby(axis=1, level=0, sort=False).sum().reindex(columns=index[0])
        cont['line_plot_data'] = [bar_plot_data[:-2].sum(axis=1).rename(name), {}, change]
        cont['bar_plot_data'] = [bar_plot_data[-2:], dict(groups=_itemdata_properties(bar_plot_data)), change]
        # del tempdf
        # del bar_plot_data
        # del temp_bar_data

    return cont

def calculatestar(arg):
    func, args, id, t = arg
    data = args.get('data')
    level = args.get('level', 0)
    iec2 = args.get('iec2', False)
    name = args.get('name')
    index = args.get('index')
    return t, id, func(data, level, iec2, name, index)


class MMPBSA_API():
    """ Main class that holds all the Free Energy data """

    def __init__(self):
        self.print_keys = None
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
        if not fname and not self.fname:
            raise

        if not self.fname:
            self.fname = fname if isinstance(fname, Path) else Path(fname)
            if not self.fname.exists():
                raise NoFileExists(f"cannot find {self.fname}!")
            os.chdir(self.fname.parent)

        # if self.fname.suffix == '.h5':
        #     self._get_fromH5(fname)
        # else:
        self._get_fromApp(self.fname)
        # print('API', self.app_namespace)

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

        return {'map': e_map, 'data': energy, 'summary': summ_df}

    def _model2df(self, energy, index):
        energy_df = pd.DataFrame(flatten(energy), index=index)
        s = pd.concat([energy_df.mean(), energy_df.std(ddof=0)], axis=1)
        s.columns = ['Average', 'SD']
        summary_df = s.T
        df = pd.concat([energy_df, summary_df])
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
                energy[et]['nmode'], summ_df[et]['nmode'] = self._model2df('entropy', startframe, endframe, interval,
                                                                           model_energy)

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
            elif recalc:
                d = self._recalc_iec2('c2', et, startframe, endframe, interval)
                entropy[et] = {'c2': {x: None for x in ['c2', 'sigma']}}
                entropy_df[et] = {'c2': pd.DataFrame({'c2': [d['c2data'], d['c2_std']], 'sigma': [d['sigma'], 0]},
                                                     index=['Average', 'SD'])}
            else:
                entropy[et] = {'c2': {x: None for x in ['c2', 'sigma']}}
                entropy_df[et] = {'c2': pd.DataFrame({'c2': [self.data[et]['c2']['c2data'],
                                                                self.data[et]['c2']['c2_std']],
                                                      'sigma': [self.data[et]['c2']['sigma'], 0]},
                                                     index=['Average', 'SD'])}
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
            df3 = pd.DataFrame({'ie': [float(self.data[et]['ie']['iedata'].mean()),
                                       float(self.data[et]['ie']['iedata'].std())],
                                'sigma': [self.data[et]['ie']['sigma'], 0]}, index=['Average', 'SD'])
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
        if energy_summary and entropy_summary:
            for et, ev in energy_summary.items():
                binding[et] = {}
                b_map[et] = {}
                for em, emv in ev.items():
                    binding[et][em] = {}
                    b_map[et][em] = []
                    # FIXME: check for stability
                    if not self.app_namespace.FILES.stability:
                        print(True)
                    edata = emv[('delta', 'TOTAL')]
                    edata.name = 'ΔH'
                    if et in entropy_summary:
                        for ent, etv in entropy_summary[et].items():
                            b_map[et][em].append(ent)
                            if ent in ['nmode', 'qh']:
                                entdata = etv[('delta', 'TOTAL')]
                            elif ent == 'ie':
                                entdata = etv['ie']
                            else:
                                entdata = etv['c2']

                            entdata.name = '-TΔS'
                            dgdata = pd.Series([edata[0] + entdata[0], utils.get_std(edata[1], entdata[1])],
                                               index=['Average', 'SD'], name='ΔG')
                            binding[et][em][ent] = pd.concat([edata, entdata, dgdata], axis=1)
        return {'map': b_map, 'data': binding}

    def get_decomp_energy(self, etype: tuple = None, model: tuple=None, mol: tuple = None, contribution: tuple = None,
                          res1: tuple = None, res2: tuple = None, term: tuple = None, res_threshold=0.5,
                          startframe=None, endframe=None, interval=None):

        s, e, index = self._get_frames_index('energy', startframe, endframe, interval)
        name = index.name
        index = pd.concat([index, pd.Series(['Average', 'SD'])])
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
                                                            temp_energy[t] = temp_energy[t].append([mean, std])
                                                            if t == 'tot' and mean < res_threshold:
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
                                                                    temp_energy_r2[t] = temp_energy_r2[t].append([mean, std])
                                                                    if t == 'tot' and r2 != r1:
                                                                        res1_contrib += mean
                                                            if temp_energy_r2:
                                                                temp_energy[r2] = temp_energy_r2
                                                    if float(abs(res1_contrib)) < res_threshold:
                                                        remove = True
                                                if not remove:
                                                    model_decomp_energy[m1][c][r1] = temp_energy
                                                    e_map[etkey][m][m1][c][r1] = temp_emap
                        tdf = pd.DataFrame(flatten(model_decomp_energy), index=index)
                        decomp_energy[etkey][m] = tdf.reindex(sorted(tdf.columns), axis=1)
        return {'map': e_map, 'data': decomp_energy}

    def get_ana_data(self, energy_options=None, entropy_options=None, decomp_options=None, basic_options=None,
                     performance_options=None):
        TASKs = []
        d = {}
        if energy_options is None:
            energy_options = {}
        if entropy_options is None:
            entropy_options = {}
        if decomp_options is None:
            decomp_options = {}

        energy_map, energy, energy_summary = self.get_energy(**energy_options).values()
        if energy_map:
            d['enthalpy'] = {'map': energy_map, 'keys': {}, 'summary': energy_summary}
            for level, value in energy_map.items():
                for level1, value1 in value.items():
                    for level2, value2 in value1.items():
                        TASKs.append([_setup_data, dict(data=energy[level][level1][level2][value2], level=1),
                                      (level, level1, level2), 'enthalpy'])
                        TASKs.extend([_setup_data, dict(data=energy[level][level1][(level2, level3)], level=0),
                                      (level, level1, level2,level3), 'enthalpy'] for level3 in value2)
        entropy_map, entropy, entropy_summary = self.get_entropy(**entropy_options).values()
        if entropy_map:
            d['entropy'] = {'map': entropy_map, 'keys': {}, 'summary': entropy_summary}
            for level, value in entropy_map.items():
                for level1, value1 in value.items():
                    if level1 in ['nmode', 'qh']:
                        for level2, value2 in value1.items():
                            TASKs.append([_setup_data, dict(data=entropy[level][level1][level2], level=1),
                                          (level, level1, level2), 'entropy'])
                            TASKs.extend([_setup_data, dict(data=entropy[level][level1][(level2, level3)], level=0),
                                          (level, level1, level2, level3), 'entropy'] for level3 in value2)
                    elif level1 == 'c2':
                        TASKs.append([_setup_data, dict(data=entropy[level][level1], level=1, iec2=True),
                                      (level, level1), 'entropy'])
                    elif level1 == 'ie':
                        TASKs.append([_setup_data, dict(data=entropy[level][level1], level=1.5),
                                      (level, level1), 'entropy'])

        bind_map, binding = self.get_binding(energy_summary, entropy_summary).values()
        if bind_map:
            d['binding'] = {'map': bind_map, 'keys': {}, 'summary': binding}
            for level, value in bind_map.items():
                for level1, value1 in value.items():
                    TASKs.extend([_setup_data, dict(data=binding[level][level1][level2], level=1), (level, level1, level2), 'binding'] for level2 in value1)

        decomp_map, decomp = self.get_decomp_energy(**decomp_options).values()
        if decomp_map:
            d['decomposition'] = {'map': decomp_map, 'keys': {}}
            for level, value in decomp_map.items():
                for level1, value1 in value.items():
                    for level2, value2 in value1.items():
                        for level3, value3 in value2.items():
                            item_lvl = 2 if self.app_namespace.INPUT['idecomp'] in [1, 2] else 3
                            index = [list(value3.keys())]
                            for level4, value4 in value3.items():
                                if item_lvl == 3 and len(index) == 1:
                                    index.append(list(value4.keys()))
                                item_lvl2 = 1 if self.app_namespace.INPUT['idecomp'] in [1, 2] else 2
                                TASKs.append(
                                    [_setup_data, dict(data=decomp[level][level1][level2][level3][level4],
                                                       level=item_lvl2, name=level4, index=list(value4.keys())),
                                     (level, level1, level2, level3, level4), 'decomposition'])
                                if self.app_namespace.INPUT['idecomp'] in [1, 2]:
                                    TASKs.extend([_setup_data,
                                                  dict(data=decomp[level][level1][level2][level3][level4][level5]),
                                          (level, level1, level2, level3, level4, level5), 'decomposition'] for level5 in
                                                 value4)
                                else:
                                    TASKs.extend([_setup_data, dict(data=decomp[level][level1][level2][level3][level4][level5], level=1, index=value5), (level, level1, level2, level3, level4, level5), 'decomposition'] for level5, value5 in value4.items())
                            TASKs.append([_setup_data, dict(data=decomp[level][level1][(level2, level3)], level=item_lvl,
                                                            name=level3, index=index),
                                      (level, level1, level2, level3), 'decomposition'])

        with ThreadPool() as pool:
            imap_unordered_it = pool.imap_unordered(calculatestar, TASKs)
            for t, id, result in imap_unordered_it:
                d[t]['keys'][id] = result

        return d

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
                         # 'decomp_mutant-normal': app.calc_types.decomp_mut_norm
                         }
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

        ts = 1 if not self.timestep else self.timestep

        INPUT = self.app_namespace.INPUT
        numframes = self.app_namespace.INFO['numframes']
        nmnumframes = self.app_namespace.INFO['numframes_nmode']

        frames_list = list(range(INPUT['startframe'],
                                 INPUT['startframe'] + numframes * INPUT['interval'],
                                 INPUT['interval']))
        INPUT['endframe'] = frames_list[-1]
        time_step_list = list(range(self.starttime,
                                    self.starttime + len(frames_list) * ts * INPUT['interval'],
                                    ts * INPUT['interval']))
        self.frames = dict(zip(frames_list, time_step_list))

        if INPUT['nmoderun']:
            nmframes_list = list(range(INPUT['nmstartframe'],
                                       INPUT['nmstartframe'] + nmnumframes * INPUT['nminterval'],
                                       INPUT['interval']))
            INPUT['nmendframe'] = nmframes_list[-1] if nmframes_list else None

            nm_start = (nmframes_list[0] - frames_list[0]) * INPUT['interval']
            nmtime_step_list = list(range(self.starttime + nm_start,
                                          self.starttime + nm_start + len(nmframes_list) * ts * INPUT['nminterval'],
                                          ts * INPUT['nminterval']))

            self.nmframes = dict(zip(nmframes_list, nmtime_step_list))

    def _update_eframes(self, startframe, endframe, interval):

        # get the original frames list
        self._get_frames()

        curr_frames = list(range(startframe, endframe + interval, interval))
        frames = list(self.frames.keys())
        for x in frames:
            if x not in curr_frames:
                self.frames.pop(x)

    def _update_nmframes(self, startframe, endframe, interval):

        # get the original frames list
        self._get_frames()

        curr_frames = list(range(startframe, endframe + interval, interval))
        frames = list(self.nmframes.keys())
        for x in frames:
            if x not in curr_frames:
                self.nmframes.pop(x)

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
