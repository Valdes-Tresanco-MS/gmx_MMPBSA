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
from PyQt5.QtWidgets import QTreeWidgetItem
from PyQt5.QtCore import Qt
import copy
from types import SimpleNamespace as Namespace
import numpy as np
import pandas as pd

class CorrelationItem(QTreeWidgetItem):
    def __init__(self, parent, stringlist, model=None, enthalpy=None, dgie=None, dgnmode=None, dgqh=None, col_box=None):
        super(CorrelationItem, self).__init__(parent, stringlist)

        self.model = model
        self.enthalpy = enthalpy
        self.dgie = dgie
        self.dgnmode = dgnmode
        self.dgqh = dgqh
        self.chart_title = f'Correlation Using {stringlist[0].upper()} model'
        self.chart_subtitle = ['Exp. Energy vs Enthalpy (ΔH)', 'Exp. Energy vs Pred. Energy (ΔH+IE)',
                               'Exp. Energy vs Pred. Energy (ΔH+NMODE)', 'Exp. Energy vs Pred. Energy (ΔH+QH)']
        self.item_name = stringlist[0]
        if col_box:
            for col in col_box:
                self.setCheckState(col, Qt.Unchecked)

        self.dh_sw = None
        self.dgie_sw = None
        self.dgnmode_sw = None
        self.dgqh_sw = None

    def getplotdata(self):
        # correlation_data[sys_name] = {'ΔG': {
        #     'gb': {'ie': 0, 'qh': 0, 'nmode': 0},
        #     'pb': {'ie': 0, 'qh': 0, 'nmode': 0},
        #     'rism std': {'ie': 0, 'qh': 0, 'nmode': 0},
        #     'rism gf': {'ie': 0, 'qh': 0, 'nmode': 0}},
        #     'ΔH': {'gb': 0, 'pb': 0, 'rism std': 0, 'rism gf': 0},
        #     'Exp.Energy': ki2energy(topItem.exp_ki, topItem.app.INPUT['temperature'])}
        for system in self.cdata:

            pass

class CustomItem(QTreeWidgetItem):
    def __init__(self, parent, stringlist, system=None, app=None, has_chart=True, cdata=None, level=0,
                 chart_title='Binding Free Energy', chart_subtitle='', col_box=(), remove_empty_terms=False):
        super(CustomItem, self).__init__(parent, stringlist)

        if isinstance(parent, CustomItem):
            self.syspath = parent.syspath
            self.exp_ki = parent.exp_ki
            self.sysname = parent.sysname
            self.temp = parent.temp
            self.app = parent.app
            self.frames = parent.frames
            self.start = parent.start
            self.end = parent.end
            self.interval = parent.interval
            self.idecomp = parent.idecomp
            self.remove_empty_terms = parent.remove_empty_terms
        else:
            self.syspath = system[1]
            self.exp_ki = system[2]
            self.sysname = system[0]
            self.temp = system[3]
            self.app = app
            self.interval = app.INPUT['interval']
            self.start = app.INPUT['startframe']
            self.frames = np.array([x for x in range(self.start, self.start +
                                                     self.app.numframes * self.interval,
                                                     self.interval)])
            self.end = app.INPUT['endframe']
            self.idecomp = app.INPUT['idecomp']
            self.remove_empty_terms = remove_empty_terms
        self.cdata = cdata
        self.level = level
        self.chart_title = chart_title
        self.chart_subtitle = chart_subtitle
        self.has_chart = has_chart
        self.item_name = stringlist[0]
        self.col_box = col_box

        self.lp_subw = None
        self.bp_subw = None
        self.hmp_subw = None
        self.pymol_process = None

        if col_box:
            for col in col_box:
                self.setCheckState(col, Qt.Unchecked)

    def get_data(self):
        self.gmxMMPBSA_data = self.getplotdata()
        self.gmxMMPBSA_current_data = copy.deepcopy(self.gmxMMPBSA_data)

    def reset(self):
        self.gmxMMPBSA_current_data = copy.deepcopy(self.gmxMMPBSA_data)

    def update_data(self, start, end, interval):
        self.gmxMMPBSA_current_data = self.getplotdata(start, end, interval)
        # print(self.gmxMMPBSA_current_data)

    def getplotdata(self, start=0, end=None, interval=1):
        """
        Method to get data at different levels and lengths
        :param start: first frame
        :param end: last frame
        """

        return_data = Namespace(line_plot_dat=None, bar_plot_dat=None, heatmap_plot_dat=None)

        if self.level == 0:
            return_data.line_plot_dat = pd.DataFrame(data={'frames': self.frames[start:end:interval],
                                                      'Energy': self.cdata[start:end:interval]})
        elif self.level == 1:
            dat = {}
            for p, d in self.cdata.items():
                # FIXME: ignore empty or 0 array ???
                if type(d) not in [list, np.ndarray] :
                    dat[p] = [d]
                else:
                    if self.remove_empty_terms:
                        if abs(d[start:end:interval].mean()) > 0.1:
                            dat[p] = d[start:end:interval]
                    else:
                        dat[p] = d
            return_data.bar_plot_dat = pd.DataFrame(data=dat)

        elif self.level == 2:
            bar = {}
            data = {'frames': self.frames[start:end:interval], 'Residues': [], 'Energy': []}
            for p, d in self.cdata.items():
                data['Residues'].append(p)
                for p1, d1 in d.items():
                    if 'tot' in str(p1).lower():
                        if self.remove_empty_terms:
                            if abs(d1[start:end:interval].mean()) > 0.1:
                                bar[p] = d1[start:end:interval]
                        else:
                            bar[p] = d1[start:end:interval]
                        data['Energy'].append(d1[start:end:interval])

            bar_plot_data = pd.DataFrame(data=bar)
            line_plot_data = pd.DataFrame(data={'frames': self.frames[start:end:interval],
                                                'Energy': np.array(data['Energy']).sum(axis=0)})

            heatmap_plot_data = pd.DataFrame(data=data['Energy'], index=data['Residues'], columns=data['frames'])

            return_data = Namespace(line_plot_dat=line_plot_data, bar_plot_dat=bar_plot_data,
                                    heatmap_plot_dat=heatmap_plot_data)
        elif self.level == 3:
            bar = {}
            data = {'frames': self.frames[start:end:interval], 'Residues': [], 'Energy': [], 'Per-pair Energy': []}
            per_pair_data = {'Residues': [], 'Pair': [], 'Energy': []}
            for p, d in self.cdata.items():
                data['Residues'].append(p)
                res_t = []
                for p1, d1 in d.items():
                    per_pair_data['Residues'].append(p)
                    per_pair_data['Pair'].append(p1)
                    for p2, d2 in d1.items():
                        if 'tot' in str(p2).lower():
                            per_pair_data['Energy'].append(d2[start:end:interval].mean())
                            res_t.append(d2[start:end:interval])
                if self.remove_empty_terms:
                    if abs(np.sum(res_t, axis=0).mean()) > 0.1:
                        bar[p] = np.sum(res_t, axis=0)
                else:
                    bar[p] = np.sum(res_t, axis=0)
                res_t_np = np.array(res_t)
                data['Energy'].append(res_t_np.sum(axis=0))
                data['Per-pair Energy'].append(res_t_np.mean(axis=1))

            bar_plot_data = pd.DataFrame(data=bar)
            line_plot_data = pd.DataFrame(data={'frames': self.frames[start:end:interval],
                                                'Energy': np.array(data['Energy']).sum(axis=0)})

            heatmap_plot_data = pd.DataFrame(data=data['Per-pair Energy'], index=data['Residues'], columns=data['Residues'])

            return_data = Namespace(line_plot_dat=line_plot_data, bar_plot_dat=bar_plot_data,
                                    heatmap_plot_dat=heatmap_plot_data)
        return return_data
