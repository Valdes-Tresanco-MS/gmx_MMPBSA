# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
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
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pandas as pd
from .utils import com2str, energy2pdb_pml


class SpacerItem(QToolButton):
    def __init__(self, parent=None):
        super(SpacerItem, self).__init__(parent)
        self.setDisabled(True)
        self.setContentsMargins(0, 0, 0, 0)
        self.setStyleSheet("QToolButton { /* mimic the look of the QToolButton with MenuButtonPopup */ "
                           "padding-right: 15px; /* make way for the popup button */}")


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
    def __init__(self, parent, stringlist, data=None, app=None, level=0, chart_title='Binding Free Energy',
                 chart_subtitle='', system_index=1, iec2_data=None, item_type='energy',
                 buttons=(), options_btn=False, remove_empty_terms=False):
        super(CustomItem, self).__init__(parent, stringlist)

        self.remove_empty_terms = remove_empty_terms
        self.data = data
        self.system_index = system_index
        self.app = app
        self.level = level
        self.chart_title = chart_title
        self.chart_subtitle = chart_subtitle
        self.item_name = stringlist[0]
        self.buttons = buttons
        self.options_btn = options_btn
        self.iec2_data = iec2_data
        self.item_type = item_type
        self.properties = {}

        self.lp_subw = None
        self.bp_subw = None
        self.hmp_subw = None
        self.pymol_process = None
        self.pymol_data_change = False
        self.bfactor_pml = None
        self.output_file_subw = None
        self.decomp_output_file_subw = None
        self.line_table_subw = None
        self.bar_table_subw = None
        self.heatmap_table_subw = None
        self.ie_plot_data = None

        # changes
        self.frange = []
        self.line_change = False
        self.bar_change = False
        self.heatmap_change = False


        self.changed = False

        self.tb = QToolBar()
        self.tb.setStyleSheet("QToolBar {padding: 0, 20, 0, 20;}")

        self.tb.setIconSize(QSize(16, 16))
        self.tb.setToolButtonStyle(Qt.ToolButtonIconOnly)
        self.item_charts = []

        self.btn_group = QButtonGroup()
        self.btn_group.setExclusive(False)
        self.btn_group.buttonToggled.connect(self.fn_btn_group)

        self.mark_all = QCheckBox()
        self.mark_all.setToolTip('Mark all actions at the same time')
        self.mark_all.setTristate(True)
        self.mark_all.setCheckable(True)
        self.mark_all.stateChanged.connect(self.fn_mark_all)

        self.charts_action = {
            1: self._define_line_chart_btn,
            2: self._define_bar_chart_btn,
            3: self._define_heatmap_chart_btn,
            4: self._define_vis_btn,
        }

    def changes(self, line, bar, heatmap):
        self.line_change = line
        self.bar_change = bar
        self.heatmap_change = heatmap

            bar_plot_data = pd.DataFrame(data=bar)
            line_plot_data = pd.DataFrame(data={'frames': self.frames[start:end:interval],
                                                'Energy': np.array(data['Energy']).sum(axis=0)})

            heatmap_plot_data = pd.DataFrame(data=data['Energy'], index=data['Residues'], columns=data['frames'])

            return_data = Namespace(line_plot_dat=line_plot_data, bar_plot_dat=bar_plot_data,
                                    heatmap_plot_dat=heatmap_plot_data)
        elif self.level == 2.1:
            bar = {}
            data = {'frames': self.nmode_frames, 'Component': [], 'Energy': []}
            for p, d in self.cdata.items():
                data['Component'].append(p)
                bar[p] = d[start:end:interval]
                data['Energy'].append(d[start:end:interval])

            bar_plot_data = pd.DataFrame(data=bar)
            line_plot_data = pd.DataFrame(data={'frames': self.frames[start:end:interval],
                                                'Energy': np.array(data['Energy'])})
            return_data.line_plot_dat = line_plot_data
            return_data.bar_plot_dat=bar_plot_data

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

    def setup_data(self, frange, iec2frames=0):
        if self.data is None:
            return
        if self.frange == frange:
            return
        self.frange = frange
        if self.level == 0:
            self.line_plot_data = self.data.loc[frange[0]:frange[1]:frange[2]]
            if self.item_type == 'ie':
                tempserie = self.line_plot_data[-iec2frames:]
                self.ie_plot_data = pd.concat([tempserie, pd.Series([self.iec2_data['sigma']] * iec2frames,
                                                                    index=tempserie.index, name='sigma')], axis=1)
        elif self.level == 1:
            self.bar_plot_data = self.data.loc[frange[0]:frange[1]:frange[2]]
            # IMPORTANT: can be used in nmode
        elif self.level == 2:
            tempdf = self.data.loc[frange[0]:frange[1]:frange[2], self.data.columns.get_level_values(1) == 'tot']
            self.bar_plot_data = tempdf.droplevel(level=1, axis=1)
            self.line_plot_data = self.bar_plot_data.sum(axis=1)
            self.heatmap_plot_data = self.bar_plot_data.transpose(copy=True)
            del tempdf
        elif self.level == 3:
            # Select only the "tot" column, remove the level, change first level of comlumn to rows and remove the mean
            # index
            tempdf = self.data.loc[frange[0]:frange[1]:frange[2], self.data.columns.get_level_values(2) == 'tot']
            self.heatmap_plot_data = tempdf.aggregate(["mean"]).droplevel(level=2, axis=1).stack().droplevel(level=0)
            self.bar_plot_data = tempdf.sum(axis=1, level=0)
            self.line_plot_data = self.bar_plot_data.sum(axis=1)
            del tempdf

            return_data = Namespace(line_plot_dat=line_plot_data, bar_plot_dat=bar_plot_data,
                                    heatmap_plot_dat=heatmap_plot_data)
        return return_data
