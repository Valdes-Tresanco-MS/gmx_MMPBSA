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

    def _define_line_chart_btn(self):
        line_menu = QMenu()
        line_menu.setTitle('Line charts menu')
        self.line_table_action = line_menu.addAction('Show in table')
        self.line_table_action.setCheckable(True)
        self.line_table_action.toggled.connect(self._show_line_table)

        self.line_chart_action = QToolButton()
        self.line_chart_action.setIcon(QIcon(
            '/home/mario/PycharmProjects/gmx_MMPBSA/GMXMMPBSA/analyzer/style/line-chart.svg'))
        self.line_chart_action.setText('Line Chart')
        self.line_chart_action.setCheckable(True)
        self.line_chart_action.setPopupMode(QToolButton.MenuButtonPopup)
        self.line_chart_action.setContentsMargins(0, 0, 0, 0)
        self.line_chart_action.setMenu(line_menu)
        self.btn_group.addButton(self.line_chart_action, 1)

        self.tb.addWidget(self.line_chart_action)


    def _define_bar_chart_btn(self):
        bar_menu = QMenu()
        bar_menu.setTitle('Bar charts menu')
        self.bar_table_action = bar_menu.addAction('Show in table')
        self.bar_table_action.setCheckable(True)
        self.bar_table_action.toggled.connect(self._show_bar_table)

        self.bar_chart_action = QToolButton()
        self.bar_chart_action.setIcon(
            QIcon('/home/mario/PycharmProjects/gmx_MMPBSA/GMXMMPBSA/analyzer/style/bar-chart.png'))
        self.bar_chart_action.setText('Bar Chart')
        self.bar_chart_action.setCheckable(True)
        self.bar_chart_action.setPopupMode(QToolButton.MenuButtonPopup)
        self.bar_chart_action.setContentsMargins(0, 0, 0, 0)
        self.bar_chart_action.setMenu(bar_menu)
        self.btn_group.addButton(self.bar_chart_action, 2)

        self.tb.addWidget(self.bar_chart_action)


    def _define_heatmap_chart_btn(self):
        heatmap_menu = QMenu()
        heatmap_menu.setTitle('heatmap charts menu')
        self.heatmap_table_action = heatmap_menu.addAction('Show in table')
        self.heatmap_table_action.setCheckable(True)
        self.heatmap_table_action.toggled.connect(self._show_heatmap_table)

        self.heatmap_chart_action = QToolButton()
        self.heatmap_chart_action.setIcon(
            QIcon('/home/mario/PycharmProjects/gmx_MMPBSA/GMXMMPBSA/analyzer/style/heatmap_icon.svg'))
        self.heatmap_chart_action.setText('Heatmap Chart')
        self.heatmap_chart_action.setCheckable(True)
        self.heatmap_chart_action.setPopupMode(QToolButton.MenuButtonPopup)
        self.heatmap_chart_action.setContentsMargins(0, 0, 0, 0)
        self.heatmap_chart_action.setMenu(heatmap_menu)
        self.btn_group.addButton(self.heatmap_chart_action, 3)

        self.tb.addWidget(self.heatmap_chart_action)

    def _define_vis_btn(self):
        heatmap_menu = QMenu()
        heatmap_menu.setTitle('Visualization menu')
        # self.heatmap_table_action = heatmap_menu.addAction('Show in text')
        # self.heatmap_table_action.setCheckable(True)
        # self.heatmap_table_action.toggled.connect(self._show_heatmap_table)

        self.vis_action = QToolButton()
        self.vis_action.setIcon(QIcon('/home/mario/PycharmProjects/gmx_MMPBSA/GMXMMPBSA/analyzer/style/molecule.svg'))
        self.vis_action.setText('Visualization')
        self.vis_action.setCheckable(True)
        self.vis_action.setPopupMode(QToolButton.MenuButtonPopup)
        self.vis_action.setContentsMargins(0, 0, 0, 0)
        # self.vis_action.setMenu(heatmap_menu)
        self.btn_group.addButton(self.vis_action, 4)

        self.tb.addWidget(self.vis_action)

    def _define_option_button(self):
        options_menu = QMenu()
        options_menu.setTitle('System Options')
        self.output_action = options_menu.addAction('Show Output file')
        self.output_action.setCheckable(True)
        self.output_action.toggled.connect(self._show_output_file)
        if self.app.systems[self.system_index]['namespace'].INFO['decomp_output_file']:
            self.decomp_output_action = options_menu.addAction('Show Decomp Output file')
            self.decomp_output_action.setCheckable(True)
            self.decomp_output_action.toggled.connect(self._show_decomp_output_file)

        self.options_button = QToolButton()
        self.options_button.setIcon(QIcon('/home/mario/PycharmProjects/PyQtRibbon/error_checker.png'))
        self.options_button.setText('Options')
        self.options_button.setPopupMode(QToolButton.InstantPopup)
        self.options_button.setContentsMargins(0, 0, 0, 0)
        self.options_button.setMenu(options_menu)

    def fn_btn_group(self, btn, checked):
        all_checked = all(x.isChecked() for x in self.btn_group.buttons())
        all_unchecked = not any(x.isChecked() for x in self.btn_group.buttons())
        if all_checked:
            self.mark_all.setCheckState(Qt.Checked)
        elif all_unchecked:
            self.mark_all.setCheckState(Qt.Unchecked)
        else:
            self.mark_all.setCheckState(Qt.PartiallyChecked)
        if self.btn_group.id(btn) == 1:
            self.plotting_line(checked)
        elif self.btn_group.id(btn) == 2:
            self.plotting_bar(checked)
        elif self.btn_group.id(btn) == 3:
            self.plotting_heatmap(checked)
        elif self.btn_group.id(btn) == 4:
            self.visualizing(checked)

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
