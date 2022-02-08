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

import sys
import os

from queue import Queue, Empty
from pathlib import Path

import pandas
import pandas as pd
from GMXMMPBSA.API import load_gmxmmpbsa_info, MMPBSA_API
from PyQt6.QtWidgets import *
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from GMXMMPBSA.analyzer.dialogs import InitDialog
from GMXMMPBSA.analyzer.customitem import CustomItem, CorrelationItem
from GMXMMPBSA.analyzer.utils import energy2pdb_pml, ki2energy, make_corr_DF, multiindex2dict
from GMXMMPBSA.analyzer.chartsettings import ChartSettings
import math
import numpy as np


def run(infofile):
    info = Path(infofile)
    app = QApplication(sys.argv)
    app.setApplicationName('GMX-MMPBSA Analyzer Tool')
    w = GMX_MMPBSA_ANA()
    w.gettting_data([info])
    w.show()
    sys.exit(app.exec())


class GMX_MMPBSA_ANA(QMainWindow):
    def __init__(self):
        super(GMX_MMPBSA_ANA, self).__init__()
        self.corr_data = {'mutant': {}}

        self.systems = {}
        self.current_system_index = None
        self.pymol_p_list = []

        self.items_counter = {'charts': 0, 'pymol': [], 'bars': 0, 'line': 0, 'heatmap': 0}

        self.mdi = QMdiArea(self)
        self.mdi.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setCentralWidget(self.mdi)

        self.menubar = self.menuBar()
        # self.opendirAct = QAction("Open gmx_MMPBSA Dir...", self)
        # self.opendirAct.triggered.connect(self.getInfoFile)
        self.fileMenu = self.menuBar().addMenu("&File")
        # self.fileMenu.addAction(self.opendirAct)
        self.fileMenu.addAction('Close All..', self.mdi.closeAllSubWindows)
        self.viewMenu = self.menuBar().addMenu("&View")
        self.viewMenu.addAction('Tile SubWindows', self.mdi.tileSubWindows)
        self.viewMenu.addAction('Cascade SubWindows', self.mdi.cascadeSubWindows)
        self.statusbar = self.statusBar()

        self.treeDockWidget = QDockWidget('Data', self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.treeDockWidget)

        self.correlation_DockWidget = QDockWidget('Correlations', self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.correlation_DockWidget)

        self.treeWidget = QTreeWidget(self)
        self.treeWidget.setMinimumWidth(380)
        self.treeWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        # self.treeWidget.customContextMenuRequested.connect(self.data_context_menu)
        self.treeWidget.itemSelectionChanged.connect(self.update_system_selection)

        self._make_options_panel()

        self.data_container = QSplitter(Qt.Orientation.Vertical)
        self.data_container.addWidget(self.treeWidget)
        self.data_container.addWidget(self.optionWidget)
        self.data_container.setStretchFactor(0, 10)
        # self.data_container.setStretchFactor(1, 1)

        self.treeDockWidget.setWidget(self.data_container)
        sys_label = QTreeWidgetItem(['System', 'Charts'])
        sys_label.setToolTip(0, 'System')
        sys_label.setToolTip(1, 'Charts')

        self.treeWidget.setHeaderItem(sys_label)
        self.treeWidget.setColumnHidden(4, True)
        header = self.treeWidget.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

        self.correlation_treeWidget = QTreeWidget(self)
        self.correlation_treeWidget.itemClicked.connect(self.update_table)
        self.correlation_treeWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        # self.correlation_treeWidget.customContextMenuRequested.connect(self.corr_context_menu)
        model_label = QTreeWidgetItem(['MODEL', 'ΔH', 'ΔH+IE', 'ΔH+NMODE', 'ΔH+QH'])
        model_label.setToolTip(0, 'Selected Model')
        model_label.setToolTip(1, 'Correlation plot for ΔG = ΔH+IE')
        model_label.setToolTip(2, 'Correlation plot for ΔG = ΔH+NMODE')
        model_label.setToolTip(3, 'Correlation plot for ΔG = ΔH+QH')
        model_label.setToolTip(4, 'Correlation plot for ΔH only')
        self.correlation_treeWidget.setHeaderItem(model_label)
        cheader = self.correlation_treeWidget.header()
        cheader.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)

        self.data_table_widget = QTableWidget(0, 3)
        self.data_table_widget.setHorizontalHeaderLabels(['System', 'Exp.ΔG'])
        self.data_table_energy_col = QTableWidgetItem('None')
        self.data_table_widget.setHorizontalHeaderItem(2, self.data_table_energy_col)
        self.data_table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.data_table_widget.horizontalHeader().setStretchLastSection(True)

        self.corr_container_widget = QSplitter(Qt.Orientation.Vertical)
        self.corr_container_widget.addWidget(self.correlation_treeWidget)
        self.corr_container_widget.addWidget(self.data_table_widget)
        self.correlation_DockWidget.setWidget(self.corr_container_widget)

        self.init_dialog = InitDialog(self)

    def _initialize_systems(self):

        # Checks all systems has the same number of frames
        frames = [self.systems[s]['namespace'].INFO['numframes'] for s in self.systems]
        frames_same = frames.count(frames[0]) == len(frames)
        start = [self.systems[s]['namespace'].INPUT['startframe'] for s in self.systems]
        start_same = start.count(start[0]) == len(start)
        interval = [self.systems[s]['namespace'].INPUT['interval'] for s in self.systems]
        interval_same = interval.count(interval[0]) == len(interval)

        # Fixme: hacer lo mismo para comprobar si los sistemas son homogeneos
        nmode_frames = [self.systems[s]['namespace'].INFO['numframes_nmode'] for s in self.systems]
        nmode_frames_same = nmode_frames.count(nmode_frames[0]) == len(nmode_frames)

        if all([frames_same, start_same, interval_same]):
            self.all_frb.setEnabled(True)

        # Select automatically the first item to update the option panel
        topitem = self.treeWidget.topLevelItem(0)
        topitem.setSelected(True)

        f_start = self.systems[1]['namespace'].INPUT['startframe']
        f_interval = self.systems[1]['namespace'].INPUT['interval']
        f_end = self.systems[1]['namespace'].INPUT['endframe']

        self.eframes_start_sb.setRange(f_start, f_end)
        self.eframes_start_sb.setSingleStep(f_interval)
        self.eframes_start_sb.setValue(f_start)

        self.eframes_inter_sb.setRange(1, f_end - f_start)

        self.eframes_end_sb.setRange(f_start, f_end)
        self.eframes_end_sb.setSingleStep(f_interval)
        self.eframes_end_sb.setValue(f_end)

    def _make_options_panel(self):

        self.optionWidget = QWidget(self)
        optionWidget_l = QVBoxLayout(self.optionWidget)
        optionWidget_l.setContentsMargins(0, 0, 0, 0)
        # optionWidget_l.addLayout(self.btn_l)
        hl = QFrame()
        hl.setFrameShape(QFrame.Shape.HLine)
        hl.setFrameShadow(QFrame.Shadow.Sunken)
        optionWidget_l.addWidget(hl)

        update_btn = QPushButton('Update')
        update_btn.clicked.connect(self.update_fn)
        resetchart_btn = QPushButton('Reset')
        btn_l = QHBoxLayout()
        btn_l.addWidget(QLabel('Options'))
        btn_l.addWidget(resetchart_btn)
        btn_l.addWidget(update_btn)
        optionWidget_l.addLayout(btn_l)

        selection_group = QGroupBox('Selection')
        selection_group_layout = QHBoxLayout(selection_group)
        selection_group_layout.setContentsMargins(10, 0, 10, 0)
        # self.curr_frb = QRadioButton('Current chart')
        # self.curr_frb.setEnabled(False)
        # selection_group_layout.addWidget(self.curr_frb)
        self.curr_sys_frb = QRadioButton('Selected System')
        self.curr_sys_frb.setChecked(True)
        selection_group_layout.addWidget(self.curr_sys_frb)
        self.all_frb = QRadioButton('All Systems')
        self.all_frb.setEnabled(False)
        selection_group_layout.addWidget(self.all_frb)

        optionWidget_l.addWidget(selection_group)

        optionWidget_c = QTabWidget(self)
        optionWidget_l.addWidget(optionWidget_c)

        frames_w = QTabWidget(optionWidget_c)
        frames_w.setTabPosition(QTabWidget.TabPosition.South)
        optionWidget_c.addTab(frames_w, 'Frames')

        frames_group = QWidget()

        self.eframes_start_sb = QSpinBox()
        # self.eframes_start_sb.setRange(1, 10000)
        self.eframes_start_sb.setAccelerated(True)
        self.eframes_start_sb.valueChanged.connect(self.frames_start_sb_update)
        self.eframes_inter_sb = QSpinBox()
        self.eframes_inter_sb.valueChanged.connect(self.frames_inter_sb_update)
        self.eframes_end_sb = QSpinBox()
        # self.eframes_end_sb.setRange(1, 10000)
        self.eframes_end_sb.setAccelerated(True)
        self.eframes_end_sb.valueChanged.connect(self.frames_end_sb_update)
        self.numframes_le = QLineEdit()
        self.numframes_le.setReadOnly(True)

        fg_l1 = QFormLayout()
        # fg_l1.setContentsMargins(5, 5, 5, 5)
        # fg_l1.setSpacing(3)
        fg_l1.addRow('Start', self.eframes_start_sb)
        fg_l1.addRow('End', self.eframes_end_sb)
        fg_l2 = QFormLayout()
        # fg_l2.setContentsMargins(5, 5, 5, 5)
        # fg_l2.setSpacing(3)
        fg_l2.addRow('Interval', self.eframes_inter_sb)
        fg_l2.addRow('Nr. frames', self.numframes_le)

        frames_group_l = QHBoxLayout(frames_group)
        frames_group_l.addLayout(fg_l1, 1)
        frames_group_l.addLayout(fg_l2, 1)

        frames_w.addTab(frames_group, 'Energy')
        frames_w.setTabToolTip(0, 'Frames defined for all energy calculations')

        nmframes_group = QWidget()
        self.nmframes_start_sb = QSpinBox()
        # self.eframes_start_sb.setRange(1, 10000)
        self.nmframes_start_sb.setAccelerated(True)
        self.nmframes_start_sb.valueChanged.connect(self.frames_start_sb_update)
        self.nmframes_inter_sb = QSpinBox()
        self.nmframes_inter_sb.valueChanged.connect(self.frames_inter_sb_update)
        self.nmframes_end_sb = QSpinBox()
        # self.eframes_end_sb.setRange(1, 10000)
        self.nmframes_end_sb.setAccelerated(True)
        self.nmframes_end_sb.valueChanged.connect(self.frames_end_sb_update)
        self.nmnumframes_le = QLineEdit()
        self.nmnumframes_le.setReadOnly(True)

        nmfg_l1 = QFormLayout()
        # nmfg_l1.setContentsMargins(5, 5, 5, 5)
        # nmfg_l1.setSpacing(3)
        nmfg_l1.addRow('Start', self.nmframes_start_sb)
        nmfg_l1.addRow('End', self.nmframes_end_sb)
        nmfg_l2 = QFormLayout()
        # nmfg_l2.setContentsMargins(5, 5, 5, 5)
        # nmfg_l2.setSpacing(3)
        nmfg_l2.addRow('Interval', self.nmframes_inter_sb)
        nmfg_l2.addRow('Nr. frames', self.nmnumframes_le)

        nmframes_group_l = QHBoxLayout(nmframes_group)
        nmframes_group_l.addLayout(nmfg_l1, 1)
        nmframes_group_l.addLayout(nmfg_l2, 1)

        frames_w.addTab(nmframes_group, 'NMODE')
        frames_w.setTabToolTip(2, 'Frames defined for NMODE calculation')

        ieframes_group = QWidget()

        self.iesegment_sb = QSpinBox()
        self.iesegment_sb.setRange(1, 100)
        self.iesegment_sb.setAccelerated(True)
        # self.segment_sb.valueChanged.connect(self.frames_start_sb_update)

        self.ienumframes_le = QLineEdit()
        self.ienumframes_le.setReadOnly(True)

        ie_l = QFormLayout(ieframes_group)
        # ie_l.setContentsMargins(5, 5, 5, 5)
        # ie_l.setSpacing(3)
        ie_l.addRow('Segment', self.iesegment_sb)
        ie_l.addRow('Nr. frames', self.ienumframes_le)

        frames_w.addTab(ieframes_group, 'IE')

        # Charts options
        from GMXMMPBSA.analyzer.parametertree import ParameterTree, Parameter
        # from GMXMMPBSA.analyzer.propertyeditor import p
        self.chart_options_param = Parameter().create()
        self.chart_options_w = ParameterTree()

        # properties_w.setParameters(p, showTop=False)
        # charts_options_w = QWidget()
        optionWidget_c.addTab(self.chart_options_w, 'Charts Options')

    def update_system_selection(self):
        if not self.treeWidget.selectedItems():
            return
        parent_item = self.treeWidget.selectedItems()[0]

        if self.current_system_index == parent_item.system_index:
            return
        while parent_item:
            next_item = parent_item.parent()
            if next_item:
                parent_item = next_item
            else:
                break
        if not self.all_frb.isChecked():
            current_start, current_end, current_interval = self.systems[parent_item.system_index]['current_frames']
            self.eframes_start_sb.setValue(current_start)
            self.eframes_inter_sb.setValue(current_interval)
            self.eframes_end_sb.setValue(current_end)
            self.numframes_le.setText(f"{int((current_end - current_start) // current_interval) + 1}")
            self.chart_options_param.restoreState(self.systems[parent_item.system_index]['chart_options'])
            self.chart_options_w.setParameters(self.chart_options_param, showTop=False)
        self.current_system_index = parent_item.system_index

    def update_fn(self):

        recalc_energy = []
        recalc_nmode = []
        recalc_ie = []
        repaint = []

        act_frange = [self.eframes_start_sb.value(), self.eframes_end_sb.value(), self.eframes_inter_sb.value()]
        act_sett = self.chart_options_param.saveState()

        processed_sys = []

        if self.all_frb.isChecked():
            for x in self.systems:
                processed_sys.append(x)
                if act_frange != self.systems[x]['current_frames']:
                    recalc_energy.append(x)

                self.systems[x]['current_frames'] = act_frange
                if self.systems[x]['chart_options'].is_changed(act_sett):
                    repaint.append(x)
        else:
            processed_sys.append(self.current_system_index)
            if act_frange != self.systems[self.current_system_index]['current_frames']:
                recalc_energy.append(self.current_system_index)

            if self.systems[self.current_system_index]['chart_options'].is_changed(act_sett):
                repaint.append(self.current_system_index)

        maximum = len(recalc_energy) + len(recalc_ie) + len(recalc_nmode) + len(repaint)

        if not maximum:
            return
        maximum += 1 # close and re-open current active windows

        qpd = QProgressDialog('Creating systems tree', 'Abort', 0, maximum, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(0)
        v = 0
        if recalc_energy:
            qpd.setLabelText('Recalculating energies for selected frames range')
            for e in recalc_energy:
                v += 1
                self.systems[e]['api'].update_energy_frames(*act_frange)
                self.systems[e]['data'] = self.systems[e]['api'].get_energy()
                self.systems[e]['current_frames'] = act_frange
                qpd.setValue(v)
                qpd.setValue(v)
        if recalc_nmode:
            qpd.setLabelText('Recalculating nmode for selected frames range')
            for e in recalc_nmode:
                v += 1
                # self.systems[e]['api'].update_energy_frames(*act_frange)
                self.systems[e]['current_frames'] = act_frange
                qpd.setValue(v)
        if recalc_ie:
            qpd.setLabelText('Recalculating Interaction Entropy for selected frames range')
            for e in recalc_ie:
                v += 1
                # self.systems[e]['api'].update_energy_frames(*act_frange)
                self.systems[e]['current_frames'] = act_frange
                qpd.setValue(v)
        if repaint:
            qpd.setLabelText('Updating charts options')
            for e in repaint:
                v += 1
                self.systems[e]['chart_options'].get_changes(act_sett)
                qpd.setValue(v)

        comp = []
        if recalc_energy and recalc_ie and recalc_nmode:
            comp.append('all')
        else:
            if recalc_energy:
                comp.append('energy')
            if recalc_ie:
                comp.append('ie')
            if recalc_nmode:
                comp.append('nmode')

        qpd.setLabelText('Setting data...')
        for s in processed_sys:
            parts = list(self.systems[s]['data'].keys())
            for p in parts:
                self.setting_item_data(s, p, comp=tuple(comp))
            self.systems[s]['items_summary'] = self.systems[s]['api'].get_summary()

        qpd.setLabelText('Updating opened charts')


        subwindows = self.mdi.subWindowList()
        for sub in subwindows:
            if sub.isVisible():
                sub.button.setChecked(False)
                sub.button.setChecked(True)


        # re-assign changes to native state after replot all open charts
        for s in processed_sys:
            self.systems[s]['chart_options'].changes = dict(line_action=0,
                                                            line_ie_action=0,
                                                            bar_action=0,
                                                            heatmap_action=0,
                                                            visualization_ation=0)


        qpd.setValue(maximum)

    def reset_dc(self):
        sub = self.mdi.activeSubWindow()
        sub.item.reset_data()
        sub.make_chart()

    def frames_start_sb_update(self, value):
        if value > self.eframes_end_sb.value():
            self.eframes_end_sb.setValue(value)

        # print(self.current_system_index, '####')
        if self.current_system_index:
            current_start = value
            current_interval = self.eframes_inter_sb.value()
            current_end = self.eframes_end_sb.value()
            # current_start, current_interval, current_end = self.systems[self.current_system_index]['current_frames']
            interval = self.systems[self.current_system_index]['namespace'].INPUT['interval']
            self.numframes_le.setText(f"{int((current_end - value) // current_interval) + 1}")

    def frames_end_sb_update(self, value):
        if value < self.eframes_start_sb.value():
            self.eframes_start_sb.setValue(value)
        if self.current_system_index:
            current_start = self.eframes_start_sb.value()
            current_interval = self.eframes_inter_sb.value()
            current_end = value
            # current_start, current_interval, current_end = self.systems[self.current_system_index]['current_frames']
            interval = self.systems[self.current_system_index]['namespace'].INPUT['interval']
            self.numframes_le.setText(f"{int((current_end - current_start) // current_interval) + 1}")

    def frames_inter_sb_update(self, value):
        # self.eframes_start_sb.setSingleStep(value)
        # self.eframes_end_sb.setSingleStep(value)
        if self.current_system_index:
            current_start = self.eframes_start_sb.value()
            current_interval = value
            current_end = self.eframes_end_sb.value()

            # current_start, current_interval, current_end = self.systems[self.current_system_index]['current_frames']
            interval = self.systems[self.current_system_index]['namespace'].INPUT['interval']
            self.numframes_le.setText(f"{int((current_end - current_start) // (value * interval)) + 1}")

    def gettting_data(self, info_files):

        self.init_dialog.get_files_info(info_files)
        self.init_dialog.show()

    # def showcorr(self, item: CorrelationItem, col):
    #     self.treeWidget.clearSelection()
    #     # self.update_options(item)   # FIXME: only when we able the options
    #     if col == 1:
    #         s = item.dh_sw
    #         if item.checkState(col) == Qt.Checked:
    #             item.setSelected(True)
    #             if s:
    #                 s.show()
    #             else:
    #                 sub = Charts(item=item, col=col, options={'chart_type': [Charts.SCATTER], 'hide_toolbar':
    #                     self.data_options['hide_toolbar']})
    #                 sub.make_chart()
    #                 self.mdi.addSubWindow(sub)
    #                 sub.show()
    #         else:
    #             if s:
    #                 self.mdi.activatePreviousSubWindow()
    #                 s.close()
    #     elif col == 2:
    #         s = item.dgie_sw
    #         if item.checkState(col) == Qt.Checked:
    #             item.setSelected(True)
    #             if s:  # check if any subwindow has been store
    #                 s.show()
    #             else:
    #                 sub = Charts(item=item, col=col, options={'chart_type': [Charts.SCATTER], 'hide_toolbar':
    #                     self.data_options['hide_toolbar']})
    #                 sub.make_chart()
    #                 self.mdi.addSubWindow(sub)
    #                 sub.show()
    #         else:
    #             if s:
    #                 self.mdi.activatePreviousSubWindow()
    #                 s.close()
    #     elif col == 3:
    #         s = item.dgnmode_sw
    #         if item.checkState(col) == Qt.Checked:
    #             item.setSelected(True)
    #             if s:  # check if any subwindow has been store
    #                 s.show()
    #             else:
    #                 sub = Charts(item=item, col=col, options={'chart_type': [Charts.SCATTER], 'hide_toolbar':
    #                     self.data_options['hide_toolbar']})
    #                 sub.make_chart()
    #                 self.mdi.addSubWindow(sub)
    #                 sub.show()
    #         else:
    #             if s:
    #                 self.mdi.activatePreviousSubWindow()
    #                 s.close()
    #     elif col == 4:
    #         s = item.dgqh_sw
    #         if item.checkState(col) == Qt.Checked:
    #             item.setSelected(True)
    #             if s:  # check if any subwindow has been store
    #                 s.show()
    #             else:
    #                 sub = Charts(item=item, col=col, options={'chart_type': [Charts.SCATTER], 'hide_toolbar':
    #                     self.data_options['hide_toolbar']})
    #                 sub.make_chart()
    #                 self.mdi.addSubWindow(sub)
    #                 sub.show()
    #         else:
    #             if s:
    #                 self.mdi.activatePreviousSubWindow()
    #                 s.close()

    def update_table(self, item: CorrelationItem, col):

        if col == 1:
            data = item.enthalpy
        elif col == 2:
            data = item.dgie
        elif col == 3:
            data = item.dgnmode
        elif col == 4:
            data = item.dgqh
        else:
            return
        col_label = self.correlation_treeWidget.headerItem().text(col)
        self.data_table_energy_col.setText(f"{item.text(0)}({col_label})")

        for row, v in enumerate(data[col_label]):
            titem = QTableWidgetItem(f'{v:.2f}')
            self.data_table_widget.setItem(row, 2, titem)

    def make_correlation(self):

        self.data_table_widget.setRowCount(len(self.corr_data) - 1)
        sys_with_ki = 0
        for x in self.corr_data:
            if x == 'mutant':
                continue
            if np.isnan(self.corr_data[x]['Exp.Energy']):
                continue
            item_s = QTableWidgetItem(x)
            self.data_table_widget.setItem(sys_with_ki, 0, item_s)
            item_e = QTableWidgetItem(f"{self.corr_data[x]['Exp.Energy']:.2f}")
            self.data_table_widget.setItem(sys_with_ki, 1, item_e)
            sys_with_ki += 1
        if sys_with_ki < 3:
            m = QMessageBox.critical(self, 'Unable to calculate correlation',
                                     'Three or more systems are needed to calculate the correlation.',
                                     QMessageBox.Ok)
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
            return

        df = make_corr_DF(self.corr_data) # FIXME:
        # get df for each model
        models = ['gb', 'pb', 'rism std', 'rism gf']
        columns = ['ΔH', 'ΔH+IE', 'ΔH+NMODE', 'ΔH+QH']
        hide_col = []
        c = 1
        for x in columns:
            if df[x].isnull().all():
                hide_col.append(c)
            c += 1

        for m in models:
            model_df = df[df['MODEL'].isin([m])]
            m_col_box = []
            model_data = []
            c = 1
            for col in columns:
                item_col_df = model_df[['System', col, 'Exp.Energy']]
                if not item_col_df[col].isnull().all():
                    m_col_box.append(c)
                    model_data.append(item_col_df)
                else:
                    model_data.append(None)
                c += 1
            if m_col_box:
                item = CorrelationItem(self.correlation_treeWidget, [m.upper()], model=model_df, enthalpy=model_data[0],
                                       dgie=model_data[1], dgnmode=model_data[2], dgqh=model_data[3], col_box=m_col_box)
        for x in hide_col:
            self.correlation_treeWidget.hideColumn(x)

    def process_data(self, rqueue: Queue, options):
        self.data_options = options
        self.init_dialog.close()
        maximum = rqueue.qsize()
        qpd = QProgressDialog('Creating systems tree', 'Abort', 0, maximum, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(1000)

        for sys_index, i in enumerate(range(maximum), start=1):
            qpd.setValue(i)
            if qpd.wasCanceled():
                break
            system, api_data = rqueue.get()
            # print(system)
            name, path, norm_mut, settings = system
            result, namespace = api_data

            # namespace.INPUT['exp_ki'] = exp_ki
            # namespace.INPUT['temperature'] = temp

            if settings == 'User-Default':
                config = 'User-Default'
            elif settings == 'Custom':
                config = path
            else:
                config = None

            self.systems[sys_index] = {'name': name, 'path': path,
                                       'namespace': namespace, 'data': result,
                                       'current_frames': [namespace.INPUT['startframe'],
                                                          namespace.INPUT['endframe'],
                                                          namespace.INPUT['interval']],
                                       'current_nmode_frames': [namespace.INPUT['nmstartframe'],
                                                                namespace.INPUT['nmendframe'],
                                                                namespace.INPUT['nminterval']],
                                       'current_ie_frames': math.ceil(
                                           namespace.INFO['numframes'] * (namespace.INPUT['ie_segment'] / 100)),
                                       'current_c2_frames': math.ceil(
                                           namespace.INFO['numframes'] * (namespace.INPUT['c2_segment'] / 100)),
                                       'chart_options': ChartSettings(config)}


            self.makeTree(sys_index, options)

        qpd.setLabelText(f"Processing data of {self.items_counter['charts']} items")
        qpd.setMaximum(self.items_counter['charts'])
        i = 0
        it = QTreeWidgetItemIterator(self.treeWidget)
        while it.value():
            item = it.value()
            if item:
                qpd.setValue(i)
                # item.get_data()
                i += 1
            if qpd.wasCanceled():
                break
            it += 1
        qpd.setValue(self.items_counter['charts'])

        self._initialize_systems()

        # if i != self.items_counter['charts']:
        #     self.close()

        if not options['corr_sys']:  # FIXME:
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
        else:
            self.make_correlation()

        # some late signal/slot connections
        # self.treeWidget.itemChanged.connect(self.showdata)
        self.correlation_treeWidget.itemChanged.connect(self.showcorr)

    def makeTree(self, sys_index, options):

        # make system item

        self.sys_item = CustomItem(self.treeWidget, [self.systems[sys_index]['name']], app=self, system_index=sys_index,
                                   buttons=(-1,), remove_empty_terms=options['remove_empty_terms'])
        self.normalitem = CustomItem(self.sys_item, ['Normal'])
        self.makeItems(sys_index, self.normalitem, options)
        self.normalitem.setExpanded(True)

        if self.systems[sys_index]['data'].mutant:
            self.mutantitem = CustomItem(self.sys_item, ['Mutant'])
            self.makeItems(sys_index, self.mutantitem, options, 1)
            self.mutantitem.setExpanded(True)

            self.mut_norm_item = CustomItem(self.sys_item, ['Mutant-Normal'])

            norm_data = self.systems[sys_index]['data']
            mut_data = self.systems[sys_index]['data'].mutant
            delta_data = self.systems[sys_index]['data'].delta = {}
            for x, y in zip(norm_data, mut_data):
                if x in ['ie', 'c2']:
                    continue
                elif x == 'decomp':
                    continue
                else:
                    delta_data[x] = mut_data[x] - norm_data[x]
                # if isinstance(x, pandas.DataFrame):
                #     self.systems[sys_index]['data'].delta = {x: }

            self.makeItems(sys_index, self.mut_norm_item, options, mutant=2)
            self.mut_norm_item.setExpanded(True)

        self.sys_item.setExpanded(True)

        itemiter = QTreeWidgetItemIterator(self.sys_item)
        while itemiter.value():
            item = itemiter.value()
            if item.item_type in ['energy', 'ie', 'c2']:
                frange = self.systems[item.system_index]['current_frames']
            else:
                frange = self.systems[item.system_index]['current_nmode_frames']

            if item.item_type == 'ie':
                eframes = self.systems[item.system_index]['current_ie_frames']
            elif item.item_type == 'c2':
                eframes = self.systems[item.system_index]['current_c2_frames']
            else:
                eframes = 0
            item.setup_data(frange, eframes)
            sb = item.setup_buttons()

            if sb:
                self.treeWidget.setItemWidget(item, 1, sb)
            itemiter += 1

    def _remove_empty(self, data, options, namespace):
        if options['remove_empty_terms'] and (
                'UB' in data.name
                or 'IMP' in data.name
                or 'CMAP' in data.name
                and not namespace.INFO['using_chamber']
        ):
            return True
        elif options['remove_empty_charts'] and (all(data > -0.01) and all(data < 0.01)):
            # FIXME: Do we need to clarify that they are not terms like Gsolv, Ggas and TOTAL?
            return True

    def _itemdata_properties(self, data, decomp=False):
        """
        Pre-processing the items data.
        Get the following properties:
        - scalable: if any value is gt 10^2
        - separable: if contains subcategories (DH [energetics components], Per-residue[receptor and ligand])

        Also, remove empty terms and charts according to selected options
        @param data:
        @return:
        """
        groups = {}
        if not decomp:
            sep_ggas_keys = []
            sep_gsolv_keys = []
            # remove empty charts? (BOND, ANGLE and DIHEDRAL for STP)
            ggas_keys = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'UB', 'IMP', 'CMAP']
            gsolv_keys = ['EGB', 'ESURF', 'EPB', 'ENPOLAR', 'POLAR SOLV', 'APOLAR SOLV']
            for k in data.columns:
                if k in ggas_keys:
                    sep_ggas_keys.append(k)
                elif k in gsolv_keys:
                    sep_gsolv_keys.append(k)
            if sep_ggas_keys:
                groups['GGAS'] = sep_ggas_keys
                groups['GSOLV'] = sep_gsolv_keys
                groups['TOTAL'] = ['GGAS', 'GSOLV', 'TOTAL']
        else:
            groups['Receptor'] = []
            groups['Ligand'] = []
            for k in data.columns:
                if k[0].startswith('R:') and k[0] not in groups['Receptor']:
                    groups['Receptor'].append(k[0])
                elif k[0].startswith('L:') and k[0] not in groups['Ligand']:
                    groups['Ligand'].append(k[0])
        return groups

    def makeItems(self, sys_index, topItem, options, mutant=0):
        correlation_data = self.corr_data
        mut_pre = ''
        if mutant == 1:
            mut_pre = 'Mut. '
            correlation_data = self.corr_data['mutant']

        sys_name = self.systems[sys_index]['name']
        # correlation_data[sys_name] = {'ΔG': {
        #                                     'gb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'pb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism std': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism gf': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan}},
        #                               'Exp.Energy': ki2energy(topItem.exp_ki, topItem.temp)}

        if mutant == 1:
            data = self.systems[sys_index]['data'].mutant
        elif mutant == 2:
            data = self.systems[sys_index]['data'].delta
        else:
            data = self.systems[sys_index]['data']
        namespace = self.systems[sys_index]['namespace']

        parts = options['components'] + ['delta']
        if namespace.FILES.stability:
            parts.append('complex')

        for level in data:
            if level in ['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh']:
                # self._make_eitems(topItem, level, data[level])
                titem = CustomItem(topItem, [level.upper()])
                str_dict = multiindex2dict(data[level].columns)
                for level1 in str_dict:
                    if level1 not in parts:
                        continue
                    item1 = CustomItem(
                        titem,
                        [level1.upper()],
                        app=self,
                        level=1,
                        buttons=(2,),
                        chart_title="Energetic Components",
                        chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | {level1.upper()}"
                    )
                    for level2 in str_dict[level1]:
                        if self._remove_empty(data[level][(level1, level2)], options, namespace):
                            del data[level][(level1, level2)]
                            continue
                        item2 = CustomItem(
                            item1,
                            [level2.upper()],
                            data=data[level][(level1, level2)],
                            app=self,
                            buttons=(1,),
                            chart_title="Energetic Components",
                            chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | {level1.upper()} | {level2.upper()}"
                        )
                        self.items_counter['charts'] += 1
                    if level1 in data[level]:
                        item1.data = data[level][(level1,)]
                        item1.properties['scalable'] = data[level][(level1,)].gt(300).any().any()
                        item1.properties['groups'] = self._itemdata_properties(data[level][(level1,)])
                    self.items_counter['charts'] += 1
                self.items_counter['charts'] += 1
            elif level == 'decomp':
                # omit decomp data
                if not options['decomposition']:
                    continue
                titem = CustomItem(topItem, [level.upper()])
                for level1 in data[level]:
                    # GB or PB
                    dat = data[level][level1]
                    item = CustomItem(titem, [level1.upper()])
                    str_dict = multiindex2dict(dat.columns)
                    # Complex, receptor, ligand and delta
                    for level2 in str_dict:
                        if level2 not in parts:
                            continue
                        item2 = CustomItem(item, [level2.upper()])
                        # TDC, SDC, BDC
                        item_lvl = 2 if namespace.INPUT['idecomp'] in [1, 2] else 3
                        title = '[Per-residue]' if namespace.INPUT['idecomp'] in [1, 2] else '[Per-wise]'
                        for level3 in str_dict[level2]:
                            item3 = CustomItem(
                                item2,
                                [level3.upper()],
                                data=dat[(level2, level3)],
                                app=self,
                                level=item_lvl,
                                buttons=(1, 2, 3, 4),
                                chart_title=f"Energetic Components {title}",
                                chart_subtitle=f"{mut_pre}{sys_name} | "
                                               f"{str(level).upper()} | "
                                               f"{str(level1).upper()} | "
                                               f"{str(level2).upper()} | "
                                               f"{str(level3).upper()}"
                            )
                            item3.properties['scalable'] = dat[(level2, level3)].gt(300).any().any()
                            item3.properties['groups'] = self._itemdata_properties(dat[(level2, level3)], decomp=True)
                            # # residue first level
                            btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                            item_lvl2 = 1 if namespace.INPUT['idecomp'] in [1, 2] else 2

                            for level4 in str_dict[level2][level3]:
                                item4 = CustomItem(
                                    item3,
                                    [level4.upper()],
                                    data=dat[(level2, level3, level4)],
                                    app=self,
                                    level=item_lvl2,
                                    buttons=btns,
                                    chart_title="Energetic Components",
                                    chart_subtitle=f"{mut_pre}{sys_name} | "
                                                   f"{str(level).upper()} | "
                                                   f"{str(level1).upper()} | "
                                                   f"{str(level2).upper()} | "
                                                   f"{str(level3).upper()} | "
                                                   f"{str(level4).upper()}"
                                )
                                item4.properties['scalable'] = dat[(level2, level3, level4)].gt(300).any().any()
                                # energetics terms
                                for level5 in str_dict[level2][level3][level4]:
                                    if namespace.INPUT['idecomp'] in [1, 2]:
                                        item5 = CustomItem(item4, [level5.upper()],
                                                           data=dat[(level2, level3, level4, level5)], app=self,
                                                           buttons=(1,),
                                                           chart_title="Energetic Components [Per-residue]",
                                                           chart_subtitle=f"{mut_pre}{sys_name} | "
                                                                          f"{str(level).upper()} | "
                                                                          f"{str(level1).upper()} | "
                                                                          f"{str(level2).upper()} | "
                                                                          f"{str(level3).upper()} | "
                                                                          f"{str(level4).upper()} | "
                                                                          f"{str(level5).upper()}"
                                                           )
                                        self.items_counter['charts'] += 1
                                    else:
                                        item5 = CustomItem(item4, [level5.upper()],
                                                           data=dat[(level2, level3, level4, level5)], app=self,
                                                           level=1, buttons=(2,))
                                        self.items_counter['charts'] += 1
                                        # energetics terms
                                        for level6 in str_dict[level2][level3][level4][level5]:
                                            item6 = CustomItem(
                                                item5,
                                                [level6.upper()],
                                                data=dat[(level2, level3, level4, level5, level6)],
                                                app=self,
                                                buttons=(1,),
                                                chart_title="Energetic Components [Per-wise]",
                                                chart_subtitle=f"{mut_pre}{sys_name} | "
                                                               f"{str(level).upper()} | "
                                                               f"{str(level1).upper()} | "
                                                               f"{str(level2).upper()} | "
                                                               f"{str(level3).upper()} | "
                                                               f"{str(level4).upper()} | "
                                                               f"{str(level5).upper()} | "
                                                               f"{str(level6).upper()}"
                                            )
                                            self.items_counter['charts'] += 1
                                self.items_counter['charts'] += 1
                            self.items_counter['charts'] += 1
                        self.items_counter['charts'] += 1
                    self.items_counter['charts'] += 1
                self.items_counter['charts'] += 1
            elif level == 'ie':
                titem = CustomItem(topItem, [level.upper()])
                # print(data[level])

                for level1 in data[level]:
                    dat = data[level][level1]
                    item1 = CustomItem(
                        titem,
                        [level1.upper()],
                        data=dat['data']['data'],
                        app=self,
                        level=0,
                        item_type='ie',
                        buttons=(1,),
                        iec2_data={'sigma': dat['sigma'], 'frames': dat['ieframes']},
                        chart_title="Interaction Entropy",
                        chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | {level1.upper()}"
                    )
                    # for level2 in dat:
                    # item2 = CustomItem(titem, [level2.upper()], data=dat[level2], app=self, level=1,
                    #                        buttons=(1,), iec2_frames=dat[])
