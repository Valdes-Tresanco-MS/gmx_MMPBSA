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
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from GMXMMPBSA.analyzer.dialogs import InitDialog
from GMXMMPBSA.analyzer.customitem import CustomItem, CorrelationItem
from GMXMMPBSA.analyzer.utils import energy2pdb_pml, ki2energy, make_corr_DF, multiindex2dict
from GMXMMPBSA.analyzer.chartsettings import ChartSettings
from GMXMMPBSA.analyzer.parametertree import ParameterTree, Parameter
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
        self.showMaximized()
        self.corr_data = {'mutant': {}}

        self.systems = {}
        self.current_system_index = None
        self.pymol_p_list = []

        self.items_counter = {'charts': 0, 'pymol': [], 'bars': 0, 'line': 0, 'heatmap': 0}

        self.mdi = QMdiArea(self)
        self.mdi.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setCentralWidget(self.mdi)

        self.treeDockWidget = QDockWidget('Data', self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.treeDockWidget)

        self.correlation_DockWidget = QDockWidget('Correlations', self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.correlation_DockWidget)

        self.optionDockWidget = QDockWidget('Options', self)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.optionDockWidget)

        self.treeWidget = QTreeWidget(self)
        self.treeWidget.setMinimumWidth(380)
        self.treeWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        # self.treeWidget.customContextMenuRequested.connect(self.data_context_menu)
        self.treeWidget.itemSelectionChanged.connect(self.update_system_selection)

        self._make_options_panel()

        self.treeDockWidget.setWidget(self.treeWidget)
        self.optionDockWidget.setWidget(self.optionWidget)

        sys_label = QTreeWidgetItem(['System', 'Charts'])
        sys_label.setToolTip(0, 'System')
        sys_label.setToolTip(1, 'Charts')

        self.treeWidget.setHeaderItem(sys_label)
        self.treeWidget.setColumnHidden(4, True)
        header = self.treeWidget.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.Stretch)

        self.correlation_treeWidget = QTreeWidget(self)
        self.correlation_treeWidget.itemClicked.connect(self.update_table)
        self.correlation_treeWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        # self.correlation_treeWidget.customContextMenuRequested.connect(self.corr_context_menu)
        model_label = QTreeWidgetItem(['MODEL', 'ΔGeff', 'ΔGie', 'ΔGc2', 'ΔGnm', 'ΔGqh'])
        model_label.setToolTip(0, 'Selected Model')
        model_label.setToolTip(1, '''<html>Correlation plot for ΔG<sub>effective</sub>. Energy 
        when entropy contribution is neglected. Calculated as ΔG<sub>effective</sub> = ΔH</html>''')
        model_label.setToolTip(2, '''<html>Correlation plot for ΔG<sub>ie</sub>. Energy 
        when entropy contribution is calculated using the Interaction Entropy method. Calculated as 
        ΔG<sub>ie</sub> = ΔH + -TΔS<sub>IE</sub></html>''')
        model_label.setToolTip(3, '''<html>Correlation plot for ΔG<sub>c2</sub>. Energy 
        when entropy contribution is calculated using the C2 Entropy method. Calculated as 
        ΔG<sub>c2</sub> = ΔH + -TΔS<sub>C2</sub></html>''')
        model_label.setToolTip(4, '''<html>Correlation plot for ΔG<sub>nm</sub>. Energy 
        when entropy contribution is calculated using the Normal Modes Entropy method. Calculated as 
        ΔG<sub>nm</sub> = ΔH + -TΔS<sub>NMODE</sub></html>''')
        model_label.setToolTip(5, '''<html>Correlation plot for ΔG<sub>ie</sub>. Energy 
        when entropy contribution is calculated using the Quasi-Harmonic Entropy method. Calculated as 
        ΔG<sub>qh</sub> = ΔH + -TΔS<sub>QH</sub></html>''')
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

        self.menubar = self.menuBar()
        # self.opendirAct = QAction("Open gmx_MMPBSA Dir...", self)
        # self.opendirAct.triggered.connect(self.getInfoFile)
        self.fileMenu = self.menuBar().addMenu("&File")
        # self.fileMenu.addAction(self.opendirAct)
        self.fileMenu.addAction('Close All..', self.mdi.closeAllSubWindows)
        self.viewMenu = self.menuBar().addMenu("&View")
        self.viewMenu.addAction('Show Data', self.treeDockWidget.show)
        self.viewMenu.addAction('Show Correlation', self.correlation_DockWidget.show)
        self.viewMenu.addAction('Show Options', self.optionDockWidget.show)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction('Tile SubWindows', self.mdi.tileSubWindows)
        self.viewMenu.addAction('Cascade SubWindows', self.mdi.cascadeSubWindows)
        self.aboutMenu = self.menuBar().addMenu("&About")
        self.aboutMenu.addAction('Help', self._help)
        self.aboutMenu.addAction('Documentation', self.treeDockWidget.show)
        self.aboutMenu.addAction('Report a bug', self.treeDockWidget.show)
        self.aboutMenu.addAction('Google group', self.treeDockWidget.show)
        self.aboutMenu.addAction('About gmx_MMPBSA', self._about_dialog)
        self.statusbar = self.statusBar()


        self.init_dialog = InitDialog(self)

    def _about_dialog(self):
        from GMXMMPBSA import __version__
        QMessageBox.about(self, "About gmx_MMPBSA",
                          "<h2>About gmx_MMPBSA</h2>"
                          "<b>gmx_MMPBSA</b> is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free "
                          "energy calculations with GROMACS files.<br><br>"
                          f"<b>Version:</b> {__version__}"
                          "<h2>Cite gmx_MMPBSA</h2>"
                          "Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. "
                          "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. "
                          "Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291. "
                          "<a href='https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645'>"
                          "https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645</a>")

    def _help(self):
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/analyzer/'))

    def _doc(self):
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/'))

    def _bug(self):
        QDesktopServices().openUrl(QUrl('https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/new/choose'))

    def _group(self):
        QDesktopServices().openUrl(QUrl('https://groups.google.com/g/gmx_mmpbsa'))

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

        conf_menu = QMenu()
        conf_menu.setTitle('Option configuratio')
        self.default_action = conf_menu.addAction('Set as default')
        self.import_action = conf_menu.addAction('Import')
        # self.line_table_action.toggled.connect(self._show_line_table)

        self.options_conf_btn = QToolButton()
        self.options_conf_btn.setIcon(self.style().standardIcon(QStyle.SP_DriveFDIcon))
        self.options_conf_btn.setPopupMode(QToolButton.InstantPopup)
        self.options_conf_btn.setMenu(conf_menu)

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
        optionWidget_c.setCornerWidget(self.options_conf_btn)

        # Charts options
        self.chart_options_param = Parameter().create()
        self.chart_options_w = ParameterTree()

        # properties_w.setParameters(p, showTop=False)
        # charts_options_w = QWidget()
        optionWidget_c.addTab(self.chart_options_w, 'Charts Options')


        frames_w = QWidget(optionWidget_c)
        optionWidget_c.addTab(frames_w, 'Frames')
        frames_l = QVBoxLayout(frames_w)

        frames_group = QGroupBox('Energy')
        frames_l.addWidget(frames_group)
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
        fg_l1.addRow('Start', self.eframes_start_sb)
        fg_l1.addRow('End', self.eframes_end_sb)
        fg_l2 = QFormLayout()
        fg_l2.addRow('Interval', self.eframes_inter_sb)
        fg_l2.addRow('Nr. frames', self.numframes_le)

        frames_group_l = QHBoxLayout(frames_group)
        frames_group_l.addLayout(fg_l1, 1)
        frames_group_l.addLayout(fg_l2, 1)

        frames_l.addWidget(frames_group)
        # frames_w.setTabToolTip(0, 'Frames defined for all energy calculations')

        nmframes_group = QGroupBox('nmode')
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

        frames_l.addWidget(nmframes_group)
        # frames_w.setTabToolTip(2, 'Frames defined for NMODE calculation')
        #
        ieframes_group = QGroupBox('IE')

        self.iesegment_sb = QSpinBox()
        self.iesegment_sb.setRange(1, 100)
        self.iesegment_sb.setAccelerated(True)
        # self.segment_sb.valueChanged.connect(self.frames_start_sb_update)

        self.ienumframes_le = QLineEdit()
        self.ienumframes_le.setReadOnly(True)

        ie_l = QFormLayout(ieframes_group)
        ie_l.addRow('Segment', self.iesegment_sb)
        ie_l.addRow('Nr. frames', self.ienumframes_le)

        frames_l.addWidget(ieframes_group)
        frames_l.addStretch(1)



        update_btn = QPushButton('Update')
        update_btn.clicked.connect(self.update_fn)
        resetchart_btn = QPushButton('Reset')

        btn_l = QHBoxLayout()
        btn_l.addWidget(resetchart_btn)
        btn_l.addWidget(update_btn)
        optionWidget_l.addLayout(btn_l)

        optionWidget_c.currentChanged.connect(self._control_options_config)

    def _control_options_config(self, ind):
        if ind:
            self.options_conf_btn.setEnabled(False)
        else:
            self.options_conf_btn.setEnabled(True)

    def update_system_selection(self):
        if not self.treeWidget.selectedItems():
            return
        parent_item = self.treeWidget.selectedItems()[0]

        if self.current_system_index == parent_item.system_index:
            return
        while parent_item:
            if next_item := parent_item.parent():
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
                self.systems[e]['current_nmode_frames'] = act_frange
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

        qpd.setLabelText('Updating open charts')

        subwindows = self.mdi.subWindowList()
        for s in processed_sys:
            changes = self.systems[s]['chart_options'].changes
            for sub in subwindows:
                if not sub.item_parent or sub.item_parent.system_index != s:
                    continue
                if not sub.isVisible():
                    continue
                if (changes['bar_action'] or changes['line_ie_action'] or
                        changes['line_action'] or changes['heatmap_action']):
                    sub.button.setChecked(False)
                    sub.button.setChecked(True)
            pymol_items = [[p, item] for p, item in self.pymol_p_list if p.state() == QProcess.Running]
            for p, item in pymol_items:
                p.kill()
                p.waitForFinished()
                item.vis_action.setChecked(True)

        # re-assign changes to native state after replot all open charts
        for s in processed_sys:
            self.systems[s]['chart_options'].changes = dict(line_action=0,
                                                            line_ie_action=0,
                                                            bar_action=0,
                                                            heatmap_action=0,
                                                            visualization_action=0)


        qpd.setValue(maximum)

    def reset_dc(self):
        sub = self.mdi.activeSubWindow()
        sub.item.reset_data()
        sub.make_chart()

    def frames_start_sb_update(self, value):
        if value > self.eframes_end_sb.value():
            self.eframes_end_sb.setValue(value)

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

        df = make_corr_DF(self.corr_data)  # FIXME:
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

    def read_data(self, queue: Queue, options):
        self.init_dialog.close()
        max_sixe = queue.qsize()
        qpd = QProgressDialog('Reading output files', 'Abort', 0, max_sixe, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(0)
        # qpd.setRange(0, queue.qsize())
        results = []

        for x in range(max_sixe):
            sys_name, path, exp_ki, options_file = queue.get()
            gmx_mmpbsa_api = MMPBSA_API()
            gmx_mmpbsa_api.set_config(options['timestart'], options['timestep'], options['timeunit'])
            gmx_mmpbsa_api.load_file(path)
            results.append([(sys_name, path, exp_ki, options_file), gmx_mmpbsa_api])
            qpd.setValue(x)
            queue.task_done()
            if qpd.wasCanceled():
                break
        qpd.setValue(max_sixe)

        self.process_data(results, options)

    def process_data(self, results: list, options):
        self.data_options = options
        maximum = len(results) + 3
        qpd = QProgressDialog('Creating systems tree', 'Abort', 0, maximum, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(0)

        for i, c in enumerate(range(len(results)), start=1):
            qpd.setValue(i)
            if qpd.wasCanceled():
                break
            system, api = results[c]
            name, path, exp_ki, settings = system
            energy = api.get_energy()
            namespace = api.app_namespace
            summary = api.get_summary()

            if settings == 'User-Default':
                config = 'User-Default'
            elif settings == 'Custom':
                config = path
            else:
                config = None

            self.systems[i] = {'name': name, 'path': path, 'api': api,
                               'namespace': namespace, 'data': energy,
                               'current_frames': [namespace.INPUT['startframe'],
                                                  namespace.INPUT['endframe'],
                                                  namespace.INPUT['interval']],
                               'current_nmode_frames': [namespace.INPUT['nmstartframe'],
                                                        namespace.INPUT['nmendframe'],
                                                        namespace.INPUT['nminterval']],
                               'current_ie_frames': math.ceil(
                                   namespace.INFO['numframes'] * (namespace.INPUT['ie_segment'] / 100)),
                               'chart_options': ChartSettings(config),
                               'options': options,
                               'items_data': {}, 'items_summary': summary,
                               'exp_ki': exp_ki
                               }
            self.makeTree(i)
        qpd.setLabelText('Initializing first system...')
        self._initialize_systems()
        qpd.setValue(i + 1)

        if not options['corr_sys']:  # FIXME:
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
        else:
            qpd.setLabelText('Calculating correlation...')
            self.make_correlation()
            qpd.setValue(i + 2)

        qpd.setValue(maximum)

        # some late signal/slot connections
        # self.treeWidget.itemChanged.connect(self.showdata)
        # self.correlation_treeWidget.itemChanged.connect(self.showcorr)

    def makeTree(self, sys_index):

        # make system item
        sys_item = CustomItem(self.treeWidget, [self.systems[sys_index]['name']], app=self, system_index=sys_index,
                                   buttons=(-1,))

        # FIXME: Binding
        classif_item = CustomItem(sys_item, ['ΔH/-TΔS/ΔG'])
        classif_item.setExpanded(True)

        print_keys = list(self.systems[sys_index]['exp_ki'].keys())
        decomp_print_keys = [f"decomp_{x}" for x in self.systems[sys_index]['exp_ki'].keys()]
        if 'normal' in print_keys and 'mutant' in print_keys:
            print_keys.append('mutant-normal')

        for part in print_keys:
            if part not in self.systems[sys_index]['data']:
                continue
            self.makeItems(sys_index, part, classif_item)
            self.setting_item_data(sys_index, part)

        if self.systems[sys_index]['namespace'].INPUT['idecomp']:
            classif_item = CustomItem(sys_item, ['Decomposition'])
            for part in decomp_print_keys:
                if part not in self.systems[sys_index]['data']:
                    continue
                self.makedecompItems(sys_index, part, classif_item)
                self.setting_item_data(sys_index, part)

        sys_item.setExpanded(True)

        # setup items to qtreewidget
        itemiter = QTreeWidgetItemIterator(sys_item)
        while itemiter.value():
            item = itemiter.value()
            if sb := item.setup_buttons():
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
            # FIXME: Do we need to clarify that they are not terms like GSOLV, GGAS and TOTAL?
            return True

    def _itemdata_properties(self, data, decomp=False):
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

    def _setup_data(self, data, iec2=False, level=0):

        # this variable show if the data changed or not. At first time, it is true, then when plotting become in false
        change = True

        cont = {'ie_plot_data': None, 'line_plot_data': None, 'bar_plot_data': None, 'heatmap_plot_data': None}
        if level == 0:
            options = {'ie': True} if iec2 else {}
            cont['line_plot_data'] = [data, options, change]
        elif level == 1:
            options = {'c2': True} if iec2 else {}
            options.update(dict(groups=self._itemdata_properties(data)))
            cont['bar_plot_data'] = [data, options, change]
        elif level == 2:
            tempdf = data.loc[:, data.columns.get_level_values(1) == 'tot']
            bar_plot_data = tempdf.droplevel(level=1, axis=1)
            cont['bar_plot_data'] = [bar_plot_data, dict(groups=self._itemdata_properties(bar_plot_data)), change]
            cont['line_plot_data'] = [bar_plot_data.sum(axis=1), {}, change]
            cont['heatmap_plot_data'] = [bar_plot_data.transpose(copy=True), {}, change]
            del bar_plot_data
            del tempdf
        elif level == 3:
            # Select only the "tot" column, remove the level, change first level of columns to rows and remove the mean
            # index
            tempdf = data.loc[:, data.columns.get_level_values(2) == 'tot']
            cont['heatmap_plot_data'] = [
                tempdf.aggregate(["mean"]).droplevel(level=2, axis=1).stack().droplevel(level=0), {}, change]
            bar_plot_data = tempdf.groupby(axis=1, level=0, sort=False).sum()
            cont['bar_plot_data'] = [bar_plot_data, dict(groups=self._itemdata_properties(bar_plot_data)), change]
            cont['line_plot_data'] = [bar_plot_data.sum(axis=1), {}, change]
            del tempdf
            del bar_plot_data

        return cont

    def makeItems(self, sys_index, part, classif_item, mutant=0):

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
        #                               'Exp.Energy': ki2energy(sys_item.exp_ki, sys_item.temp)}

        data = self.systems[sys_index]['data'][part]
        namespace = self.systems[sys_index]['namespace']


        top_item = CustomItem(classif_item, [part.capitalize()])
        top_item.setExpanded(True)
        parts = self.systems[sys_index]['options']['components'] + ['delta']
        if namespace.FILES.stability:
            parts.append('complex')

        for level in ['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh', 'ie', 'c2', 'binding']:
            if level in data:
                if level in ['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh']:
                    titem = CustomItem(top_item, [level.upper()])
                    titem.setExpanded(True)

                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        if level1 not in parts:
                            continue
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(2, -2),
                                           title="Energetic Components",
                                           subtitle=f"{sys_name} | {level.upper()} | {level1.upper()}",
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )
                        if level == 'qh':
                            continue
                        for level2 in str_dict[level1]:
                            item2 = CustomItem(item1,
                                               [level2.upper()],
                                               app=self,
                                               buttons=(1,),
                                               title="Energetic Components",
                                               subtitle=f"{sys_name} | {level.upper()} | "
                                                        f"{level1.upper()} | {level2.upper()}",
                                               keys_path=(part, level, (level1, level2)),
                                               part=part
                                               )
                            self.items_counter['charts'] += 1
                elif level == 'c2':
                    titem = CustomItem(top_item, [level.upper()])
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(2,),
                                           title="C2 Entropy",
                                           subtitle=f"{sys_name} | {level.upper()} | {level1.upper()}",
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )
                elif level == 'ie':
                    titem = CustomItem(top_item, [level.upper()])
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(1,),
                                           title="Interaction Entropy",
                                           subtitle=f"{sys_name} | {level.upper()} | {level1.upper()}",
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )
                elif level == 'binding':
                    titem = CustomItem(top_item, ['Binding Energy'])
                    for level1 in data[level]:
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(2,),
                                           title="Binding Energy",
                                           subtitle=f"{sys_name} | ΔG | {level1.upper()}",
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )

    def makedecompItems(self, sys_index, part, classif_item):
        sys_name = self.systems[sys_index]['name']
        data = self.systems[sys_index]['data'][part]
        namespace = self.systems[sys_index]['namespace']
        top_item = CustomItem(classif_item, [part.capitalize()])
        top_item.setExpanded(True)
        parts = self.systems[sys_index]['options']['components'] + ['delta']
        if namespace.FILES.stability:
            parts.append('complex')

        for level in ['gb', 'pb']:
            if level in data:
                titem = CustomItem(top_item, [level.upper()])
                titem.setExpanded(True)

                str_dict = multiindex2dict(data[level].columns)
                for level1 in str_dict:
                    if level1 not in parts:
                        continue
                    item1 = CustomItem(titem, [level1.upper()])
                    for level2 in str_dict[level1]:
                        title = '[Per-residue]' if namespace.INPUT['idecomp'] in [1, 2] else '[Per-wise]'
                        item2 = CustomItem(item1,
                                           [level2.upper()],
                                           app=self,
                                           buttons=(1, 2, 3, 4),
                                           title=f"Energetic Components {title}",
                                           subtitle=f"{sys_name} | {level.upper()} | "
                                                    f"{level1.upper()} | {level2.upper()}",
                                           keys_path=(part, level, (level1, level2)),
                                           part=part
                                           )
                        self.items_counter['charts'] += 1
                        # residue first level
                        btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                        for level3 in str_dict[level1][level2]:

                            item3 = CustomItem(item2,
                                               [level3.upper()],
                                               app=self,
                                               buttons=btns,
                                               # title=f"Energetic Components {title}",
                                               # subtitle=f"{sys_name} ({part}) | "
                                               #                f"{str(level).upper()} | "
                                               #                f"{str(level1).upper()} | "
                                               #                f"{str(level2).upper()} | "
                                               #                f"{str(level3).upper()}"
                                               keys_path=(part, level, (level1, level2, level3))
                                               )
                            # residue first level
                            #             btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                            for level4 in str_dict[level1][level2][level3]:
                                if namespace.INPUT['idecomp'] in [1, 2]:
                                    item5 = CustomItem(item3,
                                                       [level4.upper()],
                                                       app=self,
                                                       buttons=(1,),
                                                       title="Energetic Components [Per-residue]",
                                                       subtitle=f"{sys_name} ({part})| "
                                                                      f"{str(level).upper()} | "
                                                                      f"{str(level1).upper()} | "
                                                                      f"{str(level2).upper()} | "
                                                                      f"{str(level3).upper()} | "
                                                                      f"{str(level4).upper()}",
                                                       keys_path=(part, level, (level1, level2, level3, level4))
                                                       )
                                    self.items_counter['charts'] += 1
                                else:
                                    item4 = CustomItem(item3,
                                                       [level4.upper()],
                                                       app=self,
                                                       buttons=(2,),
                                                       keys_path=(part, level, (level1, level2, level3, level4))
                                                       )
                                    self.items_counter['charts'] += 1
                                    # energetics terms
                                    for level5 in str_dict[level1][level2][level3][level4]:
                                        item5 = CustomItem(item4,
                                                           [level5.upper()],
                                                           app=self,
                                                           buttons=(1,),
                                                           title="Energetic Components [Per-wise]",
                                                           subtitle=f"{sys_name} ({part}) | "
                                                                    f"{str(level).upper()} | "
                                                                    f"{str(level1).upper()} | "
                                                                    f"{str(level2).upper()} | "
                                                                    f"{str(level3).upper()} | "
                                                                    f"{str(level4).upper()} | "
                                                                    f"{str(level5).upper()}",
                                                           keys_path=(
                                                           part, level, (level1, level2, level3, level4, level5))
                                                           )

    def setting_item_data(self, sys_index, part, comp=('all')):
        correlation_data = self.corr_data
        mut_pre = ''
        # if mutant == 1:
        #     mut_pre = 'Mut. '
        #     correlation_data = self.corr_data['mutant']

        # correlation_data[sys_name] = {'ΔG': {
        #                                     'gb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'pb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism std': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism gf': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan}},
        #                               'Exp.Energy': ki2energy(topItem.exp_ki, topItem.temp)}

        data = self.systems[sys_index]['data'][part]
        namespace = self.systems[sys_index]['namespace']

        parts = self.systems[sys_index]['options']['components'] + ['delta']
        if namespace.FILES.stability:
            parts.append('complex')

        key_list = []
        if 'nmode' in comp:
            key_list.append('nmode')
        elif 'energy' in comp:
            key_list.extend(['gb', 'pb', 'rism gf', 'rism std'])
        elif 'all' in comp:
            key_list.extend(['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh', 'c2'])
        for level in ['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh', 'ie', 'c2', 'binding']:
            if level in data:
                if level in ['gb', 'pb', 'rism gf', 'rism std', 'nmode', 'qh']:
                    # FIXME: include the decomp
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        if level1 not in parts:
                            continue
                        if not part.startswith('decomp'):
                            self.systems[sys_index]['items_data'][(part, level, (level1,))] = self._setup_data(
                                    data[level][(level1,)], level=1)
                        for level2 in str_dict[level1]:
                            # if self._remove_empty(data[level][(level1, level2)], options, namespace):
                            #     del data[level][(level1, level2)]
                            #     continue
                            if not part.startswith('decomp'):
                                temp_dat = data[level][(level1, level2)]
                                temp_dat.name = level2
                                self.systems[sys_index]['items_data'][(part, level, (level1, level2))] = self._setup_data(
                                    temp_dat)
                                del temp_dat
                            else:
                                item_lvl = 2 if namespace.INPUT['idecomp'] in [1, 2] else 3
                                temp_dat = data[level][(level1, level2)]
                                temp_dat.name = level2
                                self.systems[sys_index]['items_data'][(part, level, (level1, level2))] = \
                                    self._setup_data(temp_dat, level=item_lvl)
                                del temp_dat
                                # residue first level
                                btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                                for level3 in str_dict[level1][level2]:
                                    item_lvl2 = 1 if namespace.INPUT['idecomp'] in [1, 2] else 2
                                    temp_dat = data[level][(level1, level2, level3)]
                                    temp_dat.name = level3
                                    self.systems[sys_index]['items_data'][(part, level, (level1, level2, level3))] = \
                                        self._setup_data(temp_dat, level=item_lvl2)
                                    del temp_dat
                                    # residue first level
                                    #             btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                                    for level4 in str_dict[level1][level2][level3]:
                                        temp_dat = data[level][(level1, level2, level3, level4)]
                                        temp_dat.name = level4
                                        if namespace.INPUT['idecomp'] in [1, 2]:
                                            self.systems[sys_index]['items_data'][
                                                (part, level, (level1, level2, level3, level4))] = \
                                                self._setup_data(temp_dat)
                                            del temp_dat
                                        else:
                                            self.systems[sys_index]['items_data'][
                                                (part, level, (level1, level2, level3, level4))] = \
                                                self._setup_data(temp_dat, level=1)
                                            del temp_dat
                                            # energetics terms
                                            for level5 in str_dict[level1][level2][level3][level4]:
                                                temp_dat = data[level][(level1, level2, level3, level4, level5)]
                                                temp_dat.name = level5
                                                self.systems[sys_index]['items_data'][
                                                    (part, level, (level1, level2, level3, level4, level5))] = \
                                                    self._setup_data(temp_dat)
                                                del temp_dat



                elif level == 'c2':
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        self.systems[sys_index]['items_data'][(part, level, (level1,))] = self._setup_data(
                            data[level][(level1,)], iec2=True, level=1)
                elif level == 'ie':
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        self.systems[sys_index]['items_data'][(part, level, (level1,))] = self._setup_data(
                            data[level][(level1,)], iec2=True, level=0)
                elif level == 'binding':
                    for level1 in data[level]:
                        self.systems[sys_index]['items_data'][(part, level, (level1,))] = self._setup_data(
                            data[level][level1], level=1)