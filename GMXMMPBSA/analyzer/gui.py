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
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from GMXMMPBSA.analyzer.dialogs import InitDialog
from GMXMMPBSA.analyzer.customitem import CustomItem, CorrelationItem
from GMXMMPBSA.analyzer.plots import Charts
from GMXMMPBSA.analyzer.utils import energy2pdb_pml, ki2energy, make_corr_DF, multiindex2dict
from GMXMMPBSA.analyzer.chartsettings import ChartSettings
# from GMXMMPBSA.analyzer.propertyeditor import properties_w
import math
import parmed
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

        # five PyMOL instances for reckless
        # FIXME: now we need to iterate over all subwindows and close the PyMOL subprocess if exists
        # self.pymol_p1 = QProcess()
        # self.pymol_p2 = QProcess()
        # self.pymol_p3 = QProcess()
        # self.pymol_p4 = QProcess()
        # self.pymol_p5 = QProcess()
        # self.pymol_p_list = [self.pymol_p1, self.pymol_p2, self.pymol_p3, self.pymol_p4, self.pymol_p5]

        self.items_counter = {'charts': 0, 'pymol': [], 'bars': 0, 'line': 0, 'heatmap': 0}

        self.mdi = QMdiArea(self)
        self.mdi.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
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
        self.addDockWidget(Qt.RightDockWidgetArea, self.treeDockWidget)

        self.correlation_DockWidget = QDockWidget('Correlations', self)
        self.addDockWidget(Qt.RightDockWidgetArea, self.correlation_DockWidget)

        self.treeWidget = QTreeWidget(self)
        self.treeWidget.setMinimumWidth(380)
        self.treeWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        # self.treeWidget.customContextMenuRequested.connect(self.data_context_menu)
        self.treeWidget.itemSelectionChanged.connect(self.update_system_selection)

        self._make_options_panel()

        self.data_container = QSplitter(Qt.Vertical)
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
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.Stretch)

        self.correlation_treeWidget = QTreeWidget(self)
        self.correlation_treeWidget.itemClicked.connect(self.update_table)
        self.correlation_treeWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.correlation_treeWidget.customContextMenuRequested.connect(self.corr_context_menu)
        model_label = QTreeWidgetItem(['MODEL', 'ΔH', 'ΔH+IE', 'ΔH+NMODE', 'ΔH+QH'])
        model_label.setToolTip(0, 'Selected Model')
        model_label.setToolTip(1, 'Correlation plot for ΔG = ΔH+IE')
        model_label.setToolTip(2, 'Correlation plot for ΔG = ΔH+NMODE')
        model_label.setToolTip(3, 'Correlation plot for ΔG = ΔH+QH')
        model_label.setToolTip(4, 'Correlation plot for ΔH only')
        self.correlation_treeWidget.setHeaderItem(model_label)
        cheader = self.correlation_treeWidget.header()
        cheader.setSectionResizeMode(QHeaderView.ResizeToContents)

        self.data_table_widget = QTableWidget(0, 3)
        self.data_table_widget.setHorizontalHeaderLabels(['System', 'Exp.ΔG'])
        self.data_table_energy_col = QTableWidgetItem('None')
        self.data_table_widget.setHorizontalHeaderItem(2, self.data_table_energy_col)
        self.data_table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.data_table_widget.horizontalHeader().setStretchLastSection(True)

        self.corr_container_widget = QSplitter(Qt.Vertical)
        self.corr_container_widget.addWidget(self.correlation_treeWidget)
        self.corr_container_widget.addWidget(self.data_table_widget)
        self.correlation_DockWidget.setWidget(self.corr_container_widget)

        # self.exportpdb = ExportDialog(self)
        # self.exportcsv = ExportDialogCSV(self)
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
        # self.current_system_index = 1

        f_start = self.systems[1]['namespace'].INPUT['startframe']
        f_interval = self.systems[1]['namespace'].INPUT['interval']
        f_end = self.systems[1]['namespace'].INPUT['endframe']

        self.eframes_start_sb.setRange(f_start, f_end)
        self.eframes_start_sb.setSingleStep(f_interval)
        self.eframes_start_sb.setValue(f_start)

        self.eframes_inter_sb.setRange(1, f_end - f_start)
        # self.eframes_inter_sb.setValue(current_interval)

        self.eframes_end_sb.setRange(f_start, f_end)
        self.eframes_end_sb.setSingleStep(f_interval)
        self.eframes_end_sb.setValue(f_end)

    def _make_options_panel(self):

        self.optionWidget = QWidget(self)
        optionWidget_l = QVBoxLayout(self.optionWidget)
        optionWidget_l.setContentsMargins(0, 0, 0, 0)
        # optionWidget_l.addLayout(self.btn_l)
        hl = QFrame()
        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        optionWidget_l.addWidget(hl)

        update_frames_btn = QPushButton('Update')
        update_frames_btn.clicked.connect(self.update_frames_fn)
        resetchart_btn = QPushButton('Reset')
        btn_l = QHBoxLayout()
        btn_l.addWidget(QLabel('Options'))
        btn_l.addWidget(resetchart_btn)
        btn_l.addWidget(update_frames_btn)
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
        frames_w.setTabPosition(QTabWidget.South)
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

        c2frames_group = QWidget()
        self.c2segment_sb = QSpinBox()
        self.c2segment_sb.setRange(1, 100)
        self.c2segment_sb.setAccelerated(True)
        # self.segment_sb.valueChanged.connect(self.frames_start_sb_update)

        self.c2numframes_le = QLineEdit()
        self.c2numframes_le.setReadOnly(True)

        c2_l = QFormLayout(c2frames_group)
        # c2_l.setContentsMargins(5, 5, 5, 5)
        # c2_l.setSpacing(3)
        c2_l.addRow('Segment', self.c2segment_sb)
        c2_l.addRow('Nr. frames', self.c2numframes_le)

        frames_w.addTab(c2frames_group, 'C2')
        frames_w.setTabToolTip(3, 'Segment defined for Interaction Entropy calculation')

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
            print(self.systems[parent_item.system_index]['chart_options'])
            self.chart_options_param.restoreState(self.systems[parent_item.system_index]['chart_options_state'])
            self.chart_options_w.setParameters(self.chart_options_param, showTop=False)
        self.current_system_index = parent_item.system_index

    def _get_changes(self, sett, osett):
        # changes
        line = False
        bar = False
        heatmap = False

        for k, ko in zip(sett, osett):
            if k == 'bar_options':
                if sett[k] != osett[ko]:
                    bar = True
            elif k == 'general_options':
                if sett[k] != osett[ko]:
                    line = True
                    bar = True
                    heatmap = True
            elif k == 'line_options':
                if sett[k] != osett[ko]:
                    line = True
            elif sett[k] != osett[ko]:
                heatmap = True
        return line, bar, heatmap

    def update_frames_fn(self):

        cs = ChartSettings()
        if self.all_frb.isChecked():
            for x in self.systems:
                self.systems[x]['current_frames'] = [self.eframes_start_sb.value(),
                                                     self.eframes_end_sb.value(),
                                                     self.eframes_inter_sb.value()]
                # FIXME: update frames for entropy
        # TODO: update frames for individual charts?
        else:
            act_frange = [self.eframes_start_sb.value(), self.eframes_end_sb.value(), self.eframes_inter_sb.value()]
            act_sett = self.chart_options_param.saveState()
            curr_frange = self.systems[self.current_system_index]['current_frames']
            curr_sett = self.systems[self.current_system_index]['chart_options_state']

            rdcs = cs.get_settings(act_sett)
            rdcs_c = cs.get_settings(curr_sett)
            if act_frange == curr_frange and rdcs == rdcs_c:
                return

            self.systems[self.current_system_index]['current_frames'] = act_frange
            # chart options
            self.systems[self.current_system_index]['chart_options'] = cs.get_settings(act_sett)
            self.systems[self.current_system_index]['chart_options_state'] = act_sett
            self.systems[self.current_system_index]['changes'] = self._get_changes(rdcs, rdcs_c)

        # FIXME: create a progress bar if the task take long time

        # IMPORTANT: recalculate IE and C2
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
            item.changes(*self.systems[item.system_index]['changes'])
            itemiter += 1

        # restore to False all changes in each system
        for x in self.systems:
            self.systems[x]['changes'] = [False, False, False]

        # update current open charts. This must be made after update the variables current_frames, etc., in the system
        subwindows = self.mdi.subWindowList()
        for sub in subwindows:
            if sub.isVisible():
                sub.button.setChecked(False)
                sub.button.setChecked(True)

    def reset_dc(self):
        sub = self.mdi.activeSubWindow()
        sub.item.reset_data()
        sub.make_chart()

    def frames_start_sb_update(self, value):
        if value > self.frames_end_sb.value():
            self.frames_end_sb.setValue(value)

    def frames_end_sb_update(self, value):
        if value < self.frames_start_sb.value():
            self.frames_start_sb.setValue(value)

    def frames_inter_sb_update(self, value):
        self.frames_start_sb.setSingleStep(value)
        self.frames_end_sb.setSingleStep(value)

    def initialize(self, info_files):

        self.init_dialog.get_files_info(info_files)
        self.init_dialog.show()

    def update_options(self, item: CustomItem):
        self.frames_start_sb.setRange(item.start, item.end+item.interval)
        self.frames_start_sb.setValue(item.start)
        self.frames_inter_sb.setValue(item.interval)
        self.frames_end_sb.setValue(item.end)

        self.frames_start_sb.setSingleStep(self.frames_inter_sb.value())

    def showdata(self, item: CustomItem, col):
        self.treeWidget.clearSelection()
        # self.update_options(item)   # FIXME: only when we able the options
        if col == 1:
            s = item.lp_subw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.LINE, Charts.ROLLING],
                                                              'hide_toolbar': self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 2:
            s = item.bp_subw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:  # check if any subwindow has been store
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.BAR], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 3:
            s = item.hmp_subw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:  # check if any subwindow has been store
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.HEATMAP], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 4:
            pymol_p = item.pymol_process
            pymol_path = [os.path.join(path, 'pymol') for path in os.environ["PATH"].split(os.pathsep)
                          if os.path.exists(os.path.join(path, 'pymol')) and
                          os.access(os.path.join(path, 'pymol'), os.X_OK)]
            if not pymol_path:
                m = QMessageBox.critical(self, 'PyMOL not found!', 'PyMOL not found!. Make sure PyMOL is in the '
                                                                   'PATH.', QMessageBox.Ok)
                item.setCheckState(4, Qt.Unchecked)
                return
            else:
                pymol = pymol_path[0]

            if hasattr(item.app.FILES, 'complex_fixed'):
                com_pdb = item.syspath.parent.joinpath(item.app.FILES.complex_fixed)
            else:
                self.statusbar.showMessage(f'{item.app.FILES.prefix + "FIXED_COM.pdb"} not exits. The modified PDB file can '
                                f'be inconsistent. Please, consider use the latest version of gmx_MMPBSA')
                com_pdb = item.syspath.parent.joinpath(item.app.FILES.prefix + 'COM.pdb')
            bfactor_pml = item.syspath.parent.joinpath('bfactor.pml')
            output_path = com_pdb.parent.joinpath(f'{item.sysname}_energy2bfactor.pdb')
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                available_instance = self.get_pymol_instance()
                if not available_instance:
                    m = QMessageBox.critical(self, 'Error trying to open multiple instances of PyMOL',
                                             'Only 5 instance of PyMOL is allowed and 5 are already running. '
                                             'If you want to view this, please close the some one.',
                                             QMessageBox.Ok)
                    item.setCheckState(4, Qt.Unchecked)
                    return

                if not pymol_p or pymol_p.state() == QProcess.Running:
                    pymol_p =  available_instance # Keep a reference to the QProcess (e.g. on self) while it's running.
                    item.pymol_process = pymol_p # store pymol instance until we finish the process
                qpd = QProgressDialog('Generate modified pdb and open it in PyMOL', 'Abort', 0, 2, self)
                qpd.setWindowModality(Qt.WindowModal)
                qpd.setMinimumDuration(1500)

                for i in range(2):
                    qpd.setValue(i)
                    if qpd.wasCanceled():
                        break
                    if i == 0:
                        com_pdb_str = parmed.read_PDB(com_pdb.as_posix())
                        res_dict = item.gmxMMPBSA_current_data.bar_plot_dat.mean().to_dict()
                        for res in com_pdb_str.residues:
                            res_notation = f'{res.chain}:{res.name}:{res.number}'
                            if res_notation in res_dict:
                                res_energy = res_dict[res_notation]
                            else:
                                res_energy = 0.00
                            for at in res.atoms:
                                at.bfactor = res_energy
                        com_pdb_str.save(output_path.as_posix(), 'pdb', True, renumber=False)
                        energy2pdb_pml(res_dict, bfactor_pml, output_path)
                qpd.setValue(2)
                pymol_p.start(pymol, [bfactor_pml.as_posix()])
                pymol_p.finished.connect(lambda : item.setCheckState(4, Qt.Unchecked))
            else:
                if pymol_p and pymol_p.state() == QProcess.Running:
                    pymol_p.terminate()
                    item.pymol_process = None

    def showcorr(self, item: CorrelationItem, col):
        self.treeWidget.clearSelection()
        # self.update_options(item)   # FIXME: only when we able the options
        if col == 1:
            s = item.dh_sw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.SCATTER], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 2:
            s = item.dgie_sw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:  # check if any subwindow has been store
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.SCATTER], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 3:
            s = item.dgnmode_sw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:  # check if any subwindow has been store
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.SCATTER], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()
        elif col == 4:
            s = item.dgqh_sw
            if item.checkState(col) == Qt.Checked:
                item.setSelected(True)
                if s:  # check if any subwindow has been store
                    s.show()
                else:
                    sub = Charts(item=item, col=col, options={'chart_type':[Charts.SCATTER], 'hide_toolbar':
                        self.data_options['hide_toolbar']})
                    sub.make_chart()
                    self.mdi.addSubWindow(sub)
                    sub.show()
            else:
                if s:
                    self.mdi.activatePreviousSubWindow()
                    s.close()

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

        row = 0
        for v in data[col_label]:
            titem = QTableWidgetItem(f'{v:.2f}')
            self.data_table_widget.setItem(row, 2, titem)
            row += 1

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
        qpd.setWindowModality(Qt.WindowModal)
        qpd.setMinimumDuration(1000)

        for i in range(maximum):
            qpd.setValue(i)
            if qpd.wasCanceled():
                break
            system, api_data = rqueue.get()
            result, app = api_data
            self.makeTree(system, result, app, options)
        qpd.setLabelText(f"Processing data of {self.items_counter['charts']} items")
        qpd.setMaximum(self.items_counter['charts'])
        i = 0
        it = QTreeWidgetItemIterator(self.treeWidget)
        while it.value():
            item = it.value()
            if item.has_chart:
                qpd.setValue(i)
                item.get_data()
                i += 1
            if qpd.wasCanceled():
                break
            it += 1
        qpd.setValue(self.items_counter['charts'])

        if i != self.items_counter['charts']:
            self.close()

        if not options['correlation']: # FIXME:
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
        else:
            self.make_correlation()

        # some late signal/slot connections
        self.treeWidget.itemChanged.connect(self.showdata)
        self.correlation_treeWidget.itemChanged.connect(self.showcorr)

    def makeTree(self, system, data, app, options):

        sys_name = system[0]
        if app.FILES.stability:
            sys_name += ' (Stability)'
        idecomp = app.INPUT['idecomp']
        # set visible the pymol checkbox
        if idecomp:
            self.treeWidget.setColumnHidden(4, False)
        # make system item
        self.sys_item = CustomItem(self.treeWidget, [sys_name], system, app, False,
                                   remove_empty_terms=options['remove_empty_terms'])
        self.normalitem = CustomItem(self.sys_item, ['Normal'], has_chart=False)
        self.makeItems(data, self.normalitem, options)
        self.normalitem.setExpanded(True)

        if data.mutant:
            self.mutantitem = CustomItem(self.sys_item, ['Mutant'], has_chart=False)
            self.makeItems(data.mutant, self.mutantitem, options, True)
            self.mutantitem.setExpanded(True)
        self.sys_item.setExpanded(True)

    def makeItems(self, data, topItem, options, mutant=False):

        correlation_data = self.corr_data
        mut_pre = ''
        if mutant:
            mut_pre = 'Mut. '
            correlation_data = self.corr_data['mutant']

        parts = options['components'] + ['delta']
        if topItem.app.FILES.stability:
            parts.append('complex')
        sys_name = topItem.sysname
        correlation_data[sys_name] = {'ΔG': {
                                            'gb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
                                            'pb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
                                            'rism std': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
                                            'rism gf': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan}},
                                      'Exp.Energy': ki2energy(topItem.exp_ki, topItem.temp)}
        all_keys = {'models_keys': ['gb', 'pb', 'rism gf', 'rism std'], 'entropy_keys': ['ie', 'nmode', 'qh'],
                    'decomp_keys': ['decomp']}
        models_keys = []
        entropy_keys = []
        decomp_keys = []
        keys = []
        for key in all_keys:
            if key == 'models_keys':
                for mk in all_keys['models_keys']:
                    if mk in data:
                        models_keys.append(mk)
            if key == 'entropy_keys':
                for ek in all_keys['entropy_keys']:
                    if ek in data:
                        entropy_keys.append(ek)
            if key == 'decomp_keys':
                for dk in all_keys['decomp_keys']:
                    if dk in data:
                        decomp_keys.append(dk)

        for level in models_keys:
            item = CustomItem(topItem, [level.upper()], has_chart=False)
            for level1 in data[level]:
                # LEVEL-1 [complex, receptor, ligand] + delta
                if level1 in parts: # only show graphs for selected parts
                    item1 = CustomItem(item, [str(level1).upper()], cdata=data[level][level1], level=1,
                                       chart_title=f"Energetic Components",
                                       chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | {level1.upper()}",
                                       col_box=[2])
                    self.items_counter['charts'] += 1
                    for level2 in data[level][level1]:
                        # LEVEL-2 [GB, PB or 3D-RISM components]
                        if level1 == 'delta' and level2 == 'DELTA TOTAL':
                            correlation_data[sys_name]['ΔG'][level]['ΔH'] = data[level][level1][level2]
                        if (level2 != 'DELTA TOTAL' and options['remove_empty_charts'] and
                                abs(data[level][level1][level2].mean()) < 0.1):
                            continue
                        item2 = CustomItem(item1, [str(level2).upper()], cdata=data[level][level1][level2],
                                               level=0, chart_title=f"Energetic Components",
                                               chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | "
                                                              f"{level1.upper()} | {level2.upper()}",
                                           col_box=[1])
                        self.items_counter['charts'] += 1

        # check if any entropy approach
        if entropy_keys:
            item = CustomItem(topItem, ['Entropy'], has_chart=False)
            itemd = CustomItem(topItem, ['ΔG Binding'], has_chart=False)

            for level in entropy_keys:
                if level == 'ie':
                    item1 = CustomItem(item, [str(level).upper()], cdata=data[level]['data'],
                                       level=0, chart_title=f"Interaction Entropy (-TΔS)",
                                       chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()}",
                                       col_box=[1])
                    # TODO: improve ieframes calculation
                    ieframes = int((data[level]['frames'][1] - data[level]['frames'][0])/item.interval) + 1
                    ent = data[level]['data'][-ieframes:]
                    item1.ie = [data[level]['frames'], data[level]['value']]
                    self.items_counter['charts'] += 1
                else:
                    iteme = CustomItem(item, [str(level).upper()], has_chart=False)
                    for level1 in data[level]:
                        # LEVEL-1 [complex, receptor, ligand] + delta
                        if level1 in parts:  # only show graphs for selected parts
                            ent_dict = {datum: data[level][level1][datum] * -1 for datum in data[level][level1]}
                            item1 = CustomItem(iteme, [str(level1).upper()], cdata=ent_dict,
                                               level=1, chart_title=f"Entropy {str(level).upper()} approximation ("
                                                                    f"-TΔS)",
                                               chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | "
                                                              f"{level1.upper()}",
                                               col_box=[2])
                            # for level2 in data[level][level1]:
                            #     item2 = CustomItem(item1, [f"{str(level).upper()} {str(level2).upper()}"],
                            #                        cdata=data[level][level1][level2],
                            #                        level=0.1, chart_title=f"{str(level).upper()} Entropy",
                            #                        chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | "
                            #                                       f"{level1.upper()} | {level2.upper()}",
                            #                        col_box=[1])
                            self.items_counter['charts'] += 1
                    ent = data[level]['delta']['Total'] * -1 # Since nmode and qh calculate the TΔS term


                for model in correlation_data[sys_name]['ΔG']:
                    if correlation_data[sys_name]['ΔG'][model]['ΔH'] is np.nan:
                        continue
                    itemd1 = CustomItem(itemd, [str(model).upper()], has_chart=False)

                    if len(correlation_data[sys_name]['ΔG'][model]['ΔH']) != len(ent):
                        dg = correlation_data[sys_name]['ΔG'][model]['ΔH'] + ent.mean()
                        if len(correlation_data[sys_name]['ΔG'][model]['ΔH']) > len(ent):
                            for t in range(abs(len(correlation_data[sys_name]['ΔG'][model]['ΔH']) - len(ent))):
                                ent = np.append(ent, np.nan)
                        # else:
                        #     dh = correlation_data[sys_name]['ΔG'][model]['ΔH']
                        #     for t in range(abs(len(correlation_data[sys_name]['ΔG'][model]['ΔH']) - len(ent))):
                        #         dh = np.append(dh, np.nan)
                    else:
                        dg = correlation_data[sys_name]['ΔG'][model]['ΔH'] + np.nanmean(ent)
                    itemd2 = CustomItem(itemd1, [f'ΔG Binding ({str(level).upper()})'],
                                        cdata={'ΔH': correlation_data[sys_name]['ΔG'][model]['ΔH'], '-TΔS': ent,
                                               'ΔG': dg},
                                        level=1.1, chart_title=f"ΔG Binding",
                                        chart_subtitle=f"{mut_pre}{sys_name} | {str(model).upper()} | {level.upper()}",
                                        col_box=[2])
                    self.items_counter['charts'] += 1
                    correlation_data[sys_name]['ΔG'][model][level] = dg

        for level in decomp_keys:
            # omit decomp data
            if not options['decomposition']:
                continue
            item = CustomItem(topItem, [str(level).upper()], has_chart=False)
            for level1 in data[level]:
                # LEVEL-1 [GB or PB]
                item1 = CustomItem(item, [str(level1).upper()], has_chart=False)
                for level2 in data[level][level1]:
                    # LEVEL-2 [complex, receptor, ligand, delta]
                    if level2 in parts:
                        item2 = CustomItem(item1, [str(level2).upper()], has_chart=False)
                        for level3 in data[level][level1][level2]:
                            #  LEVEL-3 [TDC, SDC, BDC]
                            if topItem.idecomp in [1, 2]:
                                # Per-residue
                                item_level = 2
                            else:
                                # Per-wise
                                item_level = 3
                            col_box = [1, 2, 3]
                            # only make a checkbox for TDC
                            if level3 == 'TDC' and level2 == 'delta':
                                col_box.append(4)
                            item3 = CustomItem(item2, [str(level3).upper()],
                                               cdata=data[level][level1][level2][level3],
                                               level=item_level, chart_title=f"Energetic Components [Per-residue]",
                                               chart_subtitle=f"{mut_pre}{sys_name} | {level.upper()} | "
                                                              f"{level1.upper()} | {level2.upper()} | "
                                                              f"{level3.upper()}",
                                               col_box=col_box)
                            self.items_counter['charts'] += 1
                            if 4 in col_box:
                                self.items_counter['pymol'].append(item3)
                            for level4 in data[level][level1][level2][level3]:
                                # LEVEL-4 Selected residues
                                col_box = [2]
                                if topItem.idecomp in [1, 2]:
                                    # Per-residue
                                    item_level = 1
                                else:
                                    # Per-wise
                                    item_level = 2
                                    col_box.extend([1, 3])
                                item4 = CustomItem(item3, [str(level4).upper()],
                                                   cdata=data[level][level1][level2][level3][level4],
                                                   level=item_level,
                                                   chart_title=f"Energetic Components [Per-residue]",
                                                   chart_subtitle=f"{mut_pre}{sys_name} | {str(level).upper()} | "
                                                                  f"{str(level1).upper()} | {str(level2).upper()} | "
                                                                  f"{str(level3).upper()} | {str(level4).upper()}",
                                                   col_box=col_box)
                                self.items_counter['charts'] += 1

                                for level5 in data[level][level1][level2][level3][level4]:
                                    # LEVEL-5 AA energetic terms if per-residue or per-wise residues
                                    if topItem.idecomp in [1, 2]:
                                        item5 = CustomItem(item4, [str(level5).upper()],
                                                           cdata=data[level][level1][level2][level3][level4][
                                                               level5],
                                                           level=0,
                                                           chart_title=f"Energetic Components [Per-wise]",
                                                           chart_subtitle=f"{mut_pre}{sys_name} | "
                                                                          f"{str(level).upper()} | "
                                                                          f"{str(level1).upper()} | "
                                                                          f"{str(level2).upper()} | "
                                                                          f"{str(level3).upper()} | "
                                                                          f"{str(level4).upper()} | "
                                                                          f"{str(level5).upper()}",
                                                           col_box=[1])
                                        self.items_counter['charts'] += 1
                                    else:
                                        item5 = CustomItem(item4, [str(level5).upper()],
                                                           cdata=data[level][level1][level2][level3][level4][level5],
                                                           level=1,
                                                           chart_title=f"Energetic Components [Per-wise]",
                                                           chart_subtitle=f"{mut_pre}{sys_name} | "
                                                                          f"{str(level).upper()} | "
                                                                          f"{str(level1).upper()} | "
                                                                          f"{str(level2).upper()} | "
                                                                          f"{str(level3).upper()} | "
                                                                          f"{str(level4).upper()} | "
                                                                          f"{str(level5).upper()}",
                                                           col_box=[2])
                                        self.items_counter['charts'] += 1

                                        for level6 in data[level][level1][level2][level3][level4][level5]:
                                            item6 = CustomItem(item5, [str(level6).upper()],
                                                               cdata=data[level][level1][level2][level3][level4][
                                                                       level5][level6],
                                                               level=0,
                                                               chart_title=f"Energetic Components [Per-wise]",
                                                               chart_subtitle=f"{mut_pre}{sys_name} | "
                                                                              f"{str(level).upper()} | "
                                                                              f"{str(level1).upper()} | "
                                                                              f"{str(level2).upper()} | "
                                                                              f"{str(level3).upper()} | "
                                                                              f"{str(level4).upper()} | "
                                                                              f"{str(level5).upper()} | "
                                                                              f"{str(level6).upper()}",
                                                               col_box=[1])
                                            self.items_counter['charts'] += 1

