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
try:
    from PyQt6.QtWidgets import *
    from PyQt6.QtCore import *
    from PyQt6.QtGui import *
except:
    from PyQt5.QtWidgets import *
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *


import sys
import os

from queue import Queue, Empty
from pathlib import Path

import pandas as pd
from GMXMMPBSA.API import load_gmxmmpbsa_info, MMPBSA_API
from GMXMMPBSA.analyzer.dialogs import InitDialog
from GMXMMPBSA.analyzer.customitem import CustomItem, CorrelationItem
from GMXMMPBSA.analyzer.style import save_default_config, default_config, save_user_config, user_config, toc_img, logo, \
    alert, config
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
    def __init__(self, ifiles):
        super(GMX_MMPBSA_ANA, self).__init__()
        self.showMaximized()
        self.setWindowIcon(QIcon(logo))
        self.corr_data = {'mutant': {}}

        self.ifiles = ifiles

        self.systems = {}
        self.all_systems_active = True
        self.current_system_index = None
        self.pymol_p_list = []
        self.removed_items = {}

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
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

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
        self.aboutMenu.addAction('Documentation', self._doc)
        self.aboutMenu.addAction('Report a bug', self._bug)
        self.aboutMenu.addAction('Google group', self._group)
        self.aboutMenu.addAction('About gmx_MMPBSA', self._about_dialog)
        self.statusbar = self.statusBar()

        self.init_dialog = InitDialog(self)
        self.init_dialog.rejected.connect(self.close)
        self.gettting_data()

    def closeEvent(self, a0: QCloseEvent) -> None:

        changes = [
            x
            for x in self.systems
            if self.systems[x]['chart_options'].is_default_changed()
        ]

        savequitbtn = None
        quitbtn = None
        if changes:
            msgBox = QMessageBox()
            msgBox.setWindowTitle('Quit...')
            msgBox.setText("Are you sure to quit?")
            msgBox.setInformativeText("Some graphics settings were modified. Do you want to save your changes?")
            msgBox.setDetailedText("The systems below have chances:\n - " +
                                   '\n - '.join([self.systems[x]['name'] for x in changes]))
            savequitbtn = msgBox.addButton("Save and Quit", QMessageBox.ButtonRole.AcceptRole)
            abortbtn = msgBox.addButton(QMessageBox.StandardButton.Cancel)
            quitbtn = msgBox.addButton("Quit", QMessageBox.ButtonRole.AcceptRole)
            msgBox.setDefaultButton(savequitbtn)
            msgBox.exec()
            action = msgBox.clickedButton()

        else:
            action = QMessageBox.question(self.mdi, 'Message', "Are you sure to quit?", QMessageBox.StandardButton.Yes,
                                          QMessageBox.StandardButton.No)

        if action in [QMessageBox.StandardButton.Yes, quitbtn, savequitbtn]:

            if action == savequitbtn:
                for x in changes:
                    self.systems[x]['chart_options'].write_system_config(self.systems[x]['path'])

            if pymol_items := [
                p for p, _ in self.pymol_p_list if p.state() == QProcess.ProcessState.Running
            ]:
                qpd = QProgressDialog('Closing PyMOL instances', 'Abort', 0, len(pymol_items), self)
                qpd.setWindowModality(Qt.WindowModality.WindowModal)
                qpd.setMinimumDuration(1000)
                qpd.setRange(0, len(self.pymol_p_list))
                for i, (p, _) in enumerate(self.pymol_p_list):
                    if p.state() == QProcess.ProcessState.Running:
                        qpd.setValue(i)
                        p.kill()
                        p.waitForFinished()
                qpd.setValue(len(self.pymol_p_list))
            a0.accept()
        else:
            a0.ignore()
            if self.in_init_dialog:
                self.init_dialog.show()

    def _about_dialog(self):
        from GMXMMPBSA import __version__
        QMessageBox.about(self.mdi, "About gmx_MMPBSA",
                              "<html>"
                              "<body>"
                              "<h2 style='text-align:center'>About gmx_MMPBSA</h2>"
                              "<div text-align:center;>"
                              f"<a href='https://valdes-tresanco-ms.github.io/gmx_MMPBSA/'>"
                              f"<img src={toc_img} alt='gmx_MMPBSA TOC' width='450' height='350'>"
                              "</a>"
                              "</div>"
                              "<p style='text-align:center'><b>gmx_MMPBSA</b> is a new tool based on AMBER's MMPBSA.py "
                              "aiming to perform end-state free energy calculations with GROMACS files.</p>"
                              f"<p style='text-align:center'><b style='text-align:center'>Version:</b> {__version__}</p>"
                              "<h2 style='text-align:center'>Cite gmx_MMPBSA</h2>"
                              "<p style='text-align:center'> "
                              "Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. "
                              "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. "
                              "Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291. "
                              "<a href='https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645'>"
                              "https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645</a> </p>"
                              "</body>"
                              "</html>")

    def _help(self):
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/analyzer/'))

    def _doc(self):
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/getting-started/'))

    def _bug(self):
        QDesktopServices().openUrl(QUrl('https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/new/choose'))

    def _group(self):
        QDesktopServices().openUrl(QUrl('https://groups.google.com/g/gmx_mmpbsa'))

    def _initialize_systems(self):

        if len(self.systems):
            self.all_systems_active = False
        if self.all_systems_active:
            self.all_frb.setEnabled(True)

        # Select automatically the first item to update the option panel
        topitem = self.treeWidget.topLevelItem(0)
        topitem.setSelected(True)

        # energy config
        if (
                self.systems[1]['namespace'].INPUT['gbrun'] or
                self.systems[1]['namespace'].INPUT['pbrun'] or
                self.systems[1]['namespace'].INPUT['rismrun']
        ):
            self.eframes_group.setEnabled(True)
            f_start = self.systems[1]['namespace'].INPUT['startframe']
            f_interval = self.systems[1]['namespace'].INPUT['interval']
            f_end = self.systems[1]['namespace'].INPUT['endframe']

            self.eframes_start_sb.setRange(f_start, f_end)
            self.eframes_start_sb.setSingleStep(f_interval)
            self.eframes_start_sb.setValue(f_start)

            self.eframes_inter_sb.setRange(1, f_end - f_start)
            self.eframes_inter_sb.setSingleStep(f_interval)
            self.eframes_inter_sb.setValue(f_interval)

            self.eframes_end_sb.setRange(f_start, f_end)
            self.eframes_end_sb.setSingleStep(f_interval)
            self.eframes_end_sb.setValue(f_end)
        else:
            self.eframes_group.setEnabled(False)

        # nmode config
        if self.systems[1]['namespace'].INPUT['nmoderun']:
            self.nmframes_group.setEnabled(True)
            nmf_start = self.systems[1]['namespace'].INPUT['nmstartframe']
            nmf_interval = self.systems[1]['namespace'].INPUT['nminterval']
            nmf_end = self.systems[1]['namespace'].INPUT['nmendframe']

            self.nmframes_start_sb.setRange(nmf_start, nmf_end)
            self.nmframes_start_sb.setSingleStep(nmf_interval)
            self.nmframes_start_sb.setValue(nmf_start)

            self.nmframes_inter_sb.setRange(1, nmf_end - nmf_start)
            self.nmframes_inter_sb.setSingleStep(nmf_interval)
            self.nmframes_inter_sb.setValue(nmf_interval)

            self.nmframes_end_sb.setRange(nmf_start, nmf_end)
            self.nmframes_end_sb.setSingleStep(nmf_interval)
            self.nmframes_end_sb.setValue(nmf_end)
        else:
            self.nmframes_group.setEnabled(False)

        # ie config
        if self.systems[1]['namespace'].INPUT['interaction_entropy']:
            self.ieframes_group.setEnabled(True)
            self.iesegment_sb.setValue(self.systems[1]['namespace'].INPUT['ie_segment'])
        else:
            self.ieframes_group.setEnabled(False)

    def _set_as_default(self):
        self.systems[self.current_system_index]['chart_options'].write_system_config()
        self.statusbar.showMessage(f"Setting this setting as default in "
                                   f"{self.systems[self.current_system_index]['chart_options'].filename.as_posix()}...")
        self.systems[self.current_system_index]['chart_options'].set_as_default()
        self.update_fn()
        self.chart_options_param.restoreState(self.systems[self.current_system_index]['chart_options'])
        self.parm_tree_w.setParameters(self.chart_options_param, showTop=False)

    def _reset2default(self):
        new_settings = ChartSettings()
        self.statusbar.showMessage("Restore default settings...")
        self.chart_options_param.restoreState(new_settings)
        self.parm_tree_w.setParameters(self.chart_options_param, showTop=False)
        self.systems[self.current_system_index]['chart_options'] = new_settings
        self.update_fn(repaint_arg=True)

    def _set_user_config(self):
        p = self.systems[self.current_system_index]['path']
        fn = p.joinpath('setting.json')
        if fn.exists():
            user_settings = ChartSettings(fn)
            self.statusbar.showMessage(f"Setting user configuration for this system from {fn.as_posix()}...")
            self.update_fn()
            self.chart_options_param.restoreState(self.systems[self.current_system_index]['chart_options'])
            self.parm_tree_w.setParameters(self.chart_options_param, showTop=False)
            self.systems[self.current_system_index]['chart_options'] = user_settings
        else:
            self.statusbar.showMessage('No setting file was found for this system...')

    def _save_user_config(self):
        p = self.systems[self.current_system_index]['path']
        fn = p.joinpath('settings.json')
        self.systems[self.current_system_index]['chart_options'].write_system_config(p)
        self.statusbar.showMessage(f"Saving user config for this system {fn.as_posix()}...")

    def _make_options_panel(self):

        self.optionWidget = QWidget(self)
        self.optionWidget.setMinimumWidth(300)
        optionWidget_l = QVBoxLayout(self.optionWidget)

        selection_group = QGroupBox('Selection')
        selection_group_layout = QHBoxLayout(selection_group)
        selection_group_layout.setContentsMargins(10, 0, 10, 0)
        self.curr_sys_frb = QRadioButton('Selected System')
        self.curr_sys_frb.setChecked(True)
        selection_group_layout.addWidget(self.curr_sys_frb)
        self.all_frb = QRadioButton('All Systems')
        self.all_frb.setEnabled(False)
        selection_group_layout.addWidget(self.all_frb)

        optionWidget_l.addWidget(selection_group)

        update_btn = QPushButton('Update')
        update_btn.clicked.connect(self.update_fn)

        optionWidget_c = QTabWidget(self)
        optionWidget_l.addWidget(optionWidget_c)
        optionWidget_c.setCornerWidget(update_btn)

        # Charts options
        self.chart_options_param = Parameter().create()
        self.chart_options_w = QWidget(self)
        self.chart_options_l = QVBoxLayout(self.chart_options_w)
        self.chart_options_l.setContentsMargins(0, 0, 0, 0)
        charts_opt_tb = QToolButton()
        charts_opt_tb.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        charts_opt_tb.setIcon(QIcon(config))
        self.chart_options_l.addWidget(charts_opt_tb, alignment=Qt.AlignmentFlag.AlignRight)
        charts_opt_menu = QMenu()
        charts_opt_menu.setToolTipsVisible(True)

        self.parm_tree_w = ParameterTree()
        self.chart_options_l.addWidget(self.parm_tree_w)

        self.set_as_default_action = charts_opt_menu.addAction(QIcon(save_default_config), 'Set as default',
                                                             self._set_as_default)
        self.set_as_default_action.setToolTip('Save the current setting as the global default settings.')
        self.default_action = charts_opt_menu.addAction(QIcon(default_config), 'Reset to default', self._reset2default)
        self.default_action.setToolTip('Restores the default settings set by the developers as the settings '
                                       'for the selected systems.')
        charts_opt_menu.addSeparator()
        self.set_user_config_action = charts_opt_menu.addAction(QIcon(save_user_config), 'Save User-config',
                                                              self._set_as_default)
        self.set_user_config_action.setToolTip('Force save the current setting for this system. Before gmx_MMPBSA_ana '
                                               'closes, if there are configuration changes, you can decide if you '
                                               'want to save them.')

        self.use_user_config_action = charts_opt_menu.addAction(QIcon(user_config), 'User-config', self._set_as_default)
        self.use_user_config_action.setToolTip('Uses the specific settings for this system if it was saved.')
        charts_opt_tb.setMenu(charts_opt_menu)

        optionWidget_c.addTab(self.chart_options_w, 'Charts Options')

        frames_w = QWidget(optionWidget_c)
        optionWidget_c.addTab(frames_w, 'Frames')
        frames_l = QVBoxLayout(frames_w)

        self.eframes_group = QGroupBox('Energy')
        frames_l.addWidget(self.eframes_group)
        self.eframes_start_sb = QSpinBox()
        self.eframes_start_sb.setAccelerated(True)
        self.eframes_start_sb.valueChanged.connect(self.frames_start_sb_update)
        self.eframes_inter_sb = QSpinBox()
        self.eframes_inter_sb.valueChanged.connect(self.frames_inter_sb_update)
        self.eframes_end_sb = QSpinBox()
        self.eframes_end_sb.setAccelerated(True)
        self.eframes_end_sb.valueChanged.connect(self.frames_end_sb_update)
        self.numframes_le = QLineEdit()
        self.numframes_le.setReadOnly(True)

        self.e_changed = QLabel()
        self.e_changed.setToolTip('The frames range has changed. Please, press the "Update" button to recalculate the '
                                   'energy or reset to the default frame range')
        self.e_changed.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.e_changed.setPixmap(QPixmap(alert))
        self.e_changed.hide()

        fg_l1 = QFormLayout()
        fg_l1.addRow('Start', self.eframes_start_sb)
        fg_l1.addRow('End', self.eframes_end_sb)
        fg_l2 = QFormLayout()
        fg_l2.addRow('Interval', self.eframes_inter_sb)
        fg_l2.addRow('Nr. frames', self.numframes_le)

        reset_energy_btn = QPushButton('Reset')
        reset_energy_btn.clicked.connect(self._reset_e_frames)

        frames_group_l = QGridLayout(self.eframes_group)
        frames_group_l.addLayout(fg_l1, 0, 0)
        frames_group_l.addLayout(fg_l2, 0, 2)
        frames_group_l.addWidget(self.e_changed, 1, 0, 1, 2, Qt.AlignmentFlag.AlignRight)
        frames_group_l.addWidget(reset_energy_btn, 1, 2)
        frames_group_l.setColumnStretch(0, 5)
        frames_group_l.setColumnStretch(1, 1)
        frames_group_l.setColumnStretch(2, 5)

        frames_l.addWidget(self.eframes_group)

        self.nmframes_group = QGroupBox('NMODE')
        self.nmframes_start_sb = QSpinBox()
        self.nmframes_start_sb.setAccelerated(True)
        self.nmframes_start_sb.valueChanged.connect(self.nmframes_start_sb_update)
        self.nmframes_inter_sb = QSpinBox()
        self.nmframes_inter_sb.valueChanged.connect(self.nmframes_inter_sb_update)
        self.nmframes_end_sb = QSpinBox()
        self.nmframes_end_sb.setAccelerated(True)
        self.nmframes_end_sb.valueChanged.connect(self.nmframes_end_sb_update)
        self.nmnumframes_le = QLineEdit()
        self.nmnumframes_le.setReadOnly(True)

        self.nm_changed = QLabel()
        self.nm_changed.setToolTip('The frames range has changed. Please, press the "Update" button to recalculate the '
                                   'energy or reset to the default frame range')
        self.nm_changed.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.nm_changed.setPixmap(QPixmap(alert))
        self.nm_changed.hide()

        nmfg_l1 = QFormLayout()
        nmfg_l1.addRow('Start', self.nmframes_start_sb)
        nmfg_l1.addRow('End', self.nmframes_end_sb)
        nmfg_l2 = QFormLayout()
        nmfg_l2.addRow('Interval', self.nmframes_inter_sb)
        nmfg_l2.addRow('Nr. frames', self.nmnumframes_le)

        reset_nmode_btn = QPushButton('Reset')
        reset_nmode_btn.clicked.connect(self._reset_nm_frames)

        nmframes_group_l = QGridLayout(self.nmframes_group)
        nmframes_group_l.addLayout(nmfg_l1, 0, 0)
        nmframes_group_l.addLayout(nmfg_l2, 0, 2)
        nmframes_group_l.addWidget(self.nm_changed, 1, 0, 1, 2, Qt.AlignmentFlag.AlignRight)
        nmframes_group_l.addWidget(reset_nmode_btn, 1, 2)
        nmframes_group_l.setColumnStretch(0, 5)
        nmframes_group_l.setColumnStretch(1, 1)
        nmframes_group_l.setColumnStretch(2, 5)

        frames_l.addWidget(self.nmframes_group)

        self.ieframes_group = QGroupBox('IE')
        self.ie_changed = QLabel()
        self.ie_changed.setToolTip('The frames range has changed. Please, press the "Update" button to recalculate the '
                                   'energy or reset to the default frame range')
        self.ie_changed.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.ie_changed.setPixmap(QPixmap(alert))
        self.ie_changed.hide()
        self.iesegment_sb = QSpinBox()
        self.iesegment_sb.setRange(1, 100)
        self.iesegment_sb.setAccelerated(True)
        self.iesegment_sb.setSuffix('%')
        self.iesegment_sb.valueChanged.connect(self.iesegment_sb_update)

        self.ienumframes_le = QLineEdit()
        self.ienumframes_le.setReadOnly(True)

        ieseg_l = QFormLayout()
        ieseg_l.addRow('Segment', self.iesegment_sb)
        ienf_l = QFormLayout()
        ienf_l.addRow('Nr. frames', self.ienumframes_le)

        ie_reset_btn = QPushButton('Reset')
        ie_reset_btn.clicked.connect(self._reset_iesegment)

        ieframes_group_l = QGridLayout(self.ieframes_group)
        ieframes_group_l.addLayout(ieseg_l, 0, 0)
        ieframes_group_l.addLayout(ienf_l, 0, 2)
        ieframes_group_l.addWidget(self.ie_changed, 1, 0, 1, 2, Qt.AlignmentFlag.AlignRight)
        ieframes_group_l.addWidget(ie_reset_btn, 1, 2)
        ieframes_group_l.setColumnStretch(0, 5)
        ieframes_group_l.setColumnStretch(1, 1)
        ieframes_group_l.setColumnStretch(2, 5)
        frames_l.addWidget(self.ieframes_group)
        frames_l.addStretch(1)

    def _reset_iesegment(self):
        self.iesegment_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['ie_segment'])

    def _reset_e_frames(self):
        self.eframes_start_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['startframe'])
        self.eframes_end_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['endframe'])
        self.eframes_inter_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['interval'])

    def _reset_nm_frames(self):
        self.nmframes_start_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['nmstartframe'])
        self.nmframes_end_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['nmendframe'])
        self.nmframes_inter_sb.setValue(self.systems[self.current_system_index]['namespace'].INPUT['nminterval'])

    def iesegment_sb_update(self):
        current_estart = self.eframes_start_sb.value()
        current_einterval = self.eframes_inter_sb.value()
        current_eend = self.eframes_end_sb.value()
        ieframes = math.ceil(
            (((current_eend - current_estart) // current_einterval) + 1) * (self.iesegment_sb.value() / 100)
        )
        self.ienumframes_le.setText(str(ieframes))
        if self.current_system_index:
            curr_eframes = [self.eframes_start_sb.value(), self.eframes_end_sb.value(), self.eframes_inter_sb.value()]
            if (
                    self.iesegment_sb.value() != self.systems[self.current_system_index]['current_ie_segment'] or
                    curr_eframes != self.systems[self.current_system_index]['current_frames']
            ):
                self.ie_changed.show()
            else:
                self.ie_changed.hide()

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
        # FIXME: check if all systems are selected?
        # if not self.all_frb.isChecked():
        # energy
        current_start, current_end, current_interval = self.systems[parent_item.system_index]['current_frames']
        self.eframes_start_sb.setValue(current_start)
        self.eframes_inter_sb.setValue(current_interval)
        self.eframes_end_sb.setValue(current_end)
        # nmode
        nmcurrent_start, nmcurrent_end, nmcurrent_interval = self.systems[parent_item.system_index][
            'current_nmode_frames']
        self.nmframes_start_sb.setValue(nmcurrent_start)
        self.nmframes_inter_sb.setValue(nmcurrent_end)
        self.nmframes_end_sb.setValue(nmcurrent_interval)
        # IE
        self.iesegment_sb.setValue(self.systems[parent_item.system_index]['current_ie_segment'])
        # self.systems[parent_item.system_index]['chart_options'] = self.chart_options_param.saveState()
        self.chart_options_param.restoreState(self.systems[parent_item.system_index]['chart_options'])
        self.parm_tree_w.setParameters(self.chart_options_param, showTop=False)
        self.current_system_index = parent_item.system_index

    def update_fn(self, repaint_arg=False):
        self.statusbar.showMessage('Updating...')
        recalc_energy = []
        recalc_nmode = []
        recalc_ie = []
        repaint = []

        curr_eframes = [self.eframes_start_sb.value(), self.eframes_end_sb.value(), self.eframes_inter_sb.value()]
        curr_nmframes = [self.nmframes_start_sb.value(), self.nmframes_end_sb.value(), self.nmframes_inter_sb.value()]
        curr_iesegment = self.iesegment_sb.value()
        act_sett = self.chart_options_param.saveState()

        processed_sys = []

        if self.all_frb.isChecked():
            for x in self.systems:
                processed_sys.append(x)
                energyrun = (self.systems[x]['namespace'].INPUT['gbrun'] or
                             self.systems[x]['namespace'].INPUT['pbrun'] or
                             self.systems[x]['namespace'].INPUT['rismrun'])
                if energyrun and curr_eframes != self.systems[x]['current_frames']:
                    recalc_energy.append(x)
                if (self.systems[x]['namespace'].INPUT['nmoderun'] and
                        curr_nmframes != self.systems[x]['current_nmode_frames']):
                    recalc_nmode.append(x)
                if curr_iesegment != self.systems[x]['current_ie_segment']:
                    recalc_ie.append(x)
                if self.systems[x]['chart_options'].is_changed(act_sett) or repaint_arg:
                    repaint.append(x)
        else:
            processed_sys.append(self.current_system_index)
            energyrun = (self.systems[self.current_system_index]['namespace'].INPUT['gbrun'] or
                         self.systems[self.current_system_index]['namespace'].INPUT['pbrun'] or
                         self.systems[self.current_system_index]['namespace'].INPUT['rismrun'])
            if energyrun and curr_eframes != self.systems[self.current_system_index]['current_frames']:
                recalc_energy.append(self.current_system_index)
            if (self.systems[self.current_system_index]['namespace'].INPUT['nmoderun'] and
                    curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']):
                recalc_nmode.append(self.current_system_index)
            if curr_iesegment != self.systems[self.current_system_index]['current_ie_segment']:
                recalc_ie.append(self.current_system_index)
            if self.systems[self.current_system_index]['chart_options'].is_changed(act_sett) or repaint_arg:
                repaint.append(self.current_system_index)

        maximum = len(recalc_energy) + len(recalc_ie) + len(recalc_nmode) + len(repaint)

        if not maximum:
            return
        maximum += 1  # close and re-open current active windows

        qpd = QProgressDialog('Creating systems tree', 'Abort', 0, maximum, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(0)
        v = 0
        if recalc_energy:
            qpd.setLabelText('Recalculating energies for selected frames range')
            for e in recalc_energy:
                v += 1
                self.systems[e]['api'].update_energy(*curr_eframes)
                self.systems[e]['current_frames'] = curr_eframes
                qpd.setValue(v)
        if recalc_nmode:
            qpd.setLabelText('Recalculating nmode for selected frames range')
            for e in recalc_nmode:
                v += 1
                self.systems[e]['api'].update_nmode(*curr_nmframes)
                self.systems[e]['current_nmode_frames'] = curr_nmframes
                qpd.setValue(v)
        if recalc_ie:
            qpd.setLabelText('Recalculating Interaction Entropy for selected frames range')
            for e in recalc_ie:
                v += 1
                self.systems[e]['api'].update_ie(curr_iesegment)
                self.systems[e]['current_ie_segment'] = curr_iesegment
                qpd.setValue(v)
        if recalc_energy or recalc_nmode or recalc_ie:
            slist = recalc_energy or recalc_nmode or recalc_ie
            for e in slist:
                self.systems[e]['data'] = self.systems[e]['api'].get_energy()
        if repaint:
            qpd.setLabelText('Updating charts options')
            for e in repaint:
                v += 1
                if repaint_arg:
                    self.systems[e]['chart_options'].changes = dict(
                        line_action=3,
                        line_ie_action=3,
                        bar_action=3,
                        heatmap_action=3,
                        visualization_action=3,
                        figure=3
                    )
                else:
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
                elif changes['figure']:
                    chart_sett = self.systems[s]['chart_options'].get_settings()
                    options = {'save-format': chart_sett[('General', 'figure-format', 'save-format')],
                               'dpi-save': chart_sett[('General', 'figure-format', 'dpi-save')]
                               }
                    sub.mpl_toolbar.update_options(options)
                    sub.fbtn.setChecked(chart_sett[('General', 'toolbar')])
            pymol_items = [[p, item] for p, item in self.pymol_p_list if p.state() == QProcess.ProcessState.Running]
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
                                                            visualization_action=0,
                                                            figure=0)
        self.ie_changed.hide()
        self.e_changed.hide()
        self.nm_changed.hide()
        qpd.setValue(maximum)
        self.statusbar.showMessage('Updating... Done.')


    def reset_dc(self):
        sub = self.mdi.activeSubWindow()
        sub.item.reset_data()
        sub.make_chart()

    def frames_start_sb_update(self, value):
        if value > self.eframes_end_sb.value():
            self.eframes_end_sb.setValue(value)

        if self.current_system_index:
            current_interval = self.eframes_inter_sb.value()
            current_end = self.eframes_end_sb.value()
            self.numframes_le.setText(f"{int((current_end - value) // current_interval) + 1}")

            curr_eframes = [value, self.eframes_end_sb.value(), self.eframes_inter_sb.value()]
            if curr_eframes != self.systems[self.current_system_index]['current_frames']:
                self.e_changed.show()
            else:
                self.e_changed.hide()
            self.iesegment_sb_update()

    def frames_end_sb_update(self, value):
        if value < self.eframes_start_sb.value():
            self.eframes_start_sb.setValue(value)
        if self.current_system_index:
            current_start = self.eframes_start_sb.value()
            current_interval = self.eframes_inter_sb.value()
            self.numframes_le.setText(f"{int((value - current_start) // current_interval) + 1}")

            curr_eframes = [self.eframes_start_sb.value(), value, self.eframes_inter_sb.value()]
            if curr_eframes != self.systems[self.current_system_index]['current_frames']:
                self.e_changed.show()
            else:
                self.e_changed.hide()
            self.iesegment_sb_update()

    def frames_inter_sb_update(self, value):
        if self.current_system_index:
            current_start = self.eframes_start_sb.value()
            current_end = self.eframes_end_sb.value()
            self.numframes_le.setText(f"{int((current_end - current_start) // value) + 1}")

            curr_eframes = [self.eframes_start_sb.value(), self.eframes_end_sb.value(), value]
            if curr_eframes != self.systems[self.current_system_index]['current_frames']:
                self.e_changed.show()
            else:
                self.e_changed.hide()
            self.iesegment_sb_update()

    def nmframes_start_sb_update(self, value):
        if value > self.nmframes_end_sb.value():
            self.nmframes_end_sb.setValue(value)

        if self.current_system_index:
            current_interval = self.nmframes_inter_sb.value()
            current_end = self.nmframes_end_sb.value()
            self.numframes_le.setText(f"{int((current_end - value) // current_interval) + 1}")

            curr_nmframes = [value, self.nmframes_end_sb.value(), self.nmframes_inter_sb.value()]
            if curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']:
                self.nm_changed.show()
            else:
                self.nm_changed.hide()

    def nmframes_end_sb_update(self, value):
        if value < self.nmframes_start_sb.value():
            self.nmframes_start_sb.setValue(value)
        if self.current_system_index:
            current_start = self.nmframes_start_sb.value()
            current_interval = self.nmframes_inter_sb.value()
            self.numframes_le.setText(f"{int((value - current_start) // current_interval) + 1}")

            curr_nmframes = [self.nmframes_start_sb.value(), value, self.nmframes_inter_sb.value()]
            if curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']:
                self.nm_changed.show()
            else:
                self.nm_changed.hide()

    def nmframes_inter_sb_update(self, value):
        if self.current_system_index:
            current_start = self.nmframes_start_sb.value()
            current_end = self.nmframes_end_sb.value()
            self.numframes_le.setText(f"{int((current_end - current_start) // value) + 1}")

            curr_nmframes = [self.nmframes_start_sb.value(), self.nmframes_end_sb.value(), value]
            if curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']:
                self.nm_changed.show()
            else:
                self.nm_changed.hide()

    def gettting_data(self):
        self.in_init_dialog = True
        self.init_dialog.get_files_info(self.ifiles)
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
            QMessageBox.critical(self.mdi, 'Unable to calculate correlation',
                                     'Three or more systems are needed to calculate the correlation.',
                                     QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
            return

        df = make_corr_DF(self.corr_data)  # FIXME:
        # get df for each model
        models = ['gb', 'pb', 'rism std', 'rism gf', 'rism pcplus']
        columns = ['ΔH', 'ΔH+IE', 'ΔH+NMODE', 'ΔH+QH']
        hide_col = [c for c, x in enumerate(columns, start=1) if df[x].isnull().all()]
        for m in models:
            model_df = df[df['MODEL'].isin([m])]
            m_col_box = []
            model_data = []
            for c, col in enumerate(columns, start=1):
                item_col_df = model_df[['System', col, 'Exp.Energy']]
                if not item_col_df[col].isnull().all():
                    m_col_box.append(c)
                    model_data.append(item_col_df)
                else:
                    model_data.append(None)
            if m_col_box:
                item = CorrelationItem(self.correlation_treeWidget, [m.upper()], model=model_df, enthalpy=model_data[0],
                                       dgie=model_data[1], dgnmode=model_data[2], dgqh=model_data[3], col_box=m_col_box)
        for x in hide_col:
            self.correlation_treeWidget.hideColumn(x)

    def read_data(self, queue: Queue, options):
        self.init_dialog.accept()
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

        # check if all systems have the same frames range
        frange_base = []

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
                config = path.parent
            else:
                config = None

            self.systems[i] = {'name': name, 'path': path.parent, 'api': api,
                               'namespace': namespace, 'data': energy,
                               'current_frames': [namespace.INPUT['startframe'],
                                                  namespace.INPUT['endframe'],
                                                  namespace.INPUT['interval']],
                               'current_nmode_frames': [namespace.INPUT['nmstartframe'],
                                                        namespace.INPUT['nmendframe'],
                                                        namespace.INPUT['nminterval']],
                               'current_ie_segment': namespace.INPUT['ie_segment'],
                               'chart_options': ChartSettings(config),
                               'options': options,
                               'items_data': {}, 'items_summary': summary,
                               'exp_ki': exp_ki
                               }
            self.makeTree(i)
            if not frange_base:
                frange_base = [namespace.INPUT['startframe'], namespace.INPUT['endframe'],
                               namespace.INPUT['interval'], namespace.INPUT['nmstartframe'],
                               namespace.INPUT['nmendframe'], namespace.INPUT['nminterval']]
                continue
            if frange_base != [namespace.INPUT['startframe'], namespace.INPUT['endframe'],
                               namespace.INPUT['interval'], namespace.INPUT['nmstartframe'],
                               namespace.INPUT['nmendframe'], namespace.INPUT['nminterval']]:
                self.all_systems_active = False

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
        for c in [0, 1]:
            sys_item.setBackground(c, QBrush(QColor(100, 100, 100)))
            sys_item.setForeground(c, QBrush(QColor(220, 220, 255)))
            font = sys_item.font(c)
            font.setWeight(QFont.Weight.DemiBold)
            font.setPointSize(font.pointSize() + 1)
            sys_item.setFont(c, font)

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
            self.setting_item_data(sys_index, part)
            self.makeItems(sys_index, part, classif_item)

        if self.systems[sys_index]['namespace'].INPUT['idecomp']:
            classif_item = CustomItem(sys_item, ['Decomposition'])
            for part in decomp_print_keys:
                if part not in self.systems[sys_index]['data']:
                    continue
                self.setting_item_data(sys_index, part)
                self.makedecompItems(sys_index, part, classif_item)

        sys_item.setExpanded(True)

        # setup items to qtreewidget
        itemiter = QTreeWidgetItemIterator(sys_item)
        while itemiter.value():
            item = itemiter.value()
            if sb := item.setup_buttons():
                self.treeWidget.setItemWidget(item, 1, sb)
            itemiter += 1

    def _remove_empty_terms(self, data):
        if self.data_options['remove_empty_terms'] and isinstance(
            data, pd.DataFrame
        ):
            columns = data.columns
            for col in columns:
                if -0.01 < data[col].mean() < 0.01 and col not in ['GSOLV', 'GGAS', 'TOTAL', 'tot']:
                    del data[col]
            return data
        return data

    def _remove_empty_charts(self, data):
        if isinstance(data, pd.Series):
            if self.data_options['remove_empty_charts'] and ((data > -0.01).all() and (data < 0.01).all()):
                return True
        elif self.data_options['remove_empty_charts'] and ((data > -0.01).all().all() and (data < 0.01).all().all()):
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
            temp_data = self._remove_empty_terms(data)
            options = {'c2': True} if iec2 else {}
            options.update(dict(groups=self._itemdata_properties(temp_data)))
            cont['bar_plot_data'] = [temp_data, options, change]
            del temp_data
        elif level == 2:
            tempdf = data.loc[:, data.columns.get_level_values(1) == 'tot']
            bar_plot_data = tempdf.droplevel(level=1, axis=1)
            cont['line_plot_data'] = [bar_plot_data.sum(axis=1), {}, change]
            cont['heatmap_plot_data'] = [bar_plot_data.transpose(copy=True), {}, change]
            temp_bar_data = self._remove_empty_terms(bar_plot_data)
            cont['bar_plot_data'] = [temp_bar_data, dict(groups=self._itemdata_properties(temp_bar_data)), change]
            del bar_plot_data
            del temp_bar_data
            del tempdf
        elif level == 3:
            # Select only the "tot" column, remove the level, change first level of columns to rows and remove the mean
            # index
            tempdf = data.loc[:, data.columns.get_level_values(2) == 'tot']
            cont['heatmap_plot_data'] = [
                tempdf.aggregate(["mean"]).droplevel(level=2, axis=1).stack().droplevel(level=0), {}, change]
            bar_plot_data = tempdf.groupby(axis=1, level=0, sort=False).sum()
            cont['line_plot_data'] = [bar_plot_data.sum(axis=1), {}, change]
            temp_bar_data = self._remove_empty_terms(bar_plot_data)
            cont['bar_plot_data'] = [temp_bar_data, dict(groups=self._itemdata_properties(temp_bar_data)), change]
            del tempdf
            del bar_plot_data
            del temp_bar_data

        return cont

    def makeItems(self, sys_index, part, classif_item):

        correlation_data = self.corr_data
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

        for level in ['gb', 'pb', 'rism gf', 'rism std', 'rism pcplus', 'nmode', 'qh', 'ie', 'c2', 'binding']:
            if level in data:
                if level in ['gb', 'pb', 'rism gf', 'rism std', 'rism pcplus', 'nmode', 'qh']:
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
                                           system_index=sys_index,
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )
                        if level == 'qh':
                            print('Currently not implemented. Please contact us to use your files to setup this analysis'
                                  ' in gmx_MMPBSA ')
                            continue
                        for level2 in str_dict[level1]:
                            keys_path = (part, level, (level1, level2))
                            if keys_path in self.removed_items[sys_index]:
                                continue
                            item2 = CustomItem(item1,
                                               [level2.upper()],
                                               app=self,
                                               buttons=(1,),
                                               title="Energetic Components",
                                               subtitle=f"{sys_name} | {level.upper()} | "
                                                        f"{level1.upper()} | {level2.upper()}",
                                               system_index=sys_index,
                                               keys_path=keys_path,
                                               part=part
                                               )
                            self.items_counter['charts'] += 1
                elif level == 'c2':
                    titem = CustomItem(top_item,
                                       [level.upper()],
                                       app=self,
                                       buttons=(-2,),
                                       title="C2 Entropy",
                                       subtitle=f"{sys_name} | {level.upper()}",
                                       system_index=sys_index,
                                       keys_path=(part, level),
                                       part=part
                                       )
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(2,),
                                           title="C2 Entropy",
                                           subtitle=f"{sys_name} | {level.upper()} | {level1.upper()}",
                                           system_index=sys_index,
                                           keys_path=(part, level, (level1,)),
                                           part=part
                                           )
                elif level == 'ie':
                    titem = CustomItem(top_item,
                                       [level.upper()],
                                       app=self,
                                       buttons=(-2,),
                                       title="Interaction Entropy",
                                       subtitle=f"{sys_name} | {level.upper()}",
                                       system_index=sys_index,
                                       keys_path=(part, level),
                                       part=part)
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        item1 = CustomItem(titem,
                                           [level1.upper()],
                                           app=self,
                                           buttons=(1,),
                                           title="Interaction Entropy",
                                           subtitle=f"{sys_name} | {level.upper()} | {level1.upper()}",
                                           system_index=sys_index,
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
                                           system_index=sys_index,
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
                                           system_index=sys_index,
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
                                               title=f"Energetic Components {title}",
                                               subtitle=f"{sys_name} ({part}) | "
                                                          f"{str(level).upper()} | "
                                                          f"{str(level1).upper()} | "
                                                          f"{str(level2).upper()} | "
                                                          f"{str(level3).upper()}",
                                               system_index=sys_index,
                                               keys_path=(part, level, (level1, level2, level3))
                                               )
                            # residue first level
                            #             btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                            for level4 in str_dict[level1][level2][level3]:
                                if namespace.INPUT['idecomp'] in [1, 2]:
                                    item4 = CustomItem(item3,
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
                                                       system_index=sys_index,
                                                       keys_path=(part, level, (level1, level2, level3, level4))
                                                       )
                                    self.items_counter['charts'] += 1
                                else:
                                    item4 = CustomItem(item3,
                                                       [level4.upper()],
                                                       app=self,
                                                       buttons=(2,),
                                                       title="Energetic Components [Per-wise]",
                                                       subtitle=f"{sys_name} ({part})| "
                                                                f"{str(level).upper()} | "
                                                                f"{str(level1).upper()} | "
                                                                f"{str(level2).upper()} | "
                                                                f"{str(level3).upper()} | "
                                                                f"{str(level4).upper()}",
                                                       system_index=sys_index,
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
                                                           system_index=sys_index,
                                                           keys_path=(
                                                               part, level, (level1, level2, level3, level4, level5))
                                                           )

    def setting_item_data(self, sys_index, part, comp=('all',)):
        correlation_data = self.corr_data
        # correlation_data[sys_name] = {'ΔG': {
        #                                     'gb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'pb': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism std': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan},
        #                                     'rism gf': {'ΔH': np.nan, 'ie': np.nan, 'nmode': np.nan, 'qh': np.nan}},
        #                               'Exp.Energy': ki2energy(topItem.exp_ki, topItem.temp)}

        data = self.systems[sys_index]['data'][part]
        self.removed_items[sys_index] = []
        namespace = self.systems[sys_index]['namespace']

        parts = self.systems[sys_index]['options']['components'] + ['delta']
        if namespace.FILES.stability:
            parts.append('complex')

        key_list = []
        if 'nmode' in comp:
            key_list.append('nmode')
        elif 'ie' in comp:
            key_list.append('ie')
        elif 'energy' in comp:
            # Include ie and c2 since they dependent of the ggas energy
            key_list.extend(['gb', 'pb', 'rism gf', 'rism std', 'rism pcplus', 'binding', 'ie', 'c2'])
        elif 'all' in comp:
            key_list.extend(['gb', 'pb', 'rism gf', 'rism std', 'rism pcplus', 'nmode', 'qh', 'ie', 'c2', 'binding'])
        for level in key_list:
            if level in data:
                if level in ['gb', 'pb', 'rism gf', 'rism std', 'rism pcplus', 'nmode', 'qh']:
                    str_dict = multiindex2dict(data[level].columns)
                    for level1 in str_dict:
                        if level1 not in parts:
                            continue
                        if not part.startswith('decomp'):
                            self.systems[sys_index]['items_data'][(part, level, (level1,))] = self._setup_data(data[level][(level1,)],
                                                                                                               level=1)
                        for level2 in str_dict[level1]:
                            if self._remove_empty_charts(data[level][(level1, level2)]):
                                del data[level][(level1, level2)]
                                self.removed_items[sys_index].append((part, level, (level1, level2)))
                                continue
                            if not part.startswith('decomp'):
                                temp_dat = data[level][(level1, level2)]
                                temp_dat.name = level2
                                self.systems[sys_index]['items_data'][
                                    (part, level, (level1, level2))] = self._setup_data(
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
                                for level3 in str_dict[level1][level2]:
                                    item_lvl2 = 1 if namespace.INPUT['idecomp'] in [1, 2] else 2
                                    temp_dat = data[level][(level1, level2, level3)]
                                    temp_dat.name = level3
                                    self.systems[sys_index]['items_data'][(part, level, (level1, level2, level3))] = \
                                        self._setup_data(temp_dat, level=item_lvl2)
                                    del temp_dat
                                    # residue first level
                                    for level4 in str_dict[level1][level2][level3]:
                                        if self._remove_empty_charts(
                                                data[level][(level1, level2, level3, level4)]
                                        ):
                                            del data[level][(level1, level2, level3, level4)]
                                            self.removed_items[sys_index].append(
                                                (part, level, (level1, level2, level3, level4)))
                                            continue

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
                                                if self._remove_empty_charts(
                                                        data[level][(level1, level2, level3, level4, level5)]
                                                ):
                                                    del data[level][(level1, level2, level3, level4, level5)]
                                                    self.removed_items[sys_index].append(
                                                        (part, level, (level1, level2, level3, level4, level5)))
                                                    continue
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
