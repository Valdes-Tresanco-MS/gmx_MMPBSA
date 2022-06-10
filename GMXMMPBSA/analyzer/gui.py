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

from queue import Queue
from functools import partial
import pandas as pd

from GMXMMPBSA import utils
from GMXMMPBSA.analyzer.items_delegate import KiTableDelegate
from GMXMMPBSA.analyzer.plots import Tables
from GMXMMPBSA.API import MMPBSA_API
from GMXMMPBSA.analyzer.dialogs import InitDialog, ProcessingProgressBar
from GMXMMPBSA.analyzer.customitem import CustomItem, CustomCorrItem
from GMXMMPBSA.analyzer.style import save_default_config, default_config, save_user_config, user_config, toc_img, logo, \
    alert, config
from GMXMMPBSA.analyzer.utils import ki2energy
from GMXMMPBSA.analyzer.chartsettings import ChartSettings, CorrChartSettings
from GMXMMPBSA.analyzer.parametertree import ParameterTree, Parameter
import math


class GMX_MMPBSA_ANA(QMainWindow):
    def __init__(self, ifiles):
        super(GMX_MMPBSA_ANA, self).__init__()
        self.showMaximized()
        self.setWindowIcon(QIcon(logo))
        self.correlation = {'chart_options': CorrChartSettings(), 'items_data': {}, 'data': {}}

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
        self.tabifyDockWidget(self.treeDockWidget, self.correlation_DockWidget)


        self.optionDockWidget = QDockWidget('Options', self)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.optionDockWidget)

        self.treeWidget = QTreeWidget(self)
        self.treeWidget.setIndentation(12)
        self.treeWidget.setUniformRowHeights(True)
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

        self.correlation_tableWidget = QTableWidget(self)

        self.correlation_tableWidget.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.correlation_tableWidget.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.correlation_tableWidget.itemSelectionChanged.connect(self.update_table)
        self.correlation_tableWidget.setColumnCount(6)
        self.correlation_tableWidget.verticalHeader().hide()
        header_tt = {'MODEL': 'Selected Model',
                     'ΔGeff': '''<html>Correlation plot for ΔG<sub>effective</sub>. Energy when entropy contribution is 
        neglected. Calculated as ΔG<sub>effective</sub> = ΔH</html>''',
                     'ΔGie': '''<html>Correlation plot for ΔG<sub>ie</sub>. Energy 
        when entropy contribution is calculated using the Interaction Entropy method. Calculated as 
        ΔG<sub>ie</sub> = ΔH + -TΔS<sub>IE</sub></html>''',
                     'ΔGc2': '''<html>Correlation plot for ΔG<sub>c2</sub>. Energy 
        when entropy contribution is calculated using the C2 Entropy method. Calculated as 
        ΔG<sub>c2</sub> = ΔH + -TΔS<sub>C2</sub></html>''',
                     'ΔGnm': '''<html>Correlation plot for ΔG<sub>nm</sub>. Energy 
        when entropy contribution is calculated using the Normal Modes Entropy method. Calculated as 
        ΔG<sub>nm</sub> = ΔH + -TΔS<sub>NMODE</sub></html>''',
                     'ΔGqh': '''<html>Correlation plot for ΔG<sub>ie</sub>. Energy 
        when entropy contribution is calculated using the Quasi-Harmonic Entropy method. Calculated as 
        ΔG<sub>qh</sub> = ΔH + -TΔS<sub>QH</sub></html>'''}

        for c, (i, t) in enumerate(header_tt.items()):
            hitem = QTableWidgetItem(i)
            hitem.setToolTip(t)
            self.correlation_tableWidget.setHorizontalHeaderItem(c, hitem)
            if c !=0:
                self.correlation_tableWidget.setColumnHidden(c, True)
        cheader = self.correlation_tableWidget.horizontalHeader()
        cheader.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        self.data_table_widget = QTableWidget(0, 6)
        self.data_table_widget.setHorizontalHeaderLabels(['Sys.', 'Type', 'Exp.ΔG', 'Avg.', 'SD', 'SEM'])
        self.data_table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.data_table_widget.verticalHeader().hide()

        self.corr_parm_tree_w = ParameterTree()
        self.corr_parms = Parameter().create()
        self.corr_parms.restoreState(CorrChartSettings())
        self.corr_parm_tree_w.setParameters(self.corr_parms, showTop=False)

        self.corr_sys_sel_w = QTableWidget()
        self.corr_sys_sel_w.verticalHeader().hide()
        self.corr_sys_sel_w.setColumnCount(6)
        self.corr_sys_sel_w.setHorizontalHeaderLabels(['Id', 'Sel.', 'Ref.', 'Type', '  Exp.Ki  ', 'Name'])
        self.corr_sys_sel_w.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.corr_sys_sel_w.horizontalHeader().setStretchLastSection(True)
        delegate = KiTableDelegate(self.corr_sys_sel_w)
        self.corr_sys_sel_w.setItemDelegateForColumn(4, delegate)
        self.corr_sys_sel_w.cellChanged.connect(self.corr_sys_sel_itemchange)

        self.corr_options_widget = QTabWidget()
        self.corr_options_widget.addTab(self.corr_parm_tree_w, 'Charts')
        self.corr_options_widget.addTab(self.corr_sys_sel_w, 'Selections')

        self.update_corr_btn = QPushButton(f'Update {self.corr_options_widget.tabText(0)}')
        self.update_corr_btn.clicked.connect(self.update_corr_fn)
        self.corr_options_widget.setCornerWidget(self.update_corr_btn)

        self.corr_options_widget.currentChanged.connect(partial(self._change_update_btn, tw='corr'))

        self.corr_container_widget = QWidget()
        self.corr_container_widget_layout = QVBoxLayout(self.corr_container_widget)
        # self.corr_container_widget.addWidget(self.correlation_treeWidget)
        self.corr_container_widget_layout.addWidget(self.correlation_tableWidget, 1)
        self.corr_container_widget_layout.addWidget(self.data_table_widget, 5)
        self.corr_container_widget_layout.addWidget(self.corr_options_widget, 5)
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
            msgBox.setDetailedText("The systems below has some changes:\n - " +
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

    def corr_sys_sel_itemchange(self, row, col):
        citem = self.corr_sys_sel_w.item(row, col)
        if col == 1:
            if ritemw := self.corr_sys_sel_w.cellWidget(row, 2):
                if citem.checkState() == Qt.CheckState.Unchecked:
                    if ritemw.isChecked():
                        selected = {i: abs(row - i) for i in range(self.corr_sys_sel_w.rowCount())
                                    if self.corr_sys_sel_w.item(i, 1).checkState() == Qt.CheckState.Checked}
                        if not selected:
                            citem.setCheckState(Qt.CheckState.Checked)
                            return
                        nr = list(selected.keys())[list(selected.values()).index(min(selected.values()))]
                        self.corr_sys_sel_w.cellWidget(nr, 2).setChecked(True)
                    ritemw.setEnabled(False)
                else:
                    ritemw.setEnabled(True)

    def _about_dialog(self):
        from GMXMMPBSA import __version__
        QMessageBox.about(self.mdi, "About gmx_MMPBSA",
                              "<html>"
                              "<body>"
                              "<h2 style='text-align:center'>About gmx_MMPBSA</h2>"
                              "<div text-align:center;>"
                              f"<a href='https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/'>"
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
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/analyzer/'))

    def _doc(self):
        QDesktopServices().openUrl(QUrl('https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/getting-started/'))

    def _bug(self):
        QDesktopServices().openUrl(QUrl('https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/new/choose'))

    def _group(self):
        QDesktopServices().openUrl(QUrl('https://groups.google.com/g/gmx_mmpbsa'))

    def _initialize_systems(self):

        if self.all_systems_active and len(self.systems) > 1:
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
        self.statusbar.showMessage(f"Using this setting as default in "
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

        self.optionWidget_c = QTabWidget(self)
        optionWidget_l.addWidget(self.optionWidget_c)

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
                                       'for the selected systems')
        charts_opt_menu.addSeparator()
        self.set_user_config_action = charts_opt_menu.addAction(QIcon(save_user_config), 'Save User-config',
                                                              self._set_as_default)
        self.set_user_config_action.setToolTip('Force save the current setting for this system. Before gmx_MMPBSA_ana '
                                               'closes, if there are configuration changes, you can decide if you '
                                               'want to save them')

        self.use_user_config_action = charts_opt_menu.addAction(QIcon(user_config), 'User-config', self._set_as_default)
        self.use_user_config_action.setToolTip('Use the specific settings for this system if it was saved')
        charts_opt_tb.setMenu(charts_opt_menu)

        self.optionWidget_c.addTab(self.chart_options_w, 'Charts')

        frames_w = QWidget(self.optionWidget_c)
        self.optionWidget_c.addTab(frames_w, 'Frames')
        frames_l = QVBoxLayout(frames_w)

        self.eframes_group = QGroupBox('Energy')
        frames_l.addWidget(self.eframes_group)
        self.eframes_start_sb = QSpinBox()
        self.eframes_start_sb.setAccelerated(True)
        self.eframes_start_sb.valueChanged.connect(self.frames_start_sb_update)

        self.eframes_inter_sb = QSpinBox()
        self.eframes_inter_sb.setAccelerated(True)
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
        fg_l2.addRow('Nº frames', self.numframes_le)

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
        self.nmframes_inter_sb.setAccelerated(True)
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
        nmfg_l2.addRow('Nº frames', self.nmnumframes_le)

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
        ienf_l.addRow('Nº frames', self.ienumframes_le)

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

        # correlation_w = QWidget(self.optionWidget_c)
        # self.optionWidget_c.addTab(correlation_w, 'Correlation')
        self.update_btn = QPushButton(f'Update {self.optionWidget_c.tabText(0)}')
        self.update_btn.clicked.connect(self.update_fn)

        self.optionWidget_c.setCornerWidget(self.update_btn)

        #
        self.optionWidget_c.currentChanged.connect(self._change_update_btn)

    def _change_update_btn(self, i, tw=None):
        if tw == 'corr':
            text = self.corr_options_widget.tabText(i)
            self.update_corr_btn.setText(f"Update {text}")
        else:
            text = self.optionWidget_c.tabText(i)
            self.update_btn.setText(f"Update {text}")

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
            curr_eframes = dict(
                startframe=self.eframes_start_sb.value(),
                endframe=self.eframes_end_sb.value(),
                interval=self.eframes_inter_sb.value()
            )
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
        # energy
        if self.systems[parent_item.system_index].get('current_frames'):
            current_start, current_end, current_interval = self.systems[parent_item.system_index]['current_frames'].values()
            self.eframes_start_sb.setValue(current_start)
            self.eframes_inter_sb.setValue(current_interval)
            self.eframes_end_sb.setValue(current_end)
        # nmode
        if self.systems[parent_item.system_index].get('current_nmode_frames'):
            nmcurrent_start, nmcurrent_end, nmcurrent_interval = self.systems[parent_item.system_index][
                'current_nmode_frames'].values()
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

        replot_energy = False
        sys_list = self.systems if self.all_frb.isChecked() else [self.current_system_index]
        if self.optionWidget_c.currentIndex() == 0:  # Chart Options
            self.statusbar.showMessage('Updating Charts...')
            act_sett = self.chart_options_param.saveState()

            if repaint := [x for x in sys_list if self.systems[x]['chart_options'].is_changed(act_sett) or repaint_arg]:
                for e in repaint:
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
        else:  # Frames Options
            self.statusbar.showMessage('Updating Energy for selected Frame range...')
            curr_eframes = dict(
                startframe=self.eframes_start_sb.value(),
                endframe=self.eframes_end_sb.value(),
                interval=self.eframes_inter_sb.value()
            )

            curr_nmframes = dict(
                nmstartframe=self.nmframes_start_sb.value(),
                nmendframe=self.nmframes_end_sb.value(),
                nminterval=self.nmframes_inter_sb.value()
            )
            curr_iesegment = self.iesegment_sb.value()
            iq = Queue()
            rq = Queue()

            for x in sys_list:
                energyrun = (self.systems[x]['namespace'].INPUT['gbrun'] or
                             self.systems[x]['namespace'].INPUT['pbrun'] or
                             self.systems[x]['namespace'].INPUT['rismrun'])
                if (energyrun and curr_eframes != self.systems[x]['current_frames'] or
                        self.systems[x]['namespace'].INPUT['nmoderun'] and
                        curr_nmframes != self.systems[x]['current_nmode_frames'] or
                        curr_iesegment != self.systems[x]['current_ie_segment']
                ):
                    if not self.systems[x]['namespace'].INPUT['nmoderun']:
                        curr_nmframes = {}
                    self._extracted_from_update_fn_69(x, curr_eframes, curr_nmframes, curr_iesegment, iq)
                    self.systems[x]['current_frames'] = curr_eframes
                    self.systems[x]['current_nmode_frames'] = curr_nmframes
                    self.systems[x]['current_ie_segment'] = curr_iesegment

            if iq.qsize():
                replot_energy = True
                from GMXMMPBSA.analyzer.dialogs import ProcessingProgressBar
                pbd = ProcessingProgressBar(self, iq, rq, self.options['performance']['jobs'], 'Updating...')
                pbd.rejected.connect(lambda z: print(f'rejected> exit code: {z}'))
                pbd.accepted.connect(lambda: self._update_itemdata_repaint(rq))
                pbd.exec()

        subwindows = self.mdi.subWindowList()
        for c, s in enumerate(sys_list):
            changes = self.systems[s]['chart_options'].changes
            for sub in subwindows:
                if not sub.item_parent or sub.item_parent.system_index != s:
                    if isinstance(sub, Tables):
                        sub.close()
                        sub.button.setChecked(True)
                    continue
                if not sub.isVisible():
                    if self.systems[s]['chart_options'].changes.get('line_action'):
                        sub.item_parent.lp_subw = None
                    elif self.systems[s]['chart_options'].changes.get('bar_action'):
                        sub.item_parent.bp_subw = None
                    elif self.systems[s]['chart_options'].changes.get('heatmap_action'):
                        sub.item_parent.hmp_subw = None

                if changes['figure']:
                    chart_sett = self.systems[s]['chart_options'].get_settings()
                    options = {'save-format': chart_sett[('General', 'figure-format', 'save-format')],
                               'dpi-save': chart_sett[('General', 'figure-format', 'dpi-save')]
                               }
                    sub.mpl_toolbar.update_options(options)
                    sub.fbtn.setChecked(chart_sett[('General', 'toolbar')])
                elif sub.isVisible() and any(changes.values()) or replot_energy:
                    sub.close()
                    sub.button.setChecked(True)
            pymol_items = [[p, item] for p, item in self.pymol_p_list if p.state() == QProcess.ProcessState.Running]
            for p, item in pymol_items:
                p.kill()
                p.waitForFinished()
                item.vis_action.setChecked(True)
            # re-assign changes to native state after replot all open charts
            self.systems[s]['chart_options'].changes = dict(line_action=0, line_ie_action=0, bar_action=0,
                                                            heatmap_action=0, visualization_action=0, figure=0)
        self.ie_changed.hide()
        self.e_changed.hide()
        self.nm_changed.hide()
        self.statusbar.showMessage('Updating... Done.', 2000)

    def update_corr_fn(self, repaint_arg=False):

        if self.corr_options_widget.currentIndex() == 0: # Chart Options
            self.statusbar.showMessage('Updating Charts...')
            act_sett = self.corr_parms.saveState()
            if repaint := self.correlation['chart_options'].is_changed(act_sett) or repaint_arg:
                if repaint_arg:
                    self.correlation['chart_options'].changes = True
                else:
                    self.correlation['chart_options'].get_changes(act_sett)
        else:   # calculate correlation
            self.statusbar.showMessage('Calculating correlation...')
            # Get any modification in selection systems table
            cols = {1: ('System', 'Selection'), 2: ('System', 'Reference'), 4: ('System', 'ExpΔG')}
            for c in range(self.corr_sys_sel_w.columnCount()):
                col_values = []
                for r in range(self.corr_sys_sel_w.rowCount()):
                    if c == 1:
                        col_values.append(self.corr_sys_sel_w.item(r, c).checkState() == Qt.CheckState.Checked)
                    elif c == 2:
                        col_values.append(self.corr_sys_sel_w.cellWidget(r, c).isChecked())
                    elif c == 4:
                        col_values.append(
                            ki2energy(float(self.corr_sys_sel_w.item(r, c).text()),
                                      self.systems[int(self.corr_sys_sel_w.item(r, 0).text())]['namespace'].INPUT[
                                          'temperature']))
                if col_values:
                    for m, v in self.correlation['data'].items():
                        v[cols[c]] = col_values
            self.data_table_widget.clearContents()
            item_list = self.correlation_tableWidget.selectedItems()
            self.correlation_tableWidget.clearSelection()

            for m, v in self.correlation['data'].items():
                met = v.columns.get_level_values(0).unique().to_list()[1:]
                for m1 in met:
                    self.correlation['items_data'][(m, m1)] = [v.loc[  # FIXME: user option to show the ref in reg plot ???
                                                           v[('System', 'Reference')] == False
                                                           ].loc[v[('System', 'Selection')] == True].loc[:,
                                                       v.columns.get_level_values(0).isin(['System', m1])
                                                       ].droplevel(0, axis=1), True]
            if item_list:
                item_list[0].setSelected(True)
            self.statusbar.showMessage('Calculating correlation... Done.', 2000)

        for c in range(self.correlation_tableWidget.columnCount()):
            for r in range(self.correlation_tableWidget.rowCount()):
                itemw = self.correlation_tableWidget.cellWidget(r, c)
                if itemw and itemw.reg_chart_action.isChecked():
                    itemw.reg_chart_action.setChecked(False)
                    itemw.reg_chart_action.setChecked(True)

    def _update_itemdata_repaint(self, r):
        maximum = r.qsize()
        for i, c in enumerate(range(maximum), start=1):
            sys_id, data, a = r.get()
            self.systems[sys_id]['items_data'] = {k:v1 for x, v in data.items()
                                                  if x != 'correlation' and v
                                                  for k, v1 in v['keys'].items()}
            self.systems[sys_id]['correlation_data'] = data.get('correlation')
        if self.options.get('correlation')['corr']:
            self.get_corr_data()

    # TODO Rename this here and in `update_fn`
    def _extracted_from_update_fn_69(self, sys_index, curr_eframes, curr_nmframes, curr_iesegment, iq):
        o = self.systems[sys_index]['anaoptions']
        o['energy_options'].update(curr_eframes)
        for d in [curr_eframes, curr_nmframes, dict(ie_segment=curr_iesegment)]:
            o['entropy_options'].update(d)
        o['decomp_options'].update(curr_eframes)
        d = {self.systems[sys_index]['api'].get_ana_data: {
            'object': 'function',
            'args': o}
        }
        iq.put((sys_index, d))

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

            curr_eframes = dict(startframe=value, endframe=self.eframes_end_sb.value(),
                                interval=self.eframes_inter_sb.value())
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

            curr_eframes = dict(startframe=self.eframes_start_sb.value(), endframe=value,
                                interval=self.eframes_inter_sb.value())
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

            curr_eframes = dict(startframe=self.eframes_start_sb.value(), endframe=self.eframes_end_sb.value(),
                                interval=value)
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
            self.nmnumframes_le.setText(f"{int((current_end - value) // current_interval) + 1}")
            curr_nmframes = dict(nmstartframe=value, nmendframe=self.nmframes_end_sb.value(),
                                 nminterval=self.nmframes_inter_sb.value())
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
            self.nmnumframes_le.setText(f"{int((value - current_start) // current_interval) + 1}")

            curr_nmframes = dict(nmstartframe=self.nmframes_start_sb.value(), nmendframe=value,
                                 nminterval=self.nmframes_inter_sb.value())
            if curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']:
                self.nm_changed.show()
            else:
                self.nm_changed.hide()

    def nmframes_inter_sb_update(self, value):
        if self.current_system_index:
            current_start = self.nmframes_start_sb.value()
            current_end = self.nmframes_end_sb.value()
            self.nmnumframes_le.setText(f"{int((current_end - current_start) // value) + 1}")

            curr_nmframes = dict(nmstartframe=self.nmframes_start_sb.value(), nmendframe=self.nmframes_end_sb.value(),
                                 nminterval=value)
            if curr_nmframes != self.systems[self.current_system_index]['current_nmode_frames']:
                self.nm_changed.show()
            else:
                self.nm_changed.hide()

    def gettting_data(self):
        self.in_init_dialog = True
        self.init_dialog.get_files_info(self.ifiles)
        self.init_dialog.show()

    def update_table(self):

        items = self.correlation_tableWidget.selectedItems()
        if not items:
            return
        item = items[0]
        if not item.keys_path:
            return
        df = self.correlation['items_data'][item.keys_path][0]
        data = df.loc[df['Selection'] == True][['Number', 'Type', 'ExpΔG', 'Average', 'SD', 'SEM']]
        self.data_table_widget.setRowCount(len(data.index))
        for i, c in enumerate(data):
            for j, r in enumerate(data[c]):
                it = QTableWidgetItem(f'{r:.2f}' if isinstance(r, float) else str(r))
                it.setFlags(it.flags() ^ Qt.ItemFlag.ItemIsEditable)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.data_table_widget.setItem(j, i, it)

    def make_correlation(self):

        self.get_corr_data()
        self.corr_sys_sel_w.setRowCount(len(self.systems))
        ref_group = QButtonGroup()
        ref = False
        r = 0
        num_sel_sys = 0
        for sys_id in self.systems:
            for ct, cv in self.systems[sys_id].get('correlation').items():
                idi = QTableWidgetItem(str(sys_id))
                idi.setFlags(idi.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.corr_sys_sel_w.setItem(r, 0, idi)
                ti = QTableWidgetItem()
                ti.setCheckState(Qt.CheckState.Checked if cv else Qt.CheckState.Unchecked)
                ti.setFlags(ti.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.corr_sys_sel_w.setItem(r, 1, ti)
                if self.systems[sys_id].get('reference'):
                    ref = True
                ri = QRadioButton()
                ri.setChecked(self.systems[sys_id].get('reference'))
                if not cv:
                    ri.setEnabled(False)
                ref_group.addButton(ri, r)
                self.corr_sys_sel_w.setCellWidget(r, 2, ri)

                tyi = QTableWidgetItem(ct[0].upper())
                tyi.setFlags(tyi.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.corr_sys_sel_w.setItem(r, 3, tyi)
                self.corr_sys_sel_w.setItem(r, 4, QTableWidgetItem(str(self.systems[sys_id]['exp_ki'][ct])))
                ni = QTableWidgetItem(self.systems[sys_id]['name'])
                ni.setFlags(ni.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.corr_sys_sel_w.setItem(r, 5, ni)
                if cv:
                    num_sel_sys += 1
                r += 1
        if not ref:
            self.corr_sys_sel_w.hideColumn(2)
        self.data_table_widget.setRowCount(num_sel_sys)

        self.correlation_tableWidget.setRowCount(len(self.correlation['data']))
        height = ((self.correlation_tableWidget.model().rowCount() - 1) +
                 self.correlation_tableWidget.horizontalHeader().height())
        for row in range(self.correlation_tableWidget.model().rowCount()):
            height += self.correlation_tableWidget.rowHeight(row)
        self.correlation_tableWidget.setMinimumHeight(height + 5)

        met_col = {'ΔGeff': 1, 'ΔGie': 2, 'ΔGc2': 3, 'ΔGnm': 4, 'ΔGqh': 5}
        for r, (m, v) in enumerate(self.correlation['data'].items()):
            citem = CustomCorrItem(text=True, model=m)
            citem.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.correlation_tableWidget.setItem(r, 0, citem)

            met = v.columns.get_level_values(0).unique().to_list()[1:]
            for m1 in met:
                self.correlation_tableWidget.setColumnHidden(met_col[m1], False)
                item = CustomCorrItem(app=self, model=m, keys_path=(m, m1))
                self.correlation_tableWidget.setItem(r, met_col[m1], item)
                item.define_button(r, met_col[m1])

    def get_corr_data(self):
        self.correlation['data'] = {}
        for s, sys_id in enumerate(self.systems):
            dta = {('System', 'Number'): [], ('System', 'Type'): [], ('System', 'Reference'): [],
                   ('System', 'ExpΔG'): [], ('System', 'Selection'): []}
            for ct, cv in self.systems[sys_id].get('correlation').items():
                dta[('System', 'Number')].append(sys_id)
                dta[('System', 'Type')].append(ct[0].upper())
                if ct == 'normal':
                    dta[('System', 'Reference')].append(self.systems[sys_id].get('reference'))
                else:
                    dta[('System', 'Reference')].append(False)
                dta[('System', 'ExpΔG')].append(
                    ki2energy(self.systems[sys_id]['exp_ki'][ct],
                              self.systems[sys_id]['namespace'].INPUT['temperature']))
                dta[('System', 'Selection')].append(cv)

            g = pd.DataFrame(dta, index=range(len(self.systems[sys_id].get('correlation'))))
            for et, v in self.systems[sys_id]['correlation_data'].items():
                for m, v1 in v.items():
                    cdf = pd.concat([g, v1], axis=1)
                    d = self.correlation['data'].setdefault(m, cdf)
                    if d.equals(cdf):
                        continue
                    self.correlation['data'][m] = pd.concat([d, cdf], ignore_index=True).sort_values(
                        by=[('System', 'Number'), ('System', 'Type')])

        def subtract_custom_value(x, custom_value):
            return x - custom_value

        for r, (m, v) in enumerate(self.correlation['data'].items()):
            ref = v.loc[v[('System', 'Reference')] == True].loc[:, ~v.columns.get_level_values(1).isin(
                ['Number', 'Type', 'Reference', 'Selection'])]
            if len(ref.index):
                for k in ref:
                    if k[-1] == 'SEM':
                        continue
                    elif k[-1] == 'SD':
                        v[k] = v[k].apply(utils.get_corrstd, args=(ref[k].values[0],))
                        v[(k[0], 'SEM')] = v[k]
                    else:
                        v[k] = v[k].apply(subtract_custom_value, args=(ref[k].values[0],))

            met = v.columns.get_level_values(0).unique().to_list()[1:]
            for m1 in met:
                self.correlation['items_data'][(m, m1)] = [v.loc[  # FIXME: user option to show the ref in reg plot ???
                                                               v[('System', 'Reference')] == False
                                                               ].loc[v[('System', 'Selection')] == True].loc[:,
                                                           v.columns.get_level_values(0).isin(['System', m1])
                                                           ].droplevel(0, axis=1), True]

    def read_data(self, queue: Queue, options):
        self.init_dialog.accept()
        self.in_init_dialog = False
        self.options = options
        max_sixe = queue.qsize()

        iq = Queue()
        self.rq = Queue()

        for i, _ in enumerate(range(max_sixe), start=1):
            sys_name, path, stability, exp_ki, corr, ref, options_file = queue.get()
            opts4ana = {'energy_options': {},
                        'entropy_options': {},
                        'decomp_options': {'res_threshold': options.get('decomposition').get('res_threshold')},
                        'correlation': options.get('correlation').get('corr'),
                        'performance_options': options.get('performance')}
            opts4ana['energy_options']['etype'] = opts4ana['entropy_options']['etype'] = (None if len(exp_ki) == 2 else
                                                                                          tuple(exp_ki.keys()))
            opts4ana['decomp_options']['etype'] = None if len(exp_ki) == 2 else tuple(f"decomp_{x}" for x in exp_ki)
            opts4ana['energy_options']['mol'] = opts4ana['entropy_options']['mol'] = (['complex'] if stability else
                tuple(options.get('energy').get('mols')))
            opts4ana['decomp_options']['mol'] = (['complex'] if stability else
                                                        tuple(options.get('decomposition').get('mols')))
            d = {MMPBSA_API: {'object': 'class',
                              'args': None},
                 'setting_time': {'object': 'method',
                                  'args': options.get('frames2time', {})},
                 'load_file': {'object': 'method',
                               'args': dict(fname=path)},
                 'get_ana_data': {'object': 'method',
                                  'args': opts4ana}
                 }
            iq.put((i, d))

            if options_file == 'User-Default':
                config = 'User-Default'
            elif options_file == 'Custom':
                config = path.parent
            else:
                config = None

            self.systems[i] = {'name': sys_name, 'path': path.parent,
                               'chart_options': ChartSettings(config),
                               'options': options,
                               'anaoptions': opts4ana,
                               'items_data': {},
                               'exp_ki': exp_ki,
                               'correlation': corr,
                               'reference': ref
                               }

            queue.task_done()
        pbd = ProcessingProgressBar(self, iq, self.rq, options['performance']['jobs'], 'Reading files...')
        pbd.rejected.connect(lambda x: print(f'rejected> exit code: {x}'))
        pbd.accepted.connect(lambda: self.process_data(self.rq))
        pbd.exec()

    def process_data(self, results: Queue):
        size = results.qsize()
        maximum = size * 2
        qpd = QProgressDialog('Creating systems tree', 'Abort', 0, maximum + 1, self)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(0)

        # check if all systems have the same frames range
        frange_base = []
        c = 0
        for i, _ in enumerate(range(size), start=1):
            qpd.setValue(c)
            if qpd.wasCanceled():
                break
            sys_id, api, data = results.get()
            self.systems[sys_id].update({
                'api': api,
                'namespace': api.app_namespace,
                'map': {x:v['map'] for x, v in data.items() if v and x != 'correlation'},
                'current_frames': dict(startframe=api.app_namespace.INPUT['startframe'],
                                       endframe=api.app_namespace.INPUT['endframe'],
                                       interval=api.app_namespace.INPUT['interval']),
                'current_nmode_frames': dict(nmstartframe=api.app_namespace.INPUT['nmstartframe'],
                                             nmendframe=api.app_namespace.INPUT['nmendframe'],
                                             nminterval=api.app_namespace.INPUT['nminterval']),
                'current_ie_segment': api.app_namespace.INPUT['ie_segment'],
                'items_data': {k:v1 for x, v in data.items() if v and x != 'correlation' for k, v1 in v['keys'].items()},
                'correlation_data': data.get('correlation')  # FIXME: fail when decomp
                # 'items_summary': summary,
                # 'exp_ki': exp_ki
            })
            c += 1

            if not self.all_systems_active:
                continue
            if not frange_base:
                frange_base = [api.app_namespace.INPUT['startframe'], api.app_namespace.INPUT['endframe'],
                               api.app_namespace.INPUT['interval'], api.app_namespace.INPUT['nmstartframe'],
                               api.app_namespace.INPUT['nmendframe'], api.app_namespace.INPUT['nminterval']]
                continue
            if frange_base != [api.app_namespace.INPUT['startframe'], api.app_namespace.INPUT['endframe'],
                               api.app_namespace.INPUT['interval'], api.app_namespace.INPUT['nmstartframe'],
                               api.app_namespace.INPUT['nmendframe'], api.app_namespace.INPUT['nminterval']]:
                self.all_systems_active = False

        for s in sorted(self.systems.keys()):
            qpd.setValue(c)
            if qpd.wasCanceled():
                break
            self.makeTree(s)
            c += 1


        self._initialize_systems()
        qpd.setValue(c)

        if not self.options.get('correlation')['corr']:
            self.correlation_DockWidget.setEnabled(False)
            self.correlation_DockWidget.hide()
        else:
            self.make_correlation()
        qpd.setValue(maximum + 1)

        # some late signal/slot connections
        # self.treeWidget.itemChanged.connect(self.showdata)
        # self.correlation_treeWidget.itemChanged.connect(self.showcorr)

    def makeTree(self, sys_index):

        # make system item
        sys_name = self.systems[sys_index]['name']
        sys_item = CustomItem(self.treeWidget, [sys_name], app=self, system_index=sys_index,
                              buttons=(-1,))
        sys_item.setExpanded(True)

        for c in [0, 1]:
            sys_item.setBackground(c, QBrush(QColor(100, 100, 100)))
            sys_item.setForeground(c, QBrush(QColor(220, 220, 255)))
            font = sys_item.font(c)
            font.setWeight(QFont.Weight.DemiBold)
            font.setPointSize(font.pointSize() + 1)
            sys_item.setFont(c, font)

        s = {'enthalpy': 'ΔH', 'entropy': '-TΔS','binding': 'ΔG'}

        for level, v in self.systems[sys_index]['map'].items():
            t = f' ({s[level]})' if s.get(level) else ''
            item = CustomItem(sys_item, [f'{level.capitalize()}{t}'], system_index=sys_index)
            for level1, v1 in v.items():

                item1 = CustomItem(item, [f'{level1.capitalize()}'], system_index=sys_index)
                item1.setExpanded(True)
                if level == 'enthalpy':
                    self._make_enthalpy_nmodeqh_items(sys_name, sys_index, item1, level1, v1)
                elif level == 'entropy':
                    self._make_enthalpy_nmodeqh_items(sys_name, sys_index, item1, level1, v1, 'Entropic')
                    self._make_iec2_items(sys_name, sys_index, item1, level1, v1)
                elif level == 'binding':
                    self._make_binding_items(sys_name, sys_index, item1, level1, v1)
                elif level == 'decomposition':
                    self._make_decomp_items(sys_name, sys_index, item1, level1, v1)

    def _make_decomp_items(self, sys_name, sys_index, item1, level1, v1):
        namespace = self.systems[sys_index]['namespace']
        # model
        for level2, v2 in v1.items():
            item2 = CustomItem(item1, [f'{level2.upper()}'], system_index=sys_index)
            item2.setExpanded(True)
            # mols
            for level3, v3 in v2.items():
                item3 = CustomItem(item2, [f'{level3.capitalize()}'], system_index=sys_index)
                item3.setExpanded(True)

                # TDC, BDC, SDC
                for level4, v4 in v3.items():
                    title = '[Per-residue]' if namespace.INPUT['idecomp'] in [1, 2] else '[Per-wise]'
                    item4 = CustomItem(item3,
                                       [level4.upper()],
                                       app=self,
                                       buttons=(1, 2, 3, 4),
                                       title=f"Energetic Components {title}",
                                       subtitle=f"{sys_name} | {level1.upper()} | "
                                                f"{level2.upper()} | {level3.capitalize()} | {level4.upper()}",
                                       system_index=sys_index,
                                       keys_path=(level1, level2, level3, level4),
                                       # part=part
                                       )
                    self.treeWidget.setItemWidget(item4, 1, item4.setup_buttons())
                    btns = (2,) if namespace.INPUT['idecomp'] in [1, 2] else (1, 2, 3)
                    for level5, v5 in v4.items():
                        item5 = CustomItem(item4,
                                           [level5.upper()],
                                           app=self,
                                           buttons=btns,
                                           title=f"Energetic Components {title}",
                                           subtitle=f"{sys_name} | "
                                                    f"{str(level1).upper()} | "
                                                    f"{str(level2).upper()} | "
                                                    f"{str(level3).upper()} | "
                                                    f"{str(level4).upper()} | "
                                                    f"{str(level5).upper()} | ",
                                           system_index=sys_index,
                                           keys_path=(level1, level2, level3, level4, level5)
                                           )
                        self.treeWidget.setItemWidget(item5, 1, item5.setup_buttons())
                        if namespace.INPUT['idecomp'] in [1, 2]:
                            for level6 in v5:
                                item6 = CustomItem(item5,
                                                   [level6.upper()],
                                                   app=self,
                                                   buttons=(1,),
                                                   title="Energetic Components [Per-residue]",
                                                   subtitle=f"{sys_name} | "
                                                            f"{str(level1).upper()} | "
                                                            f"{str(level2).upper()} | "
                                                            f"{str(level3).upper()} | "
                                                            f"{str(level4).upper()} | "
                                                            f"{str(level5).upper()} | "
                                                            f"{str(level6).upper()}",
                                                   system_index=sys_index,
                                                   keys_path=(level1, level2, level3, level4, level5, level6)
                                                   )
                                self.treeWidget.setItemWidget(item6, 1, item6.setup_buttons())
                        else:
                            for level6, v6 in v5.items():
                                item6 = CustomItem(item5,
                                                   [level6.upper()],
                                                   app=self,
                                                   buttons=(2,),
                                                   title="Energetic Components [Per-wise]",
                                                   subtitle=f"{sys_name}| "
                                                            f"{str(level1).upper()} | "
                                                            f"{str(level2).upper()} | "
                                                            f"{str(level3).upper()} | "
                                                            f"{str(level4).upper()} | "
                                                            f"{str(level5).upper()} | "
                                                            f"{str(level6).upper()}",
                                                   system_index=sys_index,
                                                   keys_path=(level1, level2, level3, level4, level5, level6)
                                                   )
                                self.treeWidget.setItemWidget(item6, 1, item6.setup_buttons())

    def _make_binding_items(self, sys_name, sys_index, item1, level1, v1):
        for level2, v2 in v1.items():
            for level3 in v2:
                item3 = CustomItem(item1, [f'{level2.upper()}+{level3.upper()}'],
                                   app=self,
                                   buttons=(2, -2),
                                   title="Binding Energy",
                                   subtitle=f"{sys_name} | {level1.capitalize()} | {level2.upper()}+{level3.upper()}",
                                   system_index=sys_index,
                                   keys_path=(level1, level2, level3),
                                   # part=part
                                   )
                item3.setExpanded(True)
                self.treeWidget.setItemWidget(item3, 1, item3.setup_buttons())

    def _make_enthalpy_nmodeqh_items(self, sys_name, sys_index, item1, level1, v1, comp_type='Energetic'):
        for level2, v2 in v1.items():
            if level2 in ['c2', 'ie']:
                continue
            item2 = CustomItem(item1, [f'{level2.upper()}'], system_index=sys_index)
            item2.setExpanded(True)
            for level3, v3 in v2.items():
                item3 = CustomItem(item2, [f'{level3.capitalize()}'],
                                   app=self,
                                   buttons=(2, -2),
                                   title=f"{comp_type} Components",
                                   subtitle=f"{sys_name} | {level1.capitalize()} | {level2.upper()} | "
                                            f"{level3.capitalize()}",
                                   system_index=sys_index,
                                   keys_path=(level1, level2, level3),
                                   # part=part
                                   )
                item3.setExpanded(True)
                self.treeWidget.setItemWidget(item3, 1, item3.setup_buttons())
                if level2 == 'qh':
                    continue
                for level4 in v3:
                    item4 = CustomItem(item3, [f'{level4.upper()}'],
                                       app=self,
                                       buttons=(1,),
                                       title=f"{comp_type} Components",
                                       subtitle=f"{sys_name} | {level1.capitalize()} | {level2.upper()} | "
                                                f"{level3.capitalize()} | {level4.upper()}",
                                       system_index=sys_index,
                                       keys_path=(level1, level2, level3, level4),
                                       # part=part
                                       )
                    item4.setExpanded(True)
                    self.treeWidget.setItemWidget(item4, 1, item4.setup_buttons())

    def _make_iec2_items(self, sys_name, sys_index, item1, level1, v1):
        for level2, v2 in v1.items():
            if level2 == 'ie':
                titem = CustomItem(item1,
                                   [level2.upper()],
                                   app=self,
                                   buttons=(1, 2),
                                   title="Interaction Entropy",
                                   subtitle=f"{sys_name} | {level1.capitalize()} | {level2.upper()}",
                                   system_index=sys_index,
                                   keys_path=(level1, level2),
                                   # part=part
                                   )
                self.treeWidget.setItemWidget(titem, 1, titem.setup_buttons())
            elif level2 == 'c2':
                titem = CustomItem(item1,
                                   [level2.upper()],
                                   app=self,
                                   buttons=(2,),
                                   title="C2 Entropy",
                                   subtitle=f"{sys_name} | {level1.capitalize()} | {level2.upper()}",
                                   system_index=sys_index,
                                   keys_path=(level1, level2),
                                   # part=part
                                   )
                self.treeWidget.setItemWidget(titem, 1, titem.setup_buttons())
                