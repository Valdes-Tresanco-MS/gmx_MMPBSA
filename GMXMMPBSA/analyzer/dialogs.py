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
import h5py
from PyQt5.QtWidgets import (QDialog, QSpinBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QCheckBox,
                             QGroupBox, QButtonGroup, QGridLayout, QTreeWidget, QTreeWidgetItem, QRadioButton,
                             QTreeWidgetItemIterator, QHeaderView, QProgressBar, QStatusBar, QMessageBox, QComboBox,
                             QDialogButtonBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, pyqtSlot
from queue import Queue, Empty
from GMXMMPBSA import API
from pathlib import Path

from GMXMMPBSA.analyzer.chartsettings import ChartSettings
from GMXMMPBSA.analyzer.utils import worker, ncpu


class InitDialog(QDialog):
    def __init__(self, parent=None):
        super(InitDialog, self).__init__(parent)
        self.parent = parent
        self.setWindowModality(Qt.WindowModality.WindowModal)
        self.setWindowTitle('Initialization gmx_MMPBSA_ana')
        self.setMinimumWidth(650)
        self.curr_progress = 0
        self.data = []
        self.chart_default_setting = ChartSettings()
        self.processing_label = QLabel()

        self.show_decomp_btn = QCheckBox('Show decomposition data')
        self.show_decomp_btn.setToolTip('Defines whether the decomposition graphs will be available during the '
                                        'session. Useful when defining a system for correlation analysis.')
        self.show_decomp_btn.setChecked(True)
        self.remove_empty_charts_btn = QCheckBox('Remove empty charts')
        self.remove_empty_charts_btn.setToolTip('Does not show charts in which all trajectory is zero')
        self.remove_empty_charts_btn.setChecked(True)
        self.remove_empty_charts_btn.clicked.connect(self.show_warn)

        self.remove_empty_terms_btn = QCheckBox('Remove empty terms')
        self.remove_empty_terms_btn.setToolTip('Does not show null terms in graphics')
        self.remove_empty_terms_btn.setChecked(True)
        self.remove_empty_terms_btn.clicked.connect(self.show_warn)

        self.check_l = QHBoxLayout()
        self.check_l.addWidget(self.show_decomp_btn)
        self.check_l.addWidget(self.remove_empty_charts_btn)
        self.check_l.addWidget(self.remove_empty_terms_btn)
        self.warn_label_empty = QLabel('Note: Remove empty charts and terms options will hide the charts and empty '
                                       'terms respectively. This does not change your results at all, it only changes '
                                       'the representation of the data.')
        self.warn_label_empty.setStyleSheet("border:3px solid green")
        self.warn_label_empty.setWordWrap(True)

        self.corr_sys_btn = QCheckBox('Calculate correlation between systems')
        self.corr_sys_btn.setToolTip('Make correlation between systems. Only works if you define more than 3 systems')
        self.corr_sys_btn.setChecked(False)

        # mutants correlation
        self.corr_mut_btn = QCheckBox('Calculate correlation between mutants')
        self.corr_mut_btn.setToolTip('Make correlation between mutants systems. Only works if you define more than 3 '
                                     'mutants')
        self.corr_mut_btn.setChecked(False)

        self.corr_group = QGroupBox('Correlation')
        self.corr_group_layout = QGridLayout(self.corr_group)
        self.corr_group_layout.addWidget(self.corr_sys_btn, 0, 0)
        self.corr_group_layout.addWidget(self.corr_mut_btn, 0, 1)

        self.delta_btn = QCheckBox('delta')
        self.delta_btn.setChecked(True)
        self.delta_btn.setEnabled(False)
        self.com_btn = QCheckBox('complex')
        self.com_btn.clicked.connect(self.show_warn)
        self.rec_btn = QCheckBox('receptor')
        self.rec_btn.clicked.connect(self.show_warn)
        self.lig_btn = QCheckBox('ligand')
        self.lig_btn.clicked.connect(self.show_warn)
        self.warn_label_big = QLabel('Warning: This can generate a huge amount of graphics which can crash the '
                                     'application!. So far, we have not registered this issue, but be cautious')
        self.warn_label_big.setStyleSheet("border:3px solid orange")
        self.warn_label_big.setWordWrap(True)
        self.warn_label_big.hide()

        self.comp_group = QGroupBox('Computation of Components Charts')
        self.comp_group_layout = QGridLayout(self.comp_group)
        self.comp_group_layout.addWidget(self.delta_btn, 0, 0)
        self.comp_group_layout.addWidget(self.com_btn, 0, 1)
        self.comp_group_layout.addWidget(self.rec_btn, 0, 2)
        self.comp_group_layout.addWidget(self.lig_btn, 0, 3)
        self.comp_group_layout.addWidget(self.warn_label_big, 1, 0, 1, 4)

        self.hide_tb_btn = QCheckBox('Hide ToolBar')
        self.hide_tb_btn.setChecked(True)
        self.reset_default_config = QCheckBox('Reset Default Configuration')
        self.reset_default_config.clicked.connect(self.show_warn)

        self.other_options = QHBoxLayout()
        self.chart_group = QGroupBox('Charts options')
        self.other_options.addWidget(self.chart_group)
        self.chart_group_layout = QHBoxLayout(self.chart_group)
        self.chart_group_layout.addWidget(self.hide_tb_btn)
        self.chart_group_layout.addWidget(self.reset_default_config)

        self.frame2time = QGroupBox('Convert frames to time')
        self.frame2time.setCheckable(True)
        self.frame2time.setChecked(False)
        self.frame2time_layout = QHBoxLayout(self.frame2time)
        self.time_start_label = QLabel('Start:')
        self.frame2time_layout.addWidget(self.time_start_label)
        self.time_start = QSpinBox()
        self.time_start.setRange(0, 10000)
        self.time_start.setAccelerated(True)
        self.frame2time_layout.addWidget(self.time_start)

        self.time_step_label = QLabel('Scale:')
        self.frame2time_layout.addWidget(self.time_step_label)
        self.time_step = QSpinBox()
        self.time_step.setRange(1, 1000)
        self.time_step.setValue(10)
        self.time_step.setAccelerated(True)
        self.frame2time_layout.addWidget(self.time_step)
        # self.frame2time_layout.addStretch(1)
        self.time_unit_label = QLabel('Unit:')
        self.frame2time_layout.addWidget(self.time_unit_label)
        self.time_unit = QComboBox()
        self.time_unit.addItems(['ps', 'ns'])
        self.frame2time_layout.addWidget(self.time_unit)

        self.other_options.addWidget(self.frame2time)

        self.warn_label_config = QLabel(f'The global default graphics settings will be stored in:\n'
                                        f'{self.chart_default_setting.filename.as_posix()}')
        self.warn_label_config.setStyleSheet("border:3px solid orange")
        self.warn_label_config.setWordWrap(True)
        self.warn_label_config.show()
        if self.chart_default_setting.config_created():
            self.warn_label_config.hide()



        self.sys_group = QGroupBox('Systems options')
        self.sys_group_layout = QVBoxLayout(self.sys_group)
        self.sys_group_layout.addLayout(self.check_l)
        self.sys_group_layout.addWidget(self.warn_label_empty)
        self.sys_group_layout.addWidget(self.corr_group)
        self.sys_group_layout.addWidget(self.comp_group)

        self.header_item = QTreeWidgetItem(['Folder name', 'Select', 'Name', 'Exp.Ki (nM)', 'Chart Settings', 'Path'])
        self.header_item.setToolTip(0, 'Container')
        self.header_item.setToolTip(1, 'Name')
        self.header_item.setTextAlignment(0, Qt.AlignmentFlag.AlignCenter)
        self.header_item.setTextAlignment(1, Qt.AlignmentFlag.AlignCenter)
        self.header_item.setTextAlignment(2, Qt.AlignmentFlag.AlignCenter)
        self.header_item.setTextAlignment(3, Qt.AlignmentFlag.AlignCenter)
        self.header_item.setTextAlignment(4, Qt.AlignmentFlag.AlignCenter)
        self.result_tree = QTreeWidget(self)
        self.result_tree.setHeaderItem(self.header_item)
        self.result_tree.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.result_tree.itemChanged.connect(self.update_item_info)
        self.pb = QProgressBar()
        # self.pb.setRange(0, 0)
        # self.save_btn.clicked.connect(self.save)

        btnbox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        accept_btn = btnbox.button(QDialogButtonBox.Ok)
        accept_btn.setText('Accept')
        accept_btn.setDefault(True)
        accept_btn.clicked.connect(self.get_data)

        cancel_btn = btnbox.button(QDialogButtonBox.Cancel)
        cancel_btn.setDefault(False)
        cancel_btn.clicked.connect(self.close)

        self.btn_layout = QHBoxLayout()
        self.btn_layout.addWidget(self.pb, 10)
        self.btn_layout.addStretch(1)
        self.btn_layout.addWidget(btnbox)


        self.statusbar = QStatusBar(self)

        self.content_layout = QVBoxLayout(self)
        self.content_layout.addWidget(self.processing_label)
        self.content_layout.addWidget(self.sys_group)
        self.content_layout.addLayout(self.other_options)
        self.content_layout.addWidget(self.result_tree)
        self.content_layout.addWidget(self.warn_label_config)
        self.content_layout.addWidget(self.statusbar)
        self.content_layout.addLayout(self.btn_layout)

    def show_warn(self):
        if self.com_btn.isChecked() or self.rec_btn.isChecked() or self.lig_btn.isChecked():
            self.warn_label_big.show()
        else:
            self.warn_label_big.hide()

        if self.remove_empty_terms_btn.isChecked() or self.remove_empty_charts_btn.isChecked():
            self.warn_label_empty.show()
        else:
            self.warn_label_empty.hide()

        if not self.chart_default_setting.config_created() or self.reset_default_config.isChecked():
            self.warn_label_config.show()
        else:
            self.warn_label_config.hide()



    def update_item_info(self, item: QTreeWidgetItem, col):
        if item.text(0) == 'All':
            return
        # path = item.info[1]
        # basename = item.text(2)
        # exp_ki = float(item.text(3))
        # custom_sett = float(item.checkState(4))
        # item.info = [basename, path, exp_ki, custom_sett]

    def get_files_info(self, info_files):
        self.nfiles = len(info_files)
        self.f_item = QTreeWidgetItem(['All'])
        self.f_item.setCheckState(1, Qt.CheckState.Checked)
        self.f_item.info = None
        self.f_item.setFlags(self.f_item.flags() | Qt.ItemFlag.ItemIsAutoTristate)
        self.result_tree.addTopLevelItem(self.f_item)

        self.processing_label.setText(f'Processing {self.nfiles} systems...')

        files_list = info_files
        names = []
        for c, fname in enumerate(files_list, start=1):
            basename = None
            exp_ki = None
            mut_only = False
            mutant = None
            if fname.suffix == '.h5':
                with h5py.File(fname) as fi:
                    basename = fi['INPUT']['sys_name'][()].decode('utf-8')
                    if basename in names:
                        while basename in names:
                            basename = f"{basename}-{names.count(basename) + 1}"
                        names.append(basename)
                    exp_ki = float(fi['INPUT']['exp_ki'][()])
                    mut_only = int(fi['INPUT']['mutant_only'][()])
                    mutant = fi['INFO']['mut_str'][()].decode('utf-8')
            else:
                with open(fname) as fi:
                    for line in fi:
                        if line.startswith("INPUT['sys_name']"):
                            basename = str(line.split()[2]).strip('"\'')
                            if basename in names:
                                while basename in names:
                                    basename = f"{basename}-{names.count(basename) + 1}"
                                names.append(basename)
                        if line.startswith("INPUT['exp_ki']"):
                            exp_ki = float(line.split()[2])
                        if line.startswith("INPUT['mutant_only']"):
                            mut_only = int(line.split()[2])
                        if line.startswith("mut_str"):
                            mutant = line.strip('\n').split('=')[1].strip(" '")
            # check for custom settings
            custom_settings = fname.parent.joinpath('settings.json').exists()

            cb = QComboBox()
            if custom_settings:
                cb.addItem('Custom')
            cb.addItem('Default')

            if not basename:
                basename = f'System-{c}'
            if not exp_ki:
                exp_ki = 0.0

            item = QTreeWidgetItem([f'{fname.parent.name}', '', f'{basename}', '', '', f'{fname.parent.absolute()}'])
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsAutoTristate)

            self._set_item_properties(item)
            self.f_item.addChild(item)

            if not mut_only:
                witem = QTreeWidgetItem(['', '', 'normal', f'{exp_ki}', '', ''])
                # witem.info = [basename, Path(fname)]
                self._set_item_properties(witem)
                item.addChild(witem)

            if mutant:
                mitem = QTreeWidgetItem(['', '', f'{mutant}', f'{exp_ki}', '', ''])
                # mitem.info = [basename, Path(fname)]
                self._set_item_properties(mitem)
                item.addChild(mitem)

            item.info = [basename, Path(fname)]

            item.setExpanded(True)
            self.result_tree.setItemWidget(item, 4, cb)
        self.f_item.setExpanded(True)

    def _set_item_properties(self, item: QTreeWidgetItem):
        item.setCheckState(1, Qt.CheckState.Checked)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
        item.setTextAlignment(2, Qt.AlignmentFlag.AlignCenter)
        item.setTextAlignment(3, Qt.AlignmentFlag.AlignRight)
        item.setTextAlignment(4, Qt.AlignmentFlag.AlignCenter)

    def get_data(self):

        self.pb.setValue(0)
        queue = Queue()
        self.result_queue = Queue()
        self.systems_list = []
        self.options = {'corr_sys': self.corr_sys_btn.isChecked(),
                        'corr_mut': self.corr_mut_btn.isChecked(),
                        'decomposition': self.show_decomp_btn.isChecked(),
                        'components': [x.text() for x in [self.com_btn, self.rec_btn, self.lig_btn] if x.isChecked()],
                        'remove_empty_charts': self.remove_empty_charts_btn.isChecked(),
                        'remove_empty_terms': self.remove_empty_terms_btn.isChecked(),
                        'timestep': self.time_step.value() if self.frame2time.isChecked() else 0,
                        'timeunit': self.time_unit.currentText(),
                        'timestart': self.time_start.value()
                        # 'default_chart_options': self.default_settings_btn.isChecked()
                        }
        counter = 0
        for c in range(self.f_item.childCount()):
            child = self.f_item.child(c)
            if not child.checkState(1):
                continue
            t = {}
            if child.childCount():
                for c1 in range(child.childCount()):
                    if child.child(c1).checkState(1) == Qt.CheckState.Checked:
                        if child.child(c1).text(2) == 'normal':
                            t['normal'] = float(child.child(c1).text(3))
                        else:
                            t['mutant'] = float(child.child(c1).text(3))
            child.info += [t, self.result_tree.itemWidget(child, 4).currentText()]
            queue.put(child.info)
            counter += 1
        if not counter:
            QMessageBox.critical(self, 'Error processing systems', 'You must select at least one system.',
                                 QMessageBox.Ok)
            return
        self.pb.setRange(0, counter)
        # Store the global default graphics settings if not exits
        if not self.chart_default_setting.config_created() or self.reset_default_config.isChecked():
            self.chart_default_setting.write_system_config()

        self.parent.read_data(queue, self.options)
