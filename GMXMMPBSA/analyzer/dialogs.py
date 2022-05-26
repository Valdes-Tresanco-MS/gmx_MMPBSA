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
    from PyQt6.QtGui import *


import h5py
from queue import Queue, Empty
from pathlib import Path
import multiprocessing
from GMXMMPBSA.analyzer.chartsettings import ChartSettings


class InitDialog(QDialog):
    def __init__(self, parent):
        super(InitDialog, self).__init__(parent)
        self.parent = parent
        self.setWindowModality(Qt.WindowModality.WindowModal)
        self.setWindowTitle('Initialization gmx_MMPBSA_ana')
        self.setMinimumWidth(650)
        self.curr_progress = 0
        self.data = []
        self.chart_default_setting = ChartSettings()
        self.processing_label = QLabel()



        self.energy_group = QGroupBox('Energy Options')
        self.energy_group_layout = QGridLayout(self.energy_group)

        self.energy_load_com = QCheckBox('Load COM data')
        self.energy_load_com.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load COM data'))
        self.energy_group_layout.addWidget(self.energy_load_com, 0, 0)
        self.energy_load_rec = QCheckBox('Load REC data')
        self.energy_load_rec.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load REC data'))
        self.energy_group_layout.addWidget(self.energy_load_rec, 0, 1)
        self.energy_load_lig = QCheckBox('Load LIG data')
        self.energy_load_lig.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load LIG data'))
        self.energy_group_layout.addWidget(self.energy_load_lig, 0, 2)

        self.remove_empty_terms_btn = QCheckBox('Remove empty terms')
        self.remove_empty_terms_btn.setToolTip('Does not show null terms in graphics')
        self.remove_empty_terms_btn.setChecked(True)
        self.remove_empty_terms_btn.toggled.connect(lambda x: self._update_ram_indicator(x, 'Remove empty terms'))
        self.energy_group_layout.addWidget(self.remove_empty_terms_btn, 1, 0)

        self.empty_threshold = QDoubleSpinBox()
        self.empty_threshold.setRange(0.0001, 1.0)
        self.empty_threshold.setValue(0.0100)
        self.empty_threshold.setSingleStep(0.0001)
        self.empty_threshold.setDecimals(4)
        self.empty_threshold.setAccelerated(True)
        self.empty_threshold.setSuffix(' kcal/mol')

        self.empty_threshold_layout = QFormLayout()
        self.empty_threshold_layout.addRow('Term threshold:', self.empty_threshold)
        self.energy_group_layout.addLayout(self.empty_threshold_layout, 1, 1)

        self.decomp_group = QGroupBox('Decomposition Options')
        self.decomp_group.setCheckable(True)
        self.decomp_group.setChecked(True)
        self.decomp_group.toggled.connect(lambda x: self._update_ram_indicator(x, 'Decomposition Options'))

        self.decomp_group_layout = QGridLayout(self.decomp_group)

        self.decomp_load_com = QCheckBox('Load COM data')
        self.decomp_load_com.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load decomp COM data'))
        self.decomp_group_layout.addWidget(self.decomp_load_com, 0, 0)
        self.decomp_load_rec = QCheckBox('Load REC data')
        self.decomp_load_rec.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load decomp REC data'))
        self.decomp_group_layout.addWidget(self.decomp_load_rec, 0, 1)
        self.decomp_load_lig = QCheckBox('Load LIG data')
        self.decomp_load_lig.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load decomp LIG data'))
        self.decomp_group_layout.addWidget(self.decomp_load_lig, 0, 2)

        self.decomp_non_contrib_res = QCheckBox('Remove non-contributing residues')
        self.decomp_non_contrib_res.setChecked(True)
        self.decomp_non_contrib_res.toggled.connect(lambda x:
                                                    self._update_ram_indicator(x, 'Remove non-contributing residues'))
        self.decomp_group_layout.addWidget(self.decomp_non_contrib_res, 1, 0)
        self.decomp_res_threshold = QDoubleSpinBox()
        self.decomp_res_threshold.setRange(0.01, 2.0)
        self.decomp_res_threshold.setSingleStep(0.01)
        self.decomp_res_threshold.setDecimals(2)
        self.decomp_res_threshold.setValue(0.5)
        self.decomp_res_threshold.setSuffix(' kcal/mol')
        self.decomp_res_threshold_layout = QFormLayout()
        self.decomp_res_threshold_layout.addRow('Per-residue threshold:', self.decomp_res_threshold)
        self.decomp_group_layout.addLayout(self.decomp_res_threshold_layout, 1, 1)

        self.corr_group = QGroupBox('Correlation')
        self.corr_group_layout = QGridLayout(self.corr_group)
        self.corr_sys_btn = QCheckBox('Between systems')
        self.corr_sys_btn.toggled.connect(lambda x: self._update_ram_indicator(x, 'correlation'))

        self.corr_sys_btn.setToolTip('Make correlation between systems. Only works when defining more than 3 systems')
        self.corr_group_layout.addWidget(self.corr_sys_btn, 0, 0)

        self.energy_type = QComboBox()
        self.energy_type.addItems(['ΔG', 'ΔΔG'])

        self.energy_corr_layout = QFormLayout()
        self.energy_corr_layout.addRow('Energy:', self.energy_type)
        self.corr_group_layout.addLayout(self.energy_corr_layout, 0, 1)

        self.chart_group = QGroupBox('Charts options')
        self.chart_group_layout = QGridLayout(self.chart_group)
        self.hide_tb_btn = QCheckBox('Hide ToolBar')
        self.hide_tb_btn.setChecked(True)
        self.chart_group_layout.addWidget(self.hide_tb_btn, 0, 0)

        self.reset_default_config = QCheckBox('Reset Configuration')
        self.reset_default_config.setToolTip('Reset the current configuration')
        self.reset_default_config.clicked.connect(self._reset_config)
        self.chart_group_layout.addWidget(self.reset_default_config, 0, 1)

        self.frame2time = QGroupBox('Convert frames to time')
        self.frame2time.setCheckable(True)
        self.frame2time.setChecked(False)
        self.frame2time_layout = QFormLayout(self.frame2time)
        self.time_start = QSpinBox()
        self.time_start.setRange(0, 1000000)
        self.time_start.setAccelerated(True)
        self.frame2time_layout.addRow('Start:', self.time_start)
        self.time_step = QSpinBox()
        self.time_step.setRange(1, 10000)
        self.time_step.setValue(10)
        self.time_step.setAccelerated(True)
        self.frame2time_layout.addRow('Scale:', self.time_step)
        self.time_unit = QComboBox()
        self.time_unit.addItems(['ps', 'ns'])
        self.frame2time_layout.addRow('Unit:', self.time_unit)

        self.performance_group = QGroupBox('Performance Options')
        self.performance_group_layout = QGridLayout(self.performance_group)

        self.performance_slider = QSlider(Qt.Orientation.Horizontal)
        self.performance_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.performance_slider.setTickInterval(1)
        self.performance_slider.setRange(1, 4)

        self.basic_sett_group = QGroupBox('Basic settings')
        self.basic_sett_group.toggled.connect(self.configure_basic_settings)
        self.slider_layout = QGridLayout(self.basic_sett_group)
        self.basic_sett_group.setCheckable(True)
        self.basic_sett_group.setChecked(True)
        self.slider_layout.addWidget(self.performance_slider, 0, 0, 1, 3)
        self.slider_layout.addWidget(QLabel('High'), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.slider_layout.addWidget(QLabel('Low'), 1, 0)
        self.performance_group_layout.addWidget(self.basic_sett_group, 0, 0, 1, 1)

        self.advanced_group = QGroupBox('Advanced settings')
        self.performance_group_layout.addWidget(self.advanced_group, 0, 1, 1, 2)
        self.advanced_group.setCheckable(True)
        self.advanced_group.setChecked(False)
        self.advanced_group.toggled.connect(lambda x: self.basic_sett_group.setChecked(not x))
        self.advanced_group_layout = QGridLayout(self.advanced_group)

        self.energy_memory = QCheckBox('Load energy in memory')
        self.energy_memory.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load energy in memory'))
        self.advanced_group_layout.addWidget(self.energy_memory, 0, 0)

        self.decomp_memory = QCheckBox('Load decomp in memory')
        self.decomp_memory.toggled.connect(lambda x: self._update_ram_indicator(x, 'Load decomp in memory'))
        self.advanced_group_layout.addWidget(self.decomp_memory, 0, 1)

        self.multiprocessing = QCheckBox('Multiprocessing')
        self.advanced_group_layout.addWidget(self.multiprocessing, 1, 1)
        self.precomp_charts = QCheckBox('Pre-compute charts')
        self.advanced_group_layout.addWidget(self.precomp_charts, 1, 0)
        self.remove_charts = QCheckBox('Clean up charts data')
        self.advanced_group_layout.addWidget(self.remove_charts, 0, 2)

        self.njobs = QSpinBox()
        self.njobs.setRange(1, multiprocessing.cpu_count())
        self.njobs.setValue(2)
        self.njobs_layout = QFormLayout()
        self.njobs_layout.addRow('Jobs:', self.njobs)
        self.advanced_group_layout.addLayout(self.njobs_layout, 1, 2)

        self.ram_pb = CustomPB(self)

        self.performance_slider.valueChanged.connect(self.slider_changed)
        self.performance_slider.setValue(2)

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
        self.result_tree.header().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.result_tree.itemChanged.connect(self.update_item_info)
        self.pb = QProgressBar()
        # self.pb.setRange(0, 0)
        # self.save_btn.clicked.connect(self.save)

        btnbox = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        accept_btn = btnbox.button(QDialogButtonBox.StandardButton.Ok)
        accept_btn.setText('Accept')
        accept_btn.setDefault(True)
        accept_btn.clicked.connect(self.get_data)

        cancel_btn = btnbox.button(QDialogButtonBox.StandardButton.Cancel)
        cancel_btn.setDefault(False)
        cancel_btn.clicked.connect(self.reject)


        self.statusbar = QStatusBar(self)

        self.content_layout = QGridLayout(self)
        self.content_layout.addWidget(self.energy_group, 0, 0, 1, 3)
        self.content_layout.addWidget(self.decomp_group, 1, 0, 1, 3)
        self.content_layout.addWidget(self.corr_group, 2, 0, 1, 2)
        self.content_layout.addWidget(self.chart_group, 3, 0, 1, 2)
        self.content_layout.addWidget(self.frame2time, 2, 2, 2, 1)
        self.content_layout.addWidget(self.performance_group, 4, 0, 1, 3)
        self.content_layout.addWidget(self.ram_pb, 5, 0, 1, 3)

        self.content_layout.addWidget(self.result_tree, 6, 0, 1, 3)
        self.content_layout.addWidget(self.statusbar, 7, 0, 1, 2)
        self.content_layout.addWidget(btnbox, 7, 2, 1, 1, Qt.AlignmentFlag.AlignRight)

    def _update_ram_indicator(self, check, obj_name):

        ram_vaues = {'Remove empty terms': [False, 1], 'Load COM data': [True, 3], 'Load REC data': [True, 3],
                     'Load LIG data': [True, 3], 'Load energy in memory': [True, 3], 'Decomposition Options': [True, 5],
                     'Load decomp COM data': [True, 5], 'Load decomp REC data': [True, 5],
                     'Load decomp LIG data': [True, 5], 'Remove non-contributing residues': [False, 2],
                     'correlation': [True, 1], 'Load decomp in memory': [True, 5]}

        dgroup = ['Decomposition Options', 'Load decomp COM data', 'Load decomp REC data', 'Load decomp LIG data',
                  'Remove non-contributing residues']
        egroup = ['Load COM data', 'Load REC data', 'Load LIG data', 'Remove empty terms']
        load_emem = 3 + sum(
            3 for x in [self.energy_load_com, self.energy_load_rec, self.energy_load_lig] if x.isChecked())
        if not self.remove_empty_terms_btn.isChecked():
            load_emem += 1

        load_dmem = 5 + sum(5 for x in [self.decomp_load_com, self.decomp_load_rec, self.decomp_load_lig]
                            if x.isChecked())
        if load_dmem and not self.decomp_non_contrib_res.isChecked():
            load_dmem += 2

        mult = int(ram_vaues[obj_name][0] == check) or -1
        if obj_name == 'Load energy in memory':
            self.ram_pb.setValue(self.ram_pb.value() + load_emem * mult)
        elif obj_name == 'Load decomp in memory' and self.decomp_group.isChecked():
            self.ram_pb.setValue(self.ram_pb.value() + load_dmem * mult)
        elif obj_name in egroup and self.energy_memory.isChecked():
            self.ram_pb.setValue(self.ram_pb.value() + ram_vaues[obj_name][1] * mult)
        elif obj_name == 'Decomposition Options' and self.decomp_memory.isChecked():
            self.ram_pb.setValue(self.ram_pb.value() + load_dmem * mult)
        elif obj_name in dgroup and self.decomp_memory.isChecked():
            self.ram_pb.setValue(self.ram_pb.value() + ram_vaues[obj_name][1] * mult)
        elif obj_name == 'correlation':
            self.ram_pb.setValue(self.ram_pb.value() + ram_vaues[obj_name][1] * mult)

    def _change_decomp_opts_color(self, check):
        if check:
            self.decomp_group.setStyleSheet('''
                                                    QGroupBox::title {
                                                       background: qlineargradient( x1:0 y1:0, x2:1 y2:0, stop:0 #ffa200, stop:1 #fedb6f);
                                            }''')
        else:
            self.decomp_group.setStyleSheet('')

    def configure_basic_settings(self, checked):
        self.advanced_group.setChecked(not checked)
        if checked:
            self.slider_changed(self.performance_slider.value())

    def slider_changed(self, value):
        self.energy_memory.setChecked(False)
        self.decomp_memory.setChecked(False)
        self.precomp_charts.setChecked(False)
        self.remove_charts.setChecked(False)
        self.multiprocessing.setChecked(False)
        self.njobs.setValue(1)

        if value == 1:
            self.remove_charts.setChecked(True)
        elif value == 2:
            self.energy_memory.setChecked(True)
            self.multiprocessing.setChecked(True)
            self.njobs.setValue(multiprocessing.cpu_count() // 4)
        elif value == 3:
            self.energy_memory.setChecked(True)
            self.precomp_charts.setChecked(True)
            self.multiprocessing.setChecked(True)
            self.njobs.setValue(multiprocessing.cpu_count() // 2)
        else:
            self.energy_memory.setChecked(True)
            self.decomp_memory.setChecked(True)
            self.precomp_charts.setChecked(True)
            self.multiprocessing.setChecked(True)
            self.njobs.setValue(multiprocessing.cpu_count())

    def _reset_config(self, check):
        if check:
            self.reset_default_config.setStyleSheet(
                '''
                QCheckBox:checked {
                    background: qlineargradient( x1:0 y1:0, x2:1 y2:0, stop:0 #ffa200, stop:1 #fedb6f);
                }''')
            self.reset_default_config.setToolTip(f'''
            <html>
                <spam>Reset the current configuration<spam>
            
                <h3><strong>Warning!</strong></h3>
                <p>The global default graphics settings will be stored in:</p>
                <p><strong>{self.chart_default_setting.filename.as_posix()}</strong></p>
               
            </html>
            ''')
        else:
            self.reset_default_config.setStyleSheet("")
            self.reset_default_config.setToolTip('Reset the current configuration')

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
            stability = False
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
                        if line.startswith("FILES.stability"):
                            stability = eval(line.split()[2])
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

            item.info = [basename, Path(fname), stability]

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
        emols = ['delta']
        emols += ['ligand'] if self.energy_load_lig.isChecked() else []
        emols += ['receptor'] if self.energy_load_rec.isChecked() else []
        emols += ['complex'] if self.energy_load_com.isChecked() else []
        emols.reverse()
        dmols = ['delta']
        dmols += ['ligand'] if self.decomp_load_lig.isChecked() else []
        dmols += ['receptor'] if self.decomp_load_rec.isChecked() else []
        dmols += ['complex'] if self.decomp_load_com.isChecked() else []
        dmols.reverse()
        self.options = {
            'energy': {
                'mols': emols,
                'remove_empty_terms': self.remove_empty_terms_btn.isChecked(),
                'threshold': self.empty_threshold.value()
            },
            'decomposition': {
                'mols': dmols,
                'remove_noncontrib_res': self.decomp_non_contrib_res.isChecked(),
                'res_threshold': self.decomp_res_threshold.value()
            },
            'correlation': {
                'normal': self.corr_sys_btn.isChecked(),
                'mutant': self.corr_mut_btn.isChecked()
            },
            'performance': {
            }
            # 'default_chart_options': self.default_settings_btn.isChecked()
        }

        if self.frame2time.isChecked():
            self.options['frames2time'] = {
                'timestart': self.time_start.value(),
                'timestep': self.time_step.value(),
                'timeunit': self.time_unit.currentText()
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
        self.accept()

        self.parent.read_data(queue, self.options)


class worker(QThread):
    job_finished = pyqtSignal()

    def __init__(self, queue: Queue, result_q: Queue, jobs: int = 1):
        super(worker, self).__init__()
        self.queue = queue
        self.result_queue = result_q
        self.jobs = jobs

    @staticmethod
    def run_process(input_args: dict):
        ident, struct_dict = input_args
        obj_class = None
        obj_function = None
        obj = None
        for k, v in struct_dict.items():
            print(k, v)
            if v['object'] == 'class':
                obj_class = k(**v['args']) if v['args'] else k()
            elif v['object'] == 'function':
                obj_function = k(**v['args']) if v['args'] else k()
            else:
                obj_method = getattr(obj_class, k)
                obj = obj_method(**v['args']) if v['args'] else obj_method()

        return ident, obj_function or obj_class, obj

    def run(self):
        size = self.queue.qsize()
        TASKS = []
        for _ in range(size):
            i, d = self.queue.get_nowait()
            TASKS.append((i, d))
        self.jobs = min(self.jobs, len(TASKS))
        with multiprocessing.Pool(self.jobs) as pool:
            imap_unordered_it = pool.imap_unordered(self.run_process, TASKS)
            for result in imap_unordered_it:
                self.job_finished.emit()
                self.result_queue.put(result)


class ProcessingProgressBar(QDialog):
    def __init__(self, parent, input_queue, output_queue, jobs, text):
        super(ProcessingProgressBar, self).__init__(parent)
        self.setWindowFlags(Qt.WindowType.Dialog|Qt.WindowType.FramelessWindowHint)
        self.output_queue = output_queue
        self.content_layout = QVBoxLayout(self)
        self.counter = 0
        self.action_title = QLabel(text)
        self.content_layout.addWidget(self.action_title)

        self.pb = QProgressBar(self)
        size = input_queue.qsize()
        if size in [1, jobs]:
            self.pb.setRange(0, 0)
        else:
            self.pb.setRange(0, size)
        self.content_layout.addWidget(self.pb)

        self.w = worker(input_queue, output_queue, jobs)
        self.w.job_finished.connect(self.update_progress)
        self.w.finished.connect(self.all_finished)
        self.w.start()

    def update_progress(self):
        self.counter += 1
        self.pb.setValue(self.counter)

    def all_finished(self):
        self.accept()
        self.close()


class CustomPB(QProgressBar):
    def __init__(self, parent=None):
        super(CustomPB, self).__init__(parent)
        self.setRange(0, 37)
        self.setValue(1)
        self.setFormat("Maximum RAM -- %p%")

    def paintEvent(self, event:QPaintEvent):

        painter = QStylePainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

        option = QStyleOptionProgressBar()
        self.initStyleOption(option)

        opt1 = option
        opt1.rect = self.style().subElementRect(QStyle.SubElement.SE_ProgressBarContents, opt1, self)

        opt2 = option
        rect = QRect(opt1.rect.topLeft().x(),
                     opt1.rect.topLeft().y(),
                     int(opt1.rect.width() * opt1.progress / opt1.maximum),
                     int(opt1.rect.height()))

        eventRect = event.rect()
        linearGradient = QLinearGradient(eventRect.topLeft().x(), eventRect.topLeft().y(),
                                         eventRect.topRight().x(), eventRect.topRight().y())
        linearGradient.setColorAt(0.0, QColor('#ffdde1'))
        linearGradient.setColorAt(1.0, QColor('#6f0000'))

        painter.fillRect(rect, linearGradient)

        if option.textVisible:
            opt2.rect = self.style().subElementRect(QStyle.SubElement.SE_ProgressBarLabel, opt2, self)
            painter.drawControl(QStyle.ControlElement.CE_ProgressBarLabel, opt2)