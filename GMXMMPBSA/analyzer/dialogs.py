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

from PyQt5.QtWidgets import (QDialog, QSpinBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QCheckBox,
                             QGroupBox, QButtonGroup, QGridLayout, QTreeWidget, QTreeWidgetItem,
                             QTreeWidgetItemIterator, QHeaderView, QProgressBar)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, pyqtSlot
from queue import Queue, Empty
from GMXMMPBSA.API import load_gmxmmpbsa_info
from pathlib import Path

class worker(QThread):
    job_finished = pyqtSignal()
    def __init__(self):
        super(worker, self).__init__()
        # self.parent = parent
    def define_dat(self, function, queue: Queue, result_q: Queue):
        self.fn = function
        self.queue = queue
        self.result_queue = result_q
        self.running = True

    def run(self):
        while self.running:
            if self.queue.qsize() == 0:
                return
            items = self.queue.get()
            try:
                results = self.fn(items)
                if len(results) == 1:
                    self.result_queue.put(results)
                else:
                    self.result_queue.put([x for x in results])
            finally:
                self.queue.task_done()
                self.job_finished.emit()

class InitDialog(QDialog):
    def __init__(self,parent=None):
        super(InitDialog, self).__init__(parent)
        self.parent = parent
        self.setWindowModality(Qt.WindowModal)
        self.setWindowTitle('Initialization gmx_MMPBSA_ana')
        self.curr_progress = 0
        self.data = []

        self.processing_label = QLabel()
        self.result_label = QLabel(f'...')

        self.corr_btn = QCheckBox('Calculate correlation between systems')
        self.corr_btn.setChecked(False)
        self.show_decomp_btn = QCheckBox('Show decomposition data')
        self.show_decomp_btn.setChecked(True)
        self.remove_empty_btn = QCheckBox('Remove empty terms')
        # self.remove_empty_btn.setChecked(True)

        self.check_l = QHBoxLayout()
        self.check_l.addWidget(self.corr_btn)
        self.check_l.addWidget(self.show_decomp_btn)
        self.check_l.addWidget(self.remove_empty_btn)

        self.delta_btn = QCheckBox('delta')
        self.delta_btn.setChecked(True)
        self.delta_btn.setEnabled(False)
        self.com_btn = QCheckBox('complex')
        self.rec_btn = QCheckBox('receptor')
        self.lig_btn = QCheckBox('ligand')
        self.warn_label = QLabel('This can generate a huge amount of graphics which can crash the application.')

        self.comp_group = QGroupBox('Computation of Components Charts')
        self.comp_group_layout = QGridLayout(self.comp_group)
        self.comp_group_layout.addWidget(self.delta_btn, 0, 0)
        self.comp_group_layout.addWidget(self.com_btn, 0, 1)
        self.comp_group_layout.addWidget(self.rec_btn, 0, 2)
        self.comp_group_layout.addWidget(self.lig_btn, 0, 3)
        self.comp_group_layout.addWidget(self.warn_label, 1, 0, 1, 3)

        self.sys_group = QGroupBox('Systems options')
        self.sys_group_layout = QVBoxLayout(self.sys_group)
        self.sys_group_layout.addLayout(self.check_l)
        self.sys_group_layout.addWidget(self.comp_group)

        self.hide_tb_btn = QCheckBox('Hide ToolBar')
        self.hide_tb_btn.setChecked(True)

        self.chart_group = QGroupBox('Charts options')
        self.chart_group_layout = QHBoxLayout(self.chart_group)
        self.chart_group_layout.addWidget(self.hide_tb_btn)

        self.header_item = QTreeWidgetItem(['Container', 'Name', 'Select', 'Exp. Ki (nM)'])
        self.header_item.setToolTip(0, 'Container')
        self.header_item.setToolTip(1, 'Name')
        self.result_tree = QTreeWidget(self)
        self.result_tree.setHeaderItem(self.header_item)
        self.result_tree.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.result_tree.itemChanged.connect(self.update_item_info)
        self.pb = QProgressBar()
        self.pb.setRange(0, 0)
        # self.save_btn.clicked.connect(self.save)

        self.accept_btn = QPushButton('Accept')
        self.accept_btn.clicked.connect(self.get_data)

        self.cancel_btn = QPushButton('Cancel')
        self.cancel_btn.clicked.connect(self.parent.close)

        self.btn_layout = QHBoxLayout()
        self.btn_layout.addWidget(self.pb, 5)
        self.btn_layout.addStretch(1)
        self.btn_layout.addWidget(self.cancel_btn, alignment=Qt.AlignRight)
        self.btn_layout.addWidget(self.accept_btn, alignment=Qt.AlignRight)

        self.content_layout = QVBoxLayout(self)
        self.content_layout.addWidget(self.processing_label)
        self.content_layout.addWidget(self.result_label)
        # self.content_layout.addLayout(self.check_l)
        self.content_layout.addWidget(self.sys_group)
        self.content_layout.addWidget(self.chart_group)
        self.content_layout.addWidget(self.result_tree)
        self.content_layout.addLayout(self.btn_layout)

        self.worker = worker()
        self.worker.job_finished.connect(self.fin)
        self.worker.finished.connect(self.res)

    def update_item_info(self, item: QTreeWidgetItem, col):
        path = item.info[1]
        basename = item.text(1)
        exp_ki = float(item.text(3))
        item.info = [basename, path, exp_ki]

    def get_files_info(self, info_files):
        self.nfiles = len(info_files)
        self.f_item = QTreeWidgetItem([f'All', 'All'])
        self.f_item.setCheckState(2, Qt.Checked)
        self.f_item.info = None
        self.f_item.setFlags(self.f_item.flags() | Qt.ItemIsAutoTristate)
        self.result_tree.addTopLevelItem(self.f_item)

        self.processing_label.setText(f'Processing {self.nfiles} files...')

        files_list = info_files
        self.pb.setRange(0, self.nfiles)

        c = 1
        for fname in files_list:
            # fname = files_list[x]
            basename = None
            exp_ki = None
            with open(fname) as fi:
                for line in fi:
                    if line.startswith("INPUT['sys_name']"):
                        basename = str(line.split()[2]).strip('"\'')
                    if line.startswith("INPUT['exp_ki']"):
                        exp_ki = line.split()[2]
            if not basename:
                basename = f'System_{c}'

            if not exp_ki:
                exp_ki = 0
            item = QTreeWidgetItem([f'{fname}', f'{basename}','', f'{exp_ki}'])
            item.info = [basename, Path(fname), float(exp_ki)]
            item.setCheckState(2, Qt.Checked)
            item.setFlags(item.flags() |  Qt.ItemIsEditable)
            self.f_item.addChild(item)
            c += 1
        self.f_item.setExpanded(True)

        it = QTreeWidgetItemIterator(self.f_item)
        while it.value():
            print(it.value().text(0), it.value().checkState(2) == Qt.Checked)
            it += 1

    def get_data(self):

        self.pb.setValue(0)
        queue = Queue()
        self.result_queue = Queue()
        self.systems_list = []
        self.options = {'correlation': self.corr_btn.isChecked(), 'decomposition': self.show_decomp_btn.isChecked(),
                   'components': [x.text() for x in [self.com_btn, self.rec_btn, self.lig_btn] if x.isChecked()],
                        'remove_empty': self.remove_empty_btn.isChecked()}
        print('options', self.options)
        it = QTreeWidgetItemIterator(self.f_item)
        while it.value():
            item = it.value()
            if item.checkState(2) == Qt.Checked and item.info:
                queue.put(item.info[1])
                self.systems_list.append(item.info)
            it += 1

        # for basename, fname, exp_ki in files_list:
        #     print(fname)
        #     queue.put(fname)

        self.worker.define_dat(load_gmxmmpbsa_info, queue, self.result_queue)
        self.worker.start()

    def fin(self):
        # print(re)
        print('termino')
        self.curr_progress += 1
        self.pb.setValue(self.curr_progress)
        # if self.pb.value() == self.pb.maximum():
        #     self.close()

    def res(self):
        # print(re)
        # self.parent.make_graph()
        print('fin')
        self.worker.deleteLater()
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        self.parent.process_data(self.systems_list, self.result_queue, self.options)
        # self.close()


class ExportDialog(QDialog):
    def __init__(self, parent=None):
        super(ExportDialog, self).__init__(parent)
        self.parent = parent
        self.setWindowTitle('Save Energy to b-factor')
        self.current_frame_s = QSpinBox()
        self.current_frame_l = QLabel('Frame:')

        self.output_label = QLabel('Output:')
        self.output_text = QLineEdit('complex_e2bfactor')
        self.output_text.setPlaceholderText('complex_e2bfactor')
        self.save_btn = QPushButton('Save')
        self.save_btn.clicked.connect(self.save)
        self.close_btn = QPushButton('Close')
        self.close_btn.clicked.connect(self.close)

        self.out_layout = QHBoxLayout()
        self.out_layout.addWidget(self.output_label)
        self.out_layout.addWidget(self.output_text)

        self.frame_layout = QHBoxLayout()
        self.frame_layout.addWidget(self.current_frame_l)
        self.frame_layout.addWidget(self.current_frame_s)
        # self.frame_layout.addWidget(self.interval_l, alignment=Qt.AlignRight)
        # self.frame_layout.addWidget(self.interval_s)

        self.btn_layout = QHBoxLayout()
        self.btn_layout.addStretch(1)
        self.btn_layout.addWidget(self.save_btn, alignment=Qt.AlignRight)
        self.btn_layout.addWidget(self.close_btn, alignment=Qt.AlignRight)

        self.layout = QVBoxLayout(self)
        self.layout.addLayout(self.out_layout)
        # self.layout.addWidget(self.location_btn)
        self.layout.addLayout(self.frame_layout)
        self.layout.addLayout(self.btn_layout)

        self.getdata()

    def getdata(self):
        ext_prog = find_progs(self.parent.app.INPUT)
        self.ntraj = Trajectory(self.parent.app.FILES.complex_prmtop, self.parent.app.FILES.complex_trajs,
                                ext_prog['cpptraj'])
        self.ntraj.Setup(1, 999999, 1)
        last = self.parent.app.INPUT['startframe'] + (self.ntraj.processed_frames - 1) * self.parent.app.INPUT[
            'interval']
        self.current_frame_s.setRange(self.parent.app.INPUT['startframe'], last)
        self.current_frame_s.setSingleStep(self.parent.app.INPUT['interval'])
        self.current_frame_s.setValue(self.parent.app.INPUT['startframe'])

    @pyqtSlot()
    def save(self):
        self.ntraj.Outtraj('_temp_.pdb', frames=str(self.current_frame_s.value()),
                           filetype='pdb')
        self.ntraj.Run('cpptraj.out')
        pdb = PDB()
        pdb.parse('_temp_.pdb')
        with open(str(self.output_text.text()) + '.pdb', 'w') as pdbout:
            for atm in pdb.allAtoms:
                if atm['id'] == 'TER':
                    pdbout.write(pdb.getOutLine(atm, ter=True))
                else:
                    if atm['resnum'] in self.parent.decomp and atm['resname'] in std_aa:
                        atm['b_factor'] = self.parent.decomp[atm['resnum']]
                    else:
                        atm['b_factor'] = 0.0
                    pdbout.write(pdb.getOutLine(atm, ter=False))
        # os.remove('_temp_.pdb')
        self.parent.statusbar.showMessage('Saving {} file... Done'.format(self.output_text.text() + '.pdb'))
        self.close()


class ExportDialogCSV(QDialog):
    def __init__(self, parent=None):
        super(ExportDialogCSV, self).__init__(parent)
        self.parent = parent
        self.setWindowTitle('Export Energy to CSV')

        self.output_label = QLabel('Output:')
        self.output_text = QLineEdit('TOTAL_ENERGY')
        self.output_text.setPlaceholderText('TOTAL_ENERGY')
        self.save_btn = QPushButton('Save')
        self.save_btn.clicked.connect(self.save)
        self.close_btn = QPushButton('Close')
        self.close_btn.clicked.connect(self.close)

        self.out_layout = QHBoxLayout()
        self.out_layout.addWidget(self.output_label)
        self.out_layout.addWidget(self.output_text)

        self.btn_layout = QHBoxLayout()
        self.btn_layout.addStretch(1)
        self.btn_layout.addWidget(self.save_btn, alignment=Qt.AlignRight)
        self.btn_layout.addWidget(self.close_btn, alignment=Qt.AlignRight)

        self.layout = QVBoxLayout(self)
        self.layout.addLayout(self.out_layout)
        self.layout.addLayout(self.btn_layout)

    @pyqtSlot()
    def save(self):
        with open(str(self.output_text.text()) + '.csv', 'w') as out_file:
            if self.parent.gb_data:
                out_file.write('GB data\n')
                out_file.write(','.join(['FRAME'] + [x[0].upper() for x in self.parent.gb_data]) + '\n')
                self.parent.writeData(out_file, self.parent.gb_data)
            if self.parent.mut_gb_data:
                out_file.write('Mutant GB data\n')
                out_file.write(','.join(['FRAME'] + [x[0].upper() for x in self.parent.mut_gb_data]) + '\n')
                self.parent.writeData(out_file, self.parent.mut_gb_data)
            if self.parent.pb_data:
                out_file.write('PB data\n')
                out_file.write(','.join(['FRAME'] + [x[0].upper() for x in self.parent.pb_data]) + '\n')
                self.parent.writeData(out_file, self.parent.pb_data)
            if self.parent.mut_pb_data:
                out_file.write('Mutant PB data\n')
                out_file.write(','.join(['FRAME'] + [x[0].upper() for x in self.parent.mut_pb_data]) + '\n')
                self.parent.writeData(out_file, self.parent.mut_pb_data)
        self.parent.statusbar.showMessage('Exporting PB/GB energy to csv file... Done.')
        self.close()