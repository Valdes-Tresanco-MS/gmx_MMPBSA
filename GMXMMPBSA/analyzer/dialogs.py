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