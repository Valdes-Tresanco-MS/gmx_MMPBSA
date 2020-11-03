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
import sys
import os
from pathlib import Path
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from GMXMMPBSA import API
from GMXMMPBSA.make_trajs import Trajectory
from GMXMMPBSA.utils import PDB
from GMXMMPBSA.findprogs import find_progs
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import matplotlib.backend_bases
from matplotlib.backends import qt_compat
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from math import ceil

plt.style.use('seaborn')

std_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'ASH', 'CYS', 'CYX', 'GLN', 'GLU', 'GLH', 'GLY', 'HID', 'HIE', 'HIP', 'ILE',
          'LEU', 'LYS', 'LYN', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIS']


class NavigationToolbar(NavigationToolbar2QT):
    """
    overwrite NavigatorToolbar to get control over save action.
    Now we can set custom dpi value
    """

    def _init_toolbar(self):
        pass

    def save_figure(self, *args):
        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(filetypes.items())
        default_filetype = self.canvas.get_default_filetype()

        startpath = os.path.expanduser(
            matplotlib.rcParams['savefig.directory'])
        start = os.path.join(startpath, self.canvas.get_default_filename())
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(['*.%s' % ext for ext in exts])
            filter = '%s (%s)' % (name, exts_list)
            if default_filetype in exts:
                selectedFilter = filter
            filters.append(filter)
        filters = ';;'.join(filters)

        fname, filter = qt_compat._getSaveFileName(
            self.canvas.parent(), "Choose a filename to save to", start,
            filters, selectedFilter)
        if fname:
            # Save dir for next time, unless empty str (i.e., use cwd).
            if startpath != "":
                matplotlib.rcParams['savefig.directory'] = (
                    os.path.dirname(fname))
            try:
                self.canvas.figure.savefig(fname, dpi=300)
            except Exception as e:
                QMessageBox.critical(
                    self, "Error saving file", str(e),
                    QMessageBox.Ok, QMessageBox.NoButton)


class Charts(QMdiSubWindow):
    def __init__(self, *args):
        super(Charts, self).__init__(*args)
        self.item = None
        self.col = None

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.item.setCheckState(self.col, Qt.Unchecked)


class CustomItem(QTreeWidgetItem):
    def __init__(self, *args):
        super(CustomItem, self).__init__(*args)
        self.subwindows = {1: None, 2: None}
        self.dataperframe = None
        self.datamean = None
        self.name = None
        self.keyslist = []


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


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
                                ext_prog['cpptraj'].full_path)
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


class GMX_MMPBSA_GUI(QMainWindow):
    def __init__(self, info_file=None):
        super(GMX_MMPBSA_GUI, self).__init__()

        self.data = None
        self.mutant_data = None
        self.app = None
        self.infofile = Path(info_file)

        self.gb_data = None
        self.mut_gb_data = None
        self.pb_data = None
        self.mut_pb_data = None
        self.decomp = {}
        self.mdi = QMdiArea(self)
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

        self.treeDockWidget = QDockWidget(self)
        self.treeDockWidget.setWindowTitle('Data')
        self.addDockWidget(Qt.RightDockWidgetArea, self.treeDockWidget)
        self.optionDockWidget = QDockWidget(self)
        self.optionDockWidget.setWindowTitle('Options')
        self.addDockWidget(Qt.RightDockWidgetArea, self.optionDockWidget)

        self.getData()
        if 'gb' in self.data or 'pb' in self.data:
            self.fileMenu.addAction('Export GB/PB energy (csv)', self.exportCSV)
        if 'decomp' in self.data:
            self.fileMenu.addAction('Energy to bfactor', self.energy2bfactor)

        self.exportpdb = ExportDialog(self)
        self.exportcsv = ExportDialogCSV(self)

    def getInfoFile(self):
        info_file, _ = QFileDialog.getOpenFileName(self, "Output Directory", QDir.currentPath(), '*_info',
                                                   options=QFileDialog.DontResolveSymlinks)
        if info_file:
            self.infofile = Path(info_file)
            self.getData()

    def getData(self):
        os.chdir(self.infofile.parent)
        self.data, self.app = API.load_gmxmmpbsa_info(self.infofile.as_posix())
        self.makeTree()

    def showdata(self, item: CustomItem, col):
        if col == 1:
            s = item.subwindows[1]
            if item.checkState(col) == Qt.Checked:
                if s:  # check if any subwindow has been store
                    if not s.isVisible():
                        s.show()
                else:

                    sub = Charts()
                    sub.setObjectName(item.dataperframe['name'])
                    sub.item = item
                    sub.col = col
                    item.subwindows[col] = sub

                    x = item.dataperframe['frames']
                    x = np.array(x)
                    y = item.dataperframe['data']

                    sub.setMinimumSize(500, 300)
                    mainwidgetmdi = QMainWindow()

                    self._main = QWidget()
                    mainwidgetmdi.setCentralWidget(self._main)
                    layout = QVBoxLayout(self._main)
                    sc = MplCanvas(self, width=5, height=4, dpi=100)
                    sc.axes.set_title(item.dataperframe['name'])
                    sc.axes.plot(x, y, color='black', linewidth=0.5)
                    sc.axes.set_xlabel(item.dataperframe['xaxis'])
                    sc.axes.set_ylabel(item.dataperframe['yaxis'])

                    # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
                    toolbar = NavigationToolbar(sc, self)
                    layout.addWidget(sc)
                    layout.addWidget(toolbar)

                    sub.setWidget(mainwidgetmdi)
                    sub.setWindowTitle(item.dataperframe['name'])
                    self.mdi.addSubWindow(sub)
                    sub.show()
                    sc.draw()
                    sc.figure.tight_layout()
                    sc.draw()
            else:
                if s:  # check if any subwindow has been store
                    if s.isVisible():
                        s.close()

        elif col == 2:
            s = item.subwindows[2]
            if item.checkState(col) == Qt.Checked:
                if s:  # check if any subwindow has been store
                    if not s.isVisible():
                        s.show()
                else:
                    sub = Charts()
                    sub.setObjectName(item.datamean['name'])
                    sub.item = item
                    sub.col = col
                    item.subwindows[col] = sub

                    y = item.datamean['data']

                    sub.setMinimumSize(500, 300)
                    mainwidgetmdi = QMainWindow()

                    self._main = QWidget()
                    mainwidgetmdi.setCentralWidget(self._main)
                    layout = QVBoxLayout(self._main)
                    sc = MplCanvas(self, width=5, height=4, dpi=100)
                    sc.axes.set_title(item.datamean['name'])
                    data_len = len(item.datamean['xaxis'])

                    sc.axes.set_xlim(0, data_len + 1)
                    sc.axes.set_xticks([x for x in range(1, data_len + 1)])
                    sc.axes.set_xticklabels(item.datamean['xaxis'])
                    sc.axes.set_ylabel(item.datamean['yaxis'])
                    for label in sc.axes.get_xticklabels():
                        label.set_rotation(40)
                        label.set_horizontalalignment('right')

                    bottom = 0
                    for i in range(0, data_len):
                        sc.axes.bar(i + 1, y[i])
                        bottom += y[i]

                    # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
                    toolbar = NavigationToolbar(sc, self)
                    layout.addWidget(sc)
                    layout.addWidget(toolbar)

                    sub.setWidget(mainwidgetmdi)
                    sub.setWindowTitle(item.datamean['name'])
                    self.mdi.addSubWindow(sub)
                    sub.show()
                    sc.draw()
                    sc.figure.tight_layout()
                    sc.draw()
            else:
                if s:  # check if any subwindow has been store
                    if s.isVisible():
                        s.close()

    def writeData(self, outfile, data):
        for i in range(len(self.frames)):
            outfile.write(','.join([str(self.frames[i])] + [str(x[1][i]) for x in data]) + '\n')

    def energy2bfactor(self):
        self.exportpdb.show()

    def exportCSV(self):
        self.exportcsv.show()

    def makeItems(self, data, topItem, mutant=False):
        mut_pre = ''
        if mutant:
            mut_pre = 'Mutant '

        for level in data:
            item = CustomItem([level.upper()])
            topItem.addChild(item)
            dh = None
            if level in ['gb', 'pb', 'rism gf', 'rism std']:
                cd = []
                dat = []
                for level1 in data[level]['delta']:
                    item1 = CustomItem([str(level1).upper()])
                    item1.setCheckState(1, Qt.Unchecked)
                    end = len(data[level]['delta'][level1]) * self.app.INPUT['interval'] + self.app.INPUT['startframe']
                    self.frames = [x for x in range(self.app.INPUT['startframe'], end, self.app.INPUT['interval'])]
                    item1.dataperframe = {
                        'name': mut_pre + '{} {} Energy (Per-Frame)'.format(level.upper(), level1.upper()),
                        'xaxis': 'frames', 'yaxis': 'Energy (kcal/mol)', 'frames': self.frames,
                        'data': data[level]['delta'][level1]}
                    cd.append([level1, data[level]['delta'][level1].mean()])
                    dat.append([level1, data[level]['delta'][level1]])
                    if level1 == 'DELTA TOTAL':
                        item1.setCheckState(2, Qt.Unchecked)
                        item1.datamean = {
                            'name': mut_pre + '{} {} Energy (Mean)'.format(level.upper(), level1.upper()),
                            'xaxis': [x[0].upper() for x in cd], 'yaxis': 'Energy (kcal/mol)',
                            'data': np.array([x[1] for x in cd])}
                        dh = data[level]['delta'][level1].mean()
                    item.addChild(item1)
                    if level1 == 'DELTA TOTAL' and 'nmode' in data:
                        itemd = CustomItem(['DELTA G Binding (NM)'])
                        itemd.setCheckState(2, Qt.Unchecked)
                        entropy = data['nmode']['delta']['Total'][0]
                        itemd.datamean = {
                            'name': mut_pre + '{} Binding Energy (with Entropy (nmode))'.format(level.upper()),
                            'xaxis': ['ΔH', '-TΔS', 'ΔG'], 'yaxis': 'Energy (kcal/mol)',
                            'data': np.array([dh, entropy * -1, dh - entropy])}
                        item.addChild(itemd)
                    if level1 == 'DELTA TOTAL' and 'qh' in data:
                        itemd = CustomItem(['DELTA G Binding (QH)'])
                        itemd.setCheckState(2, Qt.Unchecked)
                        entropy = data['qh']['delta']['Total'][0]
                        itemd.datamean = {
                            'name': mut_pre + '{} Binding Energy (with Entropy (qh))'.format(level.upper()),
                            'xaxis': ['ΔH', '-TΔS', 'ΔG'], 'yaxis': 'Energy (kcal/mol)',
                            'data': np.array([dh, entropy * -1, dh - entropy])}
                        item.addChild(itemd)
                    if level1 == '-TDS':
                        itemd = CustomItem(['DELTA G Binding (IE)'])
                        itemd.setCheckState(2, Qt.Unchecked)
                        s = 1
                        if len(data[level]['delta']['-TDS']) * self.app.INPUT['entropy_seg'] / 100 > 1:
                            s = ceil(len(data[level]['delta']['-TDS']) * (1 - self.app.INPUT['entropy_seg'] / 100))
                        ts = data[level]['delta']['-TDS'][s:].mean()
                        itemd.datamean = {
                            'name': mut_pre + '{} Total Energy (with Entropy)'.format(level.upper()),
                            'xaxis': [['ΔH', '-TΔS', 'ΔG']], 'yaxis': 'Energy (kcal/mol)',
                            'data': np.array([dh, ts, dh + ts])}
                        item.addChild(itemd)
                # To export data in csv
                if level == 'gb':
                    if mutant:
                        self.mut_gb_data = dat
                    else:
                        self.gb_data = dat
                else:
                    if mutant:
                        self.mut_pb_data = dat
                    else:
                        self.pb_data = dat

            elif level == 'decomp':
                for level1 in data[level]:
                    item1 = CustomItem([str(level1).upper()])
                    item.addChild(item1)
                    for level2 in data[level][level1]['delta']:
                        item2 = CustomItem([str(level2).upper()])
                        item1.addChild(item2)
                        item2.setCheckState(1, Qt.Unchecked)
                        item2.setCheckState(2, Qt.Unchecked)
                        lvl2_data = np.zeros(int(self.app.numframes / self.app.INPUT['interval']))
                        lvl2_meandata = []
                        for level3 in data[level][level1]['delta'][level2]:
                            item3 = CustomItem([str(level3).upper()])
                            item2.addChild(item3)
                            lvl3_data = np.zeros(int(self.app.numframes / self.app.INPUT['interval']))
                            lvl3_meandata = []
                            lvl4_data = []
                            lvl4_meandata = []
                            for level4 in data[level][level1]['delta'][level2][level3]:
                                item4 = CustomItem([str(level4).upper()])
                                item3.addChild(item4)
                                lvl5_meandata = []
                                if self.app.INPUT['idecomp'] in [3, 4]:
                                    for level5 in data[level][level1]['delta'][level2][level3][level4]:
                                        item5 = CustomItem([str(level5).upper()])
                                        item4.addChild(item5)
                                        item5.setCheckState(1, Qt.Unchecked)
                                        item5.dataperframe = {
                                            'name': mut_pre + 'Docomposition ({}) {} Residue {} Pair {} {} Energy '
                                                              '(Per-Frame)'.format(level1.upper(), level2.upper(),
                                                                                   level3, level4, level5.upper()),
                                            'xaxis': 'frames', 'yaxis': 'Energy (kcal/mol)', 'frames': self.frames,
                                            'data': data[level][level1]['delta'][level2][level3][level4][level5]}
                                        lvl5_meandata.append(
                                            [level5,
                                             data[level][level1]['delta'][level2][level3][level4][level5].mean()])
                                        if level5 == 'tot':
                                            item5.setCheckState(2, Qt.Unchecked)
                                            lvl3_data += data[level][level1]['delta'][level2][level3][level4][level5]
                                            lvl3_meandata.append(
                                                [level4,
                                                 data[level][level1]['delta'][level2][level3][level4][level5].mean()])

                                            item5.datamean = {
                                                'name': mut_pre + 'Decomposition ({}) {} Residue {} Pair {} '
                                                                  'Total Energy (Mean)'.format(level1.upper(),
                                                                                               level2.upper(),
                                                                                               level3, level4),
                                                'xaxis': [x[0].upper() for x in lvl5_meandata],
                                                'yaxis': 'Energy (kcal/mol)',
                                                'data': np.array([x[1] for x in lvl5_meandata])}
                                else:
                                    lvl3_data += data[level][level1]['delta'][level2][level3][level4]
                                    lvl3_meandata.append(
                                        [level4, data[level][level1]['delta'][level2][level3][level4].mean()])
                                    item4.setCheckState(1, Qt.Unchecked)
                                    item4.dataperframe = {
                                        'name': mut_pre + 'Docomposition ({}) {} Residue {} {} Energy '
                                                          '(Per-Frame)'.format(level1.upper(), level2.upper(), level3,
                                                                               level4.upper()),
                                        'xaxis': 'frames', 'yaxis': 'Energy (kcal/mol)', 'frames': self.frames,
                                        'data': data[level][level1]['delta'][level2][level3][level4]}
                                    lvl4_meandata.append(
                                        [level4, data[level][level1]['delta'][level2][level3][level4].mean()])
                                    if level4 == 'tot':
                                        item4.setCheckState(2, Qt.Unchecked)
                                        item4.datamean = {
                                            'name': mut_pre + 'Decomposition ({}) {} Residue {} Total Energy ('
                                                              'Mean)'.format(level1.upper(), level2.upper(), level3),
                                            'xaxis': [x[0].upper() for x in lvl4_meandata],
                                            'yaxis': 'Energy (kcal/mol)',
                                            'data': np.array([x[1] for x in lvl4_meandata]), }

                            if 'tot' in data[level][level1]['delta'][level2][level3]:  # Per-residue
                                lvl2_data += lvl3_data
                                lvl2_meandata.append([level3, lvl3_data.mean()])
                                if level2 == 'TDC':
                                    self.decomp[level3] = lvl3_data.mean()
                            else:  # Per-wise
                                item3.setCheckState(1, Qt.Unchecked)
                                item3.dataperframe = {
                                    'name': mut_pre + 'Decomposition ({})  {} Residue {} Energy (Per-Frame)'.format(
                                        level1.upper(), level2.upper(), level3), 'xaxis': 'frames',
                                    'yaxis': 'Energy (kcal/mol)', 'frames': self.frames,
                                    'data': lvl3_data}
                                if level2 == 'TDC':
                                    self.decomp[level3] = sum([x[1] for x in lvl3_meandata])
                                lvl3_meandata.extend([['TOTAL', sum([x[1] for x in lvl3_meandata])]])
                                lvl2_data += lvl3_data
                                lvl2_meandata.append([level3, lvl3_data.mean()])

                                item3.datamean = {
                                    'name': mut_pre + 'Decomposition ({}) {} Residue {} Energy (Mean) from '
                                                      'Pairs'.format(level1.upper(), level2.upper(), level3),
                                    'xaxis': [x[0] for x in lvl3_meandata], 'yaxis': 'Energy (kcal/mol)',
                                    'data': np.array([x[1] for x in lvl3_meandata])}
                                item3.setCheckState(2, Qt.Unchecked)
                        item2.dataperframe = {
                            'name': mut_pre + 'Decomposition ({})  {} Energy (Per-Frame)'.format(
                                level1.upper(), level2.upper()), 'xaxis': 'frames', 'yaxis': 'Energy (kcal/mol)',
                            'frames': self.frames, 'data': lvl2_data}
                        lvl2_meandata.extend([['TOTAL', sum([x[1] for x in lvl2_meandata])]])
                        item2.datamean = {
                            'name': mut_pre + 'Decomposition ({}) {} Energy (Mean)'.format(
                                level1.upper(), level2.upper()), 'xaxis': [x[0] for x in lvl2_meandata],
                            'yaxis': 'Energy (kcal/mol)', 'data': np.array([x[1] for x in lvl2_meandata])}

            elif level in ['nmode', 'qh']:
                item.setCheckState(2, Qt.Unchecked)
                item.datamean = {
                    'name': mut_pre + 'Entropy ({})'.format(level.upper()),
                    'xaxis': [x[0].upper() for x in data[level]['delta'].items()], 'yaxis': 'Energy (kcal/mol)',
                    'data': np.array([x[1][0] for x in data[level]['delta'].items()])}

    def makeTree(self):
        self.treeWidget = QTreeWidget(self)
        self.treeWidget.itemChanged.connect(self.showdata)
        self.treeDockWidget.setWidget(self.treeWidget)
        self.treeWidget.setHeaderLabels(['Labels', 'Per-Frame', 'Mean'])

        self.normalitem = CustomItem(['Normal System'])
        self.normalitem.setExpanded(True)
        self.treeWidget.addTopLevelItem(self.normalitem)
        self.makeItems(self.data, self.normalitem)
        self.normalitem.setExpanded(True)

        if self.data.mutant:
            self.mutantitem = CustomItem(['Mutant System'])
            self.treeWidget.addTopLevelItem(self.mutantitem)
            self.makeItems(self.data.mutant, self.mutantitem, True)
            self.mutantitem.setExpanded(True)

        self.treeWidget.header().setSectionResizeMode(QHeaderView.ResizeToContents)


def run(infofile):
    app = QApplication(sys.argv)
    app.setApplicationName('gmx_MMPBSA')
    w = GMX_MMPBSA_GUI(infofile)
    w.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setApplicationName('gmx_MMPBSA')

    w = GMX_MMPBSA_GUI('/home/mario/Drive/scripts/gmx_MMPBSA/test/qh/_GMXMMPBSA_info')
    w.show()
    sys.exit(app.exec())
