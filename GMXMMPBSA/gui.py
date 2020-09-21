# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/GMX-MMGBSA                  #
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
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import matplotlib.backend_bases
from matplotlib.figure import Figure
import numpy as np

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

class GMX_MMPBSA_GUI(QMainWindow):
    def __init__(self, info_file=None):
        super(GMX_MMPBSA_GUI, self).__init__()

        self.data = None
        self.app = None
        self.infofile = Path(info_file)

        self.mdi = QMdiArea(self)
        self.setCentralWidget(self.mdi)

        self.menubar = self.menuBar()
        self.opendirAct = QAction("Open GMX-MMPBSA Dir...", self)
        self.opendirAct.triggered.connect(self.getInfoFile)
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.opendirAct)
        self.fileMenu.addAction('Close All..',self.mdi.closeAllSubWindows)
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

        if self.infofile:
            self.getData()

    def getInfoFile(self):
        info_file, _ = QFileDialog.getOpenFileName(self, "Output Directory", QDir.currentPath(), '*_info',
                                                options=QFileDialog.DontResolveSymlinks)
        if info_file:
            self.infofile = Path(info_file)
            self.getData()

    def getData(self):
        os.chdir(self.infofile.parent)
        print(os.getcwd())

        data, self.app = API.load_gmxmmpbsa_info(self.infofile.as_posix())
        print(self.app.INPUT)
        self.data = {}
        for calc_type in data:
            self.data[calc_type] = {}

            if calc_type == 'decomp':
                # gb or pb
                for decomp_calc_type in data[calc_type]:
                    self.data[calc_type][decomp_calc_type] = {}
                    # complex, if stability: receptor, ligand
                    com = data[calc_type][decomp_calc_type]['complex']
                    # for res in com:
                    rec = lig = None
                    if 'receptor' in data[calc_type][decomp_calc_type]:
                        rec = data[calc_type][decomp_calc_type]['receptor']
                        lig = data[calc_type][decomp_calc_type]['ligand']

                    for p in com:
                        self.data[calc_type][decomp_calc_type][p] = {}
                        for res in com[p]:
                            # print(res)
                            if res not in self.data[calc_type][decomp_calc_type][p]:
                                self.data[calc_type][decomp_calc_type][p][res] = {}
                            if rec and res in rec[p]:
                                for resp, resp1 in zip(com[p][res], rec[p][res]):
                                    # d = com[p][res][resp] - rec[p][res][resp1]
                                    self.data[calc_type][decomp_calc_type][p][res][resp] = {}
                                    for para in com[p][res][resp]:
                                        d = com[p][res][resp][para] - rec[p][res][resp1][para]
                                        self.data[calc_type][decomp_calc_type][p][res][resp][para] = d
                                for resp in com[p][res]:
                                    self.data[calc_type][decomp_calc_type][p][res][resp] = com[p][res][resp]
                            elif lig and res in lig[p]:
                                for resp, resp1 in zip(list(com[p][res])[(len(com[p][res]) - len(lig[p][res])):],
                                                       lig[p][res]):
                                    # d = com[p][res][resp] - rec[p][res][resp1]
                                    self.data[calc_type][decomp_calc_type][p][res][resp] = {}
                                    for para in com[p][res][resp]:
                                        # if para == 'tot':
                                        #     print('############',res, resp, resp1, (com[p][res][resp][para]-
                                        #           lig[p][res][resp1][para]).mean())
                                        #     break

                                        d = com[p][res][resp][para] - lig[p][res][resp1][para]
                                        self.data[calc_type][decomp_calc_type][p][res][resp][para] = d
                                for resp in com[p][res]:
                                    self.data[calc_type][decomp_calc_type][p][res][resp] = com[p][res][resp]

            else:
                com = data[calc_type]['complex']

                if not self.app.stability:
                    rec = data[calc_type]['receptor']
                    lig = data[calc_type]['ligand']

                    for k in com:
                        self.data[calc_type][k] = com[k] - (rec[k] + lig[k])
                else:
                    self.data[calc_type] = com

        self.makeTree()

    def showdata(self, item: CustomItem, col):
        print(item, col, item.subwindows)
        if col == 1:
            s = item.subwindows[1]
            if item.checkState(col) == Qt.Checked:
                if s: # check if any subwindow has been store
                    if not s.isVisible():
                        s.show()
                else:

                    sub = Charts()
                    sub.setObjectName(item.dataperframe['name'])
                    sub.item = item
                    sub.col = col
                    item.subwindows[col] = sub
                    sm = sub.systemMenu()

                    x = [x1 for x1 in range(0, len(item.dataperframe['data']), item.dataperframe['interval'])]
                    x = np.array(x)
                    y = item.dataperframe['data']

                    sub.setMinimumSize(500, 300)
                    mainwidgetmdi = QMainWindow()

                    self._main = QWidget()
                    mainwidgetmdi.setCentralWidget(self._main)
                    layout = QVBoxLayout(self._main)
                    sc = MplCanvas(self, width=5, height=4, dpi=100)
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
                    # plt.clf()
                    # plt.margins(0.01)
                    sc.draw()
                    sc.figure.tight_layout()
                    sc.draw()

                # ax = plt.gca()
                # ax.set_title(item.dataperframe['name'])
                # # sub.figure.suptitle(xvg_file.title)
                #
                # print(x, y)
                # plt.plot(x, y, color='black', linewidth=0.3)
                # sub.show()
                # sub.canvas.draw()
                # sub.figure.tight_layout()
            else:
                if s:  # check if any subwindow has been store
                    if s.isVisible():
                        s.close()

        elif col == 2:
            s = item.subwindows[2]
            if item.checkState(col) == Qt.Checked:
                if s: # check if any subwindow has been store
                    if not s.isVisible():
                        s.show()
                else:
                    sub = Charts()
                    sub.setObjectName(item.datamean['name'])
                    sub.item = item
                    sub.col = col
                    item.subwindows[col] = sub
                    sm = sub.systemMenu()

                    # x = item.datamean['data']
                    # print(x)
                    y = item.datamean['data']

                    sub.setMinimumSize(500, 300)
                    mainwidgetmdi = QMainWindow()

                    self._main = QWidget()
                    mainwidgetmdi.setCentralWidget(self._main)
                    layout = QVBoxLayout(self._main)
                    sc = MplCanvas(self, width=5, height=4, dpi=100)

                    data_len = len(item.datamean['xaxis'])
                    # print(data_len, item.datamean['xaxis'], y)

                    sc.axes.set_xlim(0, data_len+1)
                    sc.axes.set_xticks([x for x in range(1, data_len+1)])
                    sc.axes.set_xticklabels(item.datamean['xaxis'])
                    sc.axes.set_ylabel(item.datamean['yaxis'])
                    for label in sc.axes.get_xticklabels():
                        label.set_rotation(40)
                        label.set_horizontalalignment('right')


                    bottom = 0
                    for i in range(0, data_len):
                        # print(i)
                        sc.axes.bar(i+1, y[i])
                        bottom += y[i]

                    # sc.axes.plot(x, y, color='black', linewidth=0.3)

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
                    # elif item.checkState(col) == Qt.Unchecked and col == 1:
                    #     try:
                    #         item.subwindows[col].close()
                    #     except:
                    #         pass
                    # elif item.checkState(col) == Qt.Unchecked and col == 2:
                    #     try:
                    #         item.subwindows[col].close()
                    #     except:
                    #         pass
            else:
                if s:  # check if any subwindow has been store
                    if s.isVisible():
                        s.close()
        else:
            pass
            # windows = self.mdi.subWindowList()
            # for x in windows:
            #     if hasattr(item, 'name'):
            #         if x.objectName() == btn.name:
            #             self.mdi.removeSubWindow(x)


    def makeTree(self):
        self.treeWidget = QTreeWidget(self)
        self.treeWidget.itemChanged.connect(self.showdata)
        self.treeDockWidget.setWidget(self.treeWidget)
        self.treeWidget.setHeaderLabels(['Labels', 'Per-Frame', 'Mean'])

        self.normalitem = CustomItem(['Normal System'])
        self.normalitem.setExpanded(True)
        self.treeWidget.addTopLevelItem(self.normalitem)


        for level in self.data:
            item = CustomItem([level.upper()])
            # item.setCheckState(0, Qt.Unchecked)
            self.normalitem.addChild(item)
            if level in ['gb', 'pb']:
                # if level == 'gb':
                #     item.
                cd = []
                for level1 in self.data[level]:
                    item1= CustomItem([str(level1).upper()])
                    item1.setCheckState(1,Qt.Unchecked)
                    item1.dataperframe = {'name': '{} {} Energy (Per-Frame)'.format(level.upper(), level1),
                                          'xaxis': 'frames', 'yaxis': 'Energy (kcal/mol)', 'interval':
                                              self.app.INPUT['interval'], 'data': self.data[level][level1]}
                    print(level1, self.data[level][level1].mean())
                    cd.append([level1, self.data[level][level1].mean()])
                    if level1 == 'TOTAL':
                        item1.setCheckState(2,Qt.Unchecked)
                        item1.datamean = {'name': '{} {} Energy (Mean)'.format(level.upper(), level1),
                                          'xaxis': [x[0] for x in cd], 'yaxis': 'Energy (kcal/mol)',
                                          'data': np.array([x[1] for x in cd]), }

                    item.addChild(item1)
                    # for level2 in self.data[level][level1]:
                    #     item2 = CustomItem([str(level2).upper()])
                    #     item1.addChild(item2)


            elif level == 'decomp':
                for level1 in self.data[level]:
                    item1 = CustomItem([str(level1).upper()])
                    item.addChild(item1)
                    for level2 in self.data[level][level1]:
                        item2 = CustomItem([str(level2).upper()])
                        item1.addChild(item2)
                        for level3 in self.data[level][level1][level2]:
                            item3 = CustomItem([str(level3).upper()])
                            item2.addChild(item3)
                            for level4 in self.data[level][level1][level2][level3]:
                                item4 = CustomItem([str(level4).upper()])
                                item3.addChild(item4)
                                print(level4)
                                for level5 in self.data[level][level1][level2][level3][level4]:
                                    item5 = CustomItem([str(level5).upper()])
                                    item4.addChild(item5)
                                    # if self.app.INPUT['idecomp'] in (3, 4):
                                    #
                                    #     for level6 in self.data[level][level1][level2][level3][level4][level5]:
                                    #         item6 = CustomItem([str(level6).upper()])
                                    #         item5.addChild(item6)

        # if self.data.mutant:
        #     self.mutantitem = CustomItem(['Mutant System'])
        #     self.treeWidget.addTopLevelItem(self.mutantitem)
        self.normalitem.setExpanded(True)
        self.treeWidget.header().setSectionResizeMode(QHeaderView.ResizeToContents)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setApplicationName('GMX-MMPBSA')


    w = GMX_MMPBSA_GUI(sys.argv[1])
    w.show()
    sys.exit(app.exec())