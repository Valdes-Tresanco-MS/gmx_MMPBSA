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

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import matplotlib.backend_bases
from matplotlib.backends import qt_compat
from matplotlib.figure import Figure
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import seaborn as sns
import matplotlib.pyplot as plt
import math
plt.style.use('seaborn')
plt.rcParams["figure.autolayout"] = True


import os

class NavigationToolbar(NavigationToolbar2QT):
    """
    overwrite NavigatorToolbar to get control over save action.
    Now we can set custom dpi value
    """
    def __init__(self, canvas, parent, coordinates=True):
        super(NavigationToolbar, self).__init__(canvas, parent, coordinates=True)
        # NavigationToolbar2QT.__init__(canvas, parent, coordinates=True)
        l = self.layout()
        # l = QHBoxLayout()
        l.setContentsMargins(0, 0, 0, 0)

    def _init_toolbar(self):
        pass

    def save_figure(self, *args):

        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(filetypes.items())
        default_filetype = self.canvas.get_default_filetype()

        startpath = os.path.expanduser(matplotlib.rcParams['savefig.directory'])
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
                self.canvas.figure.savefig(fname, dpi=600)
            except Exception as e:
                QMessageBox.critical(
                    self, "Error saving file", str(e),
                    QMessageBox.Ok, QMessageBox.NoButton)


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)


class Charts(QMdiSubWindow):
    LINE = 1
    ROLLING = 2
    JOINPLOT = 3
    BAR = 4
    SBAR = 5
    HEATMAP = 6
    SCATTER = 7
    def __init__(self, *args, item, col, options:dict=None):
        super(Charts, self).__init__(*args)
        self.setMinimumSize(400, 400)

        self.mainwidgetmdi = QMainWindow() # must be QMainWindow to handle toolbar

        self.mpl_canvas = MplCanvas(self, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        self.mpl_toolbar = NavigationToolbar(self.mpl_canvas, self)
        self.mainwidgetmdi.setCentralWidget(self.mpl_canvas)
        self.fbtn = QPushButton(self.style().standardIcon(QStyle.SP_FileDialogDetailedView), '', self.mpl_canvas)
        self.fbtn.setToolTip('Show or Hide the Navigation Toolbar')
        self.fbtn.toggled.connect(self.mpl_toolbar.setVisible)
        self.fbtn.setCheckable(True)
        self.fbtn.setChecked(True)

        self.line_plot = None
        self.movav_plot = None
        self.bar_plot = None
        self.heatmap_plot = None
        self.reg_plot = None
        # self.main_w_layout.addWidget(self.mpl_toolbar)
        self.mainwidgetmdi.addToolBar(Qt.BottomToolBarArea, self.mpl_toolbar)
        self.setWidget(self.mainwidgetmdi)
        self.item = item
        self.col = col
        self.options = options

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.item.setCheckState(self.col, Qt.Unchecked)

    def make_chart(self):
        """

        :param graph_type: 1: line plot, 2: bar plot, 3: heatmap plot
        :return:
        """
        if Charts.SCATTER in self.options['chart_type']:
            if self.reg_plot:
                self.reg_plot.cla()
                self.reg_plot.clear()
            if self.col == 1:
                self.item.dh_sw = self
                self.data = self.item.enthalpy
                self.reg_plot = sns.regplot(data=self.data, x='Exp.Energy', y='ΔH', ax=self.mpl_canvas.axes)
            elif self.col == 2:
                self.item.dgie_sw = self
                self.data = self.item.dgie
                self.reg_plot = sns.regplot(data=self.data, x='Exp.Energy', y='ΔH+IE', ax=self.mpl_canvas.axes)
            elif self.col == 3:
                self.item.dgnmode_sw = self
                self.data = self.item.dgnmode
                self.reg_plot = sns.regplot(data=self.data, x='Exp.Energy', y='ΔH+NMODE', ax=self.mpl_canvas.axes)
            else:
                self.item.dgqh_sw = self
                self.data = self.item.dgqh
                self.reg_plot = sns.regplot(data=self.data, x='Exp.Energy', y='ΔH+QH', ax=self.mpl_canvas.axes)


            self.mpl_canvas.axes.set_xlabel('Exp. Energy (kcal/mol)')
            self.mpl_canvas.axes.set_ylabel('Pred. Energy (kcal/mol)')
            chart_subtitle = self.item.chart_subtitle[self.col - 1]

            self.mpl_canvas.axes.set_title(self.item.chart_title + '\n' + chart_subtitle)
            self.mpl_canvas.draw()
            self.mpl_canvas.figure.tight_layout()
            self.mpl_canvas.draw()
            self.setWindowTitle(self.item.chart_title)
        else:
            self.data = self.item.gmxMMPBSA_current_data
            title = ''
            if Charts.LINE in self.options['chart_type']:
                self.item.lp_subw = self
                # self.mpl_canvas.axes.cla()
                if self.line_plot:
                    self.line_plot.cla()
                    self.line_plot.clear()

                label = 'ΔG'
                if 'IE' in self.item.item_name:
                    label = '-TΔS'

                self.line_plot = sns.lineplot(data=self.data.line_plot_dat, x='frames', color='black', linewidth=0.7,
                             y='Energy', label=label, ax=self.mpl_canvas.axes)
                title = self.item.chart_title + '(P.f)'
                if 'IE' in self.item.item_name:
                    self.mpl_canvas.axes.set_ylabel('Interaction Entropy (kcal/mol)')
                    self.mpl_canvas.axes.axvline(self.item.ie[0][0], ls='--', color='r')
                    self.mpl_canvas.axes.hlines(self.item.ie[1], xmin=self.item.ie[0][0], xmax=self.item.ie[0][1], ls='--',
                                                color='g')
                else:
                    self.mpl_canvas.axes.set_ylabel('Energy (kcal/mol)')

            if Charts.ROLLING in self.options['chart_type']:
                if 'IE' not in self.item.item_name:
                    self.data.line_plot_dat['movav'] = self.data.line_plot_dat['Energy'].rolling(50).mean()
                    self.movav_plot = sns.lineplot(data=self.data.line_plot_dat, x='frames', color='red', linewidth=0.8,
                                                    y='movav', label='Mov. Av.', ax=self.mpl_canvas.axes)

                # self.setWindowTitle(title)
            if Charts.BAR in self.options['chart_type']:
                self.item.bp_subw = self
                if self.bar_plot:
                    self.bar_plot.cla()
                    self.bar_plot.clear()
                self.bar_plot  = sns.barplot(data=self.data.bar_plot_dat, ci="sd", errwidth=1,
                                             ax=self.mpl_canvas.axes)
                for label in self.mpl_canvas.axes.get_xticklabels():
                    label.set_rotation(40)
                    label.set_horizontalalignment('right')
                self.mpl_canvas.axes.set_xlabel('')
                self.mpl_canvas.axes.set_ylabel('Energy (kcal/mol)')
                title = self.item.chart_title + '(Av.)'

            if Charts.SBAR in self.options['chart_type']:
                # TODO: stacked bars
                pass

            if Charts.HEATMAP in self.options['chart_type']:
                self.item.hmp_subw = self
                if self.heatmap_plot:
                    self.mpl_canvas.axes.collections[-1].colorbar.remove()
                    self.heatmap_plot.cla()
                    self.heatmap_plot.clear()
                self.heatmap_plot = sns.heatmap(self.data.heatmap_plot_dat, ax=self.mpl_canvas.axes, center=0,
                                                yticklabels=self.data.heatmap_plot_dat.index.tolist(),
                                                cmap='seismic', cbar_kws={'label': 'Energy (kcal/mol)'})
                for label in self.heatmap_plot.axes.get_yticklabels():
                    label.set_rotation(0)
                    label.set_horizontalalignment('right')
                title = self.item.chart_title + '(P.f)' # Fixme: no frames from correlation

            if self.item.item_name != 'IE':
                self.mpl_canvas.axes.invert_yaxis()


            self.mpl_canvas.axes.set_title(title + '\n' + self.item.chart_subtitle)
            self.mpl_canvas.draw()
            self.mpl_canvas.figure.tight_layout()
            self.mpl_canvas.draw()
            self.setWindowTitle(title)
