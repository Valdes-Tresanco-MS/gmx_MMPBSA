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
import matplotlib as mpl
import matplotlib.backend_bases
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas
import pandas as pd
import seaborn as sns
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib import gridspec
from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
import itertools

# sns.set_theme()
from GMXMMPBSA.analyzer.chartsettings import Palettes

plt.rcParams["figure.autolayout"] = True

import os


def rgb2rgbf(color):
    return [x / 255 for x in color]


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
            exts_list = " ".join('*.%s' % ext for ext in exts)
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


class ChartsBase(QMdiSubWindow):
    def __init__(self, button: QToolButton, options: dict = None):
        super(ChartsBase, self).__init__()
        self.setMinimumSize(400, 400)
        self.options = options['general_options']
        self.mainwidgetmdi = QMainWindow()  # must be QMainWindow to handle the toolbar
        sns.set_theme(style=self.options['theme'])
        self.plot = None
        self.frange = []  # Frames range with which it was created
        self.button = button
        self.setWidget(self.mainwidgetmdi)

    def set_cw(self, fig=None):
        # we create the figure canvas here because the fig parameter must be defined and depende of what kind of
        # chart want to make
        fig = Figure(dpi=self.options['dpi']['plot']) if not fig else fig
        self.figure_canvas = FigureCanvas(fig)
        self.fig = self.figure_canvas.figure
        self.mainwidgetmdi.setCentralWidget(self.figure_canvas)
        self.mpl_toolbar.setVisible(self.options['toolbar'])
        # similar to figure canvas
        self.mpl_toolbar = NavigationToolbar(self.figure_canvas, self)
        # self.mpl_toolbar.setVisible(True)
        self.mainwidgetmdi.addToolBar(Qt.BottomToolBarArea, self.mpl_toolbar)

        self.fbtn = QPushButton(self.style().standardIcon(QStyle.SP_FileDialogDetailedView), '', self.figure_canvas)
        self.fbtn.setToolTip('Show or Hide the Navigation Toolbar')
        self.fbtn.toggled.connect(self.mpl_toolbar.setVisible)
        self.fbtn.setCheckable(True)
        self.fbtn.setChecked(False)

    def draw(self, title, tight_layout=True):
        if tight_layout:
            self.fig.tight_layout()
            self.figure_canvas.draw()
        self.setWindowTitle(title)

    def setup_text(self, ax, options, key='', title='', xlabel='', ylabel='Energy (kcal/mol)'):
        ax.set_title(title, fontdict={'fontsize': options['general_options']['fontsize']['suptitle']})
        ax.set_xlabel(xlabel, fontdict={'fontsize': options['general_options']['fontsize']['y-label']})
        ax.set_ylabel(ylabel, fontdict={'fontsize': options['general_options']['fontsize']['y-label']})
        for label in ax.get_xticklabels():
            label.set_rotation(options[key]['axes']['x-rotation'])
            if options[key]['axes']['x-rotation'] < 0:
                label.set_horizontalalignment('left')
            else:
                label.set_horizontalalignment('right')
            label.set_fontsize(options['general_options']['fontsize']['x-ticks'])
        for label in ax.get_yticklabels():
            label.set_rotation(options[key]['axes']['y-rotation'])
            if options[key]['axes']['y-rotation'] < 0:
                label.set_horizontalalignment('left')
            else:
                label.set_horizontalalignment('right')
            label.set_fontsize(options['general_options']['fontsize']['y-ticks'])

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.button.setChecked(False)

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

                label = 'ΔH'
                if 'IE' in self.item.item_name:
                    label = '-TΔS'
                if 'NMODE' in self.item.item_name or 'QH' in self.item.item_name:
                    self.line_plot = sns.lineplot(data=self.data.line_plot_dat, x='frames', color='black',
                                                  linewidth=0.7,
                                                  y='Entropy', label=label, ax=self.mpl_canvas.axes)
                else:
                    self.line_plot = sns.lineplot(data=self.data.line_plot_dat, x='frames', color='black', linewidth=0.7,
                             y='Energy', label=label, ax=self.mpl_canvas.axes)

                self.line_plot.xaxis.set_major_locator(mticker.MaxNLocator(nbins=10))

                title = self.item.chart_title + '(P.f)'
                if 'IE' in self.item.item_name:
                    self.mpl_canvas.axes.set_ylabel('Interaction Entropy (kcal/mol)')
                    self.mpl_canvas.axes.axvline(self.item.ie[0][0], ls='--', color='r')
                    self.mpl_canvas.axes.hlines(self.item.ie[1], xmin=self.item.ie[0][0], xmax=self.item.ie[0][1], ls='--',
                                                color='g')
                else:
                    self.mpl_canvas.axes.set_ylabel('Energy (kcal/mol)')

            if Charts.ROLLING in self.options['chart_type']:
                if ('IE' not in self.item.item_name and 'NMODE' not in self.item.item_name and 'QH' not in
                        self.item.item_name):
                    if len(self.data.line_plot_dat['Energy']) > 50:
                        self.data.line_plot_dat['movav'] = self.data.line_plot_dat['Energy'].rolling(
                            int(0.1 * len(self.data.line_plot_dat['frames']))).mean()
                        self.movav_plot = sns.lineplot(data=self.data.line_plot_dat, x='frames', color='red', linewidth=0.8,
                                                        y='movav', label='Mov. Av.', ax=self.mpl_canvas.axes)

                # self.setWindowTitle(title)
            if Charts.BAR in self.options['chart_type']:
                self.item.bp_subw = self
                if self.bar_plot:
                    self.bar_plot.cla()
                    self.bar_plot.clear()
                self.bar_plot = sns.barplot(data=self.data.bar_plot_dat, ci="sd", errwidth=1, ax=self.mpl_canvas.axes)
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

                xticklabels = self.data.heatmap_plot_dat.columns.tolist()
                window = int(len(xticklabels) / 10)
                if type(xticklabels[0]) == str:
                    self.heatmap_plot = sns.heatmap(self.data.heatmap_plot_dat, ax=self.mpl_canvas.axes, center=0,
                                                yticklabels=self.data.heatmap_plot_dat.index.tolist(),
                                                xticklabels=window,
                                                cmap='seismic', cbar_kws={'label': 'Energy (kcal/mol)'})
                    title = self.item.chart_title.replace('[Per-residue]', '[Per-wise]')
                else:
                    self.heatmap_plot = sns.heatmap(self.data.heatmap_plot_dat, ax=self.mpl_canvas.axes, center=0,
                                                yticklabels=self.data.heatmap_plot_dat.index.tolist(),
                                                xticklabels=window,
                                                cmap='seismic', cbar_kws={'label': 'Energy (kcal/mol)'})
                    title = self.item.chart_title + '(P.f)'  # Fixme: no frames from correlation

            # if Charts.RELPLOT in self.options['chart_type']:
            #     self.item.hmp_subw = self
            #     # Draw each cell as a scatter point with varying size and color
            #     self.relplot = sns.relplot(
            #         data=self.data.heatmap_plot_dat,
            #         x="Residues", y="Pair", hue="Energy", size="Energy",
            #         palette="seismic",
            #         hue_norm=(-100, 100),
            #         edgecolor=".7",
            #         height=10,
            #         sizes=(50, 200),
            #         size_norm=(0, 1),
            #     )
            #
            #     # Tweak the figure to finalize
            #     self.relplot.set(xlabel="", ylabel="", aspect="equal")
            #     self.relplot.despine(left=True, bottom=True)
            #     self.relplot.ax.margins(.02)
            #     for label in self.relplot.ax.get_xticklabels():
            #         label.set_rotation(90)
            #     for artist in self.relplot.legend.legendHandles:
            #         artist.set_edgecolor(".7")
            #
            #     self.mpl_canvas.figure = self.relplot.fig

            if self.item.item_name != 'IE':
                self.mpl_canvas.axes.invert_yaxis()


            self.mpl_canvas.axes.set_title(title + '\n' + self.item.chart_subtitle)
            self.cursor = Cursor(self.mpl_canvas.axes, useblit=True, color='black', linewidth=0.5, ls='--')
            self.mpl_canvas.draw()
            self.mpl_canvas.figure.tight_layout()
            self.mpl_canvas.draw()
            self.setWindowTitle(title)
