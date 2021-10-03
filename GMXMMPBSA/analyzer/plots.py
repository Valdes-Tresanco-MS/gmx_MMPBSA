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


class LineChart(ChartsBase):
    def __init__(self, data: pandas.DataFrame, button: QToolButton, options: dict = None, data2=None):
        super(LineChart, self).__init__(button, options)
        # figure canvas definition
        self.set_cw()
        self.fig.set_size_inches(options['line_options']['figure']['width'],
                                 options['line_options']['figure']['height'])
        axes = self.fig.subplots(1, 1)
        line_plot_ax = sns.lineplot(data=data, color=rgb2rgbf(options['line_options']['line-color']),
                                    linewidth=options['line_options']['line-width'], ax=axes)
        line_plot_ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=options['line_options']['axes']['num-xticks'],
                                                                 integer=True))
        line_plot_ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=options['line_options']['axes']['num-yticks']))
        self.setup_text(line_plot_ax, options, key='line_options', xlabel='Frames')
        if data2 is not None:
            bar_options = options['line_options']['Interaction Entropy']['bar-plot']
            axins4 = inset_axes(axes,
                                width=f"{bar_options['width']}%",
                                height=f"{bar_options['height']}%",
                                bbox_to_anchor=(
                                    options['line_options']['Interaction Entropy']['bbox_to_anchor']['x-pos'],
                                    options['line_options']['Interaction Entropy']['bbox_to_anchor']['y-pos'],
                                    0.5, 0.5),
                                bbox_transform=axes.transAxes, loc=4)
            # plot the ie segment
            ie_color = rgb2rgbf(options['line_options']['Interaction Entropy']['ie-color'])
            sns.lineplot(data=data2['data'], color=ie_color, ax=axes)
            r_sigma = rgb2rgbf(options['line_options']['Interaction Entropy']['sigma-color']['reliable'])
            nr_sigma = rgb2rgbf(options['line_options']['Interaction Entropy']['sigma-color']['non-reliable'])
            colors = [ie_color, r_sigma if np.all(data2.loc[:, ['sigma']].mean() < 3.6) else nr_sigma]
            ax1 = sns.barplot(data=data2, ax=axins4, palette=colors)
            ax1.tick_params(labelsize=bar_options['axes-fontsize'])
            ax1.bar_label(ax1.containers[0], size=bar_options['bar-label-fontsize'], fmt='%.2f',
                          padding=bar_options['bar-label-padding'])
            ax1.set(xticklabels=[f"ie\n(last\n {len(data2['data'])} frames)", "σ(Int.\nEnergy)"])
            sns.despine(ax=ax1)
            ax1.spines['bottom'].set_color('darkgray')
            ax1.spines['left'].set_color('darkgray')
            ax1.tick_params(color='black', direction='out', length=6, width=2)

        self.cursor = Cursor(axes, useblit=True, color='black', linewidth=0.5, ls='--')
        self.fig.suptitle(f"{options['chart_title']}\n{options['chart_subtitle']}",
                          fontsize=options['general_options']['fontsize']['title'])
        self.draw(options['chart_subtitle'])


class BarChart(ChartsBase):
    def __init__(self, data: pandas.DataFrame, button: QToolButton, options: dict = None):
        super(BarChart, self).__init__(button, options)

        # figure canvas definition
        self.set_cw()
        self.fig.set_size_inches(options['bar_options']['figure']['width'],
                                 options['bar_options']['figure']['height'])
        if 'groups' in options and options['bar_options']['subplot-components']:
            axes = self.fig.subplots(1, len(options['groups']), sharey=True,
                                     gridspec_kw={'width_ratios': [len(x) for x in options['groups'].values()]})
            palette = (sns.color_palette(options['bar_options']['palette'])
                       if options['bar_options']['use-palette'] else None)
            s = 0
            for c, g in enumerate(options['groups']):
                bar_plot_ax = sns.barplot(data=data[options['groups'][g]], ci="sd",
                                          palette=palette[s: s + len(options['groups'][g])] if palette else palette,
                                          color=rgb2rgbf(options['bar_options']['color']),
                                          errwidth=1, ax=axes[c])
                s += len(options['groups'][g])
                if options['bar_options']['scale-big-values'] and options['scalable']:
                    bar_plot_ax.set_yscale('symlog')
                ylabel = '' if c != 0 else 'Energy (kcal/mol)'
                self.setup_text(bar_plot_ax, options, key='bar_options', title=g, ylabel=ylabel)
                setattr(self, f'cursor{c}', Cursor(bar_plot_ax, useblit=True, color='black', linewidth=0.5, ls='--'))

        else:
            axes = self.fig.subplots(1, 1)
            bar_plot_ax = sns.barplot(data=data, ci="sd", errwidth=1, ax=axes)
            self.setup_text(bar_plot_ax, options, key='bar_options')
            self.cursor = Cursor(bar_plot_ax, useblit=True, color='black', linewidth=0.5, ls='--')

        self.fig.suptitle(f"{options['chart_title']}\n{options['chart_subtitle']}",
                          fontsize=options['general_options']['fontsize']['title'])
        self.draw(options['chart_subtitle'])


class HeatmapChart(ChartsBase):
    def __init__(self, data: pandas.DataFrame, button: QToolButton, options: dict = None):
        super(HeatmapChart, self).__init__(button, options)

        heatmap_type = int(all(data.columns == data.index)) if data.columns.size == data.index.size else 2
        fig_width = (options['heatmap_options']['figure']['width-per-wise'] if heatmap_type == 1 else
                     options['heatmap_options']['figure']['width-per-residue'])
        fig_height = options['heatmap_options']['figure']['height']
        x_rotation = (options['heatmap_options']['Per-wise']['x-rotation'] if heatmap_type == 1 else
                      options['heatmap_options']['Per-residue']['x-rotation'])
        cmap = (Palettes.get_palette(options['heatmap_options']['Per-wise']['palette']) if heatmap_type == 1 else
                Palettes.get_palette(options['heatmap_options']['Per-residue']['palette']))

        if options['heatmap_options']['highlight-components']:
            mheatmap = MHeatmap(data=data,
                                figsize=(fig_width, fig_height),
                                dpi=self.options['dpi']['plot'],
                                heatmap_type=heatmap_type,
                                rec_color=rgb2rgbf(options['heatmap_options']['receptor-color']),
                                lig_color=rgb2rgbf(options['heatmap_options']['ligand-color']),
                                show_legend=options['heatmap_options']['legend'],
                                remove_molid=options['heatmap_options']['remove-molid'],
                                leg_fontsize=self.options['fontsize']['legend'],
                                xticks_fontsize=self.options['fontsize']['x-ticks'],
                                xlabel_fontsize=self.options['fontsize']['x-label'],
                                x_rotation=x_rotation,
                                num_xticks=options['heatmap_options']['Per-residue']['num-xticks'],
                                yticks_fontsize=self.options['fontsize']['y-ticks'],
                                y_rotation=options['heatmap_options']['y-rotation'],
                                colorbar_label_fontsize=self.options['fontsize']['colorbar-label'],
                                colorbar_ticks_fontsize=self.options['fontsize']['colorbar-ticks'],
                                cmap=cmap,
                                annot=options['heatmap_options']['Per-wise']['annotation']
                                )

            # figure canvas definition
            self.set_cw(mheatmap.fig)

            self.draw(options['chart_subtitle'], tight_layout=True)
        else:
            # figure canvas definition
            self.set_cw()
            axes = self.fig.subplots(1, 1)
            nxticks = 1 if self.heatmap_type == 1 else self.data.columns.size // self.num_xticks
            # heatmap_ax = sns.heatmap(data, ax=axes, center=0,
            #                          xticklabels=nxticks, cbar_ax=self.ax_cbar,
            #                          cbar_kws=colorbar_kws,
            #                          cmap=self.cmap, center=0, annot=self.annot, fmt=".2f"
            #                          # yticklabels=data.index.tolist(),
            #                          # xticklabels=window,
            #                          cmap='seismic', cbar_kws={'label': 'Energy (kcal/mol)'})
            self.cursor = Cursor(axes, useblit=True, color='black', linewidth=0.5, ls='--')

            self.draw(options['chart_subtitle'])
        self.fig.suptitle(f"{options['chart_title']}\n{options['chart_subtitle']}",
                          fontsize=options['general_options']['fontsize']['title'])
    @staticmethod
    def get_mol(x: str, rec_color, lig_color):
        if x.startswith('R:'):
            return rgb2rgbf(rec_color)
        else:
            return rgb2rgbf(lig_color)

    def make_chart(self):
        """
        :param graph_type: 1: line plot, 2: bar plot, 3: heatmap plot
        :return:
        """
        pass


