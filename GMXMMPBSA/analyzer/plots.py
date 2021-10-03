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


class MHeatmap:
    def __init__(self, data: pd.DataFrame, figsize=None, dpi=100, heatmap_type=1, rec_color=None, lig_color=None,
                 show_legend=True, remove_molid=True, leg_fontsize=10, xticks_fontsize=10, xlabel_fontsize=11,
                 x_rotation=0, num_xticks=10, yticks_fontsize=10, y_rotation=0, colorbar_ticks_fontsize=9,
                 colorbar_label_fontsize=10, cmap='seismic', annot=False):
        """
        This class is based on the seaborn clustermap. We created it to have more control over the components and
        especially for the management and support of the tight_layout
        @param data: pandas.Dataframe
        @param figsize: Managed according to the plot decomp scheme: per-wise [symmetrical matrix] and per-residue [
                        frames in xaxis]
        @param dpi: figure dpi
        @param heatmap_type: 1 per-wise, 2 per-residue
        @param rec_color: color for receptor residues identification
        @param lig_color: color for ligand residues identification
        @param show_legend: show the receptor and ligand identification colorbar
        @param remove_molid: remove the identify label (e.g. R:C:VAL:33 --> C:VAL:33)
        @param leg_fontsize: Legend font-size
        @param xticks_fontsize: xticks font-size. Is always defined
        @param xlabel_fontsize: xlabel font-size. Only for per-residue plot
        @param x_rotation: xlabels rotation
        @param num_xticks: number of xlabels to show
        @param yticks_fontsize: yticks font-size. Is always defined
        @param y_rotation: ylabels rotation
        @param colorbar_ticks_fontsize: colorbar ticks font-size
        @param colorbar_label_fontsize: colorbar label font-size
        @param cmap: color map
        @param annot: show annotations
        """
        super(MHeatmap, self).__init__()
        ratios = {'legend': 0.025, 'colorbar': 0.025, 'ws1': 0.03, 'color_row_col': 0.015, 'ws2': 0.002}
        pos = {'legend': [0, 4], 'colorbar': [-1, 0], 'color_row': [-1, 2], 'color_col': [-3, 4], 'heatmap': [-1, 4]}
        self.data = data
        self.heatmap_type = heatmap_type
        self.remove_molid = remove_molid
        self.xticks_fontsize = xticks_fontsize
        self.xlabel_fontsize = xlabel_fontsize
        self.x_rotation = x_rotation
        self.num_xticks = num_xticks
        self.yticks_fontsize = yticks_fontsize
        self.y_rotation = y_rotation
        self.cmap = cmap
        self.annot = annot

        leg_ratios = [ratios['legend'], ratios['ws1']]
        row_color_ratios = [ratios['color_row_col'], ratios['ws2']]
        h_ratios = []
        if show_legend:
            h_ratios += leg_ratios
        if heatmap_type == 1:  # per-wise heatmap
            h_ratios += row_color_ratios
        h_ratios.append(1 - sum(h_ratios))
        self.gs = gridspec.GridSpec(len(h_ratios), 5,
                                    hspace=0,
                                    wspace=0,
                                    height_ratios=h_ratios,
                                    width_ratios=[0.025, 0.02, 0.015, 0.002, 0.938])
        self.fig = Figure(figsize=figsize, dpi=dpi)
        self.row_col_colors = self.data.index.map(lambda x: self._index_color(x, rec_color, lig_color))

        if show_legend:
            self.ax_legend = self.fig.add_subplot(self.gs[pos['legend'][0], pos['legend'][1]], label='Legend')
            color = [rec_color, lig_color]
            for c, label in enumerate(['Receptor', 'Ligand']):
                self.ax_legend.bar(0, 0, color=color[c], label=label, linewidth=0)
            self.ax_legend.legend(loc="center", ncol=2, prop={'size': leg_fontsize})
            self.ax_legend.set_facecolor('white')
            self.ax_legend.get_xaxis().set_visible(False)
            self.ax_legend.get_yaxis().set_visible(False)
            sns.despine(ax=self.ax_legend, left=True, bottom=True)

        if heatmap_type == 1:
            self.ax_col_colors = self.fig.add_subplot(self.gs[pos['color_col'][0], pos['color_col'][1]],
                                                      label='Column-molid')

        self.ax_cbar = self.fig.add_subplot(self.gs[-1, 0], label='Colorbar')
        self.ax_cbar.yaxis.tick_left()
        self.ax_cbar.yaxis.set_label_text('Energy (kcal/mol)', fontdict={'size': colorbar_label_fontsize})
        # self.ax_cbar.yaxis.set_label_position('left')
        for label in self.ax_cbar.get_yticklabels():
            label.set_fontsize(colorbar_ticks_fontsize)

        divider = make_axes_locatable(self.ax_cbar)
        leg_axestop = divider.append_axes("top", size="20%", pad=0.2)
        leg_axesbot = divider.append_axes("bottom", size="20%", pad=0.2)
        leg_axestop.axes.set_visible(False)
        leg_axesbot.axes.set_visible(False)

        self.ax_row_colors = self.fig.add_subplot(self.gs[-1, -3], label='Row-molid')
        self.ax_heatmap = self.fig.add_subplot(self.gs[-1, -1], label='Heatmap')
        self.ax_heatmap.yaxis.tick_right()

        self.plot_colors()
        self.plot_matrix({'ticklocation': 'left', 'label': 'Energy (kcal/mol)'})

    @staticmethod
    def _index_color(x, rec_color, lig_color):
        if x.startswith('R'):
            return tuple(rec_color)
        else:
            return tuple(lig_color)

    @staticmethod
    def color_list_to_matrix_and_cmap(colors, axis=0):
        """Turns a list of colors into a numpy matrix and matplotlib colormap

        These arguments can now be plotted using heatmap(matrix, cmap)
        and the provided colors will be plotted.

        Parameters
        ----------
        colors : list of matplotlib colors
            Colors to label the rows or columns of a dataframe.
        ind : list of ints
            Ordering of the rows or columns, to reorder the original colors
            by the clustered dendrogram order
        axis : int
            Which axis this is labeling

        Returns
        -------
        matrix : numpy.array
            A numpy array of integer values, where each corresponds to a color
            from the originally provided list of colors
        cmap : matplotlib.colors.ListedColormap

        """
        # check for nested lists/color palettes.
        # Will fail if matplotlib color is list not tuple
        if any(issubclass(type(x), list) for x in colors):
            all_colors = set(itertools.chain(*colors))
            n = len(colors)
            m = len(colors[0])
        else:
            all_colors = set(colors)
            n = 1
            m = len(colors)
            colors = [colors]
        color_to_value = {col: i for i, col in enumerate(all_colors)}

        matrix = np.array([color_to_value[c]
                           for color in colors for c in color])

        shape = (n, m)
        matrix = matrix.reshape(shape)
        if axis == 0:
            # row-side:
            matrix = matrix.T

        cmap = mpl.colors.ListedColormap(all_colors)
        return matrix, cmap

    def plot_colors(self, **kws):
        # Plot the row colors
        matrix, cmap = self.color_list_to_matrix_and_cmap(self.row_col_colors, axis=0)
        sns.heatmap(matrix, cmap=cmap, cbar=False, ax=self.ax_row_colors, xticklabels=False, yticklabels=False, **kws)
        # Plot the column colors
        if self.heatmap_type == 1:
            matrix, cmap = self.color_list_to_matrix_and_cmap(self.row_col_colors, axis=1)
            sns.heatmap(matrix, cmap=cmap, cbar=False, ax=self.ax_col_colors, xticklabels=False, yticklabels=False,
                        **kws)

    def plot_matrix(self, colorbar_kws, **kws):

        nxticks = 1 if self.heatmap_type == 1 else self.data.columns.size // self.num_xticks
        h = sns.heatmap(self.data, ax=self.ax_heatmap, xticklabels=nxticks, cbar_ax=self.ax_cbar, cbar_kws=colorbar_kws,
                        cmap=self.cmap, center=0, annot=self.annot, fmt=".2f", **kws)
        ytl = self.ax_heatmap.get_yticklabels()
        for ylabel in ytl:
            ylabel.set_fontsize(self.yticks_fontsize)
            ylabel.set_rotation(self.y_rotation)
            if self.remove_molid:
                ylabel.set_text(ylabel.get_text()[2:])
        xtl = self.ax_heatmap.get_xticklabels()
        for xlabel in xtl:
            xlabel.set_fontsize(self.xticks_fontsize)
            xlabel.set_rotation(self.x_rotation)
            if self.remove_molid and self.heatmap_type == 1:
                xlabel.set_text(xlabel.get_text()[2:])

        if self.remove_molid:
            self.ax_heatmap.set_yticklabels(ytl)
            if self.heatmap_type == 1:
                self.ax_heatmap.set_xticklabels(xtl)
        self.ax_heatmap.yaxis.set_ticks_position('right')
        self.ax_heatmap.yaxis.set_label_position('right')
        if self.heatmap_type == 2:
            self.ax_heatmap.set_xlabel('Frames', fontdict={'fontsize': self.xlabel_fontsize})


