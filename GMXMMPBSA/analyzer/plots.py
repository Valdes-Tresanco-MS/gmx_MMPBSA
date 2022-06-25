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
    from PyQt6.QtCore import *
    from PyQt6.QtGui import *
    from PyQt6.QtWidgets import *
except:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *


import itertools

from matplotlib.colors import ListedColormap
import matplotlib.backend_bases
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.widgets import Cursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import linregress, pearsonr, spearmanr

# sns.set_theme()
from GMXMMPBSA.analyzer.chartsettings import Palettes
from GMXMMPBSA.analyzer.style import logo
from GMXMMPBSA.analyzer.utils import bar_label

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
        self.save_format = 'png'
        self.save_dpi = 300
        self.filename = 'image'

    def update_options(self, options):
        if 'save-format' in options:
            self.save_format = options['save-format']
        if 'dpi-save' in options:
            self.save_dpi = options['dpi-save']
        if 'filename' in options:
            self.filename = options['filename'].replace(' | ', '_').replace(' ', '_')

    def _init_toolbar(self):
        pass

    def save_figure(self, *args):

        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(filetypes.items())
        default_filetype = self.canvas.get_default_filetype()

        startpath = os.path.expanduser(matplotlib.rcParams['savefig.directory'])
        filename = f'{self.filename}.{self.save_format}'
        start = os.path.join(startpath, filename)
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(f'*.{ext}' for ext in exts)
            filter = f'{name} ({exts_list})'
            if self.save_format in exts:
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
                self.canvas.figure.savefig(fname, dpi=self.save_dpi)
            except Exception as e:
                QMessageBox.critical(
                    self, "Error saving file", str(e),
                    QMessageBox.Ok, QMessageBox.NoButton)


class ChartsBase(QMdiSubWindow):
    def __init__(self, button: QToolButton, options: dict = None, item_parent=None):
        super(ChartsBase, self).__init__()
        self.setMinimumSize(400, 400)
        self.options = options
        self.item_parent = item_parent
        self.setWindowIcon(QIcon(logo))

        self.mainwidgetmdi = QMainWindow()  # must be QMainWindow to handle the toolbar
        sns.set_theme(style=self.options[('General', 'theme')])
        self.plot = None
        self.frange = []  # Frames range with which it was created
        self.button = button
        self.setWidget(self.mainwidgetmdi)

    def set_cw(self, fig=None):
        # we create the figure canvas here because the fig parameter must be defined and depend on what kind of
        # chart want to make
        fig = fig or Figure(dpi=self.options[('General', 'figure-format', 'dpi-plot')])
        self.figure_canvas = FigureCanvas(fig)
        self.fig = self.figure_canvas.figure
        self.mainwidgetmdi.setCentralWidget(self.figure_canvas)
        # similar to figure canvas
        self.mpl_toolbar = NavigationToolbar(self.figure_canvas, self)
        self.mpl_toolbar.setVisible(self.options['General', 'toolbar'])
        self.mpl_toolbar.update_options({'save-format': self.options[('General', 'figure-format', 'save-format')],
                                         'dpi-save': self.options[('General', 'figure-format', 'dpi-save')],
                                         'filename': self.options['subtitle']})
        self.mainwidgetmdi.addToolBar(Qt.ToolBarArea.BottomToolBarArea, self.mpl_toolbar)

        self.fbtn = QPushButton(self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogDetailedView), '',
                                self.figure_canvas)
        self.fbtn.setToolTip('Show or Hide the Navigation Toolbar')
        self.fbtn.toggled.connect(self.mpl_toolbar.setVisible)
        self.fbtn.setCheckable(True)
        self.fbtn.setChecked(False)

    def draw(self):
        self.fig.tight_layout()
        self.figure_canvas.draw()
        QGuiApplication.restoreOverrideCursor()

    def setup_text(self, ax, options, key='', title='', xlabel='', ylabel='Energy (kcal/mol)'):
        key_list = [key] if isinstance(key, str) else key
        ax.set_title(title, fontdict={'fontsize': options[tuple(key_list + ['fontsize', 'suptitle'])]})
        if 'Line Plot' in key_list:
            ax.legend(loc="lower right", prop={'size': options[tuple(key_list + ['fontsize', 'legend'])]})
        if xlabel:
            ax.set_xlabel(xlabel, fontdict={'fontsize': options[tuple(key_list + ['fontsize', 'x-label'])]})
        ax.set_ylabel(ylabel, fontdict={'fontsize': options[tuple(key_list + ['fontsize', 'y-label'])]})
        for label in ax.get_xticklabels():
            if options.get(tuple(key_list + ['axes', 'x-rotation'])):
                label.set_rotation(options[tuple(key_list + ['axes', 'x-rotation'])])
                if options[tuple(key_list + ['axes', 'x-rotation'])] < 0:
                    label.set_horizontalalignment('left')
                else:
                    label.set_horizontalalignment('right')
            label.set_fontsize(options[tuple(key_list + ['fontsize', 'x-ticks'])])
        for label in ax.get_yticklabels():
            if options.get(tuple(key_list + ['axes', 'y-rotation'])):
                label.set_rotation(options[tuple(key_list + ['axes', 'y-rotation'])])
                if options[tuple(key_list + ['axes', 'y-rotation'])] < 0:
                    label.set_horizontalalignment('left')
                else:
                    label.set_horizontalalignment('right')
            label.set_fontsize(options[tuple(key_list + ['fontsize', 'y-ticks'])])

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.button.setChecked(False)


class LineChart(ChartsBase):
    def __init__(self, data, button: QToolButton, options: dict = None, item_parent=None):
        super(LineChart, self).__init__(button, options, item_parent)

        self.data = pd.read_parquet(data) if isinstance(data, str) else data
        self.data = self.data.iloc[:,0] if self.data.columns.size == 1 else self.data
        # figure canvas definition
        self.set_cw()
        self.fig.set_size_inches(options[('Line Plot', 'figure', 'width')],
                                 options[('Line Plot', 'figure', 'height')])

        self.axes = self.fig.subplots(1, 1)
        if options.get('iec2'):
            self.line_plot_ax = sns.lineplot(data=self.data['AccIntEnergy'],
                                             color=rgb2rgbf(options[('Line Plot', 'line-color')]),
                                             linewidth=options[('Line Plot', 'line-width')],
                                             ax=self.axes,
                                             label='IE(all)')
            # plot the ie segment
            ie_color = rgb2rgbf(options[('Bar Plot', 'IE/C2 Entropy', 'ie-color')])
            sns.lineplot(data=self.data['ie'], color=ie_color, ax=self.axes, label='IE(selected)')
        else:
            self.line_plot_ax = sns.lineplot(data=self.data, color=rgb2rgbf(options[('Line Plot', 'line-color')]),
                                             linewidth=options[('Line Plot', 'line-width')],
                                             label=self.data.name,
                                             ax=self.axes)
            if options[('Line Plot', 'Rolling average', 'show')]:
                min_periods = 1 if options[('Line Plot', 'Rolling average', 'first_obs')] == 'start' else None
                window = options[('Line Plot', 'Rolling average', 'window')]
                moving_avg = sns.lineplot(data=self.data.rolling(window, min_periods=min_periods).mean(),
                                          color=rgb2rgbf(options[('Line Plot', 'Rolling average', 'color')]),
                                          linewidth=options[('Line Plot', 'Rolling average', 'width')],
                                          label='Mov. Av.',
                                          ls=options[('Line Plot', 'Rolling average', 'style')],
                                          ax=self.axes)
        self.cursor = Cursor(self.axes, useblit=True, color='black', linewidth=0.5, ls='--')

        self.setWindowTitle(options['subtitle'])
        self.update_config(options)
        self.draw()

    def check(self):
        pass

    def update_config(self, options):
        self.fig.suptitle(f"{options['title']}\n{options['subtitle']}",
                          fontsize=options[('Line Plot', 'fontsize', 'title')])
        self.line_plot_ax.xaxis.set_major_locator(
            mticker.MaxNLocator(nbins=options[('Line Plot', 'axes', 'num-xticks')],
                                integer=True))
        self.line_plot_ax.yaxis.set_major_locator(
            mticker.MaxNLocator(nbins=options[('Line Plot', 'axes', 'num-yticks')]))

        self.setup_text(self.line_plot_ax, options, key='Line Plot', xlabel=self.data.index.name)

        self.draw()


class BarChart(ChartsBase):
    def __init__(self, data, button: QToolButton, options: dict = None, item_parent=None):
        super(BarChart, self).__init__(button, options, item_parent)

        self.data = pd.read_parquet(data) if isinstance(data, str) else data
        self.bar_labels = []
        # figure canvas definition
        self.set_cw()
        self.fig.set_size_inches(options[('Bar Plot', 'figure', 'width')],
                                 options[('Bar Plot', 'figure', 'height')])
        self.bar_frames = False
        palette = (sns.color_palette(options[('Bar Plot', 'palette')], n_colors=self.data.columns.size)
                   if options[('Bar Plot', 'use-palette')] else None)
        if options.get('groups') and options[('Bar Plot', 'subplot-components')]:
            self.axes = self.fig.subplots(1, len(options['groups']), sharey=True,
                                          gridspec_kw={'width_ratios': [len(x) for x in options['groups'].values()]})
            if options[('Bar Plot', 'axes', 'y-inverted')]:
                self.axes[0].invert_yaxis()
            s = 0
            for c, g in enumerate(options['groups']):
                df = self.data[options['groups'][g]]
                bar_plot_ax = sns.barplot(x=df.columns,
                                          y=df.loc['Average'],
                                          yerr=df.loc[options[('Bar Plot', 'error-line', 'representation')]],
                                          palette=palette[s: s + len(options['groups'][g])] if palette else palette,
                                          color=rgb2rgbf(options[('Bar Plot', 'color')]),
                                          error_kw=dict(
                                              ecolor=rgb2rgbf(options[('Bar Plot', 'error-line', 'color')]),
                                              capsize=options[('Bar Plot', 'error-line', 'cap-size')],
                                              elinewidth=options[('Bar Plot', 'error-line', 'width')]
                                          ),
                                          ax=self.axes[c]
                                          )
                s += len(options['groups'][g])
                if options[('Bar Plot', 'scale-yaxis')]: # and options['scalable']:
                    bar_plot_ax.set_yscale('symlog')
                if options[('Bar Plot', 'bar-label', 'show')]:
                    bl = bar_label(bar_plot_ax, bar_plot_ax.containers[1],
                                   size=options[('Bar Plot', 'bar-label', 'fontsize')],
                                   fmt='%.2f',
                                   padding=options[('Bar Plot', 'bar-label', 'padding')],
                                   label_type=options[('Bar Plot', 'bar-label', 'label_type')])
                    self.bar_labels.append(bl)
                ylabel = '' if c != 0 else 'Energy (kcal/mol)'
                self.setup_text(bar_plot_ax, options, key='Bar Plot', title=g, ylabel=ylabel)
                setattr(self, f'cursor{c}', Cursor(bar_plot_ax, useblit=True, color='black', linewidth=0.5, ls='--'))
                bar_plot_ax.set_xticklabels(self._set_xticks(bar_plot_ax, options[('Bar Plot', 'remove-molid')]))
        else:
            self.axes = self.fig.subplots(1, 1)
            if options.get('iec2'):
                ie_color = rgb2rgbf(options[('Bar Plot', 'IE/C2 Entropy', 'ie-color')])
                r_sigma = rgb2rgbf(options[('Bar Plot', 'IE/C2 Entropy', 'sigma-color', 'reliable')])
                nr_sigma = rgb2rgbf(options[('Bar Plot', 'IE/C2 Entropy', 'sigma-color', 'non-reliable')])
                palette = [ie_color, r_sigma if self.data['sigma'].loc[['Average']].values[0] < 3.6 else nr_sigma]

            bar_plot_ax = sns.barplot(x=self.data.columns,
                                      y=self.data.loc['Average'],
                                      yerr=self.data.loc[options[('Bar Plot', 'error-line', 'representation')]],
                                      error_kw=dict(
                                          ecolor=rgb2rgbf(options[('Bar Plot', 'error-line', 'color')]),
                                          capsize=options[('Bar Plot', 'error-line', 'cap-size')],
                                          elinewidth=options[('Bar Plot', 'error-line', 'width')]
                                      ),
                                      ax=self.axes,
                                      palette=palette,
                                      color=rgb2rgbf(options[('Bar Plot', 'color')]),
                                      )
            if options[('Bar Plot', 'axes', 'y-inverted')] and not options.get('iec2'):
                bar_plot_ax.invert_yaxis()
            if options[('Bar Plot', 'bar-label', 'show')]:
                bl = bar_label(bar_plot_ax, bar_plot_ax.containers[1],
                               size=options[('Bar Plot', 'bar-label', 'fontsize')],
                               fmt='%.2f',
                               padding=options[('Bar Plot', 'bar-label', 'padding')],
                               label_type=options[('Bar Plot', 'bar-label', 'label_type')])
                self.bar_labels.append(bl)
            self.setup_text(bar_plot_ax, options, key='Bar Plot')
            self.cursor = Cursor(bar_plot_ax, useblit=True, color='black', linewidth=0.5, ls='--')
            bar_plot_ax.set_xticklabels(self._set_xticks(bar_plot_ax, options[('Bar Plot', 'remove-molid')]))
            self.bar_frames = 'frames' in self.data
        self.setWindowTitle(options['subtitle'])
        self.update_config(options)
        self.draw()

    @staticmethod
    def _set_xticks(axe, remove_molid):
        xlabels = []
        for xlabel in axe.get_xticklabels():
            if remove_molid and xlabel.get_text()[:2] in ['R:', 'L:']:
                xlabels.append(xlabel.get_text()[2:])
            else:
                xlabels.append(xlabel.get_text())
        return xlabels

    def update_config(self, options):
        self.fig.suptitle(f"{options['title']}\n{options['subtitle']}",
                          fontsize=options[('Bar Plot', 'fontsize', 'title')])
        if isinstance(self.axes, np.ndarray):
            for c, g in enumerate(options['groups']):
                bar_plot_ax = self.axes[c]
                if options[('Bar Plot', 'bar-label', 'show')]:
                    for blabels in self.bar_labels:
                        for label in blabels:
                            label.set_fontsize(options[('Bar Plot', 'bar-label', 'fontsize')])
                ylabel = '' if c != 0 else 'Energy (kcal/mol)'
                self.setup_text(bar_plot_ax, options, key='Bar Plot', title=g, ylabel=ylabel)
        else:
            bar_plot_ax = self.axes
            if options[('Bar Plot', 'scale-yaxis')]:
                bar_plot_ax.set_yscale('symlog')
            else:
                bar_plot_ax.set_yscale('linear')
            if options[('Bar Plot', 'bar-label', 'show')]:
                for blabels in self.bar_labels:
                    for label in blabels:
                        label.set_fontsize(options[('Bar Plot', 'bar-label', 'fontsize')])
            self.setup_text(self.axes, options, key='Bar Plot', xlabel='Frames' if self.bar_frames else '')
        self.draw()


class HeatmapChart(ChartsBase):
    def __init__(self, data, button: QToolButton, options: dict = None, item_parent=None):
        super(HeatmapChart, self).__init__(button, options, item_parent)

        self.data = pd.read_parquet(data) if isinstance(data, str) else data

        self.heatmap_type = 2 if data.columns.is_numeric() else 1

        fig_width = (options[('Heatmap Plot', 'figure', 'width-per-wise')] if self.heatmap_type == 1 else
                     options[('Heatmap Plot', 'figure', 'width-per-residue')])
        fig_height = options[('Heatmap Plot', 'figure', 'height')]
        x_rotation = (options[('Heatmap Plot', 'Per-wise', 'x-rotation')] if self.heatmap_type == 1 else
                      options[('Heatmap Plot', 'Per-residue', 'x-rotation')])
        cmap = (Palettes.get_colormap(options[('Heatmap Plot', 'Per-wise', 'palette')]) if self.heatmap_type == 1 else
                Palettes.get_colormap(options[('Heatmap Plot', 'Per-residue', 'palette')]))
        nxticks = (1 if self.heatmap_type == 1 or
                        self.data.columns.size < options[('Heatmap Plot', 'Per-residue', 'num-xticks')]
                   else self.data.columns.size // options[('Heatmap Plot', 'Per-residue', 'num-xticks')])

        annotation = options[('Heatmap Plot', 'Per-wise', 'annotation')] if self.heatmap_type == 1 else False

        if options[('Heatmap Plot', 'highlight-components')]:
            mheatmap = MHeatmap(data=self.data,
                                figsize=(fig_width, fig_height),
                                dpi=options[('General', 'figure-format', 'dpi-plot')],
                                heatmap_type=self.heatmap_type,
                                rec_color=rgb2rgbf(options[('Heatmap Plot', 'receptor-color')]),
                                lig_color=rgb2rgbf(options[('Heatmap Plot', 'ligand-color')]),
                                show_legend=options[('Heatmap Plot', 'legend')],
                                remove_molid=options[('Heatmap Plot', 'remove-molid')],
                                leg_fontsize=options[('Heatmap Plot', 'fontsize', 'legend')],
                                xticks_fontsize=options[('Heatmap Plot', 'fontsize', 'x-ticks')],
                                xlabel_fontsize=options[('Heatmap Plot', 'fontsize', 'x-label')],
                                x_rotation=x_rotation,
                                num_xticks=nxticks,
                                yticks_fontsize=options[('Heatmap Plot', 'fontsize', 'y-ticks')],
                                y_rotation=options[('Heatmap Plot', 'y-rotation')],
                                colorbar_label_fontsize=options[('Heatmap Plot', 'fontsize', 'colorbar-label')],
                                colorbar_ticks_fontsize=options[('Heatmap Plot', 'fontsize', 'colorbar-ticks')],
                                # colorbar_num_ticks=options[('Heatmap Plot', 'fontsize', 'colorbar-ticks')],
                                cmap=cmap,
                                annot=annotation,
                                annot_fs=options[('Heatmap Plot', 'fontsize', 'annotation')]
                                )

            # figure canvas definition
            self.set_cw(mheatmap.fig)
            self.axes = mheatmap.ax_heatmap
        else:
            # figure canvas definition
            self.set_cw()
            self.axes = self.fig.subplots(1, 1)
            heatmap_ax = sns.heatmap(self.data, ax=self.axes, center=0,
                                     xticklabels=nxticks, #cbar_ax=self.ax_cbar,
                                     yticklabels=1, #cbar_ax=self.ax_cbar,
                                     # cbar_kws=colorbar_kws,
                                     cmap=cmap, annot=annotation, fmt=".2f",
                                     annot_kws={'size': options[('Heatmap Plot', 'fontsize', 'annotation')]},
                                     # yticklabels=data.index.tolist(),
                                     # xticklabels=window,
                                     cbar_kws={'label': 'Energy (kcal/mol)',
                                               'ticks': mticker.MaxNLocator(
                                                   nbins=options[('Heatmap Plot', 'Per-residue', 'num-xticks')])})
            ytl = heatmap_ax.get_yticklabels()
            for ylabel in ytl:
                ylabel.set_fontsize(options[('Heatmap Plot', 'fontsize', 'y-ticks')])
                ylabel.set_rotation(options[('Heatmap Plot', 'y-rotation')])
                if options[('Heatmap Plot', 'remove-molid')]:
                    ylabel.set_text(ylabel.get_text()[2:])
            xtl = heatmap_ax.get_xticklabels()
            for xlabel in xtl:
                xlabel.set_fontsize(options[('Heatmap Plot', 'fontsize', 'x-ticks')])
                xlabel.set_rotation(x_rotation)
                if options[('Heatmap Plot', 'remove-molid')] and self.heatmap_type == 1:
                    xlabel.set_text(xlabel.get_text()[2:])
                # xlabel.set_text(xlabel.get_text()[2:])

            if options[('Heatmap Plot', 'remove-molid')]:
                heatmap_ax.set_yticklabels(ytl)
                if self.heatmap_type == 1:
                    heatmap_ax.set_xticklabels(xtl)
            # heatmap_ax.yaxis.set_ticks_position('right')
            # heatmap_ax.yaxis.set_label_position('right')
            if self.heatmap_type == 2:
                heatmap_ax.set_xlabel(self.data.columns.name,
                                      fontdict={'fontsize': options[('Heatmap Plot', 'fontsize', 'x-label')]})

        self.cursor = Cursor(self.axes, useblit=True, color='black', linewidth=0.5, ls='--')

        self.setWindowTitle(options['subtitle'])
        self.update_config(options)
        self.draw()

    def update_config(self, options):
        self.fig.suptitle(f"{options['title']}\n{options['subtitle']}",
                          fontsize=options[('Heatmap Plot', 'fontsize', 'title')])
        xlabel = self.data.columns.name if self.heatmap_type == 2 else ''
        self.setup_text_hm(self.axes, options, xlabel=xlabel)
        self.draw()

    def setup_text_hm(self, ax, options, title='', xlabel=''):
        ax.set_title(title, fontdict={'fontsize': options[('Heatmap Plot', 'fontsize', 'suptitle')]})
        if xlabel:
            ax.set_xlabel(xlabel, fontdict={'fontsize': options[('Heatmap Plot', 'fontsize', 'x-label')]})
        key = ['Per-residue'] if self.heatmap_type == 2 else ['Per-wise']
        xrot = options[tuple(['Heatmap Plot'] + key + ['x-rotation'])]
        for label in ax.get_xticklabels():
            label.set_rotation(xrot)
            if xrot < 0:
                label.set_horizontalalignment('left')
            else:
                label.set_horizontalalignment('right')
            label.set_fontsize(options[('Heatmap Plot', 'fontsize', 'x-ticks')])
        yrot = options[('Heatmap Plot', 'y-rotation')]
        for label in ax.get_yticklabels():
            label.set_rotation(yrot)
            label.set_fontsize(options[('Heatmap Plot', 'fontsize', 'y-ticks')])


class MHeatmap:
    def __init__(self, data: pd.DataFrame, figsize=None, dpi=100, heatmap_type=1, rec_color=None, lig_color=None,
                 show_legend=True, remove_molid=True, leg_fontsize=10, xticks_fontsize=10, xlabel_fontsize=11,
                 x_rotation=0, num_xticks=10, yticks_fontsize=10, y_rotation=0, colorbar_ticks_fontsize=9,
                 colorbar_label_fontsize=10, colorbar_num_ticks='auto', cmap='seismic', annot=False, annot_fs=8):
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
        @param colorbar_num_ticks: num of ticks
        @param cmap: color map
        @param annot: show annotations
        @param annot_fs: annotations fontsize
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
        self.annot_fs = annot_fs

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
        self.plot_matrix({'ticklocation': 'left', 'label': 'Energy (kcal/mol)',
                          'ticks': mticker.MaxNLocator(nbins=colorbar_num_ticks)})

    @staticmethod
    def _index_color(x, rec_color, lig_color):
        return tuple(rec_color) if x.startswith('R') else tuple(lig_color)

    @staticmethod
    def color_list_to_matrix_and_cmap(colors, axis=0):
        """Turns a list of colors into a numpy matrix and matplotlib colormap

        These arguments can now be plotted using heatmap(matrix, cmap)
        and the provided colors will be plotted.

        Parameters
        ----------
        colors : list of matplotlib colors
            Colors to label the rows or columns of a dataframe.
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

        cmap = ListedColormap(all_colors)
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
        h = sns.heatmap(self.data, ax=self.ax_heatmap, xticklabels=self.num_xticks, yticklabels=1,
                        cbar_ax=self.ax_cbar, cbar_kws=colorbar_kws, cmap=self.cmap, center=0, annot=self.annot,
                        fmt=".2f", annot_kws={'size': self.annot_fs},**kws)
        yticks = [x[2:] for x in self.data.index] if self.remove_molid else self.data.index
        self.ax_heatmap.set_yticklabels(yticks, fontdict=dict(fontsize=self.yticks_fontsize),
                                        rotation=self.y_rotation)

        if self.remove_molid and self.heatmap_type == 1:
            xtl = self.ax_heatmap.get_xticklabels()
            for xlabel in xtl:
                xlabel.set_fontsize(self.xticks_fontsize)
                xlabel.set_rotation(self.x_rotation)
                if self.remove_molid and self.heatmap_type == 1:
                    xlabel.set_text(xlabel.get_text()[2:])
            self.ax_heatmap.set_xticklabels(xtl)
        self.ax_heatmap.yaxis.set_ticks_position('right')
        self.ax_heatmap.yaxis.set_label_position('right')
        if self.heatmap_type == 2:
            self.ax_heatmap.set_xlabel('Frames', fontdict={'fontsize': self.xlabel_fontsize})


class RegChart(ChartsBase):
    def __init__(self, data: pd.DataFrame, button: QToolButton, options: dict = None, item_parent=None):
        super(RegChart, self).__init__(button, options, item_parent)

        ci = options[('Regression Plot', 'Conf. Interval')] or None
        dist_type = {'violin': sns.violinplot, 'hist': sns.histplot, 'box': sns.boxplot, 'kde': sns.kdeplot}

        line_kws = {'lw': options[('Regression Plot', 'line-width')],
                    'color': rgb2rgbf(options[('Regression Plot', 'line-color')])}

        reg_options = dict(ci=ci,
                           scatter_kws={'s': options[('Regression Plot', 'Scatter', 'marker-size')],
                                        'color': rgb2rgbf(options[('Regression Plot', 'Scatter', 'color')])},
                           line_kws=line_kws)

        # get correlation coefficients
        pearson, ppvalue = pearsonr(data['ExpΔG'], data['Average'])
        slope, intercept, r_value, p_value, std_err = linregress(data['ExpΔG'], data['Average'])
        spearman, spvalue = spearmanr(data['ExpΔG'], data['Average'])

        if options[('Regression Plot', 'Distribution', 'show')]:
            args = {}
            if options[('Regression Plot', 'Distribution', 'type')] == 'hist':
                args = {'kde': options[('Regression Plot', 'Distribution', 'Histogram', 'kde')],
                        'fill': options[('Regression Plot', 'Distribution', 'Histogram', 'fill-bars')],
                        'element': options[('Regression Plot', 'Distribution', 'Histogram', 'element')]}

            regplot = sns.JointGrid(data=data, x="ExpΔG", y='Average')
            regplot.plot_joint(sns.regplot, **reg_options)
            regplot.plot_marginals(dist_type[options[('Regression Plot', 'Distribution', 'type')]],
                                   color=rgb2rgbf(options[('Regression Plot', 'Scatter', 'color')]), **args)
            # figure canvas definition
            self.set_cw(regplot.fig)
            self.axes = regplot.ax_joint
        else:
            # figure canvas definition
            self.set_cw()
            self.axes = self.fig.subplots(1, 1)
            sns.regplot(data=data, x="ExpΔG", y='Average', ax=self.axes, **reg_options)

        if options[('Regression Plot', 'Scatter', 'error-line', 'show')]:
            error = options[('Regression Plot', 'Scatter', 'error-line', 'representation')]
            self.axes.errorbar(x=data["ExpΔG"], y=data['Average'], yerr=data[error],
                               fmt='none',
                               # zorder=1,
                               elinewidth=options[('Regression Plot', 'Scatter', 'error-line', 'width')],
                               capsize=options[('Regression Plot', 'Scatter', 'error-line', 'cap-size')],
                               color=rgb2rgbf(options[('Regression Plot', 'Scatter', 'error-line', 'color')]))

        pearson_leg = mpatches.Patch(color='white', label=f'Pearson  = {pearson:.2f}  p-value = {ppvalue:.3f}')
        spearman_leg = mpatches.Patch(color='gray', label=f'Spearman = {spearman:.2f}  p-value = {spvalue:.3f}')
        sign = '+' if intercept > 0 else '-'
        equation = Line2D([0], [0], **line_kws, label=f'y = {slope:.2f}x {sign} {abs(intercept):.2f}')
        handles = []
        if options[('Regression Plot', 'pearson')]:
            handles.append(pearson_leg)
        if options[('Regression Plot', 'spearman')]:
            handles.append(spearman_leg)
        if options[('Regression Plot', 'equation')]:
            handles.append(equation)
        if handles:
            self.axes.legend(handles=handles, handlelength=0, handletextpad=0, fancybox=True, frameon=True,
                             prop={'weight': 'bold', 'size': options[('Regression Plot', 'fontsize', 'legend')],
                                   'family': 'monospace'})
        self.cursor = Cursor(self.axes, useblit=True, color='black', linewidth=0.5, ls='--')
        self.setWindowTitle(options['subtitle'])
        self.update_config(options)


    def update_config(self, options):
        self.fig.suptitle(f"{options['title']}\n{options['subtitle']}",
                          fontsize=options[('Regression Plot', 'fontsize', 'title')])
        self.setup_text(self.axes, options, key='Regression Plot', xlabel=r'$ΔG_{Experimental} (kcal/mol)$',
                        ylabel=r'$ΔG_{Calculated} (kcal/mol)$')
        self.draw()


class OutputFiles(QMdiSubWindow):
    def __init__(self, text, button):
        super(OutputFiles, self).__init__()
        self.setMinimumSize(400, 400)
        self.textedit = QTextEdit(self)
        self.textedit.setReadOnly(True)
        self.setWidget(self.textedit)

        font = QFont("Monospace")
        font.setStyleHint(QFont.TypeWriter)
        self.textedit.setFont(font)
        self.textedit.setPlainText(''.join(text))

        self.button = button

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.button.setChecked(False)


class Tables(QMdiSubWindow):
    def __init__(self, df: pd.DataFrame, button: QToolButton, options: dict = None, summary=False):
        super(Tables, self).__init__()
        self.setMinimumSize(400, 400)
        self.container = QWidget()
        self.setWidget(self.container)
        self.container_layout = QVBoxLayout(self.container)
        self.container_layout.setContentsMargins(0, 0, 0, 0)
        self.item_parent = None
        self.options = options
        self.setWindowTitle(self.options['table_name'])
        self.table_name = self.options['table_name'].replace(' | ', '_')
        self.button = button

        self.table = QTableWidget(self)
        self.container_layout.addWidget(self.table)
        self._df = df.round(2)

        if summary:
            self.df_list = self._df.values.tolist()[1:]
            labels = list(self._df.columns)
            rows = len(self._df.index) - 1
            cols = len(self._df.columns)
        else:
            temp_df_list = [x.split(',') for x in self._df.to_csv().split('\n')]
            self.df_list = temp_df_list[1:]
            labels = temp_df_list[0]
            rows = len(self.df_list) - 1
            cols = len(self.df_list[0])

        self.table.setColumnCount(cols)
        self.table.setRowCount(rows)
        self.table.setHorizontalHeaderLabels(labels)

        for r, row in enumerate(self.df_list):
            for c, col in enumerate(row):
                text = f'{col:.2f}' if isinstance(col, float) else str(col)
                item = QTableWidgetItem(text)
                if c == 0:
                    if not summary:
                        item.setTextAlignment(Qt.AlignmentFlag.AlignRight)
                    if col in (['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'GGAS', 'UB', 'IMP',
                               'CMAP', 'ESCF'] +
                               ['ΔBOND', 'ΔANGLE', 'ΔDIHED', 'ΔVDWAALS', 'ΔEEL', 'Δ1-4 VDW', 'Δ1-4 EEL', 'ΔGGAS',
                                'ΔUB', 'ΔIMP', 'ΔCMAP', 'ΔESCF'] +
                               ['ΔΔBOND', 'ΔΔANGLE', 'ΔΔDIHED', 'ΔΔVDWAALS', 'ΔΔEEL', 'ΔΔ1-4 VDW', 'ΔΔ1-4 EEL',
                                'ΔΔGGAS', 'ΔΔUB', 'ΔΔIMP', 'ΔΔCMAP', 'ΔΔESCF']):
                        item.setBackground(QColor('#e6b8b8'))
                    elif col in (['EGB', 'ESURF',  'EPB', 'ENPOLAR', 'EDISPER', 'POLAR SOLV', 'APOLAR SOLV', 'ERISM',
                                 'GSOLV'] +
                                 ['ΔEGB', 'ΔESURF',  'ΔEPB', 'ΔENPOLAR', 'ΔEDISPER', 'ΔPOLAR SOLV', 'ΔAPOLAR SOLV',
                                  'ΔERISM', 'ΔGSOLV'] +
                                 ['ΔΔEGB', 'ΔΔESURF',  'ΔΔEPB', 'ΔΔENPOLAR', 'ΔΔEDISPER', 'ΔΔPOLAR SOLV',
                                  'ΔΔAPOLAR SOLV', 'ΔΔERISM', 'ΔΔGSOLV']):
                        item.setBackground(QColor('#97c5cc'))
                    else:
                        item.setBackground(QColor('white'))
                if c > 0:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight)
                self.table.setItem(r, c, item)

        self.table.installEventFilter(self)

        self.save_btn = QPushButton('Save')
        self.save_btn.clicked.connect(self.saveToFile)
        self.save_format = QComboBox()
        self.save_format.addItem('*.csv')
        # self.save_format.addItem('*.xlsx')
        self.save_layout = QHBoxLayout()
        self.save_layout.addStretch(10)
        self.save_layout.addWidget(self.save_format)
        self.save_layout.addWidget(self.save_btn)
        self.container_layout.addLayout(self.save_layout)

        h_header = self.table.horizontalHeader()
        h_header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        h_header.setStretchLastSection(True)
        v_header = self.table.verticalHeader()
        v_header.hide()

    def saveToFile(self):
        filter = self.save_format.currentText()
        fileName, _ = QFileDialog.getSaveFileName(self, "Save table content to file",
                                                  self.table_name + filter[1:],
                                                  "CSV (*.csv);;Excel (*.xlsx);;All Files (*)",
                                                  "CSV (*.csv)" if filter == '*.csv' else "Excel (*.xlsx)")
        if not fileName:
            return
        if filter == '*.csv':
            self._df.to_csv(fileName)
        else:
            self._df.to_excel(fileName)

    def eventFilter(self, source, event):
        if (event.type() == QEvent.Type.KeyPress and event.matches(QKeySequence.StandardKey.Copy)):
            self._copySelection()
            return True
        return super(Tables, self).eventFilter(source, event)

    def _copySelection(self):
        selection = self.table.selectedIndexes()
        if not selection:
            return
        rows = sorted(index.row() for index in selection)
        columns = sorted(index.column() for index in selection)
        rowcount = rows[-1] - rows[0] + 1
        colcount = columns[-1] - columns[0] + 1
        table = [[''] * colcount for _ in range(rowcount)]
        header = [[''] * colcount]
        for index in selection:
            row = index.row() - rows[0]
            column = index.column() - columns[0]
            table[row][column] = index.data() or ''
            if not header[0][column]:
                header[0][column] = self.table.horizontalHeaderItem(column).text()
        temp = ['\t'.join(x) + '\n' for x in header] + ['\t'.join(x) + '\n' for x in table]
        text = ''.join(temp)
        QApplication.clipboard().setText(text)

    def closeEvent(self, closeEvent: QCloseEvent) -> None:
        self.button.setChecked(False)
