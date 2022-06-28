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
except Exception:
    from PyQt5.QtWidgets import *
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *

from .utils import com2str, energy2pdb_pml
from .style import (result_file_icon, bar_plot_icon, line_plot_icon, heatmap_plot_icon, pymol_icon, result_table_icon)


class SpacerItem(QToolButton):
    def __init__(self, parent=None):
        super(SpacerItem, self).__init__(parent)
        self.setDisabled(True)
        self.setContentsMargins(0, 0, 0, 0)
        self.setStyleSheet("QToolButton { /* mimic the look of the QToolButton with MenuButtonPopup */ "
                           "padding-right: 15px; /* make way for the popup button */}")


class TableActionBtn(QWidget):
    def __init__(self, parent=None):
        super(TableActionBtn, self).__init__(parent)
        self.c_widget_layout = QHBoxLayout(self)
        self.c_widget_layout.setContentsMargins(0, 0, 0, 0)
        self.reg_chart_action = QToolButton(self)
        self.c_widget_layout.addWidget(self.reg_chart_action)
        self.reg_chart_action.setIcon(QIcon(line_plot_icon))
        self.reg_chart_action.setText('Regression Chart')
        self.reg_chart_action.setCheckable(True)
        self.reg_chart_action.setContentsMargins(0, 0, 0, 0)


class CustomCorrItem(QTableWidgetItem):
    def __init__(self, text=False, app=None, model=None, keys_path=None):
        super(CustomCorrItem, self).__init__()

        self.app = app
        self.model = model
        self.keys_path = keys_path
        self.reg_sw = None
        if text:
            self.setText(model.upper())
        self.title = f'Correlation Using {model.upper()} model'
        self.subtitle = "$ΔG_{Experimental} vs ΔG_{Calculated}$"
        self.c_widget = TableActionBtn()
        self.c_widget.reg_chart_action.toggled.connect(self.plotting_reg)

    def define_button(self, row, col):
        self.tableWidget().setCellWidget(row, col, self.c_widget)

    def plotting_reg(self, state):
        from GMXMMPBSA.analyzer.plots import RegChart
        self.app.correlation_tableWidget.clearSelection()

        options = self.app.correlation['chart_options'].get_settings()

        options.update({'title': self.title, 'subtitle': self.subtitle})
        changes = self.app.correlation['chart_options'].changes
        plot_data = self.app.correlation['items_data'][self.keys_path][0]
        datachange = self.app.correlation['items_data'][self.keys_path][1]

        if state:
            self.setSelected(True)
            if not self.reg_sw or datachange or changes:
                if len(plot_data.index) < 4 or plot_data['Average'].count() < 4 or plot_data['ExpΔG'].count() < 4:
                    QMessageBox.critical(self.tableWidget(),
                                         'Unable to calculate correlation',
                                         'More than three valid systems are needed to calculate the correlation.'
                                         'Please check that the selected systems contain both experimental and '
                                         'calculated valid ΔG.',
                                         QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    self.c_widget.reg_chart_action.setChecked(False)
                    return
                QGuiApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
                self.reg_sw = RegChart(plot_data.dropna(), button=self.c_widget.reg_chart_action, options=options,
                                       item_parent=self)
                self.app.correlation['items_data'][self.keys_path][1] = False
                self.app.mdi.addSubWindow(self.reg_sw)
            self.reg_sw.show()
        elif self.reg_sw:
            self.app.mdi.activatePreviousSubWindow()
            self.reg_sw.close()


class CustomItem(QTreeWidgetItem):
    def __init__(self, parent, stringlist, app=None, title='Binding Free Energy', subtitle='', system_index=1,
                 buttons=(), keys_path=None, part=''):
        super(CustomItem, self).__init__(parent, stringlist)

        self.system_index = system_index
        self.app = app
        self.title = title
        self.subtitle = f'({part.capitalize()}) {subtitle}' if part else subtitle
        self.item_name = stringlist[0]
        self.buttons = buttons
        self.properties = {}
        self.part = part

        self.keys_path = keys_path

        self.lp_subw = None
        self.bp_subw = None
        self.hmp_subw = None
        self.pymol_process = None
        self.pymol_data_change = False
        self.pymol_current_palette = None
        self.bfactor_pml = None
        self.output_file_subw = None
        self.decomp_output_file_subw = None
        self.line_table_subw = None
        self.bar_table_subw = None
        self.heatmap_table_subw = None

        # changes
        self.frange = []
        self.line_change = False
        self.bar_change = False
        self.heatmap_change = False

        self.changed = False

        self.tb = QToolBar()
        self.tb.setStyleSheet("QToolBar {padding: 0, 20, 0, 20;}")

        self.tb.setIconSize(QSize(16, 16))
        self.tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.item_charts = []

        self.btn_group = QButtonGroup()
        self.btn_group.setExclusive(False)
        self.btn_group.buttonToggled.connect(self.fn_btn_group)

        self.mark_all = QCheckBox()
        self.mark_all.setToolTip('Mark all actions at the same time')
        self.mark_all.setTristate(True)
        self.mark_all.setCheckable(True)
        self.mark_all.stateChanged.connect(self.fn_mark_all)

        self.charts_action = {
            1: self._define_line_chart_btn,
            2: self._define_bar_chart_btn,
            3: self._define_heatmap_chart_btn,
            4: self._define_vis_btn,
        }

    def changes(self, line, bar, heatmap):
        self.line_change = line
        self.bar_change = bar
        self.heatmap_change = heatmap

    def _define_line_chart_btn(self):
        line_menu = QMenu()
        line_menu.setTitle('Line charts menu')
        self.line_table_action = line_menu.addAction('Show in table')
        self.line_table_action.setCheckable(True)
        self.line_table_action.toggled.connect(self._show_line_table)

        self.line_chart_action = QToolButton()
        self.line_chart_action.setIcon(QIcon(line_plot_icon))
        self.line_chart_action.setText('Line Chart')
        self.line_chart_action.setCheckable(True)
        self.line_chart_action.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        self.line_chart_action.setContentsMargins(0, 0, 0, 0)
        self.line_chart_action.setMenu(line_menu)
        self.btn_group.addButton(self.line_chart_action, 1)

        return self.line_chart_action

    def _show_line_table(self, state):
        from GMXMMPBSA.analyzer.plots import Tables
        line_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['line_plot_data'][0]
        options = {'table_name': self.subtitle}
        self.app.treeWidget.clearSelection()
        if state:
            self.setSelected(True)
            self.line_table_subw = Tables(line_plot_data, self.line_table_action, options)
            self.app.mdi.addSubWindow(self.line_table_subw)
            self.line_table_subw.show()
        elif self.line_table_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.line_table_subw.close()

    def _define_bar_chart_btn(self):
        bar_menu = QMenu()
        bar_menu.setTitle('Bar charts menu')
        self.bar_table_action = bar_menu.addAction('Show in table')
        self.bar_table_action.setCheckable(True)
        self.bar_table_action.toggled.connect(self._show_bar_table)

        self.bar_chart_action = QToolButton()
        self.bar_chart_action.setIcon(QIcon(bar_plot_icon))
        self.bar_chart_action.setText('Bar Chart')
        self.bar_chart_action.setCheckable(True)
        self.bar_chart_action.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        self.bar_chart_action.setContentsMargins(0, 0, 0, 0)
        self.bar_chart_action.setMenu(bar_menu)
        self.btn_group.addButton(self.bar_chart_action, 2)

        return self.bar_chart_action

    def _show_bar_table(self, state):
        from GMXMMPBSA.analyzer.plots import Tables
        self.app.treeWidget.clearSelection()
        options = {'table_name': self.subtitle}
        bar_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][0]
        if state:
            self.setSelected(True)
            if not self.bar_table_subw:
                self.bar_table_subw = Tables(bar_plot_data, self.bar_table_action, options)
                self.app.mdi.addSubWindow(self.bar_table_subw)
            self.bar_table_subw.show()
        elif self.bar_table_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.bar_table_subw.close()

    def _define_heatmap_chart_btn(self):
        heatmap_menu = QMenu()
        heatmap_menu.setTitle('heatmap charts menu')
        self.heatmap_table_action = heatmap_menu.addAction('Show in table')
        self.heatmap_table_action.setCheckable(True)
        self.heatmap_table_action.toggled.connect(self._show_heatmap_table)

        self.heatmap_chart_action = QToolButton()
        self.heatmap_chart_action.setIcon(QIcon(heatmap_plot_icon))
        self.heatmap_chart_action.setText('Heatmap Chart')
        self.heatmap_chart_action.setCheckable(True)
        self.heatmap_chart_action.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        self.heatmap_chart_action.setContentsMargins(0, 0, 0, 0)
        self.heatmap_chart_action.setMenu(heatmap_menu)
        self.btn_group.addButton(self.heatmap_chart_action, 3)

        return self.heatmap_chart_action

    def _show_heatmap_table(self, state):
        from GMXMMPBSA.analyzer.plots import Tables
        self.app.treeWidget.clearSelection()
        options = {'table_name': self.subtitle}
        heatmap_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['heatmap_plot_data'][0]
        if state:
            self.setSelected(True)
            if not self.heatmap_table_subw:
                self.heatmap_table_subw = Tables(heatmap_plot_data.T, self.heatmap_table_action, options)
                self.app.mdi.addSubWindow(self.heatmap_table_subw)
            self.heatmap_table_subw.show()
        elif self.heatmap_table_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.heatmap_table_subw.close()

    def _define_vis_btn(self):
        heatmap_menu = QMenu()
        heatmap_menu.setTitle('Visualization menu')
        # self.heatmap_table_action = heatmap_menu.addAction('Show in text')
        # self.heatmap_table_action.setCheckable(True)
        # self.heatmap_table_action.toggled.connect(self._show_heatmap_table)

        self.vis_action = QToolButton()
        self.vis_action.setIcon(QIcon(pymol_icon))
        self.vis_action.setText('Visualization')
        self.vis_action.setCheckable(True)
        self.vis_action.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        self.vis_action.setContentsMargins(0, 0, 0, 0)
        # self.vis_action.setMenu(heatmap_menu)
        self.btn_group.addButton(self.vis_action, 4)

        return self.vis_action

    def _define_option_button(self):
        options_menu = QMenu()
        options_menu.setTitle('System Options')
        self.output_action = options_menu.addAction('Show Output file')
        self.output_action.setCheckable(True)
        self.output_action.toggled.connect(self._show_output_file)
        if self.app.systems[self.system_index]['namespace'].INPUT['decomprun']:
            self.decomp_output_action = options_menu.addAction('Show Decomp Output file')
            self.decomp_output_action.setCheckable(True)
            self.decomp_output_action.toggled.connect(self._show_decomp_output_file)

        self.options_button = QToolButton()
        self.options_button.setIcon(QIcon(result_file_icon))
        self.options_button.setText('Results files')
        self.options_button.setToolTip('Results files')
        self.options_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.options_button.setContentsMargins(0, 0, 0, 0)
        self.options_button.setMenu(options_menu)

        return self.options_button

    def _show_output_file(self, state):
        from GMXMMPBSA.analyzer.plots import OutputFiles
        self.app.treeWidget.clearSelection()
        if state:
            self.setSelected(True)
            if not self.output_file_subw:
                self.output_file_subw = OutputFiles(
                    self.app.systems[self.system_index]['namespace'].INFO['output_file'],
                    self.output_action)
                self.app.mdi.addSubWindow(self.output_file_subw)
            self.output_file_subw.show()
        elif self.output_file_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.output_file_subw.close()

    def _show_decomp_output_file(self, state):
        from GMXMMPBSA.analyzer.plots import OutputFiles
        self.app.treeWidget.clearSelection()
        if state:
            self.setSelected(True)
            if not self.decomp_output_file_subw:
                self.decomp_output_file_subw = OutputFiles(
                    self.app.systems[self.system_index]['namespace'].INFO['decomp_output_file'],
                    self.decomp_output_action)
                self.app.mdi.addSubWindow(self.decomp_output_file_subw)
            self.decomp_output_file_subw.show()
        elif self.decomp_output_file_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.decomp_output_file_subw.close()

    def _define_result_table_btn(self):

        self.result_table_action = QToolButton()
        self.result_table_action.setIcon(QIcon(result_table_icon))
        self.result_table_action.setText('Result Table')
        self.result_table_action.setCheckable(True)
        self.result_table_action.setContentsMargins(0, 0, 0, 0)
        self.btn_group.addButton(self.result_table_action, 5)

        return self.result_table_action

    def fn_mark_all(self, state):
        if state == Qt.CheckState.PartiallyChecked:
            pass
        elif state == Qt.CheckState.Checked:
            for x in self.btn_group.buttons():
                x.setChecked(True)
        else:
            for x in self.btn_group.buttons():
                x.setChecked(False)

    def fn_btn_group(self, btn, checked):
        all_checked = all(x.isChecked() for x in self.btn_group.buttons())
        all_unchecked = not any(x.isChecked() for x in self.btn_group.buttons())
        oldState = self.mark_all.blockSignals(True)
        if all_checked:
            self.mark_all.setCheckState(Qt.CheckState.Checked)
        elif all_unchecked:
            self.mark_all.setCheckState(Qt.CheckState.Unchecked)
        else:
            self.mark_all.setCheckState(Qt.CheckState.PartiallyChecked)
        self.mark_all.blockSignals(oldState)

        if self.btn_group.id(btn) == 1:
            self.plotting_line(checked)
        elif self.btn_group.id(btn) == 2:
            self.plotting_bar(checked)
        elif self.btn_group.id(btn) == 3:
            self.plotting_heatmap(checked)
        elif self.btn_group.id(btn) == 4:
            self.visualizing(checked)
        elif self.btn_group.id(btn) == 5:
            self.result_table(checked)

    def plotting_line(self, state):
        from GMXMMPBSA.analyzer.plots import LineChart
        self.app.treeWidget.clearSelection()

        options = self.app.systems[self.system_index]['chart_options'].get_settings()
        options.update({'title': self.title, 'subtitle': self.subtitle})
        changes = self.app.systems[self.system_index]['chart_options'].changes
        line_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['line_plot_data'][0]
        datachange = self.app.systems[self.system_index]['items_data'][self.keys_path]['line_plot_data'][2]
        options.update(self.app.systems[self.system_index]['items_data'][self.keys_path]['line_plot_data'][1])
        if state:
            self.setSelected(True)
            line_change3 = (changes['line_action'] == 3 or changes['line_ie_action'] == 3)
            line_change1 = (changes['line_action'] == 1 or changes['line_ie_action'] == 1)

            if not self.lp_subw or datachange or line_change3:
                QGuiApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
                self.lp_subw = LineChart(line_plot_data, self.line_chart_action, options=options, item_parent=self)
                self.app.systems[self.system_index]['items_data'][self.keys_path]['line_plot_data'][2] = False
                self.app.mdi.addSubWindow(self.lp_subw)
            elif line_change1:
                self.lp_subw.update_config(options)
            self.lp_subw.show()
        elif self.lp_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.lp_subw.close()

    def plotting_bar(self, state):
        from GMXMMPBSA.analyzer.plots import BarChart
        self.app.treeWidget.clearSelection()
        options = self.app.systems[self.system_index]['chart_options'].get_settings()
        options.update({'title': self.title, 'subtitle': self.subtitle})
        changes = self.app.systems[self.system_index]['chart_options'].changes

        bar_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][0]
        datachange = self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][2]
        options.update(self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][1])

        if state:
            self.setSelected(True)
            if not self.bp_subw or datachange or changes['bar_action'] == 3:
                QGuiApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
                self.bp_subw = BarChart(bar_plot_data, self.bar_chart_action, options=options, item_parent=self)
                self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][2] = False
                self.app.mdi.addSubWindow(self.bp_subw)
            elif changes['bar_action'] == 1:
                self.bp_subw.update_config(options)
            self.bp_subw.show()
        elif self.bp_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.bp_subw.close()

    def plotting_heatmap(self, state):
        from GMXMMPBSA.analyzer.plots import HeatmapChart
        self.app.treeWidget.clearSelection()
        options = self.app.systems[self.system_index]['chart_options'].get_settings()
        options.update({'title': self.title, 'subtitle': self.subtitle})
        changes = self.app.systems[self.system_index]['chart_options'].changes
        heatmap_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['heatmap_plot_data'][0]
        datachange = self.app.systems[self.system_index]['items_data'][self.keys_path]['heatmap_plot_data'][2]
        if state:
            self.setSelected(True)
            if not self.hmp_subw or datachange or changes['heatmap_action'] == 3:
                QGuiApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
                self.hmp_subw = HeatmapChart(heatmap_plot_data, self.heatmap_chart_action, options=options,
                                             item_parent=self)
                self.app.systems[self.system_index]['items_data'][self.keys_path]['heatmap_plot_data'][2] = False
                self.app.mdi.addSubWindow(self.hmp_subw)
            elif changes['heatmap_action'] == 1:
                self.hmp_subw.update_config(options)
            self.hmp_subw.show()
        elif self.hmp_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.hmp_subw.close()

    def visualizing(self, checked):
        from GMXMMPBSA.analyzer.chartsettings import Palettes
        self.app.treeWidget.clearSelection()
        QGuiApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
        import os
        changes = self.app.systems[self.system_index]['chart_options'].changes
        options = self.app.systems[self.system_index]['chart_options'].get_settings()
        if options[('Visualization', 'palette')] == 'auto':
            palette = options[('Heatmap Plot', 'Per-residue', 'palette')]
        else:
            palette = options[('Visualization', 'palette')]

        pymol_options = {'colors': Palettes.get_palette(palette).colors}
        for o in options:
            if o[0] != 'Visualization' or o == ('Visualization', 'palette') or o[-1] == 'default':
                continue
            if o[1] == 'background':
                pymol_options['bg_rgb'] = options[o]
            elif o[1] == 'cartoon_side_chain_helper':
                pymol_options[o[1]] = int(options[o])
            else:
                pymol_options[o[1]] = options[o]

        self.pymol_data_change = bool(changes['visualization_action'] or self.pymol_current_palette != palette)
        self.pymol_current_palette = palette

        if checked:
            self.setSelected(True)
            if pymol_path := [
                os.path.join(path, 'pymol')
                for path in os.environ["PATH"].split(os.pathsep)
                if os.path.exists(os.path.join(path, 'pymol'))
                   and os.access(os.path.join(path, 'pymol'), os.X_OK)
            ]:
                pymol = pymol_path[0]

            else:
                QMessageBox.critical(self.app, 'PyMOL not found!', 'PyMOL not found!. Make sure PyMOL is in the PATH.',
                                     QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                self.vis_action.setChecked(False)
                return
            if not self.pymol_process:
                self.pymol_process = QProcess()
                self.app.pymol_p_list.append([self.pymol_process, self])
            elif self.pymol_process.state() == QProcess.ProcessState.Running:
                QMessageBox.critical(self.app, 'This PyMOL instance already running!',
                                     'This PyMOL instance already running! Please, close it to open a new PyMOL '
                                     'instance', QMessageBox.StandardButton.Ok)
                return
            self.bfactor_pml = self._e2pdb(pymol_options)

            self.pymol_process.start(pymol, [self.bfactor_pml.as_posix()])
            self.pymol_process.finished.connect(self.pymol_finished)
            self.pymol_process.started.connect(QGuiApplication.restoreOverrideCursor)
            self.pymol_data_change = False
        elif self.pymol_process.state() == QProcess.ProcessState.Running:
            self.pymol_process.terminate()
            self.pymol_process.waitForFinished(3000)

    def pymol_finished(self):
        self.vis_action.setChecked(False)

    def _e2pdb(self, options):
        com_pdb = self.app.systems[self.system_index]['namespace'].INFO['COM_PDB']
        bfactor_pml = self.app.systems[self.system_index]['path'].parent.joinpath('bfactor.pml')
        output_path = self.app.systems[self.system_index]['path'].parent.joinpath(
            f"{self.app.systems[self.system_index]['name']}_{self.part}_energy2bfactor.pdb")
        qpd = QProgressDialog('Generating modified pdb and open it in PyMOL...', 'Abort', 0, 2, self.app)
        qpd.setWindowModality(Qt.WindowModality.WindowModal)
        qpd.setMinimumDuration(1500)

        for i in range(2):
            qpd.setValue(i)
            if qpd.wasCanceled():
                break
            if i == 0:
                com_pdb_str = com2str(com_pdb)
                bar_plot_data = self.app.systems[self.system_index]['items_data'][self.keys_path]['bar_plot_data'][0]
                temp_dict = bar_plot_data.iloc[0].to_dict()
                res_dict = {k[2:]: value for k, value in temp_dict.items()}
                for res in com_pdb_str.residues:
                    res_notation = f'{res.chain}:{res.name}:{res.number}'
                    if res.insertion_code:
                        res_notation += res.insertion_code

                    res_energy = res_dict.get(res_notation, 0.00)
                    for at in res.atoms:
                        at.bfactor = res_energy
                com_pdb_str.save(output_path.as_posix(), 'pdb', True, renumber=False)
                energy2pdb_pml(res_dict, options, bfactor_pml, output_path)
        qpd.setValue(2)

        return bfactor_pml

    def result_table(self, state):
        from GMXMMPBSA.analyzer.plots import Tables
        self.app.treeWidget.clearSelection()
        table_data = self.app.systems[self.system_index]['items_summary'][self.keys_path]
        options = {'table_name': f'Summary | {self.subtitle}'}
        if state:
            self.setSelected(True)
            self.result_table_subw = Tables(table_data, self.result_table_action, options, True)
            self.app.mdi.addSubWindow(self.result_table_subw)
            self.result_table_subw.show()
        elif self.result_table_subw:
            self.app.mdi.activatePreviousSubWindow()
            self.result_table_subw.close()

    def setup_buttons(self):
        if not self.buttons:
            return
        for b in self.charts_action:
            if b not in self.buttons:
                self.tb.addWidget(SpacerItem())
            else:
                self.tb.addWidget(self.charts_action[b]())

        if -1 in self.buttons:
            self._define_option_button()
        elif -2 in self.buttons:
            self._define_result_table_btn()
        elif len(self.buttons) > 1:
            self.tb.addWidget(self.mark_all)

        return self.tb
