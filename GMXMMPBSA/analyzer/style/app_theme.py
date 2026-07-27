# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdés-Tresanco and Mario E. Valdés-Tresanco   #
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
    from PyQt6.QtCore import QSize, Qt
    from PyQt6.QtGui import QColor, QPalette
    from PyQt6.QtWidgets import QAbstractItemView, QApplication, QHeaderView
except ImportError:
    from PyQt5.QtCore import QSize, Qt
    from PyQt5.QtGui import QColor, QPalette
    from PyQt5.QtWidgets import QAbstractItemView, QApplication, QHeaderView

from pathlib import Path


STYLE_DIR = Path(__file__).parent
ARROW_DOWN = STYLE_DIR.joinpath('arrow_down.svg').absolute().as_posix()
ARROW_UP = STYLE_DIR.joinpath('arrow_up.svg').absolute().as_posix()


COLORS = {
    'window': '#e9eef4',
    'surface': '#fbfcfe',
    'surface_alt': '#f0f4f8',
    'border': '#cbd5e1',
    'border_strong': '#aeb9c8',
    'text': '#202733',
    'muted': '#526173',
    'disabled': '#9aa4b2',
    'primary': '#2f6fbd',
    'primary_hover': '#245a9f',
    'primary_soft': '#d9e8fb',
    'primary_selected': '#aecdf1',
    'header': '#dfe7f0',
    'menu_hover': '#e3edf9',
}


def _palette_role(name):
    color_role = getattr(QPalette, 'ColorRole', QPalette)
    return getattr(color_role, name)


APP_QSS = f"""
QMainWindow, QDialog {{
    background: {COLORS['window']};
    color: {COLORS['text']};
}}

QMdiArea {{
    background: #a7abae;
    border: 1px solid {COLORS['border']};
}}

QMdiSubWindow {{
    background: {COLORS['surface']};
}}

QMenuBar {{
    background: {COLORS['header']};
    color: {COLORS['text']};
    border-bottom: 1px solid {COLORS['border']};
    spacing: 2px;
}}

QMenuBar::item {{
    background: transparent;
    padding: 6px 10px;
}}

QMenuBar::item:selected {{
    background: {COLORS['menu_hover']};
    color: {COLORS['text']};
}}

QMenu {{
    background: {COLORS['surface']};
    color: {COLORS['text']};
    border: 1px solid {COLORS['border_strong']};
    padding: 4px 0;
}}

QMenu::item {{
    padding: 6px 28px 6px 24px;
}}

QMenu::item:selected {{
    background: {COLORS['menu_hover']};
    color: {COLORS['text']};
}}

QStatusBar {{
    background: {COLORS['header']};
    color: {COLORS['muted']};
    border-top: 1px solid {COLORS['border']};
}}

QDockWidget {{
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
    color: {COLORS['text']};
}}

QDockWidget::title {{
    background: {COLORS['header']};
    border: 1px solid {COLORS['border']};
    border-bottom: 0;
    padding: 6px 8px;
    text-align: left;
}}

QTabWidget::pane {{
    background: {COLORS['surface_alt']};
    border: 1px solid {COLORS['border']};
}}

QTabBar::tab {{
    background: {COLORS['header']};
    color: {COLORS['muted']};
    border: 1px solid {COLORS['border']};
    border-bottom: none;
    padding: 7px 12px;
    margin-right: 2px;
}}

QTabBar::tab:selected {{
    background: {COLORS['surface']};
    color: {COLORS['text']};
    border-top: 2px solid {COLORS['primary']};
}}

QTabBar::tab:hover:!selected {{
    background: {COLORS['surface_alt']};
    color: {COLORS['text']};
}}

QGroupBox {{
    background: {COLORS['surface_alt']};
    border: 1px solid {COLORS['border']};
    border-radius: 6px;
    margin-top: 10px;
    padding: 10px 8px 8px 8px;
    font-weight: 600;
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
    color: {COLORS['muted']};
}}

QWidget#AnalyzerOptionPanel, QWidget#CorrelationPanel {{
    background: {COLORS['window']};
}}

QWidget#ChartCanvasPanel, QWidget#AnalyzerTablePanel {{
    background: {COLORS['surface']};
}}

QTextEdit#AnalyzerTextOutput {{
    background: #f8fafc;
    color: {COLORS['text']};
    border: 1px solid {COLORS['border']};
    selection-background-color: {COLORS['primary_selected']};
    selection-color: {COLORS['text']};
}}

QToolBar#ChartControlBar {{
    background: {COLORS['header']};
    border: 0;
    border-bottom: 1px solid {COLORS['border']};
    spacing: 4px;
    padding: 5px 8px;
}}

QTreeWidget, QTableWidget, QTreeView, QTableView {{
    background: {COLORS['surface']};
    alternate-background-color: {COLORS['surface_alt']};
    color: {COLORS['text']};
    border: 1px solid {COLORS['border']};
    gridline-color: {COLORS['border']};
    selection-background-color: {COLORS['primary_selected']};
    selection-color: {COLORS['text']};
    outline: 0;
}}

QTreeWidget::item, QTableWidget::item {{
    padding: 4px 6px;
}}

QTreeWidget::item:selected, QTableWidget::item:selected {{
    background: {COLORS['primary_selected']};
    color: {COLORS['text']};
}}

QTreeWidget::item:hover, QTableWidget::item:hover {{
    background: {COLORS['menu_hover']};
}}

QHeaderView::section {{
    background: {COLORS['header']};
    color: {COLORS['text']};
    border: 0;
    border-right: 1px solid {COLORS['border']};
    border-bottom: 1px solid {COLORS['border_strong']};
    padding: 6px 8px;
    font-weight: 600;
}}

QPushButton, QToolButton {{
    background: {COLORS['surface']};
    color: {COLORS['text']};
    border: 1px solid {COLORS['border_strong']};
    border-radius: 5px;
    padding: 5px 10px;
}}

QPushButton:hover, QToolButton:hover {{
    background: {COLORS['surface_alt']};
    border-color: {COLORS['primary']};
}}

QPushButton:pressed, QToolButton:pressed {{
    background: {COLORS['primary_soft']};
}}

QToolButton:checked, QPushButton:checked {{
    background: {COLORS['primary_soft']};
    border-color: {COLORS['primary']};
    color: {COLORS['text']};
}}

QPushButton:disabled, QToolButton:disabled {{
    background: {COLORS['header']};
    color: {COLORS['disabled']};
    border-color: {COLORS['border']};
}}

QToolButton#ToolbarSpacer {{
    background: transparent;
    border: 0;
    padding: 0;
}}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background: {COLORS['surface']};
    color: {COLORS['text']};
    border: 1px solid {COLORS['border_strong']};
    border-radius: 5px;
    padding: 4px 6px;
    min-height: 24px;
}}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border-color: {COLORS['primary']};
}}

QLineEdit:read-only {{
    background: {COLORS['surface_alt']};
    color: {COLORS['muted']};
}}

QComboBox, QSpinBox, QDoubleSpinBox {{
    padding-right: 24px;
}}

QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 22px;
    border-left: 1px solid {COLORS['border']};
    border-top-right-radius: 5px;
    border-bottom-right-radius: 5px;
    background: transparent;
}}

QComboBox::down-arrow {{
    image: url({ARROW_DOWN});
    width: 8px;
    height: 8px;
}}

QSpinBox::up-button, QDoubleSpinBox::up-button {{
    subcontrol-origin: border;
    subcontrol-position: top right;
    width: 22px;
    height: 13px;
    border-left: 1px solid {COLORS['border']};
    border-bottom: 0;
    margin: 1px 1px 0 0;
    background: transparent;
}}

QSpinBox::down-button, QDoubleSpinBox::down-button {{
    subcontrol-origin: border;
    subcontrol-position: bottom right;
    width: 22px;
    height: 13px;
    border-left: 1px solid {COLORS['border']};
    border-top: 0;
    margin: 0 1px 1px 0;
    background: transparent;
}}

QSpinBox::up-button:hover, QSpinBox::down-button:hover,
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover,
QComboBox::drop-down:hover {{
    background: {COLORS['menu_hover']};
}}

QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
    image: url({ARROW_UP});
    width: 7px;
    height: 7px;
}}

QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
    image: url({ARROW_DOWN});
    width: 7px;
    height: 7px;
}}

QScrollBar:vertical, QScrollBar:horizontal {{
    background: {COLORS['surface_alt']};
    border: 0;
    margin: 0;
}}

QScrollBar:vertical {{
    width: 12px;
}}

QScrollBar:horizontal {{
    height: 12px;
}}

QScrollBar::handle {{
    background: {COLORS['border_strong']};
    border-radius: 5px;
    min-height: 24px;
    min-width: 24px;
}}

QScrollBar::handle:hover {{
    background: {COLORS['muted']};
}}

QScrollBar::add-line, QScrollBar::sub-line {{
    width: 0;
    height: 0;
}}
"""


def apply_app_theme(widget):
    app = QApplication.instance()
    if app is not None:
        palette = app.palette()
        palette.setColor(_palette_role('Window'), QColor(COLORS['window']))
        palette.setColor(_palette_role('Base'), QColor(COLORS['surface']))
        palette.setColor(_palette_role('AlternateBase'), QColor(COLORS['surface_alt']))
        palette.setColor(_palette_role('Text'), QColor(COLORS['text']))
        palette.setColor(_palette_role('WindowText'), QColor(COLORS['text']))
        palette.setColor(_palette_role('Button'), QColor(COLORS['surface']))
        palette.setColor(_palette_role('ButtonText'), QColor(COLORS['text']))
        palette.setColor(_palette_role('Highlight'), QColor(COLORS['primary_selected']))
        palette.setColor(_palette_role('HighlightedText'), QColor(COLORS['text']))
        app.setPalette(palette)
        app.setStyleSheet(APP_QSS)
    else:
        widget.setStyleSheet(APP_QSS)


def polish_table(table, stretch=True, select_rows=True, single_selection=False, hide_vertical_header=True):
    table.setAlternatingRowColors(True)
    table.setShowGrid(False)
    table.setWordWrap(False)
    table.setIconSize(QSize(18, 18))
    table.verticalHeader().setDefaultSectionSize(28)
    table.verticalHeader().setMinimumSectionSize(24)
    if hide_vertical_header:
        table.verticalHeader().hide()
    if select_rows:
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    if single_selection:
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    header = table.horizontalHeader()
    header.setHighlightSections(False)
    header.setMinimumSectionSize(48)
    header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    if stretch:
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)


def polish_tree(tree):
    tree.setAlternatingRowColors(True)
    tree.setUniformRowHeights(True)
    tree.setIndentation(16)
    tree.setIconSize(QSize(18, 18))
    tree.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    tree.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    tree.header().setHighlightSections(False)
    tree.header().setMinimumSectionSize(60)


def polish_tool_button(button, icon_size=18):
    button.setIconSize(QSize(icon_size, icon_size))
    button.setMinimumSize(QSize(28, 28))
    button.setAutoRaise(False)


def polish_toolbar(toolbar, icon_size=18):
    toolbar.setIconSize(QSize(icon_size, icon_size))
    toolbar.setContentsMargins(0, 0, 0, 0)
