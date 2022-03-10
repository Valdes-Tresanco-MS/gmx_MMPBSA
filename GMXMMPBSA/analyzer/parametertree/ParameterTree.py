# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2020  pyqtgraph                                               #
#                                                                              #
#   This module is based on pyqtgraph                                          #
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
except:
    from PyQt5.QtWidgets import *
    from PyQt5.QtCore import *

from .ParameterItem import ParameterItem


class ParameterTree(QTreeWidget):
    """Widget used to display or control data from a hierarchy of Parameters"""
    
    def __init__(self, parent=None, showHeader=True):
        """
        ============== ========================================================
        **Arguments:**
        parent         (QWidget) An optional parent widget
        showHeader     (bool) If True, then the QTreeView header is displayed.
        ============== ========================================================
        """
        super(ParameterTree, self).__init__(parent)
        self.setVerticalScrollMode(QTreeWidget.ScrollMode.ScrollPerPixel)
        self.setAnimated(True)
        self.setColumnCount(2)
        self.setHeaderLabels(["Parameter", "Value"])
        self.setAlternatingRowColors(True)
        self.paramSet = None
        self.header().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.setHeaderHidden(not showHeader)
        self.lastSel = None
        self.setRootIsDecorated(False)
        self.itemClicked.connect(self.alt_expanded)

    def alt_expanded(self, item:QTreeWidgetItem, col):
        topitem =self.topLevelItem(0)
        nchild = topitem.childCount()
        childs = [topitem.child(x) for x in range(nchild)]
        if item in childs and item.isExpanded() or item not in childs:
            return
        for child in childs:
            child.setExpanded(False)
        item.setExpanded(True)
        self.scrollToItem(item)

    def setParameters(self, param, showTop=True):
        """
        Set the top-level :class:`Parameter <pyqtgraph.parametertree.Parameter>`
        to be displayed in this ParameterTree.

        If *showTop* is False, then the top-level parameter is hidden and only 
        its children will be visible. This is a convenience method equivalent 
        to::
        
            tree.clear()
            tree.addParameters(param, showTop)
        """
        self.clear()
        self.addParameters(param, showTop=showTop)
        
    def addParameters(self, param, root=None, depth=0, showTop=True):
        """
        Adds one top-level :class:`Parameter <pyqtgraph.parametertree.Parameter>`
        to the view. 
        
        ============== ==========================================================
        **Arguments:** 
        param          The :class:`Parameter <pyqtgraph.parametertree.Parameter>` 
                       to add.
        root           The item within the tree to which *param* should be added.
                       By default, *param* is added as a top-level item.
        showTop        If False, then *param* will be hidden, and only its 
                       children will be visible in the tree.
        ============== ==========================================================
        """
        item = param.makeTreeItem(depth=depth)
        if root is None:
            root = self.invisibleRootItem()
            ## Hide top-level item
            if not showTop:
                item.setText(0, '')
                item.setSizeHint(0, QSize(1,1))
                item.setSizeHint(1, QSize(1,1))
                depth -= 1
        root.addChild(item)
        item.treeWidgetChanged()
            
        for ch in param:
            self.addParameters(ch, root=item, depth=depth+1)

    def selectionChanged(self, *args):
        sel = self.selectedItems()
        if len(sel) != 1:
            sel = None
        if self.lastSel is not None and isinstance(self.lastSel, ParameterItem):
            self.lastSel.selected(False)
        if sel is None:
            self.lastSel = None
            return
        self.lastSel = sel[0]
        if hasattr(sel[0], 'selected'):
            sel[0].selected(True)
        return super().selectionChanged(*args)
