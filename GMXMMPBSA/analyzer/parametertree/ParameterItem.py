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
except:
    from PyQt5.QtWidgets import *


class ParameterItem(QTreeWidgetItem):
    """
    Abstract ParameterTree item. 
    Used to represent the state of a Parameter from within a ParameterTree.
    
    - Sets first column of item to name
    - generates context menu if item is renamable or removable
    - handles child added / removed events
    - provides virtual functions for handling changes from parameter
    
    For more ParameterItem types, see ParameterTree.parameterTypes module.
    """
    
    def __init__(self, param, depth=0):
        super(ParameterItem, self).__init__([param.name(), ''])

        self.param = param
        self.param.registerItem(self)  ## let parameter know this item is connected to it (for debugging)
        self.depth = depth
        
        param.sigValueChanged.connect(self.valueChanged)
        param.sigChildAdded.connect(self.childAdded)
        param.sigChildRemoved.connect(self.childRemoved)
        param.sigLimitsChanged.connect(self.limitsChanged)
        param.sigDefaultChanged.connect(self.defaultChanged)
        param.sigOptionsChanged.connect(self.optsChanged)
        param.sigParentChanged.connect(self.parentChanged)

        ## flag used internally during name editing
        self.ignoreNameColumnChange = False

    def valueChanged(self, param, val):
        ## called when the parameter's value has changed
        pass
    
    def isFocusable(self):
        """Return True if this item should be included in the tab-focus order"""
        return False
        
    def setFocus(self):
        """Give input focus to this item.
        Can be reimplemented to display editor widgets, etc.
        """
        pass

    def treeWidgetChanged(self):
        """Called when this item is added or removed from a tree.
        Expansion, visibility, and column widgets must all be configured AFTER 
        the item is added to a tree, not during __init__.
        """
        self.setHidden(not self.param.opts.get('visible', True))
        self.setExpanded(self.param.opts.get('expanded', True))
        
    def childAdded(self, param, child, pos):
        item = child.makeTreeItem(depth=self.depth+1)
        self.insertChild(pos, item)
        item.treeWidgetChanged()
        
        for i, ch in enumerate(child):
            item.childAdded(child, ch, i)
        
    def childRemoved(self, param, child):
        for i in range(self.childCount()):
            item = self.child(i)
            if item.param is child:
                self.takeChild(i)
                break
                
    def parentChanged(self, param, parent):
        ## called when the parameter's parent has changed.
        pass

    def limitsChanged(self, param, limits):
        """Called when the parameter's limits have changed"""
        pass
    
    def defaultChanged(self, param, default):
        """Called when the parameter's default value has changed"""
        pass

    def optsChanged(self, param, opts):
        """Called when any options are changed that are not
        name, value, default, or limits"""
        if 'visible' in opts:
            self.setHidden(not opts['visible'])

        if 'expanded' in opts and self.isExpanded() != opts['expanded']:
            self.setExpanded(opts['expanded'])

    def contextMenuTriggered(self, name):
        def trigger():
            self.param.contextMenu(name)
        return trigger

    def editName(self):
        self.treeWidget().editItem(self, 0)
        
    def selected(self, sel):
        """Called when this item has been selected (sel=True) OR deselected (sel=False)"""
        pass

    ## for python 3 support, we need to redefine hash and eq methods.
    def __hash__(self):
        return id(self)

    def __eq__(self, x):
        return x is self
