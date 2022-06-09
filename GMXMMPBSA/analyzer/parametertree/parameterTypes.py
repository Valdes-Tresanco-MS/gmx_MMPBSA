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
    from PyQt6.QtCore import *
    from PyQt6.QtWidgets import *
    from PyQt6.QtGui import *
except:
    from PyQt5.QtCore import *
    from PyQt5.QtWidgets import *
    from PyQt5.QtGui import *
import os
from .Parameter import Parameter, registerParameterType
from .ParameterItem import ParameterItem
from collections import OrderedDict
from GMXMMPBSA.analyzer.style import restore_icon


class ColorButton(QPushButton):
    """
    **Bases:** QtGui.QPushButton

    Button displaying a color and allowing the user to select a new color.

    ====================== ============================================================
    **Signals:**
    sigColorChanging(self) emitted whenever a new color is picked in the color dialog
    sigColorChanged(self)  emitted when the selected color is accepted (user clicks OK)
    ====================== ============================================================
    """
    sigColorChanging = pyqtSignal(object)  ## emitted whenever a new color is picked in the color dialog
    sigColorChanged = pyqtSignal(object)  ## emitted when the selected color is accepted (user clicks OK)

    def __init__(self, parent=None, color=(128, 128, 128)):
        super(ColorButton, self).__init__(parent)
        self.setColor(color)
        self.colorDialog = QColorDialog()
        self.colorDialog.setOption(QColorDialog.ColorDialogOption.DontUseNativeDialog, True)
        self.colorDialog.currentColorChanged.connect(self.dialogColorChanged)
        self.colorDialog.rejected.connect(self.colorRejected)
        self.colorDialog.colorSelected.connect(self.colorSelected)
        self.clicked.connect(self.selectColor)
        self.setMinimumHeight(15)
        self.setMinimumWidth(15)

    def paintEvent(self, ev):
        p = QPainter(self)
        rect = self.rect().adjusted(0, 5, -6, -5)
        # draw white base, then texture for indicating transparency, then actual color
        # p.setBrush(QColor('white'))
        # p.drawRect(rect)
        # p.setBrush(QBrush(Qt.BrushStyle.DiagCrossPattern))
        # p.drawRect(rect)
        p.setBrush(self._color)
        p.drawRect(rect)
        p.end()

    def setColor(self, color, finished=True):
        """Sets the button's color and emits both sigColorChanged and sigColorChanging."""
        if isinstance(color, (tuple, list)):
            self._color = QColor(*color)
        elif isinstance(color, QColor):
            self._color = color
        self.update()
        if finished:
            self.sigColorChanged.emit(self)
        else:
            self.sigColorChanging.emit(self)

    def selectColor(self):
        self.origColor = self.color()
        self.colorDialog.setCurrentColor(QColor(*self.color()))
        self.colorDialog.open()

    def dialogColorChanged(self, color):
        if color.isValid():
            self.setColor(color, finished=False)

    def colorRejected(self):
        self.setColor(self.origColor, finished=False)

    def colorSelected(self, color):
        self.setColor(self._color, finished=True)

    def saveState(self):
        return self._color.getRgb()

    def restoreState(self, state):
        self.setColor(state)

    def color(self):
        # if mode == 'qcolor':
        #     return self._color.getRgb()
        return self._color.getRgb()
        # elif mode == 'byte':
        #     return self._color.getRgb()
        # elif mode == 'float':
        #     return self._color.getRgbF()

    def widgetGroupInterface(self):
        return (self.sigColorChanged, ColorButton.saveState, ColorButton.restoreState)


class WidgetParameterItem(ParameterItem):
    """
    ParameterTree item with:
    
    * label in second column for displaying value
    * simple widget for editing value (displayed instead of label when item is selected)
    * button that resets value to default
    
    ==========================  =============================================================
    **Registered Types:**
    int                         Displays a :class:`SpinBox <pyqtgraph.SpinBox>` in integer
                                mode.
    float                       Displays a :class:`SpinBox <pyqtgraph.SpinBox>`.
    bool                        Displays a QCheckBox
    str                         Displays a QLineEdit
    color                       Displays a :class:`ColorButton <pyqtgraph.ColorButton>`
    colormap                    Displays a :class:`GradientWidget <pyqtgraph.GradientWidget>`
    ==========================  =============================================================
    
    This class can be subclassed by overriding makeWidget() to provide a custom widget.
    """

    def __init__(self, param, depth):
        super(WidgetParameterItem, self).__init__(param, depth)

        self.asSubItem = False  # place in a child item's column 0 instead of column 1
        self.hideWidget = True  # hide edit widget, replace with label when not selected
        # set this to False to keep the editor widget always visible

        # build widget with a display label and default button
        self.widget = self.makeWidget()

        if self.asSubItem:
            self.subItem = QTreeWidgetItem()
            self.subItem.depth = self.depth + 1
            self.subItem.setFlags(Qt.ItemFlag.NoItemFlags)
            self.addChild(self.subItem)

        self.defaultBtn = QPushButton()
        self.defaultBtn.setAutoDefault(False)
        self.defaultBtn.setFixedWidth(20)
        self.defaultBtn.setFixedHeight(20)
        modDir = os.path.dirname(__file__)
        self.defaultBtn.setIcon(QIcon(restore_icon))
        self.defaultBtn.clicked.connect(self.defaultClicked)

        self.displayLabel = QLabel()

        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        if not self.asSubItem:
            layout.addWidget(self.widget, 1)
        layout.addWidget(self.displayLabel, 1)
        layout.addStretch(0)
        layout.addWidget(self.defaultBtn)
        self.layoutWidget = QWidget()
        self.layoutWidget.setLayout(layout)

        if self.widget.sigChanged is not None:
            self.widget.sigChanged.connect(self.widgetValueChanged)

        if hasattr(self.widget, 'sigChanging'):
            self.widget.sigChanging.connect(self.widgetValueChanging)

        ## update value shown in widget. 
        # if self.param.opts.get('value', None) is not None:
        self.valueChanged(self, self.param.opts['value'], force=True)
        # else:
            ## no starting value was given; use whatever the widget has
            # self.widgetValueChanged()

        self.updateDefaultBtn()

        self.optsChanged(self.param, self.param.opts)

    def makeWidget(self):
        """
        Return a single widget whose position in the tree is determined by the
        value of self.asSubItem. If True, it will be placed in the second tree
        column, and if False, the first tree column of a child item.

        The widget must be given three attributes:
        
        ==========  ============================================================
        sigChanged  a signal that is emitted when the widget's value is changed
        value       a function that returns the value
        setValue    a function that sets the value
        ==========  ============================================================
            
        This is a good function to override in subclasses.
        """
        opts = self.param.opts
        t = opts['type']
        if t == 'int':
            defs = {
                'value': 0,
                'min': None,
                'max': None,
                'step': 1.0,
                'dec': False,
                'siPrefix': False,
                'suffix': '',
                'decimals': 3,
                'int': True,
                'minStep': 1.0,
            }

            for k in defs:
                if k in opts:
                    defs[k] = opts[k]
            # if 'limits' in opts:
            #     defs['min'], defs['max'] = opts['limits']
            w = QSpinBox()
            if 'limits' in opts:
                w.setRange(opts['limits'][0], opts['limits'][1])
            w.setSingleStep(opts.get('step', 1))
            w.setAccelerated(opts.get('accelerated', False))
            w.setSuffix(opts.get('suffix', ''))
            w.setPrefix(opts.get('prefix', ''))
            w.sigChanged = w.valueChanged
        elif t == 'float':
            w = QDoubleSpinBox()
            if 'limits' in opts:
                w.setRange(opts['limits'][0], opts['limits'][1])
            if 'step' in opts:
                w.setSingleStep(opts['step'])
            if 'accelerated' in opts:
                w.setAccelerated(opts['accelerated'])
            if 'suffix' in opts:
                w.setSuffix(opts['suffix'])
            if 'prefix' in opts:
                w.setPrefix(opts['prefix'])

            w.sigChanged = w.valueChanged
        elif t == 'bool':
            w = QCheckBox()
            w.sigChanged = w.toggled
            w.value = w.isChecked
            w.setValue = w.setChecked
            self.hideWidget = False
        elif t == 'str':
            w = QLineEdit()
            w.setStyleSheet('border: 0px')
            w.sigChanged = w.editingFinished
            w.value = lambda: w.text()
            w.setValue = lambda v: w.setText(v)
            w.sigChanging = w.textChanged
        elif t == 'color':
            w = ColorButton()
            w.sigChanged = w.sigColorChanged
            w.sigChanging = w.sigColorChanging
            w.value = w.color
            w.setValue = w.setColor
            self.hideWidget = False
            w.setFlat(True)
        elif t == 'colormap':
            w = None
        else:
            raise Exception("Unknown type '%s'" % t)
        return w

    def setFocus(self):
        self.showEditor()

    def isFocusable(self):
        return self.param.opts['visible'] and self.param.opts['enabled']

    def valueChanged(self, param, val, force=False):
        ## called when the parameter's value has changed
        if force or val != self.widget.value():
            try:
                self.widget.sigChanged.disconnect(self.widgetValueChanged)
                self.param.sigValueChanged.disconnect(self.valueChanged)
                self.widget.setValue(val)
                self.param.setValue(self.widget.value())
            finally:
                self.widget.sigChanged.connect(self.widgetValueChanged)
                self.param.sigValueChanged.connect(self.valueChanged)
        self.updateDisplayLabel()  ## always make sure label is updated, even if values match!
        self.updateDefaultBtn()

    def updateDefaultBtn(self):
        ## enable/disable default btn 
        self.defaultBtn.setEnabled(not self.param.valueIsDefault())
        # hide / show
        self.defaultBtn.setVisible(self.param.hasDefault())

    def updateDisplayLabel(self, value=None):
        """Update the display label to reflect the value of the parameter."""
        if value is None:
            value = self.param.value()
        if isinstance(self.widget, (QSpinBox, QDoubleSpinBox)):
            text = f"{self.widget.prefix()}{self.widget.value()}{self.widget.suffix()}"
        elif isinstance(self.widget, QComboBox):
            text = self.widget.currentText()
        else:
            text = value
        self.displayLabel.setText(f"{text}")

    def widgetValueChanged(self):
        ## called when the widget's value has been changed by the user
        val = self.widget.value()
        newVal = self.param.setValue(val)

    def widgetValueChanging(self, *args):
        """
        Called when the widget's value is changing, but not finalized.
        For example: editing text before pressing enter or changing focus.
        """
        self.param.sigValueChanging.emit(self.param, self.widget.value())

    def selected(self, sel):
        """Called when this item has been selected (sel=True) OR deselected (sel=False)"""
        # ParameterItem.selected(self, sel)

        if self.widget is None:
            return
        if sel:
            self.showEditor()
        elif self.hideWidget:
            self.hideEditor()

    def showEditor(self):
        self.widget.show()
        self.displayLabel.hide()
        self.widget.setFocus(Qt.FocusReason.OtherFocusReason)

    def hideEditor(self):
        self.widget.hide()
        self.displayLabel.show()

    def limitsChanged(self, param, limits):
        """Called when the parameter's limits have changed"""
        t = self.param.opts['type']
        if t in ['int', 'float']:
            self.param.setOpts(bounds=limits)
        else:
            return  ## don't know what to do with any other types..

    def defaultChanged(self, param, value):
        self.updateDefaultBtn()

    def treeWidgetChanged(self):
        """Called when this item is added or removed from a tree."""
        ## add all widgets for this item into the tree
        if self.widget is not None:
            tree = self.treeWidget()
            if tree is None:
                return
            if self.asSubItem:
                self.subItem.setFirstColumnSpanned(True)
                tree.setItemWidget(self.subItem, 0, self.widget)
            tree.setItemWidget(self, 1, self.layoutWidget)
            self.displayLabel.hide()
            self.selected(False)

    def defaultClicked(self):
        self.param.setToDefault()

    def optsChanged(self, param, opts):
        """Called when any options are changed that are not
        name, value, default, or limits"""

        if 'enabled' in opts:
            self.updateDefaultBtn()
            self.widget.setEnabled(opts['enabled'])

        if 'readonly' in opts:
            self.updateDefaultBtn()
            if hasattr(self.widget, 'setReadOnly'):
                self.widget.setReadOnly(opts['readonly'])
            else:
                self.widget.setEnabled(self.param.opts['enabled'] and not opts['readonly'])

        if 'tip' in opts:
            self.setToolTip(0, opts['tip'])
            self.setToolTip(1, opts['tip'])


class SimpleParameter(Parameter):
    """Parameter representing a single value.

    This parameter is backed by :class:`WidgetParameterItem` to represent the
    following parameter names:

    - 'int'
    - 'float'
    - 'bool'
    - 'str'
    - 'color'
    - 'colormap'
    """
    itemClass = WidgetParameterItem

    def __init__(self, *args, **kargs):
        """Initialize the parameter.

        This is normally called implicitly through :meth:`Parameter.create`.
        The keyword arguments available to :meth:`Parameter.__init__` are
        applicable.
        """
        super(SimpleParameter, self).__init__(**kargs)

        ## override a few methods for color parameters
        if self.opts['type'] == 'color':
            self.value = self.colorValue
            self.saveState = self.saveColorState

    def colorValue(self):
        return Parameter.value(self)

    def saveColorState(self, *args, **kwds):
        state = Parameter.saveState(self, *args, **kwds)
        state['value'] = self.value()
        return state

    def _interpretValue(self, v):
        fn = {
            'int': int,
            'float': float,
            'bool': bool,
            'str': str,
            'color': self._interpColor,
        }[self.opts['type']]
        return fn(v)

    def _interpColor(self, v):
        if isinstance(v, (tuple, list)):
            return list(v)
        elif isinstance(v, QColor):
            return list(v.getRgb())
        elif isinstance(v, str):
            return list(QColor(v).getRgb())

registerParameterType('int', SimpleParameter, override=True)
registerParameterType('float', SimpleParameter, override=True)
registerParameterType('bool', SimpleParameter, override=True)
registerParameterType('str', SimpleParameter, override=True)
registerParameterType('color', SimpleParameter, override=True)


# registerParameterType('colormap', SimpleParameter, override=True)


class GroupParameterItem(ParameterItem):
    """
    Group parameters are used mainly as a generic parent item that holds (and groups!) a set
    of child parameters. It also provides a simple mechanism for displaying a button or combo
    that can be used to add new parameters to the group.
    """

    def __init__(self, param, depth):
        super(GroupParameterItem, self).__init__(param, depth)
        self.updateDepth(depth)

        self.addItem = None
        if 'addText' in param.opts:
            addText = param.opts['addText']
            if 'addList' in param.opts:
                self.addWidget = QComboBox()
                self.addWidget.setSizeAdjustPolicy(QComboBox.AdjustToContents)
                self.updateAddList()
                self.addWidget.currentIndexChanged.connect(self.addChanged)
            else:
                self.addWidget = QPushButton(addText)
                self.addWidget.clicked.connect(self.addClicked)
            w = QWidget()
            l = QHBoxLayout()
            l.setContentsMargins(0, 0, 0, 0)
            w.setLayout(l)
            l.addWidget(self.addWidget)
            l.addStretch()
            self.addWidgetBox = w
            self.addItem = QTreeWidgetItem([])
            self.addItem.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.addItem.depth = self.depth + 1
            ParameterItem.addChild(self, self.addItem)
            self.addItem.setSizeHint(0, self.addWidgetBox.sizeHint())

        self.optsChanged(self.param, self.param.opts)

    def updateDepth(self, depth):
        ## Change item's appearance based on its depth in the tree
        ## This allows highest-level groups to be displayed more prominently.
        for c in [0, 1]:
            if depth == 0:
                self.setBackground(c, QBrush(QColor(100, 100, 100)))
                self.setForeground(c, QBrush(QColor(220, 220, 255)))
                font = self.font(c)
                font.setWeight(QFont.Weight.DemiBold)
                font.setPointSize(font.pointSize() + 1)
            else:
                self.setBackground(c, QBrush(QColor(220, 220, 220)))
                self.setForeground(c, QBrush(QColor(50, 50, 50)))
                font = self.font(c)
                font.setWeight(QFont.Weight.DemiBold)
            self.setFont(c, font)

    def addClicked(self):
        """Called when "add new" button is clicked
        The parameter MUST have an 'addNew' method defined.
        """
        self.param.addNew()

    def addChanged(self):
        """Called when "add new" combo is changed
        The parameter MUST have an 'addNew' method defined.
        """
        if self.addWidget.currentIndex() == 0:
            return
        typ = self.addWidget.currentText()
        self.param.addNew(typ)
        self.addWidget.setCurrentIndex(0)

    def treeWidgetChanged(self):
        ParameterItem.treeWidgetChanged(self)
        tw = self.treeWidget()
        if tw is None:
            return
        self.setFirstColumnSpanned(True)
        if self.addItem is not None:
            tw.setItemWidget(self.addItem, 0, self.addWidgetBox)
            self.addItem.setFirstColumnSpanned(True)

    def addChild(self, child):  ## make sure added childs are actually inserted before add btn
        if self.addItem is not None:
            ParameterItem.insertChild(self, self.childCount() - 1, child)
        else:
            ParameterItem.addChild(self, child)

    def optsChanged(self, param, opts):
        ParameterItem.optsChanged(self, param, opts)

        if 'addList' in opts:
            self.updateAddList()

        if hasattr(self, 'addWidget') and 'enabled' in opts:
            self.addWidget.setEnabled(opts['enabled'])

    def updateAddList(self):
        self.addWidget.blockSignals(True)
        try:
            self.addWidget.clear()
            self.addWidget.addItem(self.param.opts['addText'])
            for t in self.param.opts['addList']:
                self.addWidget.addItem(t)
        finally:
            self.addWidget.blockSignals(False)


class GroupParameter(Parameter):
    """
    Group parameters are used mainly as a generic parent item that holds (and groups!) a set
    of child parameters. 
    
    It also provides a simple mechanism for displaying a button or combo
    that can be used to add new parameters to the group. To enable this, the group 
    must be initialized with the 'addText' option (the text will be displayed on
    a button which, when clicked, will cause addNew() to be called). If the 'addList'
    option is specified as well, then a dropdown-list of addable items will be displayed
    instead of a button.
    """
    itemClass = GroupParameterItem
    sigAddNew = pyqtSignal(object, object)  # self, type

    def addNew(self, typ=None):
        """
        This method is called when the user has requested to add a new item to the group.
        By default, it emits ``sigAddNew(self, typ)``.
        """
        self.sigAddNew.emit(self, typ)

    def setAddList(self, vals):
        """Change the list of options available for the user to add to the group."""
        self.setOpts(addList=vals)


registerParameterType('group', GroupParameter, override=True)


class ListParameterItem(WidgetParameterItem):
    """
    WidgetParameterItem subclass providing comboBox that lets the user select from a list of options.
    """
    def __init__(self, param, depth):
        self.targetValue = None
        super(ListParameterItem, self).__init__(param, depth)

    def makeWidget(self):
        opts = self.param.opts
        w = QComboBox()
        w.setMaximumHeight(20)  ## set to match height of spin box and line edit
        w.sigChanged = w.currentIndexChanged
        w.value = self.value
        w.setValue = self.setValue
        self.widget = w  ## needs to be set before limits are changed
        self.limitsChanged(self.param, self.param.opts['limits'])
        if len(self.forward) > 0:
            self.setValue(self.param.value())
        return w

    def value(self):
        key = self.widget.currentText()
        return self.forward.get(key, None)

    def setValue(self, val):
        self.targetValue = val
        if val not in self.reverse[0]:
            self.widget.setCurrentIndex(0)
        else:
            key = self.reverse[1][self.reverse[0].index(val)]
            ind = self.widget.findText(key)
            self.widget.setCurrentIndex(ind)

    def limitsChanged(self, param, limits):
        # set up forward / reverse mappings for name:value
        if len(limits) == 0:
            limits = ['']  ## Can never have an empty list--there is always at least a singhe blank item.
        self.forward, self.reverse = ListParameter.mapping(limits)
        try:
            self.widget.blockSignals(True)
            val = self.targetValue  # asUnicode(self.widget.currentText())

            self.widget.clear()
            for c, k in enumerate(self.forward):
                if str(k).startswith('---'):
                    self.widget.insertSeparator(c)
                else:
                    self.widget.addItem(k)
                if k == val:
                    self.widget.setCurrentIndex(self.widget.count() - 1)
                    self.updateDisplayLabel()
        finally:
            self.widget.blockSignals(False)


class ListParameter(Parameter):
    """Parameter with a list of acceptable values.

    By default, this parameter is represented by a :class:`ListParameterItem`,
    displaying a combo box to select a value from the list.

    In addition to the generic :class:`~pyqtgraph.parametertree.Parameter`
    options, this parameter type accepts a ``limits`` argument specifying the
    list of allowed values.  ``values`` is an alias and may be used instead.

    The values may generally be of any data type, as long as they can be
    represented as a string. If the string representation provided is
    undesirable, the values may be given as a dictionary mapping the desired
    string representation to the value.
    """
    itemClass = ListParameterItem
    def __init__(self, **opts):
        self.forward = OrderedDict()  ## {name: value, ...}
        self.reverse = ([], [])  ## ([value, ...], [name, ...])
        # Parameter uses 'limits' option to define the set of allowed values
        if 'values' in opts:
            opts['limits'] = opts['values']
        if opts.get('limits', None) is None:
            opts['limits'] = []
        super(ListParameter, self).__init__(**opts)
        self.setLimits(opts['limits'])

    def setLimits(self, limits):
        """Change the list of allowed values."""
        self.forward, self.reverse = self.mapping(limits)
        if len(self.reverse[0]) > 0 and self.value() not in self.reverse[0]:
            self.setValue(self.reverse[0][0])

    @staticmethod
    def mapping(limits):
        # Return forward and reverse mapping objects given a limit specification
        forward = OrderedDict()  ## {name: value, ...}
        reverse = ([], [])  ## ([value, ...], [name, ...])
        if isinstance(limits, dict):
            for k, v in limits.items():
                forward[k] = v
                reverse[0].append(v)
                reverse[1].append(k)
        else:
            for v in limits:
                n = v
                forward[n] = v
                reverse[0].append(v)
                reverse[1].append(n)
        return forward, reverse


registerParameterType('list', ListParameter, override=True)


class ActionParameterItem(ParameterItem):
    """ParameterItem displaying a clickable button."""

    def __init__(self, param, depth):
        super(ActionParameterItem, self).__init__(param, depth)
        # ParameterItem.__init__(self, param, depth)
        self.layoutWidget = QWidget()
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layoutWidget.setLayout(self.layout)
        self.button = QPushButton()
        self.layout.addWidget(self.button)
        self.layout.addStretch()
        self.button.clicked.connect(self.buttonClicked)
        self.optsChanged(self.param, self.param.opts)

    def treeWidgetChanged(self):
        tree = self.treeWidget()
        if tree is None:
            return
        self.setFirstColumnSpanned(True)
        tree.setItemWidget(self, 0, self.layoutWidget)

    def optsChanged(self, param, opts):
        # ParameterItem.optsChanged(self, param, opts)
        if 'enabled' in opts:
            self.button.setEnabled(opts['enabled'])

    def buttonClicked(self):
        self.param.activate()


class ActionParameter(Parameter):
    """Used for displaying a button within the tree.
    ``sigActivated(self)`` is emitted when the button is clicked.
    """
    itemClass = ActionParameterItem
    sigActivated = pyqtSignal(object)

    def activate(self):
        self.sigActivated.emit(self)
        self.emitStateChanged('activated', None)


registerParameterType('action', ActionParameter, override=True)


class TextParameterItem(WidgetParameterItem):
    """ParameterItem displaying a QTextEdit widget."""

    def makeWidget(self):
        self.hideWidget = False
        self.asSubItem = True
        self.textBox = w = QTextEdit()
        w.sizeHint = lambda: QSize(300, 100)
        w.value = lambda: str(w.toPlainText())
        w.setValue = w.setPlainText
        w.sigChanged = w.textChanged
        return w


class TextParameter(Parameter):
    """Editable string, displayed as large text box in the tree."""
    itemClass = TextParameterItem


registerParameterType('text', TextParameter, override=True)
