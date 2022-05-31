try:
    from PyQt6.QtWidgets import *
    from PyQt6.QtCore import *
    from PyQt6.QtGui import *
except:
    from PyQt5.QtWidgets import *
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *


def editor(parent):
    editor = QDoubleSpinBox(parent)
    editor.setFrame(False)
    editor.setMinimum(0)
    editor.setSingleStep(0.001)
    editor.setMaximum(100000)
    editor.setAccelerated(True)
    editor.setDecimals(3)
    return editor

class KiTableDelegate(QStyledItemDelegate):
    def __init__(self, parent):
        super(KiTableDelegate, self).__init__(parent)

    def createEditor(self, parent, option: QStyleOptionViewItem, index: QModelIndex):
        return editor(parent)

    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.ItemDataRole.EditRole)
        editor.setValue(float(value))

    def setModelData(self, editor, model: QAbstractItemModel, index: QModelIndex):
        editor.interpretText()
        value = editor.value()

        model.setData(index, value, Qt.ItemDataRole.EditRole)

    def updateEditorGeometry(self, editor, option: QStyleOptionViewItem, index:QModelIndex):
        editor.setGeometry(option.rect)


class KiTreeDelegate(QStyledItemDelegate):
    def __init__(self, parent):
        super(KiTreeDelegate, self).__init__(parent)

    def createEditor(self, parent, option: QStyleOptionViewItem, index: QModelIndex):
        if index.column() != 3:
            return
        return editor(parent)

    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.ItemDataRole.EditRole)
        editor.setValue(float(value))

    def setModelData(self, editor, model: QAbstractItemModel, index: QModelIndex):
        editor.interpretText()
        value = editor.value()
        model.setData(index, value, Qt.ItemDataRole.EditRole)

    def updateEditorGeometry(self, editor, option: QStyleOptionViewItem, index:QModelIndex):
        editor.setGeometry(option.rect)
