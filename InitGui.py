# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2024 - Harry van Langen <hvlanalysis@icloud.com>        *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

import os
import dummy
import FreeCAD
import FreeCADGui
from PySide2 import QtWidgets, QtGui, QtCore

global FCmw
FCmw = FreeCADGui.getMainWindow()
global dir_name
dir_name = os.path.dirname(dummy.file_path())
global double_validator
double_validator = QtGui.QDoubleValidator()
global int_validator
int_validator = QtGui.QIntValidator()


class fcVMWorkbench(Workbench):
    Icon = os.path.join(dir_name, "icons", "fcFEM.svg")
    MenuText = "fcVM workbench"
    ToolTip = "Plastic collapse analysis with the von Mises material model"

    def GetClassName(self):
        return "Gui::PythonWorkbench"

    def Initialize(self):
        from PySide2 import QtGui
        self.appendToolbar("fcVM", [])
        self.appendMenu("fcVM", [])
        self.palette_warning = QtGui.QPalette()
        self.palette_warning.setColor(QtGui.QPalette.Base, QtGui.QColor("orange"))
        self.palette_standard = QtGui.QPalette()
        self.palette_standard.setColor(QtGui.QPalette.Base, QtGui.QColor("white"))
        # self.palette.setColor(QtGui.QPalette.Text, QtGui.QColor("red"))

    def Activated(self):
        from PySide2 import QtCore
        global w

        self.doc = FreeCAD.activeDocument()
        if self.doc == None:
            self.doc = FreeCAD.newDocument("fcVM")

        self.file_name = self.doc.Label
        print("self.file_name: ", self.file_name)
        self.macro_file_path = dir_name + '/source code/' + 'fcVM.FCMacro'

        print(self.macro_file_path)

        self.disp_option = "incremental"

        class DocObserver(object):  # document Observer
            def __init__(self, workbench_instance):
                self.workbench_instance = workbench_instance

            def slotActivateDocument(self, doc):
                print("----- slotActivateDocument -----")
                if FreeCAD.activeDocument().Label[0:7] != "Unnamed":
                    self.workbench_instance.save_clicked()
                    self.workbench_instance.file_name = FreeCAD.activeDocument().Label
                    self.workbench_instance.open_file()
                print("--------------------------------")
                print()

            def slotFinishSaveDocument(self, doc, prop):
                print("----- slotfinishSaveDocument -----")
                self.workbench_instance.save_clicked()  # save under old file name
                self.workbench_instance.file_name = doc.Label
                self.workbench_instance.save_clicked()  # save under new file name
                print("----------------------------------")
                print()

        self.obs = DocObserver(self)
        FreeCAD.addDocumentObserver(self.obs)

        ui_Path = os.path.join(dir_name, "user_interface", "fcVM.ui")
        w = FreeCADGui.PySideUic.loadUi(ui_Path)

        w.startBtn.clicked.connect(self.start_clicked)
        w.quitBtn.clicked.connect(self.Deactivated)
        w.resetBtn.clicked.connect(self.reset_clicked)
        w.saveBtn.clicked.connect(self.save_clicked)
        w.totalRbtn.toggled.connect(self.btn_state)
        w.incrRbtn.toggled.connect(self.btn_state)

        w.max_iter.textChanged.connect(self.max_iter_changed)
        w.relax.textChanged.connect(self.relax_changed)
        w.scale_1.textChanged.connect(self.scale_1_changed)
        w.scale_2.textChanged.connect(self.scale_2_changed)
        w.scale_3.textChanged.connect(self.scale_3_changed)

        w.YSinput.setValidator(double_validator)
        w.GXinput.setValidator(double_validator)
        w.GYinput.setValidator(double_validator)
        w.GZinput.setValidator(double_validator)
        w.steps.setValidator(int_validator)
        w.max_iter.setValidator(int_validator)
        w.error.setValidator(double_validator)
        w.relax.setValidator(double_validator)
        w.scale_1.setValidator(double_validator)
        w.scale_2.setValidator(double_validator)
        w.scale_3.setValidator(double_validator)

        self.YSinput_default = "240.0"
        self.GXinput_default = "0.0"
        self.GYinput_default = "0.0"
        self.GZinput_default = "-10.0"
        self.steps_default = "10"
        self.max_iter_default = "20"
        self.error_default = "1.0e-3"
        self.relax_default = "1.2"
        self.scale_1_default = "2.0"
        self.scale_2_default = "1.2"
        self.scale_3_default = "1.2"
        self.disp_option_default = "incremental"

        self.open_file()

        FCmw.addDockWidget(QtCore.Qt.RightDockWidgetArea, w.dw)

    def Deactivated(self):

        FreeCAD.removeDocumentObserver(self.obs)

        try:
            if w.dw.isVisible():
                w.dw.setVisible(False)
        except Exception:
            None

    def start_clicked(self):
        from PySide2 import QtWidgets

        self.save_clicked()

        try:
            rv = FCmw.findChild(QtWidgets.QTextEdit, "Report view")
            rv.clear()
        except Exception:
            None

        fcVM_macro = open(self.macro_file_path).read()
        exec(fcVM_macro)

    def quit_clicked(self):
        self.Deactivated()

    def reset_clicked(self):
        w.max_iter.setText(self.max_iter_default)
        w.error.setText(self.error_default)
        w.relax.setText(self.relax_default)
        w.scale_1.setText(self.scale_1_default)
        w.scale_2.setText(self.scale_2_default)
        w.scale_3.setText(self.scale_3_default)
        w.relax.setPalette(self.palette_standard)
        w.scale_1.setPalette(self.palette_standard)
        w.scale_2.setPalette(self.palette_standard)
        w.scale_3.setPalette(self.palette_standard)
        w.incrRbtn.setChecked(True)

    def save_clicked(self):
        inp_file_path = dir_name + '/control files/' + self.file_name + '.inp'
        print("save file: ", inp_file_path)
        with open(inp_file_path, "w") as f:
            f.write(w.YSinput.text() + "\n")
            f.write(w.GXinput.text() + "\n")
            f.write(w.GYinput.text() + "\n")
            f.write(w.GZinput.text() + "\n")
            f.write(w.steps.text() + "\n")
            f.write(w.max_iter.text() + "\n")
            f.write(w.error.text() + "\n")
            f.write(w.relax.text() + "\n")
            f.write(w.scale_1.text() + "\n")
            f.write(w.scale_2.text() + "\n")
            f.write(w.scale_3.text() + "\n")
            f.write(self.disp_option + "\n")

        # self.doc = FreeCAD.activeDocument()
        # self.file_name = self.doc.Label

    def open_file(self):
        inp_file_path = dir_name + '/control files/' + self.file_name + '.inp'
        print("open file: ", inp_file_path)
        try:
            with open(inp_file_path, "r") as f:
                w.YSinput.setText(str(f.readline().strip()))
                w.GXinput.setText(str(f.readline().strip()))
                w.GYinput.setText(str(f.readline().strip()))
                w.GZinput.setText(str(f.readline().strip()))
                w.steps.setText(str(f.readline().strip()))
                w.max_iter.setText(str(f.readline().strip()))
                w.error.setText(str(f.readline().strip()))
                w.relax.setText(str(f.readline().strip()))
                w.scale_1.setText(str(f.readline().strip()))
                w.scale_2.setText(str(f.readline().strip()))
                w.scale_3.setText(str(f.readline().strip()))
                if str(f.readline().strip()) == "total":
                    w.totalRbtn.setChecked(True)
                else:
                    w.incrRbtn.setChecked(True)

        except FileNotFoundError:
            print("File does not exist")
            w.YSinput.setText(self.YSinput_default)
            w.GXinput.setText(self.GXinput_default)
            w.GYinput.setText(self.GYinput_default)
            w.GZinput.setText(self.GZinput_default)
            w.steps.setText(self.steps_default)
            w.max_iter.setText(self.max_iter_default)
            w.error.setText(self.error_default)
            w.relax.setText(self.relax_default)
            w.scale_1.setText(self.scale_1_default)
            w.scale_2.setText(self.scale_2_default)
            w.scale_3.setText(self.scale_3_default)
            w.incrRbtn.setChecked(True)

    def max_iter_changed(self):
        if (w.max_iter.text() != self.max_iter_default):
            w.max_iter.setPalette(self.palette_warning)
        else:
            w.max_iter.setPalette(self.palette_standard)

    def relax_changed(self):
        if (w.relax.text() != self.relax_default):
            if w.relax.text() == "":
                w.relax.setText("0.0")
            if float(w.relax.text()) > 1.5:
                w.relax.setText("1.5")
            elif float(w.relax.text()) < 1.0:
                w.relax.setText("1.0")
            w.relax.setPalette(self.palette_warning)
        else:
            w.relax.setPalette(self.palette_standard)

    def scale_1_changed(self):
        if (w.scale_1.text() != self.scale_1_default):
            if w.scale_1.text() == "":
                w.scale_1.setText("0.0")
            if float(w.scale_1.text()) > 3.0:
                w.scale_1.setText("3.0")
            elif float(w.scale_1.text()) < 1.0:
                w.scale_1.setText("1.0")
            w.scale_1.setPalette(self.palette_warning)
        else:
            w.scale_1.setPalette(self.palette_standard)

    def scale_2_changed(self):
        if (w.scale_2.text() != self.scale_2_default):
            if w.scale_2.text() == "":
                w.scale_2.setText("0.0")
            if float(w.scale_2.text()) > 2.0:
                w.scale_2.setText("2.0")
            elif float(w.scale_2.text()) < 1.0:
                w.scale_2.setText("1.0")
            w.scale_2.setPalette(self.palette_warning)
        else:
            w.scale_2.setPalette(self.palette_standard)

    def scale_3_changed(self):
        if w.scale_3.text() == "":
            w.scale_3.setText("0.0")
        if (w.scale_3.text() != self.scale_3_default):
            if float(w.scale_3.text()) > 2.0:
                w.scale_3.setText("2.0")
            elif float(w.scale_3.text()) < 1.0:
                w.scale_3.setText("1.0")
            w.scale_3.setPalette(self.palette_warning)
        else:
            w.scale_3.setPalette(self.palette_standard)

    def btn_state(self):
        if w.totalRbtn.isChecked():
            self.disp_option = "total"
        if w.incrRbtn.isChecked():
            self.disp_option = "incremental"


FreeCADGui.addWorkbench(fcVMWorkbench)
