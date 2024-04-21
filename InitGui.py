# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2024 - Harry van Langen <hvlanalysis@gmail.com>        *
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

__Title__ = "fcVM"
__Author__ = "HarryvL"
__Url__ = "https://github.com/HarryvL/fcVM-workbench"
__Version__ = "1.1.0"
__Date__ = "2024/02/25"
__Comment__ = "first release"
__Forum__ = "https://forum.freecad.org/viewtopic.php?t=85474"
__Status__ = "initial development"
__Requires__ = "freecad version 0.19 or higher"
__Communication__ = "https://forum.freecad.org/viewtopic.php?t=85474"

import os
import sys
import dummyVM
import FreeCAD
import FreeCADGui
from PySide2 import QtWidgets, QtGui, QtCore

global FCmw
FCmw = FreeCADGui.getMainWindow()
global dir_name
dir_name = os.path.dirname(dummyVM.file_path())
global double_validator
double_validator = QtGui.QDoubleValidator()
global int_validator
int_validator = QtGui.QIntValidator()
global res_show
res_show = FreeCAD.getHomePath() + "Mod/Fem/Resources/ui/ResultShow.ui"
global source_code_path
source_code_path = os.path.join(dir_name, 'source code')

sys.path.append(source_code_path)


class fcVMWorkbench(Workbench):
    Icon = os.path.join(dir_name, "icons", "fcFEM.svg")
    MenuText = "fcVM workbench"
    ToolTip = "Plastic collapse analysis with the von Mises material model"

    def __init__(self):
        import task_result_mechanical
        sys.modules["femtaskpanels.task_result_mechanical"] = sys.modules[task_result_mechanical.__name__]

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
        global fcVM_window

        import dummyVM
        self.dir_name = os.path.dirname(dummyVM.file_path())


        self.doc = FreeCAD.activeDocument()
        if self.doc == None:
            self.doc = FreeCAD.newDocument("fcVM")

        self.file_name = self.doc.Label

        self.macro_file_path = os.path.join(self.dir_name, "source code", "fcVM.FCMacro")

        self.sum_file_path = os.path.join(self.dir_name, "source code", "fcVM_sum.FCMacro")

        self.disp_option = "incremental"

        self.csr_option = "PEEQ"

        self.averaged_option = "unaveraged"

        class DocObserver(object):  # document Observer
            def __init__(self, workbench_instance):
                self.workbench_instance = workbench_instance

            def slotActivateDocument(self, doc):
                if FreeCAD.activeDocument().Label[0:7] != "Unnamed":
                    self.workbench_instance.save_clicked()
                    self.workbench_instance.file_name = FreeCAD.activeDocument().Label
                    self.workbench_instance.open_file()
                    fcVM_window.PEEQ.setText("0.000")
                    fcVM_window.CSR.setText("0.000")

            def slotFinishSaveDocument(self, doc, prop):
                self.workbench_instance.save_clicked()  # save under old file name
                self.workbench_instance.file_name = doc.Label
                self.workbench_instance.save_clicked()  # save under new file name

        self.obs = DocObserver(self)
        FreeCAD.addDocumentObserver(self.obs)

        ui_Path = os.path.join(self.dir_name, "user_interface", "fcVM.ui")

        fcVM_window = FreeCADGui.PySideUic.loadUi(ui_Path)

        fcVM_window.startBtn.clicked.connect(self.start_clicked)
        fcVM_window.quitBtn.clicked.connect(self.Deactivated)
        fcVM_window.resetBtn.clicked.connect(self.reset_clicked)
        fcVM_window.saveBtn.clicked.connect(self.save_clicked)
        fcVM_window.sumBtn.clicked.connect(self.sum_clicked)
        fcVM_window.totalRbtn.toggled.connect(self.btn_state)
        fcVM_window.incrRbtn.toggled.connect(self.btn_state)
        fcVM_window.peeqRbtn.toggled.connect(self.btn_state)
        fcVM_window.csrRbtn.toggled.connect(self.btn_state)
        fcVM_window.averagedChk.toggled.connect(self.btn_state)

        fcVM_window.max_iter.textChanged.connect(self.max_iter_changed)
        fcVM_window.relax.textChanged.connect(self.relax_changed)
        fcVM_window.scale_1.textChanged.connect(self.scale_1_changed)
        fcVM_window.scale_2.textChanged.connect(self.scale_2_changed)
        fcVM_window.scale_3.textChanged.connect(self.scale_3_changed)
        fcVM_window.Hinput.textChanged.connect(self.Hinput_changed)

        fcVM_window.YSinput.setValidator(double_validator)
        fcVM_window.Hinput.setValidator(double_validator)
        fcVM_window.USinput.setValidator(double_validator)
        fcVM_window.GXinput.setValidator(double_validator)
        fcVM_window.GYinput.setValidator(double_validator)
        fcVM_window.GZinput.setValidator(double_validator)
        fcVM_window.steps.setValidator(int_validator)
        fcVM_window.max_iter.setValidator(int_validator)
        fcVM_window.error.setValidator(double_validator)
        fcVM_window.relax.setValidator(double_validator)
        fcVM_window.scale_1.setValidator(double_validator)
        fcVM_window.scale_2.setValidator(double_validator)
        fcVM_window.scale_3.setValidator(double_validator)
        fcVM_window.target_LF.setValidator(double_validator)

        self.YSinput_default = "240.0"
        self.Hinput_default = "0.0"
        self.USinput_default = "0.25"
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
        self.target_LF_default = "2.0"
        self.disp_option_default = "incremental"
        self.csr_option_default = "PEEQ"
        self.averaged_Chk_default = "unaveraged"

        self.open_file()

        FCmw.addDockWidget(QtCore.Qt.RightDockWidgetArea, fcVM_window.dw)

    def Deactivated(self):

        try:
            if fcVM_window.dw.isVisible():
                fcVM_window.dw.setVisible(False)
        except Exception:
            None

        FreeCAD.removeDocumentObserver(self.obs)

    def start_clicked(self):
        from PySide2 import QtWidgets

        self.save_clicked()

        try:
            rv = FCmw.findChild(QtWidgets.QTextEdit, "Report view")
            rv.clear()
        except Exception:
            None

        FreeCADGui.Selection.clearSelection()

        fcVM_macro = open(self.macro_file_path).read()
        print(self.macro_file_path)
        exec(fcVM_macro)

    def quit_clicked(self):
        self.Deactivated()

    def reset_clicked(self):
        fcVM_window.max_iter.setText(self.max_iter_default)
        fcVM_window.error.setText(self.error_default)
        fcVM_window.relax.setText(self.relax_default)
        fcVM_window.scale_1.setText(self.scale_1_default)
        fcVM_window.scale_2.setText(self.scale_2_default)
        fcVM_window.scale_3.setText(self.scale_3_default)
        fcVM_window.relax.setPalette(self.palette_standard)
        fcVM_window.scale_1.setPalette(self.palette_standard)
        fcVM_window.scale_2.setPalette(self.palette_standard)
        fcVM_window.scale_3.setPalette(self.palette_standard)
        fcVM_window.incrRbtn.setChecked(True)
        fcVM_window.averagedChk.setChecked(False)

    def save_clicked(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '.inp')
        with open(inp_file_path, "w") as f:
            f.write(fcVM_window.YSinput.text() + "\n")
            f.write(fcVM_window.GXinput.text() + "\n")
            f.write(fcVM_window.GYinput.text() + "\n")
            f.write(fcVM_window.GZinput.text() + "\n")
            f.write(fcVM_window.steps.text() + "\n")
            f.write(fcVM_window.max_iter.text() + "\n")
            f.write(fcVM_window.error.text() + "\n")
            f.write(fcVM_window.relax.text() + "\n")
            f.write(fcVM_window.scale_1.text() + "\n")
            f.write(fcVM_window.scale_2.text() + "\n")
            f.write(fcVM_window.scale_3.text() + "\n")
            f.write(self.disp_option + "\n")
            f.write(fcVM_window.USinput.text() + "\n")
            f.write(fcVM_window.Hinput.text() + "\n")
            f.write(fcVM_window.target_LF.text() + "\n")
            f.write(self.csr_option + "\n")
            f.write(self.averaged_option + "\n")

    def sum_clicked(self):
        fcVM_sum = open(self.sum_file_path).read()
        exec(fcVM_sum)

        return

    def open_file(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '.inp')
        try:
            with open(inp_file_path, "r") as f:
                fcVM_window.YSinput.setText(str(f.readline().strip()))
                fcVM_window.GXinput.setText(str(f.readline().strip()))
                fcVM_window.GYinput.setText(str(f.readline().strip()))
                fcVM_window.GZinput.setText(str(f.readline().strip()))
                fcVM_window.steps.setText(str(f.readline().strip()))
                fcVM_window.max_iter.setText(str(f.readline().strip()))
                fcVM_window.error.setText(str(f.readline().strip()))
                fcVM_window.relax.setText(str(f.readline().strip()))
                fcVM_window.scale_1.setText(str(f.readline().strip()))
                fcVM_window.scale_2.setText(str(f.readline().strip()))
                fcVM_window.scale_3.setText(str(f.readline().strip()))
                if str(f.readline().strip()) == "total":
                    fcVM_window.totalRbtn.setChecked(True)
                else:
                    fcVM_window.incrRbtn.setChecked(True)
                USinp = str(f.readline().strip())
                if USinp == "":
                    fcVM_window.USinput.setText(self.USinput_default)
                else:
                    fcVM_window.USinput.setText(USinp)
                Hinp = str(f.readline().strip())
                if Hinp == "":
                    fcVM_window.Hinput.setText(self.Hinput_default)
                else:
                    fcVM_window.Hinput.setText(Hinp)
                LFinp = str(f.readline().strip())
                if LFinp == "":
                    fcVM_window.target_LF.setText(self.target_LF_default)
                else:
                    fcVM_window.target_LF.setText(LFinp)
                csrBtninp = f.readline().strip()
                if csrBtninp == "CSR":
                    fcVM_window.csrRbtn.setChecked(True)
                else:
                    fcVM_window.peeqRbtn.setChecked(True)
                avBtninp = str(f.readline().strip())
                if avBtninp == "averaged":
                    fcVM_window.averagedChk.setChecked(True)
                else:
                    fcVM_window.averagedChk.setChecked(False)


        except FileNotFoundError:
            fcVM_window.YSinput.setText(self.YSinput_default)
            fcVM_window.USinput.setText(self.USinput_default)
            fcVM_window.GXinput.setText(self.GXinput_default)
            fcVM_window.GYinput.setText(self.GYinput_default)
            fcVM_window.GZinput.setText(self.GZinput_default)
            fcVM_window.steps.setText(self.steps_default)
            fcVM_window.max_iter.setText(self.max_iter_default)
            fcVM_window.error.setText(self.error_default)
            fcVM_window.relax.setText(self.relax_default)
            fcVM_window.scale_1.setText(self.scale_1_default)
            fcVM_window.scale_2.setText(self.scale_2_default)
            fcVM_window.scale_3.setText(self.scale_3_default)
            fcVM_window.incrRbtn.setChecked(True)
            fcVM_window.USinput.setText(self.USinput_default)
            fcVM_window.Hinput.setText(self.Hinput_default)
            fcVM_window.target_LF.setText(self.target_LF_default)
            fcVM_window.averagedChk.setChecked(False)

    def max_iter_changed(self):
        if (fcVM_window.max_iter.text() != self.max_iter_default):
            fcVM_window.max_iter.setPalette(self.palette_warning)
        else:
            fcVM_window.max_iter.setPalette(self.palette_standard)

    def relax_changed(self):
        if (fcVM_window.relax.text() != self.relax_default):
            if fcVM_window.relax.text() == "":
                fcVM_window.relax.setText("0.0")
            if float(fcVM_window.relax.text()) > 1.5:
                fcVM_window.relax.setText("1.5")
            elif float(fcVM_window.relax.text()) < 1.0:
                fcVM_window.relax.setText("1.0")
            fcVM_window.relax.setPalette(self.palette_warning)
        else:
            fcVM_window.relax.setPalette(self.palette_standard)

    def scale_1_changed(self):
        if (fcVM_window.scale_1.text() != self.scale_1_default):
            if fcVM_window.scale_1.text() == "":
                fcVM_window.scale_1.setText("0.0")
            if float(fcVM_window.scale_1.text()) > 3.0:
                fcVM_window.scale_1.setText("3.0")
            elif float(fcVM_window.scale_1.text()) < 1.0:
                fcVM_window.scale_1.setText("1.0")
            fcVM_window.scale_1.setPalette(self.palette_warning)
        else:
            fcVM_window.scale_1.setPalette(self.palette_standard)

    def scale_2_changed(self):
        if (fcVM_window.scale_2.text() != self.scale_2_default):
            if fcVM_window.scale_2.text() == "":
                fcVM_window.scale_2.setText("0.0")
            if float(fcVM_window.scale_2.text()) > 2.0:
                fcVM_window.scale_2.setText("2.0")
            elif float(fcVM_window.scale_2.text()) < 1.0:
                fcVM_window.scale_2.setText("1.0")
            fcVM_window.scale_2.setPalette(self.palette_warning)
        else:
            fcVM_window.scale_2.setPalette(self.palette_standard)

    def scale_3_changed(self):
        if fcVM_window.scale_3.text() == "":
            fcVM_window.scale_3.setText("0.0")
        if (fcVM_window.scale_3.text() != self.scale_3_default):
            if float(fcVM_window.scale_3.text()) > 2.0:
                fcVM_window.scale_3.setText("2.0")
            elif float(fcVM_window.scale_3.text()) < 1.0:
                fcVM_window.scale_3.setText("1.0")
            fcVM_window.scale_3.setPalette(self.palette_warning)
        else:
            fcVM_window.scale_3.setPalette(self.palette_standard)

    def Hinput_changed(self):
        if float(fcVM_window.Hinput.text()) < 0.0:
            fcVM_window.Hinput.setText("0.0")

    def btn_state(self):
        if fcVM_window.totalRbtn.isChecked():
            self.disp_option = "total"
        if fcVM_window.incrRbtn.isChecked():
            self.disp_option = "incremental"
        if fcVM_window.peeqRbtn.isChecked():
            self.csr_option = "PEEQ"
        if fcVM_window.csrRbtn.isChecked():
            self.csr_option = "CSR"
        if fcVM_window.averagedChk.isChecked():
            self.averaged_option = "averaged"
        else:
            self.averaged_option = "unaveraged"


FreeCADGui.addWorkbench(fcVMWorkbench)
