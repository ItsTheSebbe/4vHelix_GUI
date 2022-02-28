import os
from signal import pause
import sys
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from load_files import open_ply
import numpy as np


class Ui_MainWindow(object):

    def setupUi(self, MainWindow):

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1000, 1100)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        MainWindow.setCentralWidget(self.centralwidget)

        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(10, 20, 90, 40))
        self.pushButton.setObjectName("pushButton_open_file")

        self.pushButton_plot = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_plot.setGeometry(QtCore.QRect(110, 20, 90, 40))
        self.pushButton_plot.setObjectName("pushButton_plot")

        self.pushButton_openply = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_openply.setGeometry(QtCore.QRect(230, 20, 90, 40))
        self.pushButton_openply.setObjectName("pushButton_load_ply")
        self.pushButton_openply.clicked.connect(self.load_ply)

        self.pushButton_select = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select.setGeometry(QtCore.QRect(450, 20, 90, 40))
        self.pushButton_select.setObjectName("pushButton_select_edge")

        self.pushButton_reinforce = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_reinforce.setGeometry(QtCore.QRect(740, 20, 180, 40))
        self.pushButton_reinforce.setObjectName("pushButton_reinforce_edge")

        self.pushButton_select_all = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select_all.setGeometry(QtCore.QRect(550, 20, 90, 40))
        self.pushButton_select_all.setObjectName("pushButton_select_all_edge")

        self.glViewer = gl.GLViewWidget(self.centralwidget)
        self.glViewer.setGeometry(QtCore.QRect(50, 100, 900, 900))
        self.glViewer.setObjectName("GL_viewer")
        self.setup_glViewer(self.glViewer)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "4vHelix"))
        self.pushButton.setText(_translate("MainWindow", "Open rpoly"))
        self.pushButton_plot.setText(_translate("MainWindow", "Plot rpoly"))
        self.pushButton_openply.setText(_translate("MainWindow", "Open ply"))
        self.pushButton_select.setText(_translate(
            "MainWindow", "Select edge"))
        self.pushButton_reinforce.setText(
            _translate("MainWindow", "Reinforce selected edges"))
        self.pushButton_select_all.setText(
            _translate("MainWindow", "Select all"))

    def load_ply(self):
        """
        Load and plot ply file.
        """
        if Ply_Object.exists == False:
            self.ply = Ply_Object()  # Create new instance of Ply object
            self.ply.draw_ply(self.glViewer)
        else:
            self.ply.draw_ply(self.glViewer)


    def setup_glViewer(self, glViewer):
        glViewer.setCameraPosition(distance=10)

        grid = gl.GLGridItem()
        grid.scale(2, 2, 1)
        glViewer.addItem(grid)


class Ply_Object(QtWidgets.QWidget):

    exists = False

    def draw_ply(self, glViewer):
        """
        Draw ply file to glViewer
        """

        # Get file path
        file_path = QtWidgets.QFileDialog.getOpenFileName()
        self.vertNum, self.vertices, self.faceNum, self.faces, loadFlag = open_ply(
            str(file_path[0]))
        
        # Remove previous plot if exists
        if Ply_Object.exists == True:
            self.RemoveObject(glViewer)

        # Draw new plot
        if loadFlag == True:
            self.AddObject(glViewer)
            Ply_Object.exists = True
        else: 
            print("Unable to load .ply file!")

    def AddObject(self, glViewer):
        """
        Draw new object.
        """
        self.wireframe = gl.GLMeshItem(vertexes=self.vertices, faces=self.faces,
                                       smooth=False, drawEdges=True, drawFaces=False, edgeColor=(1, 1, 1, 1))
        glViewer.addItem(self.wireframe)

    def RemoveObject(self, glViewer):
        """
        Remove old object
        """
        glViewer.removeItem(self.wireframe)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
