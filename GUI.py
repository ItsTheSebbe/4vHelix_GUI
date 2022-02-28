# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow_1.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QVBoxLayout

import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import CheckButtons
from matplotlib.lines import Line2D

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3D

from pyquaternion import Quaternion

from tacoxDNA.src.libs import cadnano_utils as cu
from tacoxDNA.src.libs import base

from load_files import move_along_vector
from load_files import open_rpoly, open_ply
from vHelix_auto_2 import GenerateJson


class Ui_MainWindow(object):

    def setupUi(self, MainWindow):

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1000, 108)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        MainWindow.setCentralWidget(self.centralwidget)

        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(10, 20, 90, 40))
        self.pushButton.setObjectName("pushButton_open_file")
        self.pushButton.clicked.connect(self.load_rpoly)

        self.pushButton_plot = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_plot.setGeometry(QtCore.QRect(110, 20, 90, 40))
        self.pushButton_plot.setObjectName("pushButton_plot")
        self.pushButton_plot.clicked.connect(self.plot_rpoly)

        self.pushButton_loadply = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_loadply.setGeometry(QtCore.QRect(230, 20, 90, 40))
        self.pushButton_loadply.setObjectName("pushButton_load_ply")
        self.pushButton_loadply.clicked.connect(self.load_ply)

        self.pushButton_plot_ply = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_plot_ply.setGeometry(QtCore.QRect(330, 20, 90, 40))
        self.pushButton_plot_ply.setObjectName("pushButton_plot_ply")
        self.pushButton_plot_ply.clicked.connect(self.plot_ply)

        self.pushButton_select = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select.setGeometry(QtCore.QRect(450, 20, 90, 40))
        self.pushButton_select.setObjectName("pushButton_select_edge")
        self.pushButton_select.clicked.connect(self.open_checkbox)

        self.pushButton_reinforce = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_reinforce.setGeometry(QtCore.QRect(740, 20, 180, 40))
        self.pushButton_reinforce.setObjectName("pushButton_reinforce_edge")
        self.pushButton_reinforce.clicked.connect(self.reinforce_selected)

        self.pushButton_select_all = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select_all.setGeometry(QtCore.QRect(550, 20, 90, 40))
        self.pushButton_select_all.setObjectName("pushButton_select_all_edge")
        self.pushButton_select_all.clicked.connect(self.select_all)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pushButton.setText(_translate("MainWindow", "Open rpoly"))
        self.pushButton_plot.setText(_translate("MainWindow", "Plot rpoly"))
        self.pushButton_loadply.setText(_translate("MainWindow", "Open ply"))
        self.pushButton_plot_ply.setText(_translate("MainWindow", "Plot ply"))
        self.pushButton_select.setText(_translate(
            "MainWindow", "Select edge"))
        self.pushButton_reinforce.setText(
            _translate("MainWindow", "Reinforce selected edges"))
        self.pushButton_select_all.setText(
            _translate("MainWindow", "Select all"))

    def load_rpoly(self):
        self.rpoly = Rpoly_Object()
        self.rpoly.load_rpoly()

    def plot_rpoly(self):
        self.rpoly.plot()
        self.expand()
        self.rpoly.move(0, 100)
        self.rpoly.show()

    def load_ply(self):
        self.ply = Ply_Object()

    def plot_ply(self):
        self.ply.plot()
        self.expand()
        self.ply.move(0, 70)
        self.ply.show()

    def expand(self):
        MainWindow.resize(1000, 1100)

    def open_checkbox(self):
        self.rpoly.create_checkboxes()

    def reinforce_selected(self):
        # print(self.rpoly.selected_edges)
        for i in range(len(self.rpoly.x_list)):
            self.rpoly.selected_edges[i] = self.rpoly.selected_edges[i] + 1
        GenerateJson(self.rpoly.fileNameNoExt, self.rpoly.selected_edges)

    def select_all(self):
        self.rpoly.select_all()


class Ply_Object(QtWidgets.QWidget):

    ply_exists = 0

    def __init__(self):
        self.load_ply()
        self.setupUi()

    def setupUi(self):
        QtWidgets.QWidget.__init__(self, MainWindow)

        self.setWindowTitle('Ply Mesh Plot')
        self.resize(1000, 1000)
        self.main_widget = QtWidgets.QWidget()
        self.fig = Figure(figsize=(15, 14.5), dpi=80)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_widget)
        self.ax = Axes3D(self.fig)

        vbox = QVBoxLayout(self)
        vbox.addWidget(self.canvas)

    def load_ply(self):

        file_path = QtWidgets.QFileDialog.getOpenFileName()

        self.number_vertices, self.vertices_list, self.number_face, self.faces_list = open_ply(
            str(file_path[0]))

        print('.ply opened!')

    def plot(self):
        for i in range(self.number_face):
            for j in range(1, self.faces_list[i][0]+1):
                if j == len(self.faces_list[0])-1:
                    point1 = self.faces_list[i, j]
                    point2 = self.faces_list[i, 1]
                else:
                    point1 = self.faces_list[i, j]
                    point2 = self.faces_list[i, j+1]

                xx = [self.vertices_list[point1, 0],
                      self.vertices_list[point2, 0]]
                yy = [self.vertices_list[point1, 1],
                      self.vertices_list[point2, 1]]
                zz = [self.vertices_list[point1, 2],
                      self.vertices_list[point2, 2]]
                self.ax.plot(xx, yy, zz, 'r')

        self.ax.scatter(
            self.vertices_list[:, 0], self.vertices_list[:, 1], self.vertices_list[:, 2])

        for i in range(self.number_vertices):
            self.ax.text(
                self.vertices_list[i, 0], self.vertices_list[i, 1], self.vertices_list[i, 2], s=i)


class Rpoly_Object(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self, MainWindow)
        self.setupUi()
        self.created_checkboxes = False

    def setupUi(self):
        self.setWindowTitle('Rpoly Mesh Plot')
        self.resize(1000, 1000)
        self.main_widget = QtWidgets.QWidget()
        self.fig = Figure(figsize=(15, 14.5), dpi=80)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_widget)
        self.ax = Axes3D(self.fig)

        vbox = QVBoxLayout(self)
        vbox.addWidget(self.canvas)

    def load_rpoly(self):

        file_path = QtWidgets.QFileDialog.getOpenFileName()

        self.rpoly_data, _, _ = open_rpoly(
            str(file_path[0]))

        fileName = os.path.basename(file_path[0])
        fileName = os.path.splitext(fileName)[0]
        self.fileNameNoExt = str(fileName)
        self.LinePicker()
        self.numEdges = len(self.x_list)-1
        print('.rpoly opened!')

    def LinePicker(self):
        self.selected_edges = []
        position_list = []
        generator = cu.StrandGenerator()

        new_position_list = []
        self.x_list, self.y_list, self.z_list = [], [], []
        scaffold_fragments = base.System([100, 100, 100])

        for n, i in enumerate(self.rpoly_data):

            position = [float(i[3]) / 0.84, float(i[4]) / 0.84,
                        float(i[5]) / 0.84]  # 0.84 scaling is ad hoc solution to get good looking models
            position_list.append(position)

            q = Quaternion(w=float(i[9]), x=float(i[6]), y=float(i[7]),
                           z=float(i[8]))  # find the helix rotation Info from file
            vec = q.rotate(
                np.array(
                    [0.0, 0.0, 1.0]))  # use it to figure out direction vec = q.rotate(np.array([0.0, 0.0, 1.0]))
            vec2 = q.rotate([0.65, -0.76, 0.0])
            n_bp = int(i[2])
            # calculate the position of start base
            new_position = move_along_vector(position, vec, n_bp)
            new_position_list.append(new_position)

            new_strands = generator.generate_or_sq(bp=n_bp, start_pos=new_position, direction=vec,
                                                   perp=vec2)  # generate strands

            (scaffold_fragments.add_strand(new_strands[0]))

            # get the sequence only for scaffold (like the one in .conf file)
            sequence = new_strands[0]._get_Marco_output()

            sequence_list = sequence.split('\n')

            base_coord_list = []

            for sequence in sequence_list:  # get coordinates
                base_coord = sequence.split(' ')
                base_coord_list.append(base_coord)

            self.x_list.append(float(base_coord_list[0][0]))
            self.y_list.append(float(base_coord_list[0][1]))
            self.z_list.append(float(base_coord_list[0][2]))

    def plot(self):
        self.lines_list = []
        self.labels_list = []

        self.selected_edges_label = QtWidgets.QLabel(MainWindow)
        self.selected_edges_label.move(50, 70)
        self.selected_edges_label.setText(
            "Selected Edges: " + str(self.selected_edges))
        self.selected_edges_label.adjustSize()
        self.selected_edges_label.show()

        for n in range(self.numEdges):
            x1, y1, z1 = self.x_list[n], self.y_list[n], self.z_list[n]
            x2, y2, z2 = self.x_list[n +
                                     1], self.y_list[n + 1], self.z_list[n + 1]
            self.line = self.ax.plot([x1, x2], [y1, y2], [z1, z2], marker='o', color='r',
                                     linewidth=3, picker=True, label=(n+1))
            self.labels_list.append(str(n))
            self.lines_list.append(self.line)

            self.ax.text(x1, y1, z1, s=str(n + 1))

        self.line = self.ax.plot([self.x_list[-1], self.x_list[0]], [self.y_list[-1], self.y_list[0]], [self.z_list[-1], self.z_list[0]], color='r',
                                 linewidth=3, picker=True, label=(n+2))
        self.labels_list.append(str(n+2))
        self.lines_list.append(self.line)
        # print(self.lines_list[0][0].get_data_3d())

        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.set_zlabel('Z axis')
        max_range = max([max(self.x_list) - min(self.x_list), max(self.y_list) - min(self.y_list),
                         max(self.z_list) - min(self.z_list)]) / 2.0
        mid_x = (max(self.x_list) + min(self.x_list)) * 0.5
        mid_y = (max(self.y_list) + min(self.y_list)) * 0.5
        mid_z = (max(self.z_list) + min(self.z_list)) * 0.5

        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

        self.fig.canvas.mpl_connect('pick_event', self.select_line)

    def select_line(self, event):
        # Check selected line
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            label = thisline.get_label()
            label = int(label)-1
            if label in self.selected_edges:
                self.selected_edges.remove(label)
            else:
                self.selected_edges.append(label)

        # Update color of line
        self.update_color()
        if self.created_checkboxes == True:
            self.check_boxes.update_checkboxes()

    def update_color(self):
        # Update selected edges label
        self.selected_edges_label.setText(
            "Selected Edges: " + str(self.selected_edges))
        self.selected_edges_label.adjustSize()

        # Update color of line
        for label in range(self.numEdges):
            if label in self.selected_edges:
                (self.lines_list[int(label)+1][0].set_color('b'))
            else:
                (self.lines_list[int(label)+1][0].set_color('r'))
        self.canvas.draw()

    def select_all(self):
        """
        Select all available edges
        """
        self.selected_edges = []
        for i in range(self.numEdges):
            self.selected_edges.append(int(i))
        # Update color of line
        print(self.selected_edges)
        self.update_color()
        if self.created_checkboxes == True:
            self.check_boxes.update_checkboxes()
        

    def create_checkboxes(self):
        """
        Create checkboxes object
        """
        self.check_boxes = check_boxes(self)
        self.check_boxes.show()
        self.created_checkboxes = True


class check_boxes(QtWidgets.QWidget):

    def __init__(self, rpoly):
        self.rpoly = rpoly
        self.x_list = rpoly.x_list
        self.selected_edges = rpoly.selected_edges
        QtWidgets.QCheckBox.__init__(self)

        # Setup buttons
        self.setupUI()
        # Update states according to
        self.update_checkboxes()
        # Add click function
        for i in range(len(self.x_list)):
            self.box[i].stateChanged.connect(self.click_on_check_box)

    def setupUI(self):
        """
        Set up of check box UI.
        """
        self.resize(200, len(self.x_list)*22)
        i = 0
        cnt = 0
        self.box = {}
        self.label_list = []
        for x in self.x_list:
            self.label = self.x_list.index(x)
            self.label_list.append(self.label)
            self.box[cnt] = QtWidgets.QCheckBox(str(self.label), self)
            self.box[cnt].move(10, i+20)
            i += 20
            cnt += 1

    def update_checkboxes(self):
        """
        Updates check boxes and corresponding line colors
        """
        cnt = 0
        for x in self.x_list:
            # Make sure click_on_check_box is not called
            self.box[cnt].blockSignals(True)

            # Check if al checkboxes are set correctly
            if self.x_list.index(x) in self.selected_edges:
                self.box[cnt].setChecked(True)
            else:
                self.box[cnt].setChecked(False)

            # Turn signals back on
            self.box[cnt].blockSignals(False)

            cnt += 1
        self.rpoly.update_color()
        # print(str(self.selected_edges))

    def click_on_check_box(self):
        """
        Add/remove selected edges
        """
        checkBox = self.sender()
        selected = int(checkBox.text())
        if selected in self.selected_edges:
            self.selected_edges.remove(selected)
        else:
            self.selected_edges.append(selected)
        self.update_checkboxes()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
