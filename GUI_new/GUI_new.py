import os
import sys
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from load_files import open_ply, open_rpoly

from pyquaternion import Quaternion
from tacoxDNA.src.libs import cadnano_utils as cu
from tacoxDNA.src.libs import base
from load_files import move_along_vector

class Ui_MainWindow(object):

    def setupUi(self, MainWindow):

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1100, 1100)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        MainWindow.setCentralWidget(self.centralwidget)

        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(10, 20, 90, 40))
        self.pushButton.setObjectName("pushButton_open_file")
        self.pushButton.clicked.connect(self.OpenRpoly)

        self.pushButton_openply = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_openply.setGeometry(QtCore.QRect(230, 20, 90, 40))
        self.pushButton_openply.setObjectName("pushButton_load_ply")
        self.pushButton_openply.clicked.connect(self.load_ply)

        self.pushButton_select_all = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select_all.setGeometry(QtCore.QRect(350, 20, 90, 40))
        self.pushButton_select_all.setObjectName("pushButton_select_all_edge")
        self.pushButton_select_all.clicked.connect(self.AddAllHighlight)

        self.pushButton_deselect_all = QtWidgets.QPushButton(
            self.centralwidget)
        self.pushButton_deselect_all.setGeometry(QtCore.QRect(450, 20, 90, 40))
        self.pushButton_deselect_all.setObjectName(
            "pushButton_deselect_all_edge")
        self.pushButton_deselect_all.clicked.connect(self.DeselectAll)

        self.pushButton_reinforce = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_reinforce.setGeometry(QtCore.QRect(650, 20, 180, 40))
        self.pushButton_reinforce.setObjectName("pushButton_reinforce_edge")
        self.pushButton_reinforce.clicked.connect(self.Reinforce)

        self.pushButton_seq_designer = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_seq_designer.setGeometry(QtCore.QRect(850, 20, 180, 40))
        self.pushButton_seq_designer.setObjectName("pushButton_seq_designer")
        self.pushButton_seq_designer.clicked.connect(self.RunSequenceDesigner)

        self.glViewer = gl.GLViewWidget(self.centralwidget)
        self.glViewer.setGeometry(QtCore.QRect(50, 100, 900, 900))
        self.glViewer.setObjectName("GL_viewer")
        self.setup_glViewer()

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "4vHelix"))
        self.pushButton.setText(_translate("MainWindow", "Open rpoly"))
        self.pushButton_openply.setText(_translate("MainWindow", "Open ply"))
        self.pushButton_seq_designer.setText(_translate(
            "MainWindow", "Run Sequence Designer"))
        self.pushButton_reinforce.setText(
            _translate("MainWindow", "Reinforce selected edges"))
        self.pushButton_select_all.setText(
            _translate("MainWindow", "Select all"))
        self.pushButton_deselect_all.setText(
            _translate("MainWindow", "Deselect all"))

    def load_ply(self):
        """
        Load and plot ply file.
        """
        
        if Ply_Object.exists == False:
            self.ply = Ply_Object(self.glViewer)  # Create new instance of Ply object
            self.ply.draw_ply()
        else:
            self.ply.draw_ply()

    def setup_glViewer(self):
        """
        Set up initial condition of glviewer
        """
        self.glViewer.setCameraPosition(distance=10)
        grid = gl.GLGridItem()
        grid.scale(2, 2, 1)
        self.glViewer.addItem(grid)

    def AddAllHighlight(self):
        """
        Select all edges
        """
        if Ply_Object.exists == True:
            self.ply.AddAllHighlight()
        else:
            print("No file selected!")

    def DeselectAll(self):
        """
        Select all edges
        """
        if Ply_Object.exists == True:
            self.ply.RemoveAllHighlight()
        else:
            print("No file selected!")

    def RunSequenceDesigner(self):
        print("Run sequence designer placeholder")

    def Reinforce(self):
        print("Reinforce placeholder")


    def OpenRpoly(self):
        if Rpoly_Object.exists == False:
            self.ply = Rpoly_Object(self.glViewer)  # Create new instance of Ply object
            self.ply.draw_rpoly()
        else:
            self.ply.draw_rpoly()


class check_boxes(QtWidgets.QWidget):
    exists = False

    def __init__(self, plotObj):
        super(check_boxes, self).__init__(None)
    
        check_boxes.exists = True
        self.plotObj = plotObj

        # Setup buttons
        self.CreateCheckboxes()
        # Update states according to
        self.update_checkboxes()

    def RemoveCheckboxes(self):
        """
        Removes all checkboxes
        """
        for i in range(self.plotObj.edgeNum):
            self.box[i].deleteLater()

    def CreateCheckboxes(self):
        """
        Create all checkboxes
        """
        vertPos = 0
        cnt = 0
        self.box = {}
        for i in range(self.plotObj.edgeNum):
            self.box[cnt] = QtWidgets.QCheckBox(str(i), win)
            self.box[cnt].move(1000, vertPos+100)
            self.box[i].stateChanged.connect(self.click_on_check_box)
            self.box[i].show()
            vertPos += 20
            cnt += 1

    def update_checkboxes(self):
        """
        Updates check boxes and corresponding line colors
        """
        for i in range(self.plotObj.edgeNum):
            # Make sure click_on_check_box is not called
            self.box[i].blockSignals(True)

            # Check if al checkboxes are set correctly
            if i in self.plotObj.selectedEdges:
                self.box[i].setChecked(True)
            else:
                self.box[i].setChecked(False)

            # Turn signals back on
            self.box[i].blockSignals(False)

        self.plotObj.UpdateHighlight()
        print("Selected edges: " + str(self.plotObj.selectedEdges))

    def click_on_check_box(self):
        """
        Add/remove selected edges
        """
        checkBox = self.sender()
        selected = int(checkBox.text())
        if selected in self.plotObj.selectedEdges:
            self.plotObj.selectedEdges.remove(selected)
        else:
            self.plotObj.selectedEdges.append(selected)
        self.update_checkboxes()


class Ply_Object(QtWidgets.QWidget):

    exists = False # For checking whether object exists


    def __init__(self, glViewer):
        self.selectedEdges = []
        self.glViewer = glViewer


    def draw_ply(self):
        """
        Draw ply file to glViewer
        """
        # Get file path
        file_path = QtWidgets.QFileDialog.getOpenFileName()
        self.vertNum, self.vertices, self.faceNum, self.faces, loadFlag = open_ply(
            str(file_path[0]))

        # Remove previous plot if exists
        if Ply_Object.exists == True:
            self.ClearScreen()
            Ply_Object.exists = False

        # Draw new plot
        if loadFlag == True:
            self.CountEdges()
            self.PlotPly()
            Ply_Object.exists = True

            self.check_boxes = check_boxes(self)

        else:
            print("Unable to load .ply file!")


    def PlotPly(self):
        """
        Draw new wireframe object
        """
        self.wireframe = gl.GLMeshItem(vertexes=self.vertices, faces=self.faces,
                                       smooth=False, drawEdges=True, drawFaces=False, edgeColor=(1, 1, 1, 1))
        self.glViewer.addItem(self.wireframe)


    def ClearScreen(self):
        """
        Remove all old object
        """

        # Remove wireframe
        self.glViewer.removeItem(self.wireframe)

        # Remove highlights
        self.RemoveAllHighlight()

        # Remove checkboxes
        self.check_boxes.RemoveCheckboxes()

    def CountEdges(self):
        """
        Finds number of unique edges in ply.
        """
        self.edges = []
        for i in range(self.faceNum):
            for j in range(len(self.faces[i, :])):

                # Check for duplicate inversions
                if j == len(self.faces[i, :])-1:
                    point1 = self.faces[i, j]
                    point2 = self.faces[i, 0]
                else:
                    point1 = self.faces[i, j]
                    point2 = self.faces[i, j+1]

                # If no duplicates found in list, append
                if ([point1, point2] not in self.edges) and ([point2, point1] not in self.edges):
                    self.edges.append([point1, point2])

        self.edgeNum = len(self.edges)  # Number of edges
        self.edges = np.array(self.edges)  # List of edges

        # List of highlighted edges
        self.highlights = np.zeros(self.edgeNum, dtype=gl.GLLinePlotItem)

    def LoadEdge(self, lineNum):
        """
        Returns two vertices of edge for given line number.
        """
        id1 = self.edges[lineNum, 0]
        id2 = self.edges[lineNum, 1]
        point1 = (self.vertices[id1, 0],
                  self.vertices[id1, 1], self.vertices[id1, 2])
        point2 = (self.vertices[id2, 0],
                  self.vertices[id2, 1], self.vertices[id2, 2])
        return np.array([point1, point2])

    def UpdateHighlight(self):
        """
        Match plotted highlights with selected edges list.
        """
        for i in range(self.edgeNum):
            if i in self.selectedEdges and self.highlights[i] == 0:
                self.AddHighlight(i)
            elif i not in self.selectedEdges and self.highlights[i] != 0:
                self.RemoveHighlight(i)

    def RemoveHighlight(self, lineNum):
        """
        Remove highlight of specific line from glViewer.
        """
        self.glViewer.removeItem(self.highlights[lineNum])
        self.highlights[lineNum] = 0

    def AddHighlight(self, lineNum):
        """
        Highlight specific line to glViewer.
        """
        # Load points
        pts = self.LoadEdge(lineNum)

        line = gl.GLLinePlotItem(
            pos=pts, width=10, antialias=False, color=(255, 0, 0, 1))
        self.highlights[lineNum] = line
        self.glViewer.addItem(self.highlights[lineNum])

    def AddAllHighlight(self):
        """
        Select all edges.
        """
        for i in range(self.edgeNum):
            if i not in self.selectedEdges:

                # Plot in viewer
                self.AddHighlight(i)

                # Add to list of selected edges
                self.selectedEdges.append(i)
        # Update checkboxes
        self.check_boxes.update_checkboxes()

    def RemoveAllHighlight(self):
        """
        Remove all selected edges.
        """
        for i in range(self.edgeNum):
            if i in self.selectedEdges:

                # Remove from viewer
                self.RemoveHighlight(i)

                # Remove from list of selected edges
                self.selectedEdges.remove(i)

        # Update checkboxes
        self.check_boxes.update_checkboxes()


class Rpoly_Object(QtWidgets.QWidget):

    exists = False # For checking whether object exists


    def __init__(self, glViewer):
        self.selectedEdges = []
        self.glViewer = glViewer


    def draw_rpoly(self):
        """
        Draw ply file to glViewer
        """
        # Get file path
        file_path = QtWidgets.QFileDialog.getOpenFileName()

        self.rpoly_data, _, _, loadFlag = open_rpoly(
            str(file_path[0]))

        print('.rpoly opened!')

        # Remove previous plot if exists
        if Ply_Object.exists == True:
            self.ClearScreen()
            Ply_Object.exists = False

        # Draw new plot
        if loadFlag == True:
            self.CreatePointList()
            # self.CountEdges()
            self.plot()
            Ply_Object.exists = True

            self.check_boxes = check_boxes(self)
        else:
            print("Unable to load .rpoly file!")



    def plot(self):
        self.highlights = np.zeros(self.edgeNum, dtype=gl.GLLinePlotItem)
        self.wireframe = np.zeros(self.edgeNum, dtype=gl.GLLinePlotItem)
        for n in range(self.edgeNum):
            x1, y1, z1 = self.x_list[n], self.y_list[n], self.z_list[n]
            x2, y2, z2 = self.x_list[n +
                                     1], self.y_list[n + 1], self.z_list[n + 1]
            p1 = (x1, y1, z1)
            p2 = (x2, y2, z2)
            pts = np.array([p1, p2])
            print(pts)

            line = gl.GLLinePlotItem(
                pos=pts, width=1, antialias=False, color=(255, 255, 255, 1))
            self.wireframe[n] = line
            self.glViewer.addItem(self.wireframe[n])


    def CreatePointList(self):
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

            self.edgeNum = len(self.x_list)-1


    def ClearScreen(self):
        """
        Remove all old object
        """

        # Remove wireframe
        self.glViewer.removeItem(self.wireframe)

        # Remove highlights
        self.RemoveAllHighlight()

        # Remove checkboxes
        self.check_boxes.RemoveCheckboxes()


    def LoadEdge(self, lineNum):
        """
        Returns two vertices of edge for given line number.
        """
        id1 = self.edges[lineNum, 0]
        id2 = self.edges[lineNum, 1]
        point1 = (self.vertices[id1, 0],
                  self.vertices[id1, 1], self.vertices[id1, 2])
        point2 = (self.vertices[id2, 0],
                  self.vertices[id2, 1], self.vertices[id2, 2])
        return np.array([point1, point2])

    def UpdateHighlight(self):
        """
        Match plotted highlights with selected edges list.
        """
        for i in range(self.edgeNum):
            if i in self.selectedEdges and self.highlights[i] == 0:
                self.AddHighlight(i)
            elif i not in self.selectedEdges and self.highlights[i] != 0:
                self.RemoveHighlight(i)

    def RemoveHighlight(self, lineNum):
        """
        Remove highlight of specific line from glViewer.
        """
        self.glViewer.removeItem(self.highlights[lineNum])
        self.highlights[lineNum] = 0

    def AddHighlight(self, lineNum):
        """
        Highlight specific line to glViewer.
        """
        # Load points
        pts = self.LoadEdge(lineNum)

        line = gl.GLLinePlotItem(
            pos=pts, width=10, antialias=False, color=(255, 0, 0, 1))
        self.highlights[lineNum] = line
        self.glViewer.addItem(self.highlights[lineNum])

    def AddAllHighlight(self):
        """
        Select all edges.
        """
        for i in range(self.edgeNum):
            if i not in self.selectedEdges:

                # Plot in viewer
                self.AddHighlight(i)

                # Add to list of selected edges
                self.selectedEdges.append(i)
        # Update checkboxes
        self.check_boxes.update_checkboxes()

    def RemoveAllHighlight(self):
        """
        Remove all selected edges.
        """
        for i in range(self.edgeNum):
            if i in self.selectedEdges:

                # Remove from viewer
                self.RemoveHighlight(i)

                # Remove from list of selected edges
                self.selectedEdges.remove(i)

        # Update checkboxes
        self.check_boxes.update_checkboxes()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(win)
    win.show()
    sys.exit(app.exec_())
