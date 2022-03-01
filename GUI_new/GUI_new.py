import os
from signal import pause
import sys
from pygments import highlight
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from load_files import open_ply
import numpy as np


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

        self.pushButton_plot = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_plot.setGeometry(QtCore.QRect(110, 20, 90, 40))
        self.pushButton_plot.setObjectName("pushButton_plot")
        self.pushButton_plot.clicked.connect(self.PlotRpoly)

        self.pushButton_openply = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_openply.setGeometry(QtCore.QRect(230, 20, 90, 40))
        self.pushButton_openply.setObjectName("pushButton_load_ply")
        self.pushButton_openply.clicked.connect(self.load_ply)

        self.pushButton_select = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select.setGeometry(QtCore.QRect(350, 20, 90, 40))
        self.pushButton_select.setObjectName("pushButton_select_edge")
        self.pushButton_select.clicked.connect(self.select)

        self.pushButton_select_all = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_select_all.setGeometry(QtCore.QRect(450, 20, 90, 40))
        self.pushButton_select_all.setObjectName("pushButton_select_all_edge")
        self.pushButton_select_all.clicked.connect(self.SelectAll)

        self.pushButton_deselect_all = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_deselect_all.setGeometry(QtCore.QRect(550, 20, 90, 40))
        self.pushButton_deselect_all.setObjectName("pushButton_deselect_all_edge")
        self.pushButton_deselect_all.clicked.connect(self.DeselectAll)

        self.pushButton_reinforce = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_reinforce.setGeometry(QtCore.QRect(740, 20, 180, 40))
        self.pushButton_reinforce.setObjectName("pushButton_reinforce_edge")
        self.pushButton_reinforce.clicked.connect(self.Reinforce)
        
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
        self.pushButton_plot.setText(_translate("MainWindow", "Plot rpoly"))
        self.pushButton_openply.setText(_translate("MainWindow", "Open ply"))
        self.pushButton_select.setText(_translate(
            "MainWindow", "Select edge"))
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
            self.ply = Ply_Object(self.centralwidget)  # Create new instance of Ply object
            self.ply.draw_ply(self.glViewer)
        else:
            self.ply.draw_ply(self.glViewer)

    def setup_glViewer(self):
        self.glViewer.setCameraPosition(distance=10)

        grid = gl.GLGridItem()
        grid.scale(2, 2, 1)
        self.glViewer.addItem(grid)
    
    def select(self):
        print("select placeholder")

    def SelectAll(self):
        self.ply.HighlightAll(self.glViewer)
    
    def DeselectAll(self):
        self.ply.RemoveHighlightAll(self.glViewer)
    
    def Reinforce(self):
        print("Reinforce placeholder")

    def PlotRpoly(self):
        print("plot rpoly placeholder")

    def OpenRpoly(self):
        print("open rpoly placeholder")


class check_boxes(QtWidgets.QWidget):
    exists = False

    def __init__(self, parent, edgeNum, edges, ply, glViewer):

        super(check_boxes, self).__init__(parent)

        check_boxes.exists = True
        self.parent = parent
        self.edgeNum = edgeNum
        self.selected_edges = edges
        self.ply = ply
        self.glViewer = glViewer
        
        # Setup buttons
        self.setupUI()
        # Update states according to
        self.update_checkboxes()
    
    def ClearCheckbox(self):
        for i in range(self.edgeNum):
            self.box[i].deleteLater()
  

    def setupUI(self):
        """
        Set up of check box UI.
        """
        vertPos = 0
        cnt = 0
        self.box = {}
        for i in range(self.edgeNum):
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
        for i in range(self.edgeNum):
            # Make sure click_on_check_box is not called
            self.box[i].blockSignals(True)

            # Check if al checkboxes are set correctly
            if i in self.selected_edges:
                self.box[i].setChecked(True)
            else:
                self.box[i].setChecked(False)

            # Turn signals back on
            self.box[i].blockSignals(False)

        self.ply.UpdateHighlight(self.glViewer)
        print("Selected edges: " + str(self.selected_edges))

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
        self.ply.UpdateHighlight(self.glViewer)
    

class Ply_Object(QtWidgets.QWidget):

    exists = False

    def __init__(self, centralwidget):
        self.selectedEdges = []
        self.centralwidget = centralwidget

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
            self.ClearScreen(glViewer)
            self.check_boxes.ClearCheckbox()
            Ply_Object.exists = False

        # Draw new plot
        if loadFlag == True:
            self.CountEdges()
            self.PlotPly(glViewer)
            Ply_Object.exists = True

            self.check_boxes = check_boxes(self.centralwidget, self.edgeNum, self.selectedEdges, self, glViewer)

        else: 
            print("Unable to load .ply file!")

    def PlotPly(self, glViewer):
        """
        Draw new object.
        """
        self.wireframe = gl.GLMeshItem(vertexes=self.vertices, faces=self.faces,
                                       smooth=False, drawEdges=True, drawFaces=False, edgeColor=(1, 1, 1, 1))
        glViewer.addItem(self.wireframe)

    def ClearScreen(self, glViewer):
        """
        Remove old object
        """
        glViewer.removeItem(self.wireframe)
        for i in range(self.edgeNum):
            if i in self.selectedEdges:
                glViewer.removeItem(self.highlights[i])
                self.selectedEdges.remove(i)
        self.highlights = np.array([], dtype=gl.GLLinePlotItem)
    
    def CountEdges(self):
        """
        Finds number of unique edges in ply
        """
        self.edges = []
        for i in range(self.faceNum):
            for j in range(len(self.faces[i,:])):
                if j == len(self.faces[i,:])-1:
                    point1 = self.faces[i,j]
                    point2 = self.faces[i,0]
                else:
                    point1 = self.faces[i,j]
                    point2 = self.faces[i,j+1]
                
                if ([point1, point2] not in self.edges) and ([point2, point1] not in self.edges):
                    self.edges.append([point1, point2])
        
        self.edgeNum = len(self.edges)
        self.edges = np.array(self.edges)
        self.highlights = np.zeros(self.edgeNum, dtype=gl.GLLinePlotItem)


    def LoadEdge(self, lineNum):
        """
        Returns two points of edge for given line number
        """
        id1 = self.edges[lineNum,0]
        id2 = self.edges[lineNum,1]
        point1 = (self.vertices[id1,0],self.vertices[id1,1], self.vertices[id1,2])
        point2 = (self.vertices[id2,0],self.vertices[id2,1], self.vertices[id2,2])   
        return np.array([point1, point2])

    def UpdateHighlight(self, glViewer):
        for i in range(self.edgeNum):
            if i in self.selectedEdges and self.highlights[i] == 0:
                self.AddHighlight(i, glViewer)
            elif i not in self.selectedEdges and self.highlights[i] != 0:
                self.RemoveHighlight(i, glViewer)

    def RemoveHighlight(self, lineNum, glViewer):
        """
        Remove highlight of specific line
        """
        glViewer.removeItem(self.highlights[lineNum])
        self.highlights[lineNum] = 0

    def AddHighlight(self, lineNum, glViewer):
        """
        Highlight specific line
        """
        # Load points
        pts = self.LoadEdge(lineNum)
       
        line = gl.GLLinePlotItem(pos=pts, width=10, antialias=False, color=(255,0,0,1))
        self.highlights[lineNum] = line
        
        glViewer.addItem(self.highlights[lineNum])
        # check_boxes.update_checkboxes()
    
    def HighlightAll(self, glViewer):
        """
        Select all
        """
        for i in range(self.edgeNum):
            if i not in self.selectedEdges:
                pts = self.LoadEdge(i)
                line = gl.GLLinePlotItem(pos=pts, width=10, antialias=False, color=(255,0,0,1))
                self.highlights[i] = line
                self.selectedEdges.append(i)
                glViewer.addItem(self.highlights[i])
        check_boxes.update_checkboxes()

        print(self.selectedEdges)
    
    def RemoveHighlightAll(self, glViewer):
        """
        Select all
        """
        for i in range(self.edgeNum):
            if i in self.selectedEdges:
                glViewer.removeItem(self.highlights[i])
                self.selectedEdges.remove(i)
                self.highlights[i] = 0

        print(self.selectedEdges)
    
      
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(win)
    win.show()
    sys.exit(app.exec_())
