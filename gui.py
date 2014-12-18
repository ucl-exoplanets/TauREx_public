

import sys

from PyQt4 import QtCore, QtGui, uic

import matplotlib
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.use('Qt4Agg')

#loading classes
sys.path.append('./classes')
sys.path.append('./library')


import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from tp_profile import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

gui_class = uic.loadUiType('gui.ui')[0]

class Qt4MplCanvas(FigureCanvas):

    def __init__(self, parent):
        # plot definition
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        t = np.arange(0.0, 3.0, 0.01)
        s = np.cos(2*np.pi*t)
        self.axes.plot(t, s)
        # initialization of the canvas
        FigureCanvas.__init__(self, self.fig)
        # set the parent widget
        self.setParent(parent)
        # we define the widget as expandable
        FigureCanvas.setSizePolicy(self,
        QtGui.QSizePolicy.Expanding,
        QtGui.QSizePolicy.Expanding)
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self)

class PlotWindow(QtGui.QMainWindow, gui_class):
    def __init__(self):
        # initialization of Qt MainWindow widget
        QtGui.QMainWindow.__init__(self)
        # set window title
        self.setWindowTitle("Matplotlib Figure in a Qt4 Window With NavigationToolbar")
        # instantiate a widget, it will be the main one
        self.main_widget = QtGui.QWidget(self)
        # create a vertical box layout widget
        vbl = QtGui.QVBoxLayout(self.main_widget)
        # instantiate our Matplotlib canvas widget
        qmc = Qt4MplCanvas(self.main_widget)
        # instantiate the navigation toolbar
        ntb = NavigationToolbar(qmc, self.main_widget)
        # pack these widget into the vertical box
        vbl.addWidget(qmc)
        vbl.addWidget(ntb)
        # set the focus on the main widget
        self.main_widget.setFocus()
        # set the central widget of MainWindow to main_widget
        self.setCentralWidget(self.main_widget)



class ApplicationWindow(QtGui.QMainWindow, gui_class):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self)

        self.setupUi(self)

        # menu bar
        self.pushButton_open_par_file.clicked.connect(self.handle_menu_open_par_file)

        self.aw = PlotWindow()
        self.aw.show()



    def handle_menu_open_par_file(self):
        fname = QtGui.QFileDialog.getOpenFileName()


        params = parameters(fname)


        #####################################################################

        #initialising data object
        dataob = data(params)

        #adding some molecules to the atmosphere
        dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
        dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)


        #initialising TP profile object
        profileob = tp_profile(dataob)

        #initialising transmission radiative transfer code object
        if params.gen_type == 'transmission':
            transob = transmission(profileob)

            if params.trans_cpp:
                MODEL = transob.cpath_integral()  # computing transmission
            else:
                MODEL = transob.path_integral()  # computing transmission

        #initialising transmission radiative transfer code object
        if params.gen_type == 'emission':
            emisob = emission(profileob)

            MODEL = emisob.path_integral()  # computing transmission

        # #
        OUT = np.zeros((len(dataob.specgrid),2))
        OUT[:,0] = dataob.specgrid
        OUT[:,1] = MODEL
        # OUT[:,2] += 1e-5 #adding errorbars. can be commented


#put the main window here
def main():
    import sys
    qApp = QtGui.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()
    sys.exit(qApp.exec_())

if __name__ == '__main__':
    main()


exit()



import sys

from PyQt4 import QtCore, QtGui, uic

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

#loading classes
sys.path.append('./classes')
sys.path.append('./library')


import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from tp_profile import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

gui_class = uic.loadUiType('gui.ui')[0]

class MatplotlibWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)


        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


    def plotDataPoints(self, x, y):
        self.axes.plot(x,y)
        self.draw()


class ApplicationWindow(QtGui.QMainWindow, gui_class):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self)

        self.setupUi(self)

        # menu bar
        self.pushButton_open_par_file.triggered.connect(self.handle_menu_open_par_file)




        self.mplwidget = MatplotlibWidget(self.mplwidget)

        #self.mplwidget = MatplotlibWidget(self.graphicsView_matplotlib)

        #self.mplwidget.setGeometry(QtCore.QRect(70, 50, 600, 500))
        #self.mplwidget.setObjectName("mplwidget")

        self.mplwidget.plotDataPoints([5, 4],[2, 3])





        #self.pushButton_plot.clicked.connect(self.handleButton)




    def handle_menu_open_par_file(self):
        fname = QtGui.QFileDialog.getOpenFileName()


        params = parameters(fname)


        #####################################################################

        #initialising data object
        dataob = data(params)

        #adding some molecules to the atmosphere
        dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
        dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)


        #initialising TP profile object
        profileob = tp_profile(dataob)

        #initialising transmission radiative transfer code object
        if params.gen_type == 'transmission':
            transob = transmission(profileob)

            if params.trans_cpp:
                MODEL = transob.cpath_integral()  # computing transmission
            else:
                MODEL = transob.path_integral()  # computing transmission

        #initialising transmission radiative transfer code object
        if params.gen_type == 'emission':
            emisob = emission(profileob)

            MODEL = emisob.path_integral()  # computing transmission

        # #
        OUT = np.zeros((len(dataob.specgrid),2))
        OUT[:,0] = dataob.specgrid
        OUT[:,1] = MODEL
        # OUT[:,2] += 1e-5 #adding errorbars. can be commented



if __name__ == '__main__':

    import sys
    app = QtGui.QApplication(sys.argv)
    window = ApplicationWindow()
    window.show()
    sys.exit(app.exec_())
