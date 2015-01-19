

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

'''
coap

'''
import parameters,emission,transmission,output,fitting,atmosphere,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

gui_class = uic.loadUiType('gui.ui')[0]

#conversion constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg


class Qt4MplCanvas(FigureCanvas):

    def __init__(self, parent):
        # plot definition
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)

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
        self.qmc = Qt4MplCanvas(self.main_widget)
        # instantiate the navigation toolbar
        ntb = NavigationToolbar(self.qmc, self.main_widget)
        # pack these widget into the vertical box
        vbl.addWidget(self.qmc)
        vbl.addWidget(ntb)
        # set the focus on the main widget
        self.main_widget.setFocus()
        # set the central widget of MainWindow to main_widget
        self.setCentralWidget(self.main_widget)



class ApplicationWindow(QtGui.QMainWindow, gui_class):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self)

        self.setupUi(self)

        self.aw = PlotWindow()
        self.aw.show()

        # initialise parameter, data and atmosphere objects
        self.params = parameters()
        self.params.planet_molec =  ['1H2-16O']#, '1H-12C-14N', '12C-1H4', '12C-16O2', '12C-16O', '14N-1H3', '28Si-16O', '48Ti-16O', '51V-16O']
        self.params.planet_mixing = [1e-4     ]#, 0           , 1e-4     , 1e-4      , 1e-4     ,  1e-4    , 0         , 0         , 0        ]
        self.params.gen_spec_res = 1000
        self.params.gen_wavemin = 0.4
        self.params.gen_wavemax = 20
        self.params.gen_manual_waverange = True

        self.dataob = data(self.params)
        self.atmosphereob = atmosphere(self.dataob)

        self.dataob = data(self.params)
        self.atmosphereob = atmosphere(self.dataob)
        self.set_params_values()

        # connect
        self.doubleSpinBox_H2O.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_CH4.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_CO2.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_NH3.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_planet_mu.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_Rp_Rstar.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_planet_surf_pressure.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_planet_surf_pressure.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_H2.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_He.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_N2.valueChanged.connect(self.event_status_changed)

        self.pushButton_plot.clicked.connect(self.plot_forwardmodel)



    def set_params_values(self, params=None):

        if not params:
            params = self.params

        self.doubleSpinBox_H2O.setValue(params.planet_mixing[0]*100)
        #self.doubleSpinBox_HCN.setValue(params.planet_mixing[1]*100)
        #self.doubleSpinBox_CH4.setValue(params.planet_mixing[2]*100)
        #self.doubleSpinBox_CO2.setValue(params.planet_mixing[3]*100)
        #self.doubleSpinBox_CO.setValue(params.planet_mixing[4]*100)
        #self.doubleSpinBox_NH3.setValue(params.planet_mixing[5]*100)
        #self.doubleSpinBox_SiO.setValue(params.planet_mixing[6]*100)
        #self.doubleSpinBox_TiO.setValue(params.planet_mixing[7]*100)
        #self.doubleSpinBox_VO.setValue(params.planet_mixing[8]*100)
        self.doubleSpinBox_H2.setValue(params.planet_H2_fraction*100.)
        self.doubleSpinBox_He.setValue(params.planet_He_fraction*100.)
        self.doubleSpinBox_N2.setValue(params.planet_N2_fraction*100.)
        self.doubleSpinBox_planet_albedo.setValue(params.planet_albedo)
        self.doubleSpinBox_planet_mu.setValue(params.planet_mu/AMU)
        self.doubleSpinBox_Rp_Rstar.setValue(params.planet_radius/params.star_radius)
        self.spinBox_planet_T.setValue(params.planet_temp)
        #self.doubleSpinBox_planet_surf_grav.setValue(params.)
        self.doubleSpinBox_planet_surf_pressure.setValue(params.tp_max_pres/1.e6)
        self.doubleSpinBox_clouds_lower.setValue(params.in_cld_pressure[0])
        self.doubleSpinBox_clouds_upper.setValue(params.in_cld_pressure[1])
        self.doubleSpinBox_star_radius.setValue(params.star_radius/RSOL)
        self.spinBox_star_T.setValue(params.star_temp)
        self.doubleSpinBox_max_wav.setValue(0.4)
        self.doubleSpinBox_min_wav.setValue(20)
        self.spinBox_resolution.setValue(1000)
        self.checkBox_rayleigh.setCheckState(1)
        self.checkBox_induced_absorption.setCheckState(params.in_include_cia)

    def event_status_changed(self):
        self.params.planet_mixing[0] = self.doubleSpinBox_H2O.value() / 100.
        #self.params.planet_mixing[1] = self.doubleSpinBox_HCN.value() / 100.
        #self.params.planet_mixing[2] = self.doubleSpinBox_CH4.value() / 100.
        #self.params.planet_mixing[3] = self.doubleSpinBox_CO2.value() / 100.
        #self.params.planet_mixing[4] = self.doubleSpinBox_CO.value() / 100.
        #self.params.planet_mixing[5] = self.doubleSpinBox_NH3.value() / 100.
        #self.params.planet_mixing[6] = self.doubleSpinBox_SiO.value() / 100.
        #self.params.planet_mixing[7] = self.doubleSpinBox_TiO.value() / 100.
        #self.params.planet_mixing[8] = self.doubleSpinBox_VO.value() / 100.

        self.params.planet_mu = self.doubleSpinBox_planet_mu.value() * AMU
        self.params.planet_radius = self.doubleSpinBox_Rp_Rstar.value()*self.doubleSpinBox_star_radius.value()*RSOL
        self.params.tp_max_pres = self.doubleSpinBox_planet_surf_pressure.value() * 1.e6

        self.params.planet_H2_fraction = self.doubleSpinBox_H2.value()/100.
        self.params.planet_He_fraction = self.doubleSpinBox_He.value()/100.
        self.params.planet_N2_fraction = self.doubleSpinBox_N2.value()/100.

        if self.checkBox_plot_autorefresh.isChecked():
            self.plot_forwardmodel(self.params)

    def plot_forwardmodel(self, params=None):

        if not params:
            params = self.params

        self.atmosphereob = atmosphere(self.dataob, params=self.params)

        if self.params.gen_type == 'transmission':
            self.forwardmodelob = transmission(self.atmosphereob, params=self.params)

        out = np.zeros((len(self.dataob.specgrid),2))
        out[:,0] = self.dataob.specgrid
        out[:,1] = self.forwardmodelob.model()

        if not self.checkBox_plot_overplot.isChecked():
            self.aw.qmc.axes.clear()
        self.aw.qmc.axes.plot(out[:,0], out[:,1], lw=0.5)
        self.aw.qmc.axes.set_xscale('log')
        self.aw.qmc.axes.set_xlim(np.min(out[:,0]), np.max(out[:,0]))
        self.aw.qmc.draw()


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
