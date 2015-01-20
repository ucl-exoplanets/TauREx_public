

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
        self.params.planet_molec =  ['1H2-16O', '1H-12C-14N', '12C-1H4', '12C-16O2', '12C-16O', '14N-1H3', '28Si-16O', '48Ti-16O', '51V-16O']
        self.params.planet_mixing = [1e-4     , 0           , 1e-4     , 1e-4      , 1e-4     ,  1e-4    , 0         , 0         , 0        ]
        self.params.gen_spec_res = 500
        self.params.gen_wavemin = 0.4
        self.params.gen_wavemax = 20

        self.params.gen_manual_waverange = True

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
        self.doubleSpinBox_planet_mass.valueChanged.connect(self.event_status_changed)
        self.spinBox_planet_T.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_H2.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_He.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_N2.valueChanged.connect(self.event_status_changed)
        self.spinBox_star_T.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_star_radius.valueChanged.connect(self.event_status_changed)

        self.pushButton_plot.clicked.connect(self.plot_forwardmodel)



    def set_params_values(self, params=None):

        if not params:
            params = self.params

        self.doubleSpinBox_H2O.setValue(params.planet_mixing[0]*100.)
        self.doubleSpinBox_HCN.setValue(params.planet_mixing[1]*100.)
        self.doubleSpinBox_CH4.setValue(params.planet_mixing[2]*100.)
        self.doubleSpinBox_CO2.setValue(params.planet_mixing[3]*100.)
        self.doubleSpinBox_CO.setValue(params.planet_mixing[4]*100.)
        self.doubleSpinBox_NH3.setValue(params.planet_mixing[5]*100.)
        self.doubleSpinBox_SiO.setValue(params.planet_mixing[6]*100.)
        self.doubleSpinBox_TiO.setValue(params.planet_mixing[7]*100.)
        self.doubleSpinBox_VO.setValue(params.planet_mixing[8]*100.)
        self.doubleSpinBox_H2.setValue(params.planet_H2_fraction*100.)
        self.doubleSpinBox_He.setValue(params.planet_He_fraction*100.)
        self.doubleSpinBox_N2.setValue(params.planet_N2_fraction*100.)
        self.doubleSpinBox_planet_albedo.setValue(params.planet_albedo)
        self.doubleSpinBox_planet_mu.setValue(params.planet_mu/AMU)
        self.doubleSpinBox_Rp_Rstar.setValue(params.planet_radius/params.star_radius)
        self.spinBox_planet_T.setValue(params.planet_temp)

        self.doubleSpinBox_planet_mass.setValue(params.planet_mass/MJUP)

        self.doubleSpinBox_planet_surf_pressure.setValue(params.tp_max_pres/1.e6)
        self.doubleSpinBox_clouds_lower.setValue(params.in_cld_pressure[0])
        self.doubleSpinBox_clouds_upper.setValue(params.in_cld_pressure[1])
        self.doubleSpinBox_star_radius.setValue(params.star_radius/RSOL)
        self.spinBox_star_T.setValue(params.star_temp)
        self.doubleSpinBox_max_wav.setValue(self.params.gen_wavemin)
        self.doubleSpinBox_min_wav.setValue(self.params.gen_wavemax)
        self.spinBox_resolution.setValue(self.params.gen_spec_res)
        self.checkBox_rayleigh.setCheckState(1)
        self.checkBox_induced_absorption.setCheckState(params.in_include_cia)

    def event_status_changed(self):
        self.params.planet_mixing[0] = self.doubleSpinBox_H2O.value() / 100.
        self.params.planet_mixing[1] = self.doubleSpinBox_HCN.value() / 100.
        self.params.planet_mixing[2] = self.doubleSpinBox_CH4.value() / 100.
        self.params.planet_mixing[3] = self.doubleSpinBox_CO2.value() / 100.
        self.params.planet_mixing[4] = self.doubleSpinBox_CO.value() / 100.
        self.params.planet_mixing[5] = self.doubleSpinBox_NH3.value() / 100.
        self.params.planet_mixing[6] = self.doubleSpinBox_SiO.value() / 100.
        self.params.planet_mixing[7] = self.doubleSpinBox_TiO.value() / 100.
        self.params.planet_mixing[8] = self.doubleSpinBox_VO.value() / 100.

        self.params.planet_mu = self.doubleSpinBox_planet_mu.value() * AMU
        self.params.planet_temp = self.spinBox_planet_T.value()
        self.params.planet_radius = self.doubleSpinBox_Rp_Rstar.value()*self.doubleSpinBox_star_radius.value()*RSOL
        self.params.tp_max_pres = self.doubleSpinBox_planet_surf_pressure.value() * 1.e6

        self.params.planet_H2_fraction = self.doubleSpinBox_H2.value()/100.
        self.params.planet_He_fraction = self.doubleSpinBox_He.value()/100.
        self.params.planet_N2_fraction = self.doubleSpinBox_N2.value()/100.

        self.params.star_radius = self.doubleSpinBox_star_radius.value() * RSOL
        self.params.planet_mass = self.doubleSpinBox_planet_mass.value() * MJUP

        print self.doubleSpinBox_planet_mass.value() * MJUP


        if self.checkBox_plot_autorefresh.isChecked():
            self.plot_forwardmodel(self.params)

    def plot_forwardmodel(self, params=None):

        if not params:
            params = self.params

        self.atmosphereob = atmosphere(self.dataob, params=self.params)
        print 'Here the mass is ', self.params.planet_mass
        if self.params.gen_type == 'transmission':
            self.forwardmodelob = transmission(self.atmosphereob, params=self.params)

        print 'Rayeigh water', self.forwardmodelob.scatterRayleigh(0.6, 'He')

        out = np.zeros((len(self.dataob.specgrid),2))
        out[:,0] = self.dataob.specgrid
        out[:,1] = self.forwardmodelob.model()

        if not self.checkBox_plot_overplot.isChecked():
            self.aw.qmc.axes.clear()
        self.aw.qmc.axes.plot(out[:,0], out[:,1]*1.e6, lw=0.5)
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
