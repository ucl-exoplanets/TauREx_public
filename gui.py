

import sys

from PyQt4 import QtCore, QtGui, uic
import numpy as np

import matplotlib
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib as mpl
from matplotlib import rc
matplotlib.use('Qt4Agg')
mpl.rcParams['axes.linewidth'] = 2 #set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 14})

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

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
        self.fig = Figure(facecolor='white')
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
        self.setWindowTitle("Spectrum")
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


        # initialise plot window
        self.aw = PlotWindow()
        self.aw.show()

        self.observations = False

        # initialise parameter, data, atmosphere and forwardmodel objects
        logging.info('Running TauREX GUI')
        self.params = self.load_parameter_object()

        self.dataob = data(self.params)
        self.atmosphereob = atmosphere(self.dataob)
        self.forwardmodel = transmission(self.atmosphereob)
        self.set_params_values()
        self.nplot = 1
        self.cm = plt.get_cmap('Accent')
        self.plot_forwardmodel()

        # connect
        self.connect_spinboxes()
        self.pushButton_load_par_file.clicked.connect(self.load_par_file)
        self.pushButton_save_plot.clicked.connect(self.save_plot)
        self.pushButton_save_ascii.clicked.connect(self.save_ascii)
        self.pushButton_add_observations.clicked.connect(self.add_observations)
        self.pushButton_plot.clicked.connect(self.plot_forwardmodel)

    def connect_spinboxes(self):

        self.doubleSpinBox_H2O.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_HCN.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_CH4.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_CO2.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_CO.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_NH3.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_SiO.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_TiO.valueChanged.connect(self.event_status_changed)
        #self.doubleSpinBox_VO.valueChanged.connect(self.event_status_changed)
        self.doubleSpinBox_NO.valueChanged.connect(self.event_status_changed)

        if not self.params.fit_couple_mu:
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
        self.checkBox_rayleigh.stateChanged.connect(self.event_status_changed)

    def diconnect_spinboxes(self):

        self.doubleSpinBox_H2O.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_HCN.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_CH4.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_CO2.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_CO.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_NH3.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_SiO.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_TiO.valueChanged.disconnect(self.event_status_changed)
        #self.doubleSpinBox_VO.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_NO.valueChanged.disconnect(self.event_status_changed)
        if not self.params.fit_couple_mu:
            self.doubleSpinBox_planet_mu.valueChanged.disconnect(self.event_status_changed)

        self.doubleSpinBox_Rp_Rstar.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_planet_surf_pressure.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_planet_mass.valueChanged.disconnect(self.event_status_changed)
        self.spinBox_planet_T.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_H2.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_He.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_N2.valueChanged.disconnect(self.event_status_changed)
        self.spinBox_star_T.valueChanged.disconnect(self.event_status_changed)
        self.doubleSpinBox_star_radius.valueChanged.disconnect(self.event_status_changed)
        self.checkBox_rayleigh.stateChanged.disconnect(self.event_status_changed)

    def load_par_file(self):

        filename = QtGui.QFileDialog.getOpenFileName(self, 'Select parameter file', 'Parfiles/')
        if filename:
            self.params = self.load_parameter_object(filename)
            self.dataob = data(self.params)
            self.atmosphereob = atmosphere(self.dataob)
            self.forwardmodel = transmission(self.atmosphereob)
            self.diconnect_spinboxes()
            self.set_params_values()
            self.connect_spinboxes()
            self.plot_forwardmodel()

    def save_plot(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save plot', 'Output/')
        self.aw.qmc.fig.savefig(str(filename), dpi=80,  bbox_inches='tight')

    def save_ascii(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save to ascii', 'Output/')

        out = np.zeros((len(self.dataob.specgrid),2))
        out[:,0] = self.dataob.specgrid
        out[:,1] = self.forwardmodel.model()
        np.savetxt(str(filename), out)


        self.aw.qmc.fig.savefig(str(filename), dpi=80,  bbox_inches='tight')

    def load_parameter_object(self, filename=None):

        params = parameters(filename)
        params.gen_run_gui = True

        # absorbing gases

        # All molecules are pre-loaded. Set mixing ratio if specified in the .par file being loaded
        # otherwise assume 0
        planet_mixing_tmp = []
        for idx, gasname in enumerate(params.all_absorbing_gases): # list all avaialble molecules
            try:
                # the molecule appears in the .par file, load the mixing ratio
                idx_from_params = params.planet_molec.index(gasname)
                planet_mixing_tmp.append(params.planet_mixing[idx_from_params])
            except:
                # the molecule doesn't appear in the .par file, assume X = 0
                planet_mixing_tmp.append(0)
        params.planet_mixing = planet_mixing_tmp
        params.planet_molec =  params.all_absorbing_gases

        # do the same for inactive gases
        planet_inactive_gases_X_tmp = []
        for idx, gasname in enumerate(params.all_inactive_gases): # list all avaialble molecules
            try:
                # the molecule appears in the .par file, load the mixing ratio
                idx_from_params = params.planet_inactive_gases.index(gasname)
                planet_inactive_gases_X_tmp.append(params.planet_inactive_gases_X[idx_from_params])
            except:
                # the molecule doesn't appear in the .par file, assume X = 0
                planet_inactive_gases_X_tmp.append(0)

        params.planet_inactive_gases_X = planet_inactive_gases_X_tmp
        params.planet_inactive_gases =  params.all_inactive_gases

        params.gen_spec_res = 500
        params.gen_wavemin = 0.4
        params.gen_wavemax = 20.0
        params.gen_manual_waverange = True

        return params

    def set_params_values(self):

        self.doubleSpinBox_H2O.setValue(self.forwardmodel.params.planet_mixing[0]*100.)
        self.doubleSpinBox_HCN.setValue(self.forwardmodel.params.planet_mixing[1]*100.)
        self.doubleSpinBox_CH4.setValue(self.forwardmodel.params.planet_mixing[2]*100.)
        self.doubleSpinBox_CO2.setValue(self.forwardmodel.params.planet_mixing[3]*100.)
        self.doubleSpinBox_CO.setValue(self.forwardmodel.params.planet_mixing[4]*100.)
        self.doubleSpinBox_NH3.setValue(self.forwardmodel.params.planet_mixing[5]*100.)
        self.doubleSpinBox_SiO.setValue(self.forwardmodel.params.planet_mixing[6]*100.)
        self.doubleSpinBox_TiO.setValue(self.forwardmodel.params.planet_mixing[7]*100.)
        #self.doubleSpinBox_VO.setValue(self.forwardmodel.params.planet_mixing[8]*100.)
        self.doubleSpinBox_NO.setValue(self.forwardmodel.params.planet_mixing[8]*100.)
        self.doubleSpinBox_He.setValue(self.forwardmodel.atmosphere.inactive_gases_X[0]*100.)
        self.doubleSpinBox_H2.setValue(self.forwardmodel.atmosphere.inactive_gases_X[1]*100.)
        self.doubleSpinBox_N2.setValue(self.forwardmodel.atmosphere.inactive_gases_X[2]*100.)
        self.doubleSpinBox_planet_albedo.setValue(self.forwardmodel.params.planet_albedo)

        if not self.params.fit_couple_mu:
            self.doubleSpinBox_planet_mu.setValue(self.forwardmodel.params.planet_mu/AMU)

        self.doubleSpinBox_Rp_Rstar.setValue(self.forwardmodel.params.planet_radius/self.forwardmodel.params.star_radius)
        self.spinBox_planet_T.setValue(self.forwardmodel.params.planet_temp)
        self.doubleSpinBox_planet_mass.setValue(self.forwardmodel.params.planet_mass/MJUP)
        self.doubleSpinBox_planet_surf_pressure.setValue(self.forwardmodel.params.tp_max_pres/1.e6)
        self.doubleSpinBox_clouds_lower.setValue(self.forwardmodel.params.in_cld_pressure[0])
        self.doubleSpinBox_clouds_upper.setValue(self.forwardmodel.params.in_cld_pressure[1])
        self.doubleSpinBox_star_radius.setValue(self.forwardmodel.params.star_radius/RSOL)
        self.spinBox_star_T.setValue(self.forwardmodel.params.star_temp)
        self.doubleSpinBox_max_wav.setValue(self.params.gen_wavemax)
        self.doubleSpinBox_min_wav.setValue(self.params.gen_wavemin)
        self.spinBox_resolution.setValue(self.params.gen_spec_res)

        self.checkBox_induced_absorption.setCheckState(self.forwardmodel.params.in_include_cia)

    def event_status_changed(self):

        self.forwardmodel.atmosphere.absorbing_gases_X[0] =  self.doubleSpinBox_H2O.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[1] = self.doubleSpinBox_HCN.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[2] = self.doubleSpinBox_CH4.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[3] = self.doubleSpinBox_CO2.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[4] = self.doubleSpinBox_CO.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[5] = self.doubleSpinBox_NH3.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[6] = self.doubleSpinBox_SiO.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[7] = self.doubleSpinBox_TiO.value() / 100.
        #self.forwardmodel.atmosphere.absorbing_gases_X[8] = self.doubleSpinBox_VO.value() / 100.
        self.forwardmodel.atmosphere.absorbing_gases_X[8] = self.doubleSpinBox_NO.value() / 100.
        self.forwardmodel.atmosphere.X = self.forwardmodel.atmosphere.set_mixing_ratios()
        self.forwardmodel.atmosphere.inactive_gases_X[0] = self.doubleSpinBox_He.value() / 100.
        self.forwardmodel.atmosphere.inactive_gases_X[1] = self.doubleSpinBox_H2.value() / 100.
        self.forwardmodel.atmosphere.inactive_gases_X[2] = self.doubleSpinBox_N2.value() / 100.

        if self.params.fit_couple_mu:
            # get couple mu
            self.forwardmodel.atmosphere.planet_mu = self.forwardmodel.atmosphere.get_coupled_planet_mu()
            self.doubleSpinBox_planet_mu.setValue(self.forwardmodel.atmosphere.planet_mu/AMU)
        else:
            self.forwardmodel.atmosphere.planet_mu = self.doubleSpinBox_planet_mu.value() * AMU

        self.forwardmodel.atmosphere.planet_temp = self.spinBox_planet_T.value()
        self.forwardmodel.atmosphere.planet_mass = self.doubleSpinBox_planet_mass.value() * MJUP

        self.forwardmodel.atmosphere.planet_radius = self.doubleSpinBox_Rp_Rstar.value()*self.doubleSpinBox_star_radius.value()*RSOL
        self.forwardmodel.atmosphere.planet_grav = self.forwardmodel.atmosphere.get_surface_gravity()
        self.forwardmodel.atmosphere.scaleheight = self.forwardmodel.atmosphere.get_scaleheight()
        self.forwardmodel.atmosphere.max_pressure = self.doubleSpinBox_planet_surf_pressure.value() * 1.e6
        self.forwardmodel.atmosphere.pta = self.forwardmodel.atmosphere.setup_pta_grid()
        self.forwardmodel.atmosphere.P = self.forwardmodel.atmosphere.pta[:,0] # pressure array
        self.forwardmodel.atmosphere.P_bar = self.forwardmodel.atmosphere.P * 1.0e-5 #convert pressure from Pa to bar
        self.forwardmodel.atmosphere.T = self.forwardmodel.atmosphere.pta[:,1] # temperature array
        self.forwardmodel.atmosphere.z = self.forwardmodel.atmosphere.pta[:,2] # altitude array
        self.forwardmodel.atmosphere.rho = self.forwardmodel.atmosphere.get_rho()
        self.forwardmodel.atmosphere.rho = self.forwardmodel.atmosphere.get_rho()
        self.forwardmodel.params.star_radius = self.doubleSpinBox_star_radius.value() * RSOL
        self.forwardmodel.params.in_include_Rayleigh = self.checkBox_rayleigh.isChecked()

        if self.checkBox_rayleigh.isChecked():
            self.forwardmodel.Rsig = self.forwardmodel.get_Rsig()
        else:
            self.forwardmodel.Rsig = zeros((self.forwardmodel.nlambda))


        if self.checkBox_plot_autorefresh.isChecked():
            self.plot_forwardmodel()

    def add_observations(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Select observations', 'Input/observations')
        self.observations = np.loadtxt(open(filename, 'rb'))
        self.plot_observations()
        self.aw.qmc.draw()

    def plot_forwardmodel(self):


        R = self.spinBox_resolution.value()

        # bin down internal model to given resolution (default = 1000)
        wavegrid, dlamb_grid = self.dataob.get_specgrid(R=int(R),lambda_min=self.params.gen_wavemin,lambda_max=self.params.gen_wavemax)
        spec_bin_grid, spec_bin_grid_idx = self.dataob.get_specbingrid(wavegrid, self.dataob.specgrid)
        model = self.forwardmodel.model()
        model_binned = [model[spec_bin_grid_idx == i].mean() for i in xrange(1,len(spec_bin_grid))]

        out = np.zeros((len(wavegrid),3))
        out[:,0] = wavegrid
        out[:,1] = model_binned
        #out[:,2] += 1e-5 #adding errorbars. can be commented

        # out = np.zeros((len(self.dataob.specgrid),2))
        # out[:,0] = self.dataob.specgrid
        # out[:,1] = self.forwardmodel.model()


        if not self.checkBox_plot_overplot.isChecked():
            self.aw.qmc.axes.clear()
            self.nplot = 1

        color = str(self.lineEdit_color.text())

        self.plot_observations()
        self.aw.qmc.axes.plot(out[:,0], out[:,1]*1.e6, lw=2, color=str(self.lineEdit_color.text()))
        self.aw.qmc.axes.set_xscale('log')
        self.aw.qmc.axes.set_xlabel('Wavelength (micron)')
        self.aw.qmc.axes.set_ylabel('Transit depth (ppm)')
        self.aw.qmc.axes.get_xaxis().set_major_formatter(FuncFormatter(tick_formatter))
        self.aw.qmc.axes.get_xaxis().set_minor_formatter(FuncFormatter(tick_formatter))
        self.aw.qmc.axes.set_xlim(np.min(self.doubleSpinBox_min_wav.value()), np.max(self.doubleSpinBox_max_wav.value()))
        #self.aw.qmc.axes.set_ylim(np.min(out[:,1]*1.e6)-100, np.max(out[:,1]*1.e6)+100)
        #self.aw.qmc.axes.set_ylim(self.spinBox_min_y.value(), self.spinBox_max_y.value())

        set_backgroundcolor(self.aw.qmc.axes, 'white')
        self.aw.qmc.draw()
        self.nplot += 1

    def plot_observations(self):

        if isinstance(self.observations, (np.ndarray, np.generic)):
            if np.shape(self.observations)[1] == 3:
                self.aw.qmc.axes.errorbar(self.observations[:,0], self.observations[:,1]*1.e6, yerr=self.observations[:,2]*1.e6,)
            elif np.shape(self.observations)[1] == 2:
                self.aw.qmc.axes.plot(self.observations[:,0], self.observations[:,1]*1.e6,)

def tick_formatter(x, p):

    if x < 1.0:
        return "%.1f" % x
    if x >= 1.0:
        return "%i" % x

# two useful functions to set background and foreground colours

def set_foregroundcolor(ax, color):
     '''For the specified axes, sets the color of the frame, major ticks,
         tick labels, axis labels, title and legend
     '''
     for tl in ax.get_xticklines() + ax.get_yticklines():
         tl.set_color(color)
     for spine in ax.spines:
         ax.spines[spine].set_edgecolor(color)
     for tick in ax.xaxis.get_major_ticks():
         tick.label1.set_color(color)
     for tick in ax.yaxis.get_major_ticks():
         tick.label1.set_color(color)
     ax.axes.xaxis.label.set_color(color)
     ax.axes.yaxis.label.set_color(color)
     ax.axes.xaxis.get_offset_text().set_color(color)
     ax.axes.yaxis.get_offset_text().set_color(color)
     ax.axes.title.set_color(color)
     lh = ax.get_legend()
     if lh != None:
         lh.get_title().set_color(color)
         lh.legendPatch.set_edgecolor('none')
         labels = lh.get_texts()
         for lab in labels:
             lab.set_color(color)
     for tl in ax.get_xticklabels():
         tl.set_color(color)
     for tl in ax.get_yticklabels():
         tl.set_color(color)


def set_backgroundcolor(ax, color):
     '''Sets the background color of the current axes (and legend).
         Use 'None' (with quotes) for transparent. To get transparent
         background on saved figures, use:
         pp.savefig("fig1.svg", transparent=True)
     '''
     ax.patch.set_facecolor(color)
     lh = ax.get_legend()
     if lh != None:
         lh.legendPatch.set_facecolor(color)


#put the main window here
def main():
    import sys
    qApp = QtGui.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()
    sys.exit(qApp.exec_())

if __name__ == '__main__':
    main()
