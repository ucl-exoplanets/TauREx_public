#loading libraries
import numpy, pylab, sys, os, optparse, time
from ConfigParser import SafeConfigParser

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,atmosphere,data,preselector
from parameters import *
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *


parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
                  action="store_true",
)
options, remainder = parser.parse_args()
params = parameters(options.param_filename)

if params.gen_type == 'transmission' or params.fit_transmission:
    dataob = data(params)
    atmosphereob = atmosphere(dataob)
    forwardmodelob = transmission(atmosphereob)
    fittingob = fitting(forwardmodelob)
    if params.mcmc_run and pymc_import:
        fittingob.MCMC = True
    if params.nest_run and multinest_import:
        fittingob.NEST = True
    outputob = output(fittingob)
    if params.verbose or params.out_save_plots:
        outputob.plot_all(save2pdf=params.out_save_plots)
    outputob.save_ascii_spectra()