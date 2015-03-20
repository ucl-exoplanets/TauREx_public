

#loading libraries
import numpy, pylab, sys, os, optparse, time
from numpy import * #nummerical array library
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

#checking for multinest library
global multinest_import
try:
    with os.devnull as sys.stdout:
        import pymultinest
        multinest_import = True
except:
    multinest_import = False

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

## Apply TP profile
def movingaverage(values,window):
    weigths = np.repeat(1.0, window)/window
    smas = np.convolve(values, weigths, 'valid')
    return smas #smas2[::-1] # as a numpy array


AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  help='Parameter filename')
options, remainder = parser.parse_args()

#initialising parameters object
params = parameters(options.param_filename)

# set model resolution to 1000
params.gen_spec_res = 1000
params.gen_manual_waverange = True

errors = [50, 100, 300]
resolutions = [300, 100, 50, 10]

#errors = [100]
#resolutions = [100]

out_path_orig = params.out_path

for error in errors:
    for res in resolutions:

        params.out_path = os.path.join(out_path_orig, 'R%i' % res, 'E%s' % error)
        params.gen_manual_waverange = True
        params.in_spectrum_file = False

        try:
            os.makedirs(params.out_path)
        except:
            pass

        if  os.path.isfile(os.path.join(params.out_path, 'NEST_out.db')):
            logging.info('Skipping resolution %i error %i ppm, already processed' % (res, error))
        else:

            logging.info('Processing resolution %i error %i ppm' % (res, error))


            #initialising data object
            dataob = data(params)

            #initialising TP profile object
            atmosphereob = atmosphere(dataob)

            if MPIrank == 0:

                logging.info('Creating custom TP profile')
                MAX_P = atmosphereob.P[0]
                MIN_P = atmosphereob.P[-1]
                smooth_window = 5 #smoothing window size as percent of total data
                Pnodes = [MAX_P, MAX_P, 0.65e4, MIN_P]
                Tnodes = [1000,1000, 600,600]
                TP = np.interp((np.log(atmosphereob.P[::-1])), np.log(Pnodes[::-1]), Tnodes[::-1])
                atmosphereob.T = TP[::-1]
                atmosphereob.update_atmosphere()

                # save TP profile to file
                out = np.zeros((len(atmosphereob.T),2))
                out[:,0] = atmosphereob.T
                out[:,1] = atmosphereob.P
                np.savetxt(os.path.join(params.out_path, 'TP_profile.dat'), out)

                print 'The mean molecular weight is', atmosphereob.planet_mu/AMU

                # save spectrum

                #initialising transmission radiative transfer code object
                forwardmodelob = transmission(atmosphereob)

                # bin down internal model to given resolution
                wavegrid, dlamb_grid = dataob.get_specgrid(R=res,lambda_min=params.gen_wavemin,lambda_max=params.gen_wavemax)
                spec_bin_grid, spec_bin_grid_idx = dataob.get_specbingrid(wavegrid, dataob.specgrid)

                model = forwardmodelob.model()
                model_binned = [np.mean(model[spec_bin_grid_idx == i]) for i in xrange(1,len(spec_bin_grid))]

                out = np.zeros((len(wavegrid),3))
                out[:,0] = wavegrid
                out[:,1] = model_binned
                out[:,2] += float(error) * 1e-6
                out[:,1] += np.random.normal(0, float(error) * 1e-6, len(wavegrid))

                # remove nan values. Sometimes the last bin has a nan...
                out = out[~np.isnan(out).any(axis=1)]

                #initiating output object with fitted data from fitting class
                outputob = output(forwardmodel=forwardmodelob)
                outputob.save_spectrum_to_file(spectrum=out, saveas='observed_spectrum.dat')

            params.in_spectrum_file = os.path.join(params.out_path, 'observed_spectrum.dat')

            MPI.COMM_WORLD.Barrier()

            logging.info('Start fitting...')

            params.gen_manual_waverange = False

            dataob = data(params)
            atmosphereob = atmosphere(dataob)
            if params.fit_fix_temp == True:
                atmosphereob.T = TP[::-1]
            forwardmodelob = transmission(atmosphereob)
            fittingob = fitting(forwardmodelob)
            fittingob.multinest_fit() # Nested sampling fit

            # save output
            if MPIrank == 0:
                outputob = output(fittingob)
                outputob.plot_all(save2pdf=params.out_save_plots)
                outputob.save_ascii_spectra()

        MPI.COMM_WORLD.Barrier()
