import sys

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,fitting,atmosphere,data
from parameters import *
from emission import *
from transmission import *
from fitting import *
from atmosphere import *
from data import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

import create_spectrum
from create_spectrum import *

from scipy import interpolate
import glob

if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-p', '--parfile',
                      dest="param_filename",
                      default="Parfiles/default.par",
    )
    parser.add_option('-T', '--temp',
                      dest="flare_temp",
                      default=412,
    )
    parser.add_option('-b', '--bin',
                      dest="bin",
                      default='resolution',
    )
    parser.add_option('-r', '--res',    # binning resolution
                      dest="resolution",
                      default=1000,
    )
    parser.add_option('-d', '--dlambda', # delta lambda for binning in micron
                      dest="dlambda",
                      default=0.005,
    )
    parser.add_option('-e', '--error',
                      dest="error",
                      default=0,
    )

    options, remainder = parser.parse_args()
    params = parameters(options.param_filename)

    flaredir = 'Input/flares/files_%sK/' % options.flare_temp

    #make output dir
    atm_dir = os.path.join(params.out_path, 'atm_files')
    spectra_dir = os.path.join(params.out_path, 'spectra')

    if not os.path.isdir(atm_dir):
        os.makedirs(atm_dir)
    if not os.path.isdir(spectra_dir):
        os.makedirs(spectra_dir)

    print 'spectra dir is ', spectra_dir

    # read tp profile. Column 0: altitude (km), column 1: pressure (mbar), column 2: temperature (K)
    tp_file = os.path.join(flaredir, 'profil_PTn_ADLeob_%sK.dat' % options.flare_temp)
    tp_in = np.loadtxt(tp_file)

    tp_interp_pres_temp = interpolate.interp1d(tp_in[:,1], tp_in[:,2]) # interpolate pressure - temperature profile
    tp_interp_pres_alt = interpolate.interp1d(tp_in[:,1], tp_in[:,0]) # interpolate pressure - altitude profile

    # define planetary radius at 10 bar
    alt_readjust = tp_interp_pres_alt(1.e4)
    logging.info('Altitude readjustment: %2.1f km ' % alt_readjust)

    # loop through flare files
    for file in [os.path.join(flaredir, 'fractions_molaires_steady_ADLeo_%sK.dat' % options.flare_temp)] + glob.glob(os.path.join(flaredir, 'fractions_molaires_ADLEO_renorm_entete_*.dat')):

        filename, file_extension = os.path.splitext(os.path.basename(file))

        try:
            time = filename.split('_')[5]
        except:
            time = 0

        atm_profile_in = np.loadtxt(file, skiprows=2)
        altitude = atm_profile_in[:,0]
        pressure = atm_profile_in[:,1]
        temperature = tp_interp_pres_temp(pressure)

        # convert pressure from mbar to pascal
        pressure *= 100

        # readjust altitude to have 0 km at 10 bar
        altitude -= alt_readjust

        # read molecules mixing ratios
        NO = atm_profile_in[:,68]
        CO2 = atm_profile_in[:,78]
        OH = atm_profile_in[:,79]
        HCN = atm_profile_in[:,81]
        NH3 = atm_profile_in[:,82]
        CH4 = atm_profile_in[:,83]
        CO = atm_profile_in[:,85]
        H2O = atm_profile_in[:,86]

        # interpolate to new grid with 100 layers, from 0 km to max z
        alt_grid = np.linspace(0, np.max(altitude), 40)
        pressure_interp = np.interp(alt_grid, altitude, pressure)
        temperature_interp = np.interp(alt_grid, altitude, temperature)
        NO_interp = np.interp(alt_grid, altitude, NO)
        CO2_interp = np.interp(alt_grid, altitude, CO2)
        OH_interp = np.interp(alt_grid, altitude, OH)
        HCN_interp = np.interp(alt_grid, altitude, HCN)
        NH3_interp = np.interp(alt_grid, altitude, NH3)
        CH4_interp = np.interp(alt_grid, altitude, CH4)
        CO_interp = np.interp(alt_grid, altitude, CO)
        H2O_interp = np.interp(alt_grid, altitude, H2O)

        # sort in increasing pressure
        idx = np.argsort(pressure_interp)

        atm_output = np.column_stack((pressure_interp[idx], temperature_interp[idx], alt_grid[idx],
                                       NO_interp[idx], CO2_interp[idx], OH_interp[idx], HCN_interp[idx], NH3_interp[idx],
                                       CH4_interp[idx], CO_interp[idx], H2O_interp[idx]))

        header = 'New P-T grid, created for GJ1214b \n'
        header += ' Model Atmosphere Pressures, Temperatures, and Gas Mixing Ratios\n'
        header += ' \n'
        header += ' sounding: Atmospheric Sounding created by the program PROFILE\n'
        header += ' number of standard atmosphere levels, nlstd =        39\n'
        header += ' latitude of sounding (degrees)              =         0.0\n'
        header += ' solar longitude of sounding (degrees)       =         0.0\n'
        header += ' number of absorbing gases:                  =         8\n'
        header += ' \n'
        header += '   p(Pas)       t(k)    alt(km) '

        atm_filename = '%s.atm' % filename
        np.savetxt(os.path.join(atm_dir, atm_filename), atm_output, header=header, fmt='%.5e',)

        params.in_atm_file = os.path.join(atm_dir, atm_filename)

        #loading object
        createob = create_spectrum(options, params=params)

        #generating spectrum
        createob.generate_spectrum()

        #saving spectrum
        spectrum_filename = '%s.dat' % filename
        createob.save_spectrum(filename=os.path.join('spectra', spectrum_filename))