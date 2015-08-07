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

AMU   = 1.660538921e-27 #atomic mass to kg

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

    logging.info('Spectra dir is %s' % spectra_dir)

    # set some fixed params
    params.tp_couple_mu = False
    params.planet_molec = ['14N-16O', '12C-16O2', '16O-1H', '1H-12C-14N' ,'14N-1H3', '12C-1H4', '12C-16O', '1H2-16O']
    params.planet_mixing = [1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5]
    params.in_use_ATMfile = True

    # define number of layers
    nlayers = 100

    # read tp profile. Column 0: altitude (km), column 1: pressure (mbar), column 2: temperature (K)
    tp_file = os.path.join(flaredir, 'profil_PTn_ADLeob_%sK.dat' % options.flare_temp)
    tp_in = np.loadtxt(tp_file)

    tp_interp_pres_temp = interpolate.interp1d(tp_in[:,1], tp_in[:,2]) # interpolate pressure - temperature profile
    tp_interp_pres_alt = interpolate.interp1d(tp_in[:,1], tp_in[:,0]) # interpolate pressure - altitude profile

    # define planetary radius at 10 bar
    alt_readjust = tp_interp_pres_alt(1.e4)
    logging.info('Altitude readjustment: %2.1f km ' % alt_readjust)


    # loop through flare files
    for file in [os.path.join(flaredir, 'fractions_molaires_retour_ADLeo_%sK.dat' % options.flare_temp)] + [os.path.join(flaredir, 'fractions_molaires_steady_ADLeo_%sK.dat' % options.flare_temp)] + glob.glob(os.path.join(flaredir, 'fractions_molaires_ADLEO_renorm_entete_*.dat')):
    #+ glob.glob(os.path.join(flaredir, 'fractions_molaires_ADLEO_renorm_entete_*.dat')):
    #for file in [os.path.join(flaredir, 'fractions_molaires_retour_ADLeo_%sK.dat' % options.flare_temp)]:

        logging.info('Processing %s' % file)
        filename, file_extension = os.path.splitext(os.path.basename(file))

        try:
            time = filename.split('_')[5]
        except:
            time = 0

        # read list of molecules + mass
        with open(file, 'r') as f:
            molecules = f.readline().split()
            weight = f.readline().split()
        f.close()

        atm_profile_in = np.loadtxt(file, skiprows=2)
        altitude = atm_profile_in[:,0]
        pressure = atm_profile_in[:,1]
        temperature = tp_interp_pres_temp(pressure)

        # convert pressure from mbar to pascal
        pressure *= nlayers

        # readjust altitude to have 0 km at 10 bar
        altitude -= alt_readjust

        # interpolate to new grid with nlayers layers
        atm_profile_interp = np.zeros((nlayers, len(molecules)))
        alt_grid = np.linspace(0, np.max(altitude), nlayers)

        pressure_interp = np.interp(alt_grid, altitude, pressure)
        temperature_interp = np.interp(alt_grid, altitude, temperature)

        # sort by increasing pressure
        idx = np.argsort(pressure_interp)
        pressure_interp = pressure_interp[idx]
        temperature_interp = temperature_interp[idx]
        alt_grid = alt_grid[idx]

        for j in xrange(len(molecules)):
            atm_profile_interp[:,j] = np.interp(alt_grid, altitude, atm_profile_in[:,j+2]) # also sorts by idx

        # reset mean molecular weight for each layer taking into account all 105 molecules
        mu = np.zeros(nlayers)
        for i in range(nlayers):
            for j in range(len(molecules)):
                mu[i] += atm_profile_interp[i,j]*float(weight[j])*AMU
            logging.debug('The mean molecular weight for layer %i is %.4f' % (i, mu[i]/AMU))
        mu = mu[::-1] # reverse array

        atm_output = np.column_stack((pressure_interp, temperature_interp, alt_grid,
                                      atm_profile_interp[:,molecules.index('NO')],
                                      atm_profile_interp[:,molecules.index('CO2')],
                                      atm_profile_interp[:,molecules.index('OH')],
                                      atm_profile_interp[:,molecules.index('HCN')],
                                      atm_profile_interp[:,molecules.index('NH3')],
                                      atm_profile_interp[:,molecules.index('CH4')],
                                      atm_profile_interp[:,molecules.index('CO')],
                                      atm_profile_interp[:,molecules.index('H2O')]))

        header = 'P-T grid - Venot Flare. Generated on %s \n'
        header += ' Model Atmosphere Pressures, Temperatures, and Gas Mixing Ratios\n'
        header += ' \n'
        header += ' sounding: Atmospheric Sounding created by the program PROFILE\n'
        header += ' number of standard atmosphere levels, nlstd =        39\n'
        header += ' latitude of sounding (degrees)              =         0.0\n'
        header += ' solar longitude of sounding (degrees)       =         0.0\n'
        header += ' number of absorbing gases:                  =         8\n'
        header += ' \n'
        header += '   p(Pas)       t(k)    alt(km)    NO     CO2     OH     HCN     NH3     CH4     CO     H2O'



        atm_filename = '%s.atm' % filename
        np.savetxt(os.path.join(atm_dir, atm_filename), atm_output, header=header, fmt='%.5e',)

        params.in_atm_file = os.path.join(atm_dir, atm_filename)

        #loading object
        createob = create_spectrum(options, params=params)

        # set mean molecular weight for each level, based on all 105 molecules

        createob.atmosphereob.planet_mu =  mu
        # if options.flare_temp == 1303:
        #     createob.atmosphereob.planet_mu = 2.748*AMU
        # elif options.flare_temp == 412:
        #     createob.atmosphereob.planet_mu = 3.15*AMU

        #generating spectrum
        createob.generate_spectrum()

        z = createob.atmosphereob.z

        # save ap profile
        ap = np.zeros((nlayers, 2))
        ap[:,0] = z
        ap[:,1] = createob.atmosphereob.P
        np.savetxt(os.path.join(params.out_path, 'ap', '%s.dat' % filename), ap)

        # save scale height profile
        sh = np.zeros((nlayers, 2))
        sh[:,0] = z
        sh[:,1] = createob.atmosphereob.scaleheight
        np.savetxt(os.path.join(params.out_path, 'scaleheight', '%s.dat' % filename), sh)

        # gravity profile
        grav = np.zeros((nlayers, 2))
        grav[:,0] = z
        grav[:,1] = createob.atmosphereob.planet_grav
        np.savetxt(os.path.join(params.out_path, 'gravity', '%s.dat' % filename), grav)

        # T profile
        temp = np.zeros((nlayers, 2))
        temp[:,0] = z
        temp[:,1] = createob.atmosphereob.T
        np.savetxt(os.path.join(params.out_path, 'temp', '%s.dat' % filename), temp)

        # density profile
        density = np.zeros((nlayers, 2))
        density[:,0] = z
        density[:,1] = createob.atmosphereob.rho
        np.savetxt(os.path.join(params.out_path, 'density', '%s.dat' % filename), density)

        # mu profile
        mup = np.zeros((nlayers, 2))
        mup[:,0] = z
        mup[:,1] = createob.atmosphereob.planet_mu
        np.savetxt(os.path.join(params.out_path, 'mu', '%s.dat' % filename), mup)

        # mol profile
        mol = np.zeros((nlayers, len(params.planet_molec)+1))
        mol[:,0] = z
        for idx, molecule in enumerate(params.planet_molec):
            mol[:,idx+1] = createob.atmosphereob.X[idx,:]

        np.savetxt(os.path.join(params.out_path, 'molecules', '%s.dat' % filename), mol)

        #saving spectrum
        spectrum_filename = '%s.dat' % filename
        createob.save_spectrum(filename=os.path.join('spectra', spectrum_filename))
