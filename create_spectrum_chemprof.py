#! /usr/bin/python 

'''
TauREx script to read-in and caluclate chemical profile and 
temperature-pressure files and to create a forward model. 
'''
import sys
import os
import argparse
import logging
import numpy as np 
import pylab as pl 
# loading classes
sys.path.append('../')
sys.path.append('../classes')
sys.path.append('../library')

from create_spectrum import *
from parameters import *
from transmission import *
from emission import *
from atmosphere import *
from data import *


class convert_chemprof(object):
    def __init__(self,params,chem_profile): 
        
        self.params = params
        self.chem_profile_file = chem_profile
        self.load_chem_model()

    def load_chem_model(self):

            logging.info('Load atmospheric profile from external files. Format: Taurex external models')

#             self.ven_altitude_int = interp1d(self.ven_pressure, self.ven_altitude)
#             self.ven_temperature_int = interp1d(self.ven_pressure, self.ven_temperature)

            # extract list of molecules from 'fractions_molaires' file
            with open(self.chem_profile_file, 'r') as f:
                self.ven_molecules = [x.upper() for x in f.readline().split()]
                self.ven_molweight = [float(x) for x in f.readline().split()] # molecular weight in AMU
            table = np.loadtxt(self.chem_profile_file, skiprows=2)
#             self.ven_molprof_altitude = table[:,0]/1000. # convert km to m
            self.ven_pressure = table[:,0] # in pascal
            self.ven_temperature = table[:,1] #in K
            self.ven_molprof_mixratios = table[:,2:] # mixing ratios referring to self.ven_molecules
#             self.ven_molprof_mixratios_int = [interp1d(self.ven_molprof_pressure, self.ven_molprof_mixratios[:,i])
#                                               for i in xrange(np.shape(self.ven_molprof_mixratios)[1])]
           
            logging.info('Atmospheric pressure boundaries from chemical model: %f-%f' %
                         (np.min(self.ven_pressure), np.max(self.ven_pressure)))
            
#             pl.figure()
#             pl.plot(self.ven_temperature, np.log(self.ven_pressure))
           
            # determine list of molecules with cross sections. exclude the others.
            self.ven_active_gases = []
            self.ven_active_gases_idx = []
            self.ven_noxsec = []

            for mol_idx, mol_val in enumerate(self.ven_molecules):
                
                molpath1 = os.path.join(self.params.in_xsec_path, '%s.db' % mol_val)
                molpath2 = os.path.join(self.params.in_xsec_path, '%s.pickle' % mol_val)
                
                if os.path.isfile(molpath1):
                    self.ven_active_gases.append(mol_val)
                    self.ven_active_gases_idx.append(mol_idx)
                elif os.path.isfile(molpath2):
                    self.ven_active_gases.append(mol_val)
                    self.ven_active_gases_idx.append(mol_idx)
                else:
                    self.ven_noxsec.append(mol_val)

            logging.info('Molecules with available cross sections: %s' %  self.ven_active_gases)
            logging.info('Molecules with available cross sections idx: %s' %  self.ven_active_gases_idx)
            logging.info('Excluded molecules: %s' %  self.ven_noxsec)

            # determine list of inactive gases
            gases = ['H2', 'HE', 'N2']
#             gases = ['H2','N2']
            self.ven_inactive_gases = []
            self.ven_inactive_gases_idx = []
            for gasname in gases:
                if gasname in self.ven_molecules:
                    self.ven_inactive_gases.append(gasname)
                    self.ven_inactive_gases_idx.append(self.ven_molecules.index(gasname))
            logging.info('Inactive gases: %s' %  self.ven_inactive_gases)
        
    
    def modify_params(self,params):
        '''
        modifies params object to include the right molecules when loading create_spectrum. 
        returns modified params object
        '''
        self.params = params
        params.atm_active_gases = self.ven_active_gases
        params.atm_active_gases_mixratios = [1e-10 for i in range(len(self.ven_active_gases))]
        params.atm_max_pres = np.max(self.ven_pressure)
        params.atm_min_pres = np.min(self.ven_pressure)
        params.mode == 'custom'
        return params
    
    
    def set_chem_profiles(self,csob):
        '''
        takes create_spectrum instance and modifies the atmosphere to 
        the venot chemical and temperature-pressure profiles
        '''
        pressure_profile = csob.atmosphereob.pressure_profile
        
        #setting chemical profiles 
        active, inactive = self.get_chemical_profile(pressure_profile)
        csob.atmosphereob.active_mixratio_profile = active
        csob.atmosphereob.inactive_mixratio_profile = inactive
        
        #setting temperature profiles
        csob.atmosphereob.temperature_profile = self.get_tp_profile(pressure_profile)
        
        #setting mu
        csob.atmosphereob.mu_profile = self.get_atmosphere_mu(csob)
        
        # re-compute altitude profile, scale height and planet gravity arrays
        csob.atmosphereob.set_altitude_gravity_scaleheight_profile()

        # re-set density profile
        csob.atmosphereob.set_density_profile()
        
#         
#         pl.figure()
# #         pl.plot(csob.atmosphereob.pressure_profile,csob.atmosphereob.temperature_profile)
#         pl.plot(csob.atmosphereob.temperature_profile,np.log(csob.atmosphereob.pressure_profile))
#         pl.show()
#         exit()
        
        
        logging.info('Read in chemical and thermal profiles. Modifying atmosphere accordingly...')
        logging.info('Mean molecular weight (1st layer): %.5f AMU' % (csob.atmosphereob.mu_profile[0]/AMU))
        logging.info('Scale height (1st layer): %.1f km' % (csob.atmosphereob.scaleheight_profile[0]/1000.))
        logging.info('Temperature (1st layer): %.1f K' % (csob.atmosphereob.temperature_profile[0]))
        logging.info('Atmospheric max pressure: %.3f bar' % (csob.atmosphereob.params.atm_max_pres/1e5))
        
        return csob
    
    
    def get_atmosphere_mu(self,csob):
        '''
        calculating the total atmospheric mu
        '''
        AMU = 1.660538921e-27
        pressure_profile = csob.atmosphereob.pressure_profile
        # get mu for each layer
        mu = np.zeros(self.params.atm_nlayers)
        for mol_idx, mol_val in enumerate(self.ven_molecules):
            ven_molprof_mixratio = np.interp(pressure_profile[::-1],
                                             self.ven_pressure[::-1],
                                             self.ven_molprof_mixratios[:,mol_idx][::-1])[::-1]
            mu[:] += ven_molprof_mixratio* self.ven_molweight[mol_idx]*AMU
            mu = np.asarray(mu, order='C')
        return mu
    
    def get_tp_profile(self,pressure):
        ''' 
        takes taurex initialised tp_profile and adapts new tp_profile
        ''' 
        temperature_profile = np.interp(pressure[::-1],
                                                 self.ven_pressure[::-1], self.ven_temperature[::-1])[::-1]
        temperature_profile = np.asarray(temperature_profile, order='C')
        return temperature_profile
        
        
    def get_chemical_profile(self,pressure_profile):
        '''
        takes taurex initialised chemical profile and adapts it with new tp_profile
        '''
        self.nlayers = self.params.atm_nlayers
        active_mixratio_profile = np.zeros((len(self.ven_active_gases), self.nlayers))
        inactive_mixratio_profile = np.zeros((len(self.ven_inactive_gases), self.nlayers))
        for mol_idx, mol_val in enumerate(self.ven_active_gases):
            ven_molprof_idx = self.ven_molecules.index(mol_val)
            active_mixratio_profile[mol_idx, :] =  np.interp(pressure_profile[::-1],
                                                                  self.ven_pressure[::-1],
                                                                  self.ven_molprof_mixratios\
                                                                      [:,ven_molprof_idx][::-1])[::-1]
            active_mixratio_profile[mol_idx, :] = np.asarray(active_mixratio_profile[mol_idx, :],
                                                                  order='C')
        for mol_idx, mol_val in enumerate(self.ven_inactive_gases):
            ven_molprof_idx = self.ven_molecules.index(mol_val)
            inactive_mixratio_profile[mol_idx, :] =  np.interp(pressure_profile[::-1],
                                                                    self.ven_pressure[::-1],
                                                                    self.ven_molprof_mixratios\
                                                                        [:,ven_molprof_idx][::-1])[::-1]
            inactive_mixratio_profile[mol_idx, :] = np.asarray(inactive_mixratio_profile[mol_idx, :],
                                                                    order='C')
        return active_mixratio_profile,inactive_mixratio_profile



# gets called when running from command line
if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                      dest='param_filename',
                      default='Parfiles/default.par',
                       help='Input parameter file'
                      )
    
    parser.add_argument('--chem_profile',
                        dest='chemprofile',
                        default='profile.dat',
                        help='Input chemical profile file')
    
    parser.add_argument('--pickle', # store spectrum in its pickled version (faster for highres spectra)
                      action='store_true',
                      dest='pickle_save_sp',
                      default=False,
                       help='Store the final output spectrum in Python pickle format. This is much faster for high resolution '
                            'spectra than using ascii files. See also --sp_filename.')

    parser.add_argument('--sp_filename',
                      dest='sp_filename',
                      default=False,
                      help='Specify a custom file path and filename for the output spectrum (note that the path is relative'
                           'to the current working folder). Remember that if using --pickle_save_sp the ouput spectrum '
                           'is stored in Python pickle format.')

    parser.add_argument('--save_instance',
                      action='store_true',
                      dest='save_instance',
                      default=False,
                      help = 'Save a dictionary in .pickle format containing the full spectrum instance, including mixing ratio and ' \
                             'temperature profiles used, and all the paramaters used to generate the spectrum.'
                             'See also the options --opacity__contriv, --contrib_func. Default file path is in the'
                             'Output folder specified in parameter file, and the default filename is '
                             'SPECTRUM_INSTANCE_out.pickle. Otherwise see --instance_filename')

    parser.add_argument('--instance_filename',
                      dest='instance_filename',
                      default=False,
                      help = 'Specify a custom file path and filename for the output spectrum instance dictionary (stored'
                             'in Python pickle format.')

    parser.add_argument('--plot',
                      dest='plot_spectrum',
                      action='store_true',
                      default=False,
                      help='Display an instantanous plot after the spectrum is computed.')

    parser.add_argument('--opacity_contrib',  # calculates the opacity contribution for each opacity source.
                      dest='opacity_contrib', # stored in SPECTRUM_INSTANCE_out.pickle, so use '--save_instance' as well
                      action='store_true',
                      default=False,
                      help = 'Calculates the opacity contribution for each opacity source. It computes one spectrum for '
                             'each opacity source, suppressing the contribution from all the other sources. Stored in the instance'
                             'dictionary (see option --save_instance) under ["data"]["opacity_contrib"][<opacity_source>]')

    parser.add_argument('--contrib_func',
                      dest='contrib_func',
                      action='store_true',
                      default=False,
                      help = 'Only valid for emission. Store the contribution function as a funciton of pressure and wavelength.'
                             'The 2d array is stored in the instance dictionary (see option --save_instance), '
                             'under ["data"]["contrib_func"].',)

    parser.add_argument('--transmittance',
                      dest='transmittance',
                      action='store_true',
                      default=False,
                      help = 'Only valid for transmission. Store the spectral transmittance as a function of pressure'
                             ' as a funciton of pressure and wavelength. The transmittance is integrated over the path '
                             'parallel to the line of sight. The 2d array is stored in the instance dictionary (see option --save_instance), '
                             'under ["data"]["transmittance"].',)

    parser.add_argument('--nthreads',  # run forward model in multithreaded mode (use Python multiprocesing for sigma array interpolation
                       dest='nthreads', # and openmp parallel version of cpp code). You need to spcify the number of cores to use,
                       default=0,
                       type=int,
                       help = 'Run forward model in multithreaded mode (using NTHREADS cores). NTHREADS should not '
                              'be larger than the number of cores available.')
    # add command line interface to parameter file

    params_tmp = parameters(mpi=False, log=False)
    params_dict = params_tmp.params_to_dict() # get all param names
    for param in params_dict:
        if type(params_dict[param]) == list:
            # manage lists
            parser.add_argument('--%s' % param,
                                action='append',
                                dest=param,
                                default=None,
                                type = type(params_dict[param][0])
                                )
        else:
            parser.add_argument('--%s' % param,
                                dest=param,
                                default=None,
                                type = type(params_dict[param])
                                )
    options = parser.parse_args()

    # Initialise parameters instance
    params = parameters(options.param_filename, mode='forward_model', mpi=False)

    # Override params object from command line input
    for param in params_dict:
        if getattr(options, param) != None:

            value = getattr(options, param)
            if param == 'planet_mass':
                value *= MJUP
            if param == 'planet_radius':
                value *= RJUP
            if param == 'star_radius':
                value *= RSOL
            if param == 'atm_mu':
                value *= AMU
            setattr(params, param, value)

    # checks
    if params.gen_type == 'transmission' and options.contrib_func:
        logging.error('Options --contrib_func is only valid in emission. Maybe you wanted to use --transmittance? ')

    if params.gen_type == 'emission' and options.transmittance:
        logging.error('Options --transmittance is only valid in transmission. Maybe you wanted to use --contrib_func ?')

    if (options.transmittance or options.contrib_func) and not options.save_instance:
        logging.warning('Options --transmittance and --contrib_func require --save_instance. This options is '
                        'switched on automatically. The instance of the spectrum will be stored in a .pickle file'
                        'in the Output folder.')
        options.save_instance = True
        
        
    #initialising convert_venot instance
    venotob = convert_chemprof(params=params,chem_profile=options.chemprofile)
    
    #partially initialising standard create_spectrum instance 
    spectrumob = create_spectrum(params=venotob.modify_params(params), nthreads=options.nthreads,full_init=False)
#     spectrumob = create_spectrum(params=params, nthreads=options.nthreads,full_init=False)
    spectrumob.init_data() #initialising data
    spectrumob.init_atmosphere(nthreads=options.nthreads) #initialising atmosphere

    #modifying taurex atmosphere with venot chemical profile and temperature-pressure profile
    spectrumob = venotob.set_chem_profiles(spectrumob)
    spectrumob.init_fmob() #initialising forward model
    
    #generating forward model spectrum
    spectrumob.generate_spectrum(save_instance=options.save_instance,
                                 instance_filename=options.instance_filename,
                                 opacity_contrib=options.opacity_contrib,
                                 contrib_func=options.contrib_func,
                                 transmittance=options.transmittance)
    
    #saving spectrum
    spectrumob.save_spectrum(sp_filename=options.sp_filename, pickled=options.pickle_save_sp)

    #plotting spectrum
    if options.plot_spectrum:
        logging.info('Plot spectrum... Close the Plot window to terminate the script.')
        sp = spectrumob.SPECTRUM_INSTANCE_out['data']['spectrum']
        plt.plot(sp[:,0], sp[:,1])
        plt.xscale('log')
        plt.show()
    
    
    