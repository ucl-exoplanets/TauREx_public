'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    TauREx create spectrum

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

# loading libraries
import  sys
import os
import argparse
import logging

# loading classes
sys.path.append('./classes')
sys.path.append('./library')

from parameters import *
from transmission import *
from emission import *
from output import *
from atmosphere import *
from data import *


class create_spectrum(object):

    def __init__(self, params=None, param_filename=None):

        if params:
            self.params = params
        elif hasattr(options, 'param_filename'):
            self.params = parameters(options.param_filename)
        elif param_filename:
            self.params = parameters(param_filename)

        # Initialise taurex instances up to forward model
        self.params.gen_manual_waverange = True
        self.params.nest_run = False
        self.params.mcmc_run = False
        self.params.downhill_run = False

        self.dataob = data(self.params)
        self.atmosphereob = atmosphere(self.dataob)
        if self.params.gen_type == 'transmission':
            self.fmob = transmission(self.atmosphereob)
        elif self.params.gen_type == 'emission':
            self.fmob = emission(self.atmosphereob)

    def generate_spectrum(self, save_db=False, db_filename=None):

        # this function returns the SPECTRUM_out dictionary
        # if filename is not specified, store in out_path/SPECTRUM_out.db
        # also stored in self.SPECTRUM_out

        # create SPECTRUM_out
        outdb = {'type': 'create_spectrum',
                 'params': self.params.params_to_dict()}
        outdata = {}

        # compute spectrum
        outdata['spectrum'] = np.zeros((self.dataob.int_nwlgrid_obs, 3))
        outdata['spectrum'][:,0] = self.dataob.int_wlgrid_obs
        outdata['spectrum'][:,1] = self.fmob.model()

        # freeze the mixing ratio profiles, disable gen_ace
        if self.fmob.params.gen_ace:
            self.fmob.params.gen_ace = False

        # compute contribution function
        outdata['contrib_func'] = self.fmob.model(return_tau=True)

        # calculate opacity contributions
        outdata['opacity_contrib'] = {}
        active_mixratio_profile = self.fmob.atmosphere.active_mixratio_profile
        atm_rayleigh = self.fmob.params.atm_rayleigh
        atm_cia = self.fmob.params.atm_cia
        atm_clouds = self.fmob.params.atm_clouds

        # opacity from molecules
        for idx, val in enumerate(self.atmosphereob.active_gases):
            mask = np.ones(len(self.atmosphereob.active_gases), dtype=bool)
            mask[idx] = 0
            active_mixratio_profile_mask = np.copy(active_mixratio_profile)
            active_mixratio_profile_mask[mask, :] = 0
            self.fmob.atmosphere.active_mixratio_profile = active_mixratio_profile_mask
            self.fmob.params.atm_rayleigh = False
            self.fmob.params.atm_cia = False
            #self.fmob.params.atm_clouds = False
            outdata['opacity_contrib'][val] = np.zeros((self.dataob.int_nwlgrid_obs, 2))
            outdata['opacity_contrib'][val][:,0] = self.dataob.int_wlgrid_obs
            outdata['opacity_contrib'][val][:,1] = self.fmob.model()

        self.fmob.atmosphere.active_mixratio_profile = np.copy(active_mixratio_profile)

        if self.params.gen_type == 'transmission':

            self.fmob.atmosphere.active_mixratio_profile[:, :] = 0

            # opacity from rayleigh
            if atm_rayleigh:
                self.fmob.params.atm_rayleigh = True
                self.fmob.params.atm_cia = False
                #self.fmob.params.atm_clouds = False
                outdata['opacity_contrib']['rayleigh'] = np.zeros((self.dataob.int_nwlgrid_obs, 2))
                outdata['opacity_contrib']['rayleigh'][:,0] = self.dataob.int_wlgrid_obs
                outdata['opacity_contrib']['rayleigh'][:,1] = self.fmob.model()

            # opacity from cia
            if atm_cia:
                self.fmob.params.atm_rayleigh = False
                self.fmob.params.atm_cia = True
                #self.fmob.params.atm_clouds = False
                outdata['opacity_contrib']['cia'] = np.zeros((self.dataob.int_nwlgrid_obs, 2))
                outdata['opacity_contrib']['cia'][:,0] = self.dataob.int_wlgrid_obs
                outdata['opacity_contrib']['cia'][:,1] = self.fmob.model()

            # opacity from clouds
            if atm_clouds:
                self.fmob.params.atm_rayleigh = False
                self.fmob.params.atm_cia = False
                self.fmob.params.atm_clouds = True
                outdata['opacity_contrib']['clouds'] = np.zeros((self.dataob.int_nwlgrid_obs, 2))
                outdata['opacity_contrib']['clouds'][:,0] = self.dataob.int_wlgrid_obs
                outdata['opacity_contrib']['clouds'][:,1] = self.fmob.model()

            self.fmob.atmosphere.active_mixratio_profile = np.copy(active_mixratio_profile)

        # tp profile
        outdata['tp_profile'] = np.zeros((self.atmosphereob.nlayers, 2))
        outdata['tp_profile'][:,0] = self.fmob.atmosphere.pressure_profile
        outdata['tp_profile'][:,1] = self.fmob.atmosphere.temperature_profile

        # mixing ratios
        outdata['active_mixratio_profile'] = np.zeros((len(self.atmosphereob.active_gases), self.atmosphereob.nlayers, 2))
        outdata['inactive_mixratio_profile'] = np.zeros((len(self.atmosphereob.inactive_gases), self.atmosphereob.nlayers, 2))
        for i in range(len(self.atmosphereob.active_gases)):
            outdata['active_mixratio_profile'][i,:,0] = self.fmob.atmosphere.pressure_profile
            outdata['active_mixratio_profile'][i,:,1] =  self.fmob.atmosphere.active_mixratio_profile[i,:]
        for i in range(len(self.atmosphereob.inactive_gases)):
            outdata['inactive_mixratio_profile'][i,:,0] = self.fmob.atmosphere.pressure_profile
            outdata['inactive_mixratio_profile'][i,:,1] =  self.fmob.atmosphere.inactive_mixratio_profile[i,:]

        self.fmob.params.atm_rayleigh = atm_rayleigh
        self.fmob.params.atm_cia = atm_cia
        self.fmob.params.atm_clouds = atm_clouds
        self.fmob.atmosphere.inactive_mixratio_profile = np.copy(active_mixratio_profile)

        # store data
        outdb['data'] = outdata

        self.SPECTRUM_out = outdb

        if save_db:
            if not db_filename:
                db_filename = os.path.join(self.params.out_path, 'SPECTRUM_out.db')

            pickle.dump(outdb, open(db_filename, 'wb'))

        return outdb

    def save_spectrum(self, sp_filename=None):

        if not hasattr(self, 'SPECTRUM_out'):
            self.generate_spectrum()

        if not sp_filename:
            sp_filename = os.path.join(self.params.out_path , 'SPECTRUM_out.dat')

        np.savetxt(sp_filename, self.SPECTRUM_out['data']['spectrum'])


# gets called when running from command line
if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                      dest='param_filename',
                      default='Parfiles/default.par'
                      )
    parser.add_argument('--save_sp',
                      action='store_true',
                      dest='save_sp',
                      default=False)
    parser.add_argument('--save_db',
                      action='store_true',
                      dest='save_db',
                      default=False)
    parser.add_argument('--sp_filename',
                      dest='sp_filename',
                      default=False)
    parser.add_argument('--db_filename',
                      dest='db_filename',
                      default=False)
    parser.add_argument('--plot',
                      dest='plot_spectrum',
                      action='store_true',
                      default=False)

    # add command line interface to parameter file

    params_tmp = parameters()
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
    params = parameters(options.param_filename)

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

    spectrumob = create_spectrum(params=params)

    spectrumob.generate_spectrum(save_db=options.save_db, db_filename=options.db_filename)

    if options.save_sp:
        spectrumob.save_spectrum(sp_filename=options.sp_filename)

    if options.plot_spectrum:
        sp = spectrumob.SPECTRUM_out['data']['spectrum']
        plt.plot(sp[:,0], sp[:,1])
        plt.xscale('log')
        plt.show()