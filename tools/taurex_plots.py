'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Plotting routines

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

import pickle
import os
import glob
import corner
import argparse

import matplotlib as mpl
from matplotlib.pylab import *
from matplotlib import cm
import matplotlib.pylab as plt
from matplotlib import rc
import matplotlib.patches


#some global matplotlib vars
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 12})

phi = 1.618

class taurex_plots(object):

    def __init__(self, dbfname, **kwargs):

        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        else:
            cmap = 'Paired'

        #norm = mpl.colors.Normalize(vmin=0, vmax=5)
        self.cmap = cm.get_cmap('Paired')

        self.db = pickle.load(open(dbfname))
        self.type = self.db['type']

    def plot_posteriors(self):

        if self.type == 'NEST':

            labels = self.db['fit_params_texlabels']

            N = len(self.db['solutions'])

            # determine ranges for all solutions
            ranges_all = np.zeros((len(self.db['solutions']), len(self.db['fit_params_names']), 3))
            for solution_idx, solution_val in enumerate(self.db['solutions']):
                for param_idx, param_val in enumerate(self.db['fit_params_names']):

                    val = solution_val['fit_params'][param_val]['value']
                    sigma = solution_val['fit_params'][param_val]['sigma']

                    if param_val == 'T_irr':
                        print val, sigma

                    ranges_all[solution_idx, param_idx, 0] = val
                    ranges_all[solution_idx, param_idx, 1] = val - 3*sigma
                    ranges_all[solution_idx, param_idx, 2] = val + 3*sigma

            ranges = []
            for param_idx, param_val in enumerate(self.db['fit_params_names']):
                min = np.min(ranges_all[:, param_idx, 1])
                max = np.max(ranges_all[:, param_idx, 2])
                if min < self.db['fit_params_bounds'][param_idx][0]:
                    min = self.db['fit_params_bounds'][param_idx][0]
                if max > self.db['fit_params_bounds'][param_idx][1]:
                    max = self.db['fit_params_bounds'][param_idx][1]

                ranges.append((min, max))

            figs = []
            for solution_idx, solution_val in enumerate(self.db['solutions']):

                tracedata = solution_val['tracedata']
                weights = solution_val['weights']

                if solution_idx > 0:
                    figure_past = figs[solution_idx - 1]
                else:
                    figure_past = None


                fig =  corner.corner(tracedata,
                                      weights=weights,
                                      labels=labels,
                                      #truths=[0.0, 0.0, 0.0],
                                      quantiles=[0.16, 0.5, 0.84],
                                      #quantiles=[0.16, 0.5],
                                      range=ranges,
                                      ret=True,
                                      fill_contours=True,
                                      color=self.cmap(float(solution_idx)/N),
                                      bins=30,
                                      fig = figure_past)
                figs.append(fig)

            plt.savefig(os.path.join(self.db['params']['out_path'], 'posteriors.pdf'))

    def plot_fitted_spectrum(self):

        # fitted model
        fig = plt.figure(figsize=(7,7/phi))
        ax = fig.add_subplot(111)
        obs = self.db['obs_spectrum']
        plt.errorbar(obs[:,0], obs[:,1], obs[:,2], lw=1, color='black', alpha=0.5, ls='none', zorder=99, label='Observed')
        N = len(self.db['solutions'])
        for solution_idx, solution_val in enumerate(self.db['solutions']):
            if N > 1:
                label = 'Fitted model (%i)' % (solution_idx + 1)
            else:
                label = 'Fitted model'
            spectra = solution_val['obs_spectrum']
            fit_highres = solution_val['fit_spectrum_xsecres']
            fit_obsres = solution_val['fit_spectrum_xsecres']
            plt.plot(fit_highres[:,0], fit_highres[:,1], color=self.cmap(float(solution_idx)/N), label=label)
            plt.plot(spectra[:,0], spectra[:,3], 'd', markersize=3, color=self.cmap(float(solution_idx)/N))
        plt.xlim(np.min(obs[:,0])-0.05*np.min(obs[:,0]), np.max(obs[:,0])+0.05*np.max(obs[:,0]))
        plt.xscale('log')
        plt.xlabel('Wavelength ($\mu$m)')
        plt.ylabel('$(R_p/R_*)^2$')
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        plt.legend(loc='auto', ncol=2, frameon=False, prop={'size':12})
        plt.savefig('fitted_model.pdf')
        plt.savefig(os.path.join(self.db['params']['out_path'], 'fitted_model.pdf'))

        # contribution
        N = len(self.db['solutions'][0]['opacity_contrib'])
        for solution_idx, solution_val in enumerate(self.db['solutions']):
            fig = plt.figure(figsize=(7,7/phi))
            ax = fig.add_subplot(111)
            plt.plot(fit_highres[:,0], fit_highres[:,1], alpha=0.7, color='black', label='Fitted model')
            plt.errorbar(obs[:,0], obs[:,1], obs[:,2], lw=1, color='black', alpha=0.5, ls='none', zorder=99, label='Observed')
            for idx, val in enumerate(self.db['solutions'][0]['opacity_contrib']):
                sp = self.db['solutions'][0]['opacity_contrib'][val]
                plt.plot(sp[:,0], sp[:,1], color=self.cmap(float(idx)/N), label=val)
            plt.xlim(np.min(obs[:,0])-0.05*np.min(obs[:,0]), np.max(obs[:,0])+0.05*np.max(obs[:,0]))
            plt.xscale('log')
            plt.xlabel('Wavelength ($\mu$m)')
            plt.ylabel('$(R_p/R_*)^2$')
            plt.tick_params(axis='x', which='minor')
            ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            legend = plt.legend(loc='upper left', ncol=1, prop={'size':12})
            legend.get_frame().set_facecolor('white')
            legend.get_frame().set_edgecolor('white')
            legend.get_frame().set_alpha(0.8)
            plt.savefig(os.path.join(self.db['params']['out_path'], 'fitted_model_contrib_sol%i.pdf' % solution_idx))

    def plot_fitted_tp(self):

        # fitted model
        fig = plt.figure(figsize=(7,7/phi))
        ax = fig.add_subplot(111)
        N = len(self.db['solutions'])
        for solution_idx, solution_val in enumerate(self.db['solutions']):
            if N > 1:
                label = 'Fitted profile (%i)' % (solution_idx + 1)
            else:
                label = 'Fitted profile'
            tp = solution_val['tp_profile']
            plt.plot(tp[:,1], tp[:,0]/1e5, color=self.cmap(float(solution_idx)/N), label=label)
            plt.fill_betweenx(tp[:,0]/1e5,  tp[:,1]-tp[:,2],  tp[:,1]+tp[:,2], color=self.cmap(float(solution_idx)/N), alpha=0.5)

        # add input spectrum param if available
        if 'SPECTRUM_db' in self.db:
            tp = self.db['SPECTRUM_db']['data']['tp_profile']
            plt.plot(tp[:,1], tp[:,0]/1e5, color='#C04F55', label='Input profile')

        plt.gca().invert_yaxis()
        plt.yscale('log')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Pressure (bar)')
        plt.tight_layout()
        legend = plt.legend(loc='upper left', ncol=1, prop={'size':12})
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_edgecolor('white')
        legend.get_frame().set_alpha(0.8)
        plt.savefig(os.path.join(self.db['params']['out_path'], 'tp_profile.pdf'))

    def plot_fitted_xprofiles(self):

        # fitted model
        for solution_idx, solution_val in enumerate(self.db['solutions']):
            fig = plt.figure(figsize=(7,7/phi))
            ax = fig.add_subplot(111)
            N = len(self.db['params']['atm_active_gases']+self.db['params']['atm_inactive_gases'])
            for mol_idx, mol_val in enumerate(self.db['params']['atm_active_gases']):
                prof = solution_val['active_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(mol_idx)/N), label=mol_val)
                #plt.fill_betweenx(prof[:,0]/1e5,  prof[:,1]-prof[:,2],  prof[:,1]+prof[:,2], color=self.cmap(float(mol_idx)/N), alpha=0.5)
            for mol_idx, mol_val in enumerate(self.db['params']['atm_inactive_gases']):
                prof = solution_val['inactive_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(mol_idx+len(self.db['params']['atm_inactive_gases']))/N), label=mol_val)
                #plt.fill_betweenx(prof[:,0]/1e5,  prof[:,1]-prof[:,2],  prof[:,1]+prof[:,2], color=self.cmap(float(mol_idx)/N), alpha=0.5)

        # add input spectrum param if available
        if 'SPECTRUM_db' in self.db:
            for mol_idx, mol_val in enumerate(self.db['params']['atm_active_gases']):
                prof = self.db['SPECTRUM_db']['data']['active_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(mol_idx)/N), ls='dashed')
            for mol_idx, mol_val in enumerate(self.db['params']['atm_inactive_gases']):
                prof = self.db['SPECTRUM_db']['data']['inactive_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(mol_idx+len(self.db['params']['atm_inactive_gases']))/N), ls='dashed')

        plt.gca().invert_yaxis()
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e-12, 1)
        plt.xlabel('Mixing ratio')
        plt.ylabel('Pressure (bar)')
        plt.tight_layout()
        legend = plt.legend(loc='upper left', ncol=1, prop={'size':12})
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_edgecolor('white')
        legend.get_frame().set_alpha(0.8)
        plt.savefig(os.path.join(self.db['params']['out_path'], 'mixratio.pdf'))



if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--db_filename',
                          dest='db_filename',
                          default=False)

    options = parser.parse_args()

    plot = taurex_plots(options.db_filename)
    plot.plot_posteriors()
    plot.plot_fitted_spectrum()
    # plot.plot_fitted_xprofiles()
    # plot.plot_fitted_tp()