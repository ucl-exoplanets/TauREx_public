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

import warnings
warnings.filterwarnings("ignore", module="matplotlib")


#some global matplotlib vars
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['text.antialiased'] = True
#rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 11})

phi = 1.618

class taurex_plots(object):

    def __init__(self, dbfname=False, multinest_dir=False, multinest_prefix='1-', title=False, prefix='', out_folder=False, **kwargs):

        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        else:
            cmap = 'Paired'

        #norm = mpl.colors.Normalize(vmin=0, vmax=5)
        self.cmap = cm.get_cmap('Paired')

        if dbfname:
            self.db = pickle.load(open(dbfname))
            try:
                self.type = self.db['type']
            except:
                self.type = 'nest'
        elif multinest_dir:
            self.db = False
            self.type = 'nest'

        if multinest_dir:
            self.multinest_dir = multinest_dir
            self.multinest_prefix = multinest_prefix
        else:
            self.multinest_dir = False

        if not dbfname and not multinest_dir:
            print 'Specify a .db filename or a multinest directory'
            exit()

        self.title = title
        if not prefix:
            self.prefix = ''
        else:
            self.prefix = prefix

        if out_folder:
            self.out_folder = out_folder
        else:
            if dbfname:
                self.out_folder = self.db['params']['out_path']
            elif multinest_dir:
                self.out_folder = multinest_dir

        try:
            os.makedirs(self.out_folder)
        except:
            print 'Cannot create output directory: %s' % self.out_folder


    def plot_posteriors(self, **kwargs):

        if self.multinest_dir:
            self._plot_posteriors_from_multinest(**kwargs)
        elif self.db:
            self._plot_posteriors_db(**kwargs)

    def _plot_posteriors_from_multinest(self):

        if self.multinest_dir:

            # todo get labels...

            tracedata = np.loadtxt(os.path.join(self.multinest_dir, '%s.txt' % self.multinest_prefix))
            weights = tracedata[:,0]
            traces = tracedata[:,2:]

            fig =  corner.corner(traces,
                                  weights=weights,
                                  # labels=labels,
                                  smooth=True,
                                  scale_hist=True,
                                  #truths=truths,
                                  quantiles=[0.16, 0.5, 0.84],
                                  show_titles=True,
                                  ret=True,
                                  fill_contours=True,
                                  color='#4F75AD',
                                  bins=30)

            if self.title:
                fig.gca().annotate(self.title, xy=(0.5, 1.0), xycoords="figure fraction",
                  xytext=(0, -5), textcoords="offset points",
                  ha="center", va="top", fontsize=14)
        filename = os.path.join(self.out_folder, '%s%s_posteriors.pdf' % (self.prefix, self.type))
        plt.savefig(filename)
        print 'Plot saved in %s' % filename

    def _plot_posteriors_db(self, plot_truths=False):

        if self.db:

            if self.type.upper() == 'NEST':

                labels = self.db['fit_params_texlabels']

                N = len(self.db['solutions'])

                # determine ranges for all solutions
                ranges_all = np.zeros((len(self.db['solutions']), len(self.db['fit_params_names']), 3))
                for solution_idx, solution_val in enumerate(self.db['solutions']):
                    for param_idx, param_val in enumerate(self.db['fit_params_names']):

                        val = solution_val['fit_params'][param_val]['value']
                        try:
                            sigma_m = solution_val['fit_params'][param_val]['sigma_m']
                            sigma_p = solution_val['fit_params'][param_val]['sigma_p']
                        except:
                            sigma_m = sigma_p = solution_val['fit_params'][param_val]['sigma']

                        ranges_all[solution_idx, param_idx, 0] = val
                        ranges_all[solution_idx, param_idx, 1] = val - 5*sigma_m
                        ranges_all[solution_idx, param_idx, 2] = val + 5*sigma_p

                ranges = []
                for param_idx, param_val in enumerate(self.db['fit_params_names']):
                    min = np.min(ranges_all[:, param_idx, 1])
                    max = np.max(ranges_all[:, param_idx, 2])
                    if min < self.db['fit_params_bounds'][param_idx][0]:
                        min = self.db['fit_params_bounds'][param_idx][0]
                    if max > self.db['fit_params_bounds'][param_idx][1]:
                        max = self.db['fit_params_bounds'][param_idx][1]

                    ranges.append((min, max))

                # build imput value lists if input SPECTRUM_db is in NEST_db
                if 'SPECTRUM_db' in self.db and plot_truths:
                    truths = self.build_truths()
                else:
                    truths = None

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
                                          smooth=True,
                                          scale_hist=True,
                                          truths=truths,
                                          quantiles=[0.16, 0.5, 0.84],
                                          show_titles=True,
                                          range=ranges,
                                          #quantiles=[0.16, 0.5],
                                          ret=True,
                                          fill_contours=True,
                                          color=self.cmap(float(solution_idx)/N),
                                          bins=30,
                                          fig = figure_past)
                    if self.title:
                        fig.gca().annotate(self.title, xy=(0.5, 1.0), xycoords="figure fraction",
                          xytext=(0, -5), textcoords="offset points",
                          ha="center", va="top", fontsize=14)

                    figs.append(fig)

                plt.savefig(os.path.join(self.out_folder, '%s%s_posteriors.pdf' % (self.prefix, self.type)))

    def plot_fitted_spectrum(self, plot_contrib=True):

        # fitted model
        fig = plt.figure(figsize=(5.3, 3.5))
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
            rangeidx = np.logical_and(fit_highres[:,0] > np.min(obs[:,0])-0.05*np.min(obs[:,0]), fit_highres[:,0] < np.max(obs[:,0])+0.05*np.max(obs[:,0]))
            plt.plot(fit_highres[:,0][rangeidx], fit_highres[:,1][rangeidx], color=self.cmap(float(solution_idx)/N), label=label)
            plt.plot(spectra[:,0], spectra[:,3], 'd', markersize=3, color=self.cmap(float(solution_idx)/N))

            # add 1 sigma spread
            if self.db['params']['out_sigma_spectrum']:
                plt.fill_between(solution_val['fit_spectrum_xsecres'][:,0][rangeidx],
                solution_val['fit_spectrum_xsecres'][:,1][rangeidx]-solution_val['fit_spectrum_xsecres'][:,2][rangeidx],
                solution_val['fit_spectrum_xsecres'][:,1][rangeidx]+solution_val['fit_spectrum_xsecres'][:,2][rangeidx],
                alpha=0.5, zorder=-1, color=self.cmap(float(solution_idx)/N), edgecolor="none")
                plt.fill_between(solution_val['fit_spectrum_xsecres'][:,0][rangeidx],
                solution_val['fit_spectrum_xsecres'][:,1][rangeidx]-2*solution_val['fit_spectrum_xsecres'][:,2][rangeidx],
                solution_val['fit_spectrum_xsecres'][:,1][rangeidx]+2*solution_val['fit_spectrum_xsecres'][:,2][rangeidx],
                alpha=0.2, zorder=-1, color=self.cmap(float(solution_idx)/N), edgecolor="none")

        plt.xlim(np.min(obs[:,0])-0.05*np.min(obs[:,0]), np.max(obs[:,0])+0.05*np.max(obs[:,0]))
        plt.xlabel('Wavelength ($\mu$m)')
        if self.db['params']['fit_emission']:
            plt.ylabel('$F_p/F_*$')
        else:
            plt.ylabel('$(R_p/R_*)^2$')
        # set log scale only if interval is greater than 5 micron
        if np.max(obs[:,0]) - np.min(obs[:,1]) > 5:
            plt.xscale('log')
            plt.tick_params(axis='x', which='minor')
            ax.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%i"))
        plt.legend(loc='best', ncol=2, frameon=False, prop={'size':11})
        if self.title:
            plt.title(self.title, fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_folder, '%s%s_spectrum.pdf'  % (self.prefix, self.type)))

        # contribution
        if plot_contrib:
            N = len(self.db['solutions'][0]['opacity_contrib'])
            for solution_idx, solution_val in enumerate(self.db['solutions']):
                fig = plt.figure(figsize=(7,3.5))
                ax = fig.add_subplot(111)
                rangeidx = np.logical_and(fit_highres[:,0] > np.min(obs[:,0])-0.05*np.min(obs[:,0]), fit_highres[:,0] < np.max(obs[:,0])+0.05*np.max(obs[:,0]))
                plt.plot(fit_highres[:,0][rangeidx], fit_highres[:,1][rangeidx], alpha=0.7, color='black', label='Fitted model')
                plt.errorbar(obs[:,0], obs[:,1], obs[:,2], lw=1, color='black', alpha=0.5, ls='none', zorder=99, label='Observed')
                for idx, val in enumerate(self.db['solutions'][solution_idx]['opacity_contrib']):
                    sp = self.db['solutions'][solution_idx]['opacity_contrib'][val]
                    plt.plot(sp[:,0][rangeidx], sp[:,1][rangeidx], color=self.cmap(float(idx)/N), label=val)
                plt.xlim(np.min(obs[:,0])-0.05*np.min(obs[:,0]), np.max(obs[:,0])+0.05*np.max(obs[:,0]))
                plt.xlabel('Wavelength ($\mu$m)')
                plt.ylabel('$(R_p/R_*)^2$')
                plt.tick_params(axis='x', which='minor')
                # set log scale only if interval is greater than 5 micron
                if np.max(obs[:,0]) - np.min(obs[:,1]) > 5:
                    plt.xscale('log')
                    plt.tick_params(axis='x', which='minor')
                    ax.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                    ax.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                plt.tight_layout()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'size':11}, frameon=False)
                if self.title:
                    plt.title(self.title, fontsize=14)
                plt.savefig(os.path.join(self.out_folder, '%s%s_spectrum_contrib_sol%i.pdf' % (self.prefix, self.type, solution_idx)))

    def plot_fitted_tp(self):

        # fitted model
        fig = plt.figure(figsize=(5,3.5))
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
        legend = plt.legend(loc='upper left', ncol=1, prop={'size':11})
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_edgecolor('white')
        legend.get_frame().set_alpha(0.8)
        if self.title:
            plt.title(self.title, fontsize=14)
        plt.savefig(os.path.join(self.out_folder, '%s%s_tp_profile.pdf' % (self.prefix, self.type)))

    def plot_fitted_xprofiles(self):

        # fitted model
        for solution_idx, solution_val in enumerate(self.db['solutions']):
            fig = plt.figure(figsize=(7,7/phi))
            ax = fig.add_subplot(111)
            N = len(self.db['params']['atm_active_gases']+self.db['params']['atm_inactive_gases'])
            cols_mol = {}
            for mol_idx, mol_val in enumerate(self.db['params']['atm_active_gases']):
                cols_mol[mol_val] = self.cmap(float(mol_idx)/N)
                prof = solution_val['active_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], label=mol_val)
                #plt.fill_betweenx(prof[:,0]/1e5,  prof[:,1]-prof[:,2],  prof[:,1]+prof[:,2], color=self.cmap(float(mol_idx)/N), alpha=0.5)
            for mol_idx, mol_val in enumerate(self.db['params']['atm_inactive_gases']):
                cols_mol[mol_val] = self.cmap(float(mol_idx+len(self.db['params']['atm_inactive_gases']))/N)
                prof = solution_val['inactive_mixratio_profile'][mol_idx]
                plt.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], label=mol_val)
                #plt.fill_betweenx(prof[:,0]/1e5,  prof[:,1]-prof[:,2],  prof[:,1]+prof[:,2], color=self.cmap(float(mol_idx)/N), alpha=0.5)

        # add input spectrum param if available
        if 'SPECTRUM_db' in self.db:
            for mol_idx, mol_val in enumerate(self.db['SPECTRUM_db']['params']['atm_active_gases']):
                prof = self.db['SPECTRUM_db']['data']['active_mixratio_profile'][mol_idx]
                if mol_val in cols_mol:
                    plt.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], ls='dashed')
                else:
                    plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(len(self.db['params']['atm_active_gases'])+len(self.db['params']['atm_inactive_gases'])+mol_idx)/N), label=mol_val, ls='dashed')

            for mol_idx, mol_val in enumerate(self.db['SPECTRUM_db']['params']['atm_inactive_gases']):
                prof = self.db['SPECTRUM_db']['data']['inactive_mixratio_profile'][mol_idx]
                if mol_val in cols_mol:
                    plt.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], ls='dashed')
                else:
                    plt.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(len(self.db['params']['atm_active_gases'])+len(self.db['params']['atm_inactive_gases'])+mol_idx)/N), label=mol_val, ls='dashed')

        plt.gca().invert_yaxis()
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e-12, 3)
        plt.xlabel('Mixing ratio')
        plt.ylabel('Pressure (bar)')
        plt.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'size':11}, frameon=False)
        if self.title:
            plt.title(self.title, fontsize=14)
        plt.savefig(os.path.join(self.out_folder, '%s%s_mixratio.pdf' % (self.prefix, self.type)))

    def build_truths(self):

        truths = []
        in_params = self.db['SPECTRUM_db']['params']
        for param in self.db['fit_params_names']:
            if param[:3] == 'log' and (param[4:] in in_params['atm_active_gases']):
                mol_idx = in_params['atm_active_gases'].index(param[4:])
                truths.append(np.log10(in_params['atm_active_gases_mixratios'][mol_idx]))
            elif param in in_params['atm_active_gases']:
                mol_idx = in_params['atm_active_gases'].index(param)
                truths.append(in_params['atm_active_gases_mixratios'][mol_idx])
            elif param == 'ace_log_metallicity':
                truths.append(in_params['atm_ace_metallicity'])
            elif param == 'ace_co':
                truths.append(in_params['fit_fit_ace_co'])
            elif param == 'T':
                truths.append(in_params['atm_tp_iso_temp'])
            elif param == 'mu':
                truths.append(in_params['atm_mu']/1.660538921e-27) # todo check, might not be right. See how mu is coupled in the fitting
            elif param == 'Radius':
                truths.append(in_params['planet_radius']/69911000.0)
            elif param == 'clouds_topP':
                if in_params['atm_clouds']:
                    truths.append(np.log(in_params['atm_cld_topP']))
                else:
                    truths.append(6) # if no clouds in input spectrum, but fitting for clouds, assume clouds at 10 bar
            else:
                truths.append(0)
        return truths


if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--db_filename',
                          dest='db_filename',
                          default=False)
    parser.add_argument('--db_dir',
                          dest='db_dir',
                          default=False)
    parser.add_argument('--multinest_dir',
                          dest='multinest_dir',
                          default=False)
    parser.add_argument('--multi_folders',
                           dest='multi_dir',
                           default=False)
    parser.add_argument('--multinest_prefix',
                          dest='multinest_prefix',
                          default='1-')
    parser.add_argument('--title',
                          dest='title',
                          default=False)
    parser.add_argument('--prefix',
                          dest='prefix',
                          default='')
    parser.add_argument('--out_folder',
                          dest='out_folder',
                          default=False)
    parser.add_argument('--plot_all',
                          dest='plot_all',
                          action='store_true',
                          default=False)
    parser.add_argument('--plot_posteriors',
                          dest='plot_posteriors',
                          action='store_true',
                          default=False)
    parser.add_argument('--plot_spectrum',
                          dest='plot_spectrum',
                          action='store_true',
                          default=False)
    parser.add_argument('--plot_tp',
                          dest='plot_tp',
                          action='store_true',
                          default=False)
    parser.add_argument('--plot_x',
                          dest='plot_x',
                          action='store_true',
                          default=False)

    options = parser.parse_args()


    if options.db_dir:
        db_filenames = glob.glob(os.path.join(options.db_dir, '*.db'))
        out_folders = [val[:-3] for val in db_filenames]
    elif options.multi_dir:
        out_folders = glob.glob(os.path.join(options.multi_dir, '*'))
        db_filenames = []
        for folder in out_folders:
            tmp = os.path.join(folder,'nest_out.db')
            if os.path.exists(tmp):
                db_filenames.append(tmp)
        print db_filenames
    else:
        db_filenames = [options.db_filename]
        out_folders = [options.out_folder]

    for db_idx, db_val in enumerate(db_filenames):

        plot = taurex_plots(dbfname=db_val,
                            multinest_dir = options.multinest_dir,
                            multinest_prefix = options.multinest_prefix,
                            title=options.title,
                            prefix=options.prefix,
                            out_folder=out_folders[db_idx])

        if options.plot_all or options.plot_posteriors:
            print 'Plotting posteriors'
            plot.plot_posteriors()

#         if options.db_dir or options.db_filename:

        if options.plot_all or options.plot_spectrum:
                print 'Plotting fitted spectrum'
                plot.plot_fitted_spectrum()
        if options.plot_all or options.plot_x:
                print 'Plotting mixing ratio profiles'
                plot.plot_fitted_xprofiles()
        if options.plot_all or options.plot_tp:
                print 'Plotting pressure-temperature profile'
                plot.plot_fitted_tp()