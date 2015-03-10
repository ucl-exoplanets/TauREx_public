import numpy as np
import logging

def analyse_traces(self, tracedata, clr_inv=False, cluster_analysis=False):

    ntraces = len(tracedata[:,0])
    ncol = len(tracedata[0,:])

    if clr_inv == True:
        # create new array of traces with the mixing ratios obtained from CLR traces
        # log-ratios are always the first n columns
        mixing_ratios_clr = tracedata[:, :self.fitting.forwardmodel.atmosphere.nallgases-1] # get traces of log-ratios
        mixing_ratios_clr_inv, coupled_mu_trace = self.inverse_clr_transform(mixing_ratios_clr)

        self.tracedata_clr_inv = np.c_[mixing_ratios_clr_inv, coupled_mu_trace,
                                       tracedata[:, self.fitting.forwardmodel.atmosphere.nallgases-1:]]
        tracedata_analyse = self.tracedata_clr_inv # set the trace to analyse and get solutions from

        # set parameter names of new traces
        self.params_names = []
        for idx, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases +
                self.fitting.forwardmodel.atmosphere.inactive_gases):
            self.params_names.append(gasname)
        # todo careful about adding coupled_mu, it might interefer with fit-params?
        self.params_names.append('coupled_mu')
        for i in range(self.fitting.forwardmodel.atmosphere.nallgases-1, len(self.fitting.fit_params_names)):
            self.params_names.append(self.fitting.fit_params_names[i])

    else:
        self.params_names = self.fitting.fit_params_names
        tracedata_analyse = tracedata # set the trace to analyse and get solutions from


    # todo: if cluster_analysis is False, but multinest is in multimode, then we should get solutions from multinest output. See end of file with commented code.
    if cluster_analysis == True:

        # clustering analysis on traces. Split up the traces into individual clusters
        logging.info('Cluster analysis with sklearn MeanShift')

        from sklearn.cluster import MeanShift
        from sklearn.preprocessing import normalize
        data_normalized = normalize(tracedata_analyse, axis=0)  # normalize the data
        estimator = MeanShift()
        estimator.fit(data_normalized)
        labels = estimator.labels_

        # save txt with labels corresponding to different clusters found in the posterior traces
        np.savetxt('Output/labels.txt', labels)
        # labels = np.loadtxt('Output/labels.txt')
        unique_labels = set(labels)
        n_clusters = len(unique_labels) - (1 if -1 in labels else 0)
        logging.info('Estimated number of clusters: %d' % n_clusters)

    else:
        # clustering analysis is not performed. Only one mode assumed
        unique_labels = [0]
        labels = np.zeros(ntraces)

    fit_out = self.get_solutions_from_traces(tracedata_analyse, multimode=cluster_analysis, labels=labels) # todo: needs to be changed, if we get solutions from multinest multimode (see above)
    fit_out = self.add_spectra_from_solutions(fit_out)

    return fit_out, tracedata_analyse, labels


def get_solutions_from_traces(self, tracedata, multimode=False, labels=[]):

    ntraces = len(tracedata[:,0])
    ncol = len(tracedata[0,:])

    if multimode == False:
        unique_labels = [0]
        labels = np.zeros(ntraces)
    else:
        unique_labels = set(labels)
        if len(labels) != ntraces:
            logging.error('Multimode is ON and the number of labels is different by the number of samples')

    # Create list of solutions. Each solution (cluster) is stored into a dictionary. If multimode = False, then
    # only only one solution is stored.

    fit_out = []
    for k in unique_labels:
        if k >= 0:

            if (k > 0 and len(tracedata[labels == k-1, 0])/len(tracedata[labels == k, 0]) > 20):
                # exlcude solutions that have a factor of 20 fewer samples than the previous solution @todo careful!
                break
            else:

                # save individual solutions into a dictionary
                dict = {'fit_params': {}}
                sort_order = 0
                for idx, param in enumerate(self.params_names):

                    try:
                        cluster = tracedata[labels == k, idx] # get cluster of points correposning to cluster k, and param idx

                        # get errors
                        hist, bin_edges = np.histogram(cluster, bins=100, density=True)
                        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
                        idxmax = np.argmax(hist)
                        value = bin_centres[idxmax]
                        left = hist[:idxmax][::-1]
                        right = hist[idxmax:]
                        left_centres = bin_centres[:idxmax][::-1]
                        right_centres = bin_centres[idxmax:]

                        try:
                            left_err = value - left_centres[(np.abs(left-np.percentile(left,68))).argmin()]
                        except:
                            left_err = 0
                        try:
                            right_err = right_centres[(np.abs(right-np.percentile(right,68))).argmin()] - value
                        except:
                            right_err = 0
                        if left_err == 0 and right_err > 0:
                            left_err = right_err
                        if right_err == 0 and left_err > 0:
                            right_err = left_err

                        mean_err = np.mean((left_err, right_err))

                        dict['fit_params'][self.params_names[idx]] = {
                            'value': value,
                            'std': np.std(cluster),
                            #'std': mean_err, @ todo try to use mean err!
                            'std_plus': right_err,
                            'std_minus': left_err,
                            'trace': np.asarray(cluster),
                            'sort_order': sort_order
                        }
                        sort_order += 1

                    except ValueError:
                        # this parameter was not fitted...
                        pass
                fit_out.append(dict)

    logging.info('Number of solutions identified: %d' % len(fit_out))

    return fit_out


def inverse_clr_transform(self, clr_tracedata):

    # Convert the traces of the log-ratios back to the gas mixing ratios in log space. This is done by calculating the inverse
    # CLR in the tracedata array, line by line. The new traces are then returned.

    logging.info('Inverse log-ratio transformation of the traces')

    ntraces = len(clr_tracedata[:,0])
    ncol = len(clr_tracedata[0,:])

    mixing_ratios_tracedata = np.exp(np.c_[clr_tracedata, -sum(clr_tracedata, axis=1)]) # add log-ratio (= -sum(other log ratios)
    mixing_ratios_tracedata /= (np.asarray([np.sum(mixing_ratios_tracedata, axis=1),]*(ncol+1)).transpose()) # clousure operation

    # calculate couple mean molecular weight trace as a byproduct
    coupled_mu_trace = np.zeros(ntraces)
    i = 0
    for sample in mixing_ratios_tracedata:
        # get mean molecular weight from gas mixing ratios
        absorbing_gases_X = sample[:len(self.fitting.forwardmodel.atmosphere.absorbing_gases)]
        inactive_gases_X = sample[len(self.fitting.forwardmodel.atmosphere.inactive_gases):]
        coupled_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu(absorbing_gases_X, inactive_gases_X)
        coupled_mu_trace[i] = coupled_mu
        i += 1

    return np.power(10, mixing_ratios_tracedata), coupled_mu_trace

def add_spectra_from_solutions(self, fitting_out):

    for idx, solution in enumerate(fitting_out):

        fit_params = [solution['fit_params'][param]['value'] for param in self.params_names]
        self.fitting.update_atmospheric_parameters(fit_params)
        model = self.fitting.forwardmodel.model()
        model_binned = [model[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid)]

        out = np.zeros((len(self.data.spectrum[:,0]), 2))
        out[:,0] = self.data.spectrum[:,0]
        out[:,1] = model_binned

        fitting_out[idx]['spectrum'] = out

    return fitting_out
