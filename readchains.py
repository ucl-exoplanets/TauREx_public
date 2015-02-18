
import numpy as np
import matplotlib.pylab as plt

nparams = 5

plt.ion()

        # clf()
        # errorbar(self.data.spectrum[:,0],self.data.spectrum[:,1], yerr=self.data.spectrum[:,2])
        # plot(self.data.spectrum[:,0], model_binned)
        # draw()
        # pause(0.01)

while True:
    plt.clf()
    nest_data = np.loadtxt('Output/multinest/1-.txt')
    nchains = len(nest_data[:,0])
    mixing_chains = []
    for line in nest_data:
        clr = []
        for n in [5,6,7]:
            clr.append(line[n])
        clr.append(-np.sum(clr)) # append last clr. Note that sum(clr) = 0!

        # transform back to simplex space
        clr_inv_tmp = []
        for i in range(nparams):
            clr_inv_tmp.append(np.exp(clr[i]))
        clr_inv_tmp = np.asarray(clr_inv_tmp)/np.sum(clr_inv_tmp)

        mixing_chains.append(clr_inv_tmp)

    mixing_chains = np.asarray(mixing_chains).transpose()

    for i in range(nparams):
        plt.subplot(3,2,i+1)
        plt.hist(mixing_chains[i], 50)
        plt.title(i)

    plt.draw()
    plt.pause(1)
