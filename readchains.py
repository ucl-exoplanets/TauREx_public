
import numpy as np
import matplotlib.pylab as plt

nparams = 3
nest_data = np.loadtxt('Output/multinest/1-.txt')
nchains = len(nest_data[:,0])

mixing_chains = []
for line in nest_data:
    clr = []
    for n in [0,1,2]:
        clr.append(line[n+2])
    clr.append(-np.sum(clr)) # append last clr. Note that sum(clr) = 0!

    # transform back to simplex space
    clr_inv_tmp = []
    for i in range(nparams):
        clr_inv_tmp.append(np.exp(clr[i]))
    clr_inv_tmp = np.asarray(clr_inv_tmp)/np.sum(clr_inv_tmp)

    mixing_chains.append(clr_inv_tmp)

mixing_chains = np.asarray(mixing_chains).transpose()

#
# plt.plot(mixing_chains[3], mixing_chains[4], 'o', alpha=0.1)
# plt.show()
#
# exit()

for i in range(nparams):
    plt.subplot(3,2,i+1)
    plt.hist(mixing_chains[i], 500)
    plt.title(i)

plt.show()
exit()

# nested = np.zeros((nparams, nchains))
# count = 0
# for i in [3,4,5,6,7]:
#     nested[count,:] = nest_data[:, i+2]
#     count = count+1
#
# nested[5] = np.negative(np.sum(nested, axis=0))
#
# nested_exp = np.exp(nested)
# row_sums = nested_exp.sum(axis=0)
# print np.shape(row_sums)
#
# mixing_chains = nested / np.array([row_sums,] * nparams)

