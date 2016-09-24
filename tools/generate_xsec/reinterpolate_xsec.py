# python3 only

# reinterpolate pressure-dependent cross sections for a given temperature into a common grid

import numpy as np
try:
    import cPickle as pickle
except:
    import pickle

import os
import sys
import imp
import glob
import heapq
import multiprocessing
from multiprocessing import Pool

## Import input dictionary
try:
    mol = imp.load_source('molecule', sys.argv[1])
    define = mol.define
except:
    print(' Cannot import dictionary file')
    exit()

mol_name = define['molecule_name']

# group xsec into a single wn range
ncores = 27

folder_xsec = os.path.join(define['working_folder'], 'xsec_combined')
output_folder = os.path.join(define['working_folder'], 'xsec_combined_taurex')

os.system('mkdir -p %s' % output_folder)

def worker(temp_idx):

    temp_val = define['temp_list'][temp_idx]
    filenms = glob.glob('%s/*T%i*' % (folder_xsec, temp_val))
    largest_file = heapq.nlargest(1, filenms, key=os.path.getsize)[0]
    lrg_xsec = pickle.load(open(largest_file, 'rb'), encoding='latin1')
    grid = lrg_xsec[:,0]
    sigma_array = np.zeros((len(define['press_list']), 1, len(grid)))

    for press_idx, press_val in enumerate(np.sort(define['press_list'])):
        print('Doing', temp_val, press_val)
        xs = np.zeros((len(grid), 2))
        xsec_file = '%s_T%i_P%.4e' % (mol_name, temp_val, press_val)
        xsec_load = pickle.load(open(os.path.join(folder_xsec, '%s.xsec.pickle' % xsec_file), 'rb'), encoding='latin1')
        xs[:,0] = grid
        xs[:,1] = np.interp(grid, xsec_load[:,0], xsec_load[:,1])

        sigma_array[press_idx, 0, :] = xs[:,1]
    sigma_out = {
        'name': mol_name,
        'p': np.sort(define['press_list']),
        't': np.asarray([np.float(temp_val)]),
        'wno': grid,
        'xsecarr': sigma_array,
    }
    pickle.dump(sigma_out, open(os.path.join(output_folder, '%s_T%i.TauREx.pickle' % (mol_name, temp_val)), 'wb'), protocol=4)

pool = Pool(processes=int(ncores))
pool_result = pool.map(worker,[temp_idx for temp_idx in range(len(define['temp_list']))])
pool.close()
pool.join()