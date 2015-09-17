#!/usr/bin/env python
'''
This script plot the covariance matrix with usual Planck colours for the given set of data.


'''


__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"

import matplotlib.pyplot as plt
import numpy as np
import pickle
# READ DATA
# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying_lambda'
run_idx = 2
# ======

dat = np.genfromtxt(base_dir + 'data/{}/run{}/output/param_cov.txt'.format(data_type, str(run_idx)))
fid = pickle.load(open(base_dir + 'data/{}/run{}/fid_values.p'.format(data_type, str(run_idx)), "rb"))
label = {}



label['massless_neutrinos'] = r'$N_{eff}$'
label['hubble'] = r'$H_{0}$'
label['scalar_amp(1)'] = r'$A_{s}$'
label['scalar_spectral_index(1)'] = r'$n_{s}$'
label['omnuh2'] = r'$\Omega$'
label['re_optical_depth'] = r'$\tau$'
label['massless_neutrinos'] =  r'$N_{eff}$'
label['w'] = r'$w$'
label['ombh2'] = r'$\Omega_{b}$'
label['omch2'] = r'$\Omega_{c}$'
label['omnuh2'] = r'$\Omega_{\nu}$'

tick_array = np.arange(np.size(fid.keys()))

plt.imshow(dat,cmap=plt.get_cmap('bwr'),interpolation='nearest',vmin=-1,vmax=1)

plt.xticks(tick_array,[label[key] for key in fid.keys() ])
plt.yticks(tick_array,[label[key] for key in fid.keys() ])


plt.gca().xaxis.set_ticks_position('top')

plt.colorbar()

plt.savefig('../../images/parcov_cmb.pdf')

plt.savefig(base_dir + 'data/{}/run{}/output/parcov_cmb.pdf'.format(data_type, str(run_idx) ), dpi=400, papertype='Letter',
            format='pdf', bbox_inches='tight')

