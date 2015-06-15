import matplotlib.pyplot as plt
import numpy as np
import pickle

run_idx = 4
dat = np.genfromtxt('../output_cmb/param_cov.txt')
fid = pickle.load(open('../data/run{}/fid_values.p'.format(run_idx), "rb"))
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
plt.show()
