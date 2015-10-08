'''

# DATA FROM CAMB

we use the _lenspotentialcls so

 l CTT CEE CBB CTE Cdd CdT CdE
 0  1  2    3   4   5   6  7

  CX are l(l+1)Cl/2pi and Cdd=[l(l+1)]^2 Clphi/2pi, CdT=[l(l+1)]^(3/2) ClphiT/2pi, CdE=[l(l+1)]^(3/2)ClphiE/2pi


'''


import numpy as np
import pickle as pickle
#import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special
import scipy.interpolate
from itertools import cycle
#from matplotlib.colors import colorConverter
from matplotlib.ticker import MaxNLocator  # needed for tick
import n_priors_code.utils as utils

plot_type ='lmax'

no_lcdm_parameters = ['massless_neutrinos', 'w', 'omnuh2', 'helium_fraction']
plot_now = ['scalar_spectral_index(1)','omnuh2']
excluded_parameters = list(set(no_lcdm_parameters) - set(plot_now))
# omnuh2

# READ DATA
# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying+Yp'
run_idx = 1
lmax = 3000
lmin = 50
N_det = 10 ** 6
# N_phi_l = np.loadtxt('data/noise/wu_cdd_noise_5.txt')
fsky = 0.75
# ======
fid = pickle.load(open(base_dir + 'data/{}/run{}/fid_values.p'.format(data_type, str(run_idx)), "rb"))
values = pickle.load(open(base_dir + 'data/{}/run{}/grid_values.p'.format(data_type, str(run_idx)), "rb"))
par_gaps = pickle.load(open(base_dir + 'data/{}/run{}/par_gaps.p'.format(data_type, str(run_idx)), "rb"))
fisher_mat = np.loadtxt(
    base_dir + 'data/{}/run{}/output/fisher_mat_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(data_type, str(run_idx), lmin, lmax, N_det, fsky))

par_gaps, values, fid, fisher_mat = utils.exclude_parameters_from_fisher(
    excluded_parameters, par_gaps, values, fid, fisher_mat)
plot_param = list(set(fid.keys()) - set(excluded_parameters))
sigma1 = np.loadtxt(base_dir + 'data/{}/run{}/output/marginalized_ell_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(
    data_type, str(run_idx), lmin, lmax, N_det, fsky))
sigma_no = np.loadtxt(base_dir + 'data/{}/run{}/output/no_marginalized_ell_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(
    data_type, str(run_idx), lmin, lmax, N_det, fsky))
full_fisher = np.load(base_dir+ 'data/{}/run{}/output/full_fisher_mat_joint_lmin={}_lmax={}_ndet={}_fsky={}.npy'.format(data_type, str(run_idx),lmin,lmax,N_det,fsky))

# ============================================
# ============================================
# ============================================


# ============================================
# SIZE OF THE PICTURE
# ============================================
WIDTH = 246.0  # the number latex spits out
# the fraction of the width you'd like the figure to occupy, 0.5 if double
# column cause
FACTOR = 0.9
fig_width_pt = WIDTH * FACTOR

inches_per_pt = 1.0 / 72.27
# golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
ratio = 0.8
fig_width_in = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * ratio   # figure height in inches
fig_dims = [fig_width_in, fig_height_in]  # fig dims as a list
# ============================================
# ============================================

# ============================================
# SET LATEX
plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
plt.rcParams['font.weight'] = 400

# ============================================

# ============================================
# FONT SIZES
font_size = 25
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = font_size / 1.2
plt.rcParams['axes.linewidth'] = font_size / 22.
plt.rcParams['axes.titlesize'] = font_size / 1.1
plt.rcParams['legend.fontsize'] = font_size / 1.7
plt.rcParams['xtick.labelsize'] = font_size / 2
plt.rcParams['ytick.labelsize'] = font_size / 2


plt.rcParams['lines.linewidth'] = font_size / 1000


# ============================================
# ============================================
# TICKS

plt.rcParams['xtick.major.width'] = 0.13
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.width'] = 0.13
plt.rcParams['xtick.minor.size'] = 2.8
plt.rcParams['ytick.major.width'] = 0.13
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.width'] = 0.13
plt.rcParams['ytick.minor.size'] = 2.8

# fg = plt.figure(figsize=fig_dims)
fg = plt.figure(figsize=(8, 8))

ax1 = plt.subplot2grid((1, 1), (0, 0))  # , rowspan=2)
#ax2 = plt.subplot2grid((1,2), (0, 1))
#ax3 = plt.subplot2grid((2,2), (1, 1))

# ============================================
# set colors and lines
lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)
ax1.set_color_cycle(['#38197A', '#FC7D00', '#318B10', '#B61911'])
colorcycler = cycle(['#38197A', '#FC7D00', '#318B10', '#B61911'])

# ============================================


# ============================================

# Simplify paths by removing "invisible" points, useful for reducing
# file size when plotting a large number of points
plt.rcParams['path.simplify'] = False
# ============================================

# ============================================

# Have the legend only plot one point instead of two, turn off the
# frame, and reduce the space between the point and the label

plt.rcParams['axes.linewidth'] = 1.0
#plt.rc("lines", markeredgewidth=10)
# ============================================

# ============================================
# LEGEND
# ============================================
# OPTION
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.numpoints'] = 3
plt.rcParams['legend.handletextpad'] = 0.3


# ============================================
# ============================================
# ============================================
# ============================================
# ============================================
# REAL PLOT HERE
# ============================================
# data5 is> data6 so he probaly does + bedfore -
# plot1 = ax1.semilogx(data1[:,0], (data5[:,1]-data6[:,1])/pargaps[3]/data0[:,1],
# next(linecycler), linewidth=1, color='r',label=r'$C^{T}$')

# plot2 = plt.semilogx(data_fid[:,0] , 10.*np.nan_to_num((dats[:,1,]- dats[:,1,] )/data_fid[:,1]),linewidth=1, color='b',label=r'$C^{T}$')
key = 'scalar_spectral_index(1)'


if plot_type=='lmin':
    sigma_test = [ np.nan_to_num(np.sqrt((np.diag(np.linalg.inv(np.sum(full_fisher[:,:,i:], axis =2)))[fid.keys().index(key)]))) for i,j in enumerate(sigma1[:, 0])]

elif plot_type=='lmax':
    sigma_test = [ np.nan_to_num(np.sqrt((np.diag(np.linalg.inv(np.sum(full_fisher[:,:,:i], axis =2)))[fid.keys().index(key)]))) for i,j in enumerate(sigma1[:, 0])]
else:
    sys.exit('plot type wrong')


label = {}

label['massless_neutrinos'] = 'N_{eff}'
label['hubble'] = 'H_{0}'
label['mnu'] = 'M'
label['scalar_amp(1)'] = 'A_{s}'
label['scalar_spectral_index(1)'] = 'n_{s}'
label['omnuh2'] = r'\Omega_{\nu}h^{2}'
label['re_optical_depth'] = r'\tau'
label['ombh2'] = '\Omega_{b}h^{2}'
label['ombch2'] = '\Omega_{m}h^{2}'
label['omch2'] = '\Omega_{c}h^{2}'
label['helium_fraction'] = 'Y_{p}'
label['w'] = 'w'

plot = plt.plot
# plot2 = plot_type(sigma1[:, 0], sigma1[:, fid.keys().index(key) + 1], linewidth=1, label='Degeneracies')
plot2 = plot(sigma1[:80, 0], sigma_test[:80], linewidth=1)


# plot2 = plot_type(sigma_no[:, 0], sigma_no[:, fid.keys().index(key) + 1],
#                   linewidth=1, label='Perfect priors')

# plot2=plt.axvline(0.01,linewidth=1, ls= '--')


legend = ax1.legend()
ax1.legend(loc=0)

# ============================================
# FINALLY SAVE
ax1.set_ylabel(r'$ \sigma(' + label[key]+ ')$')

if plot_type =='lmax': ax1.set_xlabel(r'$\ell_{\mathrm{max}}$')

if plot_type =='lmin': ax1.set_xlabel(r'$\ell_{\mathrm{min}}$')

# ax1.set_xlim((0, 2000))
# ax1.set_ylim((10 ** -9, 10 ** -3))
ax1.minorticks_on()
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ============================================
# ============================================


# ============================================

plt.savefig(base_dir + 'data/{}/run{}/output/sigma_ell_{}_snow_mass_lmin={}_lmax={}_ndet={}_fsky={}.png'.format(data_type, str(run_idx), str(key), lmin, lmax, N_det, fsky), dpi=400, papertype='Letter',
            format='png', bbox_inches='tight')

plt.close()

