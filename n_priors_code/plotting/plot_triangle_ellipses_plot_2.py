#!/usr/bin/env python
'''
This script plot the usual tringle ellispes plot for the given dataset.
Remember this is a Fisher matrix code so everything is approximated to be gaussian.

This was written to compare our results with the one given by Zhen
'''


__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"




import numpy as np
import pickle as pickle
import utils
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special
import scipy.interpolate
from itertools import cycle
#from matplotlib.colors import colorConverter
from matplotlib.ticker import MaxNLocator  # needed for tick
import sys
from matplotlib.patches import Ellipse
import n_priors_code.utils as utils
from palettable.colorbrewer.qualitative import Set1_9



# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA


no_lcdm_parameters = ['massless_neutrinos', 'w', 'omnuh2', 'helium_fraction','wa','omk','scalar_nrun(1)']
plot_now = ['omnuh2']
excluded_parameters = list(set(no_lcdm_parameters) - set(plot_now))
# omnuh2

# READ DATA
# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying_all'
run_idx = 4
lmax = 4499
lmin = 4
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



my_fisher_inv = np.linalg.inv(fisher_mat)
index = {}

index['hubble'] = 3
index['scalar_spectral_index(1)'] = 1
index['scalar_amp(1)'] = 0
index['re_optical_depth'] = 2
index['ombh2'] = 4
index['ombch2'] = 5
index['mnu'] = 6


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
font_size = 20
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = font_size / 1.2
plt.rcParams['axes.linewidth'] = font_size / 38.
plt.rcParams['axes.titlesize'] = font_size / 1.1
plt.rcParams['legend.fontsize'] = font_size / 1.7
plt.rcParams['xtick.labelsize'] = font_size / 2
plt.rcParams['ytick.labelsize'] = font_size / 2


plt.rcParams['lines.linewidth'] = font_size / 1000


# ============================================
# ============================================
# TICKS

plt.rcParams['xtick.major.width'] = 0.13/2.
plt.rcParams['xtick.major.size'] = 5/2.
plt.rcParams['xtick.minor.width'] = 0.13/2.
plt.rcParams['xtick.minor.size'] = 2.8/2.
plt.rcParams['ytick.major.width'] = 0.13/2.
plt.rcParams['ytick.major.size'] = 5/2.
plt.rcParams['ytick.minor.width'] = 0.13/2.
plt.rcParams['ytick.minor.size'] = 2.8/2.


#ax2 = plt.subplot2grid((1,2), (0, 1))
#ax3 = plt.subplot2grid((2,2), (1, 1))

# ============================================
# set colors and lines
lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)
# ax1.set_color_cycle(['#38197A', '#FC7D00', '#318B10', '#B61911'])
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
# ============================================
# ============================================
# ============================================
# ============================================
# REAL PLOT HERE
# ============================================


label = {}

label['massless_neutrinos'] = 'N_{eff}'
label['hubble'] = 'H_{0}'
label['mnu'] = 'M'
label['scalar_amp(1)'] = 'A_{s}'
label['scalar_spectral_index(1)'] = 'n_{s}'
label['omnuh2'] = r'\Omega_{\nu}'
label['re_optical_depth'] = r'~\tau'
label['ombh2'] = '\Omega_{b}h^{2}'
label['ombch2'] = '\Omega_{m}h^{2}'
label['omch2'] = '\Omega_{c}h^{2}'
label['w'] = 'w'



fg = plt.figure(figsize=(10,10))


size = np.size(par_gaps.keys())
for i, key_i in enumerate(par_gaps.keys()):

    for j, key_j in enumerate(par_gaps.keys()):

        if(j>i) :
            print(i,j)
            print key_i,key_j
            ax1 = plt.subplot2grid((size-1, size-1), (i,j-1))

            sigmax_squared= my_fisher_inv[fid.keys().index(key_j),fid.keys().index(key_j)]
            sigmay_squared= my_fisher_inv[fid.keys().index(key_i),fid.keys().index(key_i)]
            sigmaxy= my_fisher_inv[fid.keys().index(key_j),fid.keys().index(key_i)]
            print 'rho',sigmax_squared,sigmay_squared


            a_fisher = (sigmax_squared+sigmay_squared)/2. + np.sqrt((sigmax_squared-sigmay_squared)**2/4.+sigmaxy**2)
            b_fisher =  (sigmax_squared+sigmay_squared)/2. - np.sqrt((sigmax_squared-sigmay_squared)**2/4.+sigmaxy**2)
            angle_fisher = np.degrees(np.arctan(2.*sigmaxy/(sigmax_squared-sigmay_squared)))/2.

            if sigmax_squared>sigmay_squared:

                one_sigma = Ellipse((fid[key_j], fid[key_i]), 2.*np.sqrt(np.amax([a_fisher,b_fisher]))*1.52, 2.*np.sqrt(np.amin([a_fisher,b_fisher]))*1.52,
                             angle=angle_fisher, linewidth=0.8, fill=False)
            elif sigmax_squared<sigmay_squared:
                one_sigma = Ellipse((fid[key_j], fid[key_i]), 2.*np.sqrt(np.amin([a_fisher,b_fisher]))*1.52, 2.*np.sqrt(np.amax([a_fisher,b_fisher]))*1.52,
                             angle=angle_fisher, linewidth=0.8, fill=False)


            # two_sigma = Ellipse((fid[key_j], fid[key_i]), np.sqrt(a_fisher)*2.48, np.sqrt(b_fisher)*2.48,
            #              angle=angle_fisher, linewidth=0.5, fill=False)
            ax1.add_patch(one_sigma)
            xl =fid[key_j]-2.52*np.sqrt(sigmax_squared)
            xu = fid[key_j]+2.52*np.sqrt(sigmax_squared)
            yl = fid[key_i]-2.52*np.sqrt(sigmay_squared)
            yu = fid[key_i]+2.52*np.sqrt(sigmay_squared)
            ax1.set_ylim(yl,yu)
            ax1.set_xlim(xl,xu)
            plt.locator_params(nbins=4)
            if j==i+1:
                ax1.set_xlabel(r'${}$'.format(label[key_j]))
                ax1.set_ylabel(r'${}$'.format(label[key_i]))
            else:

                ax1.xaxis.set_ticklabels([])
                ax1.yaxis.set_ticklabels([])



# fg.tight_layout()

# ============================================
# FINALLY SAVE
# ============================================
# ============================================

# ============================================
# LEGEND
# ============================================
# OPTION
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.numpoints'] = 3
plt.rcParams['legend.handletextpad'] = 0.3
# ============================================
# plt.savefig('../../images/trinagle.pdf', dpi=400, papertype='Letter',
#             format='pdf', transparent=True)
plt.savefig('triangle_all_prior_{}_snow_mass_lmin={}_lmax={}_ndet={}_fsky={}_no_planck.pdf'.format(str('joint_real'), lmin, lmax, N_det, fsky), dpi=400, papertype='Letter',format='pdf', transparent=True)
plt.close()
