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


def years2sec(years):
    return years * 365 * 24. * 60. * 60.


def noise(Y, N_det, beam_arc, ell):
    import math
    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    # Y = 0.25  # 25%yeld
    # N_det = 10 ** 6  # 1 milion of detectors
    s = 350. * math.sqrt(20626. * 60. * 60.) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # s = 0.48 as in table from paper so it is ok.

    t = beam_arc / 60. / 180. * math.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    return (s * np.pi / 180. / 60.) ** 2 * np.exp(ell * (ell + 1) * t ** 2 / 8. / math.log(2))


# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA

# =============================

run_idx = 3


# READ PARAMS
dats = np.genfromtxt('../data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx))
fid = np.genfromtxt('../data/run{}/fiducial_pars.txt'.format(run_idx))
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('../data/run{}/grid_values.p'.format(run_idx), "rb"))
par_gaps = pickle.load(open('../data/run{}/par_gaps.p'.format(run_idx), "rb"))
fac = (dats[:, 0] * (dats[:, 0] + 1.) / 2. / np.pi) / (7.4311 * 10 ** 12)
noise_l = noise(0.25, 10e6, 2., dats[:, 0])
noise1 = np.loadtxt('../data/noise/wu_cdd_noise_6.txt')
noise2 = np.loadtxt('../data/noise/wu_cdd_noise_5.txt')
noise3 = np.loadtxt('../data/noise/wu_cdd_noise_4.txt')
noise_tt = np.loadtxt('../TT_lensing_noise.txt')
noise_ee = np.loadtxt('../EE_lensing_noise.txt')
noise_eb = np.loadtxt('../EB_lensing_noise.txt')
noise_te = np.loadtxt('../TE_lensing_noise.txt')
noise_tb = np.loadtxt('../TB_lensing_noise.txt')

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

plot_type = plt.semilogy
plot2 = plot_type(dats[:, 0], dats[:, 5], linewidth=1, label=r'$\phi$')
plot2 = plot_type(dats[:, 0], dats[:, 5] / dats[:, 5] * 1e-8, linewidth=1, label=r'Noise')
# plot2 = plt.loglog(noise3[:, 0], noise3[:, 1], linewidth=1, label=r'Wu 4')
# plot3 = plt.loglog(noise3[:, 0], noise2[:, 1], linewidth=1, label=r'Wu 5')
# plot4 = plt.loglog(noise3[:, 0], noise1[:, 1], linewidth=1, label=r'Wu 6')

plot4 = plot_type(noise_tt[:, 0], (noise_tt[:, 0]*(noise_tt[:, 0]+1))** 2 * noise_tt[:, 1] / 2 / np.pi, linewidth=1, label=r'quicklens TT')
plot4 = plot_type(noise_te[:, 0],  (noise_tt[:, 0]*(noise_tt[:, 0]+1))** 2 * noise_te[:, 1] / 2 / np.pi, linewidth=1, label=r'quicklens TE')
plot4 = plot_type(noise_ee[:, 0],  (noise_tt[:, 0]*(noise_tt[:, 0]+1))** 2 * noise_ee[:, 1] / 2 / np.pi, linewidth=1, label=r'quicklens EE')
plot4 = plot_type(noise_eb[:, 0],  (noise_tt[:, 0]*(noise_tt[:, 0]+1))** 2 * noise_eb[:, 1] / 2 / np.pi, linewidth=1, label=r'quicklens EB')
plot4 = plot_type(noise_tb[:, 0],  (noise_tt[:, 0]*(noise_tt[:, 0]+1))** 2 * noise_tb[:, 1] / 2 / np.pi, linewidth=1, label=r'quicklens TB')

legend = ax1.legend()
ax1.legend(loc=0)

# ============================================
# FINALLY SAVE
# print r'$ \frac{\partial C^{T}}{\partial ~'+label[key]+"$"
ax1.set_ylabel(r'C')
ax1.set_xlabel(r'$\ell$')
ax1.set_xlim((0, 2000))
ax1.set_ylim((10 ** -9, 10 ** -3))
ax1.minorticks_on()
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ============================================
# ============================================


# ============================================
plt.savefig('../../images/PS_phi_with_noise.pdf', dpi=400, papertype='Letter',
            format='pdf', transparent=True)

plt.clf()
