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

data0 = np.genfromtxt('../data/dat0.txt', dtype=float)
data1 = np.genfromtxt('../data/dat1.txt', dtype=float)
data2 = np.genfromtxt('../data/dat2.txt', dtype=float)
data3 = np.genfromtxt('../data/dat3.txt', dtype=float)
data4 = np.genfromtxt('../data/dat4.txt', dtype=float)
data5 = np.genfromtxt('../data/dat5.txt', dtype=float)
data6 = np.genfromtxt('../data/dat6.txt', dtype=float)
data7 = np.genfromtxt('../data/dat7.txt', dtype=float)
data8 = np.genfromtxt('../data/dat8.txt', dtype=float)

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
# ============================================
# ============================================
# ============================================


print "read an array of size=", data1.shape


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


plot2 = plt.loglog(dats[:, 0], dats[:, 1] / fac, linewidth=1,label=r'TT')
plot2 = plt.loglog(dats[:, 0], dats[:, 2] / fac, linewidth=1,label=r'EE')
plot2 = plt.loglog(dats[:, 0], dats[:, 3] / fac, linewidth=1,label=r'BB')
plot2 = plt.loglog(dats[:, 0], np.abs(dats[:, 4] / fac), linewidth=1,label=r'TE')
plot2 = plt.loglog(dats[:, 0], noise_l[:], linewidth=1, label=r'Noise T')
plot2 = plt.loglog(dats[:, 0], noise_l[:]*2., linewidth=1, label=r'Noise Pol')


legend = ax1.legend()
ax1.legend(loc=0)

# ============================================
# FINALLY SAVE
# print r'$ \frac{\partial C^{T}}{\partial ~'+label[key]+"$"
ax1.set_ylabel(r'$ $')
ax1.set_xlabel(r'$\ell$')
# ax1.set_xlim((0, 3))
# ax1.set_ylim((10**-16,10**-10))
ax1.minorticks_on()
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ============================================
# ============================================


# ============================================
plt.savefig('../../images/PS_with_noise.pdf', dpi=400, papertype='Letter',
            format='pdf', transparent=True)

plt.clf()
