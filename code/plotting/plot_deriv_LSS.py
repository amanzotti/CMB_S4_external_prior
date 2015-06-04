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


# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA


run_idx = 3


# load fiducial data
dats = np.genfromtxt('../data/run{}/fiducial_matterpower.dat'.format(run_idx))
lenght = dats.shape[0]
print dats.shape[0]
# load fiducial parameters used
fid = pickle.load(open('../data/run{}/fid_values.p'.format(run_idx), "rb"))
print "fid ", fid
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('../data/run{}/grid_values.p'.format(run_idx), "rb"))
par_gaps = pickle.load(open('../data/run{}/par_gaps.p'.format(run_idx), "rb"))
n_values = np.size(values.keys())
# Load data for all parameters variations
for key, value in values.iteritems():
    for i in np.arange(0, 4):
        print dats.shape
        print key, values[key][i]
        filename = '../data/run{}/'.format(run_idx)
        filename += key + '_{:.13f}'.format(values[key][i]) + '_matterpower.dat'
        newdat = np.genfromtxt(filename)
        print newdat.shape[0]
        lenght_temp = newdat.shape[0]
        if lenght_temp<lenght:
            lenght = lenght_temp
            # print dats[:lenght_temp,:,:].shape[0],newdta.shape[0]
            try:
                 dats = np.dstack((dats[:lenght_temp,:], newdat))
            except:
                 dats = np.dstack((dats[:lenght_temp,:,:], newdat))

        elif lenght_temp>lenght:
            try:
                dats = np.dstack((dats, newdat[:lenght,:,:]))
            except:
                dats = np.dstack((dats, newdat[:lenght,:]))

        elif lenght_temp == lenght:
            dats = np.dstack((dats, newdat))



# ============================================
# ============================================
# ============================================


print "read an array of size=", dats.shape


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

label = {}

label['massless_neutrinos'] = 'N_{eff}'
label['hubble'] = 'H_{0}'
label['scalar_amp(1)'] = 'A_{s}'
label['scalar_spectral_index(1)'] = 'n_{s}'
label['omnuh2'] = '\Omega'
label['re_optical_depth'] = '\tau'


for i, key in enumerate(par_gaps.keys()):
    print key, i
    if key=='re_optical_depth':
        continue
    plt.clf()
    plot2 = plt.loglog(dats[:, 0, 0], np.abs(
        np.nan_to_num(dats[:, 1, 4 * i + 2] - dats[:, 1, 4 * i + 1]) / (par_gaps[key])/dats[:, 1, 0] * fid[key]), linewidth=1, color='k',rasterized=True)

    plot2 = plt.loglog(dats[:, 0, 0], np.abs(
        np.nan_to_num(dats[:, 1, 4 * i + 4] - dats[:, 1, 4 * i + 3]) / (par_gaps[key])/dats[:, 1, 0]* fid[key]), linewidth=1, color='g',rasterized=True)

    plot2 = plt.loglog(dats[:, 0, 0], np.abs(
        np.nan_to_num(dats[:, 1, 4 * i + 3] - dats[:, 1, 4 * i + 2]) / (2 * par_gaps[key])/dats[:, 1, 0]* fid[key]), linewidth=1, color='b',rasterized=True)

    plot3 = plt.loglog(dats[:, 0, 0], np.abs(
        np.nan_to_num(dats[:, 1, 4 * i + 4] - dats[:, 1, 4 * i + 1]) / (4 * par_gaps[key])/dats[:, 1, 0]* fid[key]), linewidth=1, color='r',rasterized=True)

    plot3 = plt.loglog(dats[:, 0, 0], np.abs(
        np.nan_to_num(-dats[:, 1, 4 * i + 4] + 8. * dats[:, 1, 4 * i + 3] - 8. * dats[:, 1, 4 * i + 2] + dats[:, 1, 4 * i + 1]) / (12 * par_gaps[key])/dats[:, 1, 0] * fid[key]), linewidth=1 ,marker='o')

    ax1.set_ylim((10**-7,10**1))
    fg.tight_layout(pad=0.1)
    ax1.set_ylabel(r'$ \frac{\partial P(k)}{\partial}$')
    ax1.set_xlabel(r'$k/h$')
    legend = ax1.legend()
    ax1.legend(loc=0)
    plt.savefig('../../images/test_PK_{}.pdf'.format(str(key)), dpi=400, papertype='Letter',
                format='pdf', transparent=True)
    ax1.minorticks_on()
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
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