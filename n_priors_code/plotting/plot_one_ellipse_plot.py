'''
Plot scipt to plot the derivatives of CMB specctra respect to all the parameters of the Fisher matrix.

Different lines corresponde to different gaps and techniques.



'''




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
import n_priors_code.utils


# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA


base_dir ='/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'

fid = pickle.load(open(base_dir+'data/{}/run3/fid_values.p'.format('test_fisher'), "rb"))
values = pickle.load(open(base_dir+'data/{}/run3/grid_values.p'.format('test_fisher'), "rb"))
par_gaps = pickle.load(open(base_dir+'data/{}/run3/par_gaps.p'.format('test_fisher'), "rb"))
index = {}

index['hubble'] = 3
index['scalar_spectral_index(1)'] = 1
index['scalar_amp(1)'] = 0
index['re_optical_depth'] = 2
index['ombh2'] = 4
index['ombch2'] = 5
index['mnu'] = 6

new_order = [index[fid.keys()[0]], index[fid.keys()[1]], index[fid.keys()[2]], index[fid.keys()[3]],
             index[fid.keys()[4]], index[fid.keys()[5]], index[fid.keys()[6]]]
# import zhen fisher
fisher_zhen = np.genfromtxt(base_dir+'data/test_fisher/F_CMBS4.dat')

reordered_fisher = np.zeros_like(fisher_zhen)
reordered_fisher = fisher_zhen[:, new_order][new_order]
fisher_zhen = reordered_fisher

my_fisher_1 = np.genfromtxt(base_dir+'output/fisher_mat_test_run1.txt')
my_fisher_2 = np.genfromtxt(base_dir+'output/fisher_mat_test_run2.txt')
my_fisher_3 = np.genfromtxt(base_dir+'output/fisher_mat_test_run3.txt')
my_fisher_4 = np.genfromtxt(base_dir+'output/fisher_mat_test_run4.txt')
my_fisher_5 = np.genfromtxt(base_dir+'output/fisher_mat_test_run5.txt')


my_fisher_inv_1 = (np.linalg.inv(my_fisher_1))
my_fisher_inv_2 = (np.linalg.inv(my_fisher_2))
my_fisher_inv_3 = (np.linalg.inv(my_fisher_3))
my_fisher_inv_4 = (np.linalg.inv(my_fisher_4))
my_fisher_inv_5 = (np.linalg.inv(my_fisher_5))

fisher_zhen_inv = (np.linalg.inv(reordered_fisher))


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
label['scalar_amp(1)'] = 'A_{s}'
label['scalar_spectral_index(1)'] = 'n_{s}'
label['omnuh2'] = r'\Omega_{\nu}'
label['re_optical_depth'] = r'~\tau'
label['ombh2'] = '\Omega_{b}'
label['omch2'] = '\Omega_{c}'
label['w'] = 'w'

# parameters =['mnu','ombch2']
parameters =['ombch2','mnu']

my_fisher_inv= fisher_zhen_inv

fg = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1, 1), (0, 0))



sigmax_squared= my_fisher_inv[fid.keys().index(parameters[0]),fid.keys().index(parameters[0])]
sigmay_squared= my_fisher_inv[fid.keys().index(parameters[1]),fid.keys().index(parameters[1])]
sigmaxy= my_fisher_inv[fid.keys().index(parameters[0]),fid.keys().index(parameters[1])]
print 'rho',np.sqrt(sigmay_squared)/fid[parameters[1]]*100.


a_fisher = (sigmax_squared+sigmay_squared)/2. + np.sqrt((sigmax_squared-sigmay_squared)**2/4.+sigmaxy**2)
b_fisher =  (sigmax_squared+sigmay_squared)/2. - np.sqrt((sigmax_squared-sigmay_squared)**2/4.+sigmaxy**2)
angle_fisher = np.degrees(np.arctan(2.*sigmaxy/(sigmax_squared-sigmay_squared)))/2.
print sigmax_squared,sigmay_squared
if sigmax_squared>sigmay_squared:

    one_sigma = Ellipse((fid[parameters[0]], fid[parameters[1]]), 2.*np.sqrt(np.amax([a_fisher,b_fisher]))*1.52, 2.*np.sqrt(np.amin([a_fisher,b_fisher]))*1.52,
                 angle=angle_fisher, linewidth=0.5, fill=True)
elif sigmax_squared<sigmay_squared:
    one_sigma = Ellipse((fid[parameters[0]], fid[parameters[1]]), 2.*np.sqrt(np.amin([a_fisher,b_fisher]))*1.52, 2.*np.sqrt(np.amax([a_fisher,b_fisher]))*1.52,
                 angle=angle_fisher, linewidth=0.5, fill=True)



# two_sigma = Ellipse((fid[key_j], fid[key_i]), np.sqrt(a_fisher)*2.48, np.sqrt(b_fisher)*2.48,
#              angle=angle_fisher, linewidth=0.5, fill=False)
ax1.add_patch(one_sigma)
plt.axhline(fid[parameters[1]]-1.52*np.sqrt(sigmay_squared),linewidth=1)
plt.axhline(fid[parameters[1]]+1.52*np.sqrt(sigmay_squared),linewidth=1)

xl =fid[parameters[0]]-2.52*np.sqrt(sigmax_squared)
xu = fid[parameters[0]]+2.52*np.sqrt(sigmax_squared)
yl = fid[parameters[1]]-2.52*np.sqrt(sigmay_squared)
yu = fid[parameters[1]]+2.52*np.sqrt(sigmay_squared)
ax1.set_ylim(yl,yu)
ax1.set_xlim(xl,xu)


fg.tight_layout()

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
plt.savefig('ellipse_zhen.pdf', dpi=400, papertype='Letter',
            format='pdf', transparent=True)
plt.close()
