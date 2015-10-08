'''
Plot scipt to plot the derivatives of CMB specctra respect to all the parameter of the Fisher matrix.

Different lines corresponde to different gaps and techniques.



'''


__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"


import numpy as np
import pickle as pickle
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
from scipy.optimize import fsolve
from scipy.integrate import quad
# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA


no_lcdm_parameter = ['massless_neutrinos', 'w', 'omnuh2']
plot_now = ['omnuh2']
excluded_parameter = list(set(no_lcdm_parameter) - set(plot_now))
# omnuh2

# READ DATA
# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying_lambda'
run_idx = 2
lmax = 4499
lmin = 50
# ======
fid = pickle.load(open(base_dir + 'data/{}/run{}/fid_values.p'.format(data_type, str(run_idx)), "rb"))
values = pickle.load(open(base_dir + 'data/{}/run{}/grid_values.p'.format(data_type, str(run_idx)), "rb"))
par_gaps = pickle.load(open(base_dir + 'data/{}/run{}/par_gaps.p'.format(data_type, str(run_idx)), "rb"))
fisher_mat = np.loadtxt(
    base_dir + 'data/{}/run{}/output/fisher_mat_joint_lmin={}_lmax={}.txt'.format(data_type, str(run_idx), lmin, lmax))

par_gaps, values, fid, fisher_mat = utils.exclude_parameters_from_fisher(
    excluded_parameter, par_gaps, values, fid, fisher_mat)
# plot_param = list(set(fid.keys()) -set(excluded_parameter))
print excluded_parameter
print fid
my_fisher_inv = np.linalg.inv(fisher_mat)




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

plt.rcParams['xtick.major.width'] = 0.13 / 2.
plt.rcParams['xtick.major.size'] = 5 / 2.
plt.rcParams['xtick.minor.width'] = 0.13 / 2.
plt.rcParams['xtick.minor.size'] = 2.8 / 2.
plt.rcParams['ytick.major.width'] = 0.13 / 2.
plt.rcParams['ytick.major.size'] = 5 / 2.
plt.rcParams['ytick.minor.width'] = 0.13 / 2.
plt.rcParams['ytick.minor.size'] = 2.8 / 2.


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

# parameter =['mnu','ombch2']
parameter = ['omnuh2']

fg = plt.figure(figsize=(10, 10))
ax1 = plt.subplot2grid((1, 1), (0, 0))


sigma_squared = my_fisher_inv[fid.keys().index(parameter[0]), fid.keys().index(parameter[0])]
center = fid[parameter[0]]
x= np.linspace(center-2.52 * np.sqrt(sigma_squared), center + 2.52 * np.sqrt(sigma_squared) )
y = np.exp(-(x-center)*(x-center)/(2*sigma_squared))/np.sqrt(2*np.pi*sigma_squared)
plt.plot(x,y,linewidth=1)
xl = fid[parameter[0]] - 2.52 * np.sqrt(sigma_squared)
xu = fid[parameter[0]] + 2.52 * np.sqrt(sigma_squared)
ax1.set_xlim(xl, xu)


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
plt.savefig('../data/{}/run{}/output/likelihood_{}.pdf'.format(data_type,run_idx,str(parameter[0])), dpi=400, papertype='Letter',
                format='pdf', transparent=True)
plt.close()
