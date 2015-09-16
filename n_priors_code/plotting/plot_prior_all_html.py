'''
This plotting script is used to plot the effect of external priors on all the parameters involved.

To put all the parameters in the same plot the prior is defined as the improvement from the CMB S4 results we get.

the same is true for the improvement on the y axis parameters



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
import n_priors_code.utils as utils
from palettable.colorbrewer.qualitative import Set1_9
from bokeh.plotting import output_file, show,save
from bokeh import mpl
# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

# READ DATA

# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying_lambda'
run_idx = 2
# ======
fid = pickle.load(open(base_dir + 'data/{}/run{}/fid_values.p'.format(data_type, str(run_idx)), "rb"))
values = pickle.load(open(base_dir + 'data/{}/run{}/grid_values.p'.format(data_type, str(run_idx)), "rb"))
par_gaps = pickle.load(open(base_dir + 'data/{}/run{}/par_gaps.p'.format(data_type, str(run_idx)), "rb"))
fisher_mat = np.loadtxt(base_dir + 'data/{}/run{}/output/fisher_mat.txt'.format(data_type, str(run_idx)))


# ============================================
# PLOTTING DEFINITION SKIP TO ~155
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
font_size = 24
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = font_size / 1.3
plt.rcParams['axes.titlesize'] = font_size / 1.1
plt.rcParams['legend.fontsize'] = font_size / 1.7
plt.rcParams['xtick.labelsize'] = font_size / 2.
plt.rcParams['ytick.labelsize'] = font_size / 2.


plt.rcParams['lines.linewidth'] = font_size / 10.


# ============================================
# ============================================
# TICKS

plt.rcParams['xtick.major.width'] = 0.13/1.1
plt.rcParams['xtick.major.size'] = 5/1.1
plt.rcParams['xtick.minor.width'] = 0.13/1.1
plt.rcParams['xtick.minor.size'] = 2.8/1.1
plt.rcParams['ytick.major.width'] = 0.13/1.1
plt.rcParams['ytick.major.size'] = 5/1.1
plt.rcParams['ytick.minor.width'] = 0.13/1.1
plt.rcParams['ytick.minor.size'] = 2.8/1.1


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
plt.rcParams['path.simplify'] = True
# ============================================

# ============================================

# Have the legend only plot one point instead of two, turn off the
# frame, and reduce the space between the point and the label

plt.rcParams['axes.linewidth'] = font_size/18.5
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

# DEFINE LABELS

label = {}

label['massless_neutrinos'] = 'N eff'
label['hubble'] = 'H0'
label['mnu'] = 'Mnu'
label['scalar_amp(1)'] = 'As'
label['scalar_spectral_index(1)'] = 'ns'
label['omnuh2'] = r'Omega nu'
label['re_optical_depth'] = r'tau'
label['ombh2'] = 'Omega b h2'
label['ombch2'] = 'Omega m h2'
label['omch2'] = 'Omega c h2'
label['w'] = 'w'
fisher_inv = np.linalg.inv(fisher_mat)


# CYCLE ON PARAMETERS (KEYS HERE)
for y, key_y in enumerate(par_gaps.keys()):
    print key_y
    fg = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1, 1), (0, 0))
    ax1.set_color_cycle(Set1_9.mpl_colors)
    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)

    for i, key in enumerate(par_gaps.keys()):

        if key == key_y:
            continue
            # DO NOT TEST PRIOR ON ONE PARMS IN ITSELF; TRIVIAL

        # ERROR ON Y IN THE CMB S4
        sigma_just_CMB_y = (np.sqrt(fisher_inv[fid.keys().index(key_y), fid.keys().index(key_y)]))
        # RELATIVE ERROR ON X IN THE CMB S4 we need relative cause this is how prior are used
        sigma_just_CMB_x = (np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / fid[key])

        prior_value = np.linspace(sigma_just_CMB_x/10., sigma_just_CMB_x*4.5, 900)

        # sigma_just_CMB_y_percent =sigma_just_CMB_y_percent /fid[key_y]
        normalize_x = prior_value / sigma_just_CMB_x  # prior respect to actual error
        # compute the new sigma on y given the prior
        new_sigma = utils.return_simgax_y_prior(fid, fisher_mat, key_y, key, prior_value)
        normalize_y = new_sigma / sigma_just_CMB_y  # make the new sigma y relative.
        # plot
        plt.plot(normalize_x, new_sigma, label=r'{0}={1:.1f} %'.format(
            str(label[key]), np.abs(sigma_just_CMB_x * 100.)), rasterized=True, linestyle=next(linecycler))

    legend = ax1.legend()
    ax1.legend(loc=0)
    ax1.minorticks_on()
    ax1.set_ylim((0, np.max(new_sigma)))
    ax1.set_xlim((0.1, 3.1))
    # ax1.set_title(r'$\sigma({0})={1:.1f}\%$'.format(str(label[key_y]), np.abs(sigma_just_CMB_y / fid[key_y] * 100.)))
    ax1.set_ylabel(r'sigma('+ label[key_y] + r')')
    ax1.set_xlabel(r'prior/sigma(x) old')
    plt.title(r'Prior results on ' + label[key_y])
    # ============================================
    # FINALLY SAVE
    # ============================================
    # ============================================

    # ============================================

    output_file("subplots.html")
    save(mpl.to_bokeh())
    plt.savefig(base_dir + 'data/{}/run{}/output/prior_{}_snow_mass.pdf'.format(data_type, str(run_idx), str(key_y)), dpi=400, papertype='Letter',
                format='pdf', transparent=True, bbox_inches='tight')
    plt.close()
