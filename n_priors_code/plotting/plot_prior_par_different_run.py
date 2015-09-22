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
# ============================================
# ============================================
# ============================================
# TO DO
#
# ============================================
# ============================================
# ============================================

no_lcdm_parameters =  ['massless_neutrinos', 'w', 'omnuh2']
plot_now = ['omnuh2']
excluded_parameters = list(set(no_lcdm_parameters) -set(plot_now))
# omnuh2
priors_to_plot = ['re_optical_depth']

# we need to cycle over this process



# READ DATA
# DEFINE YOUR FOLDER HERE
base_dir = '/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'
data_type = 'varying_lambda'
run_idx = 2
lmax = 4499
lmin = 100

lmin = [10,20,50,100]

fg = plt.figure(figsize=fig_dims)

ax1 = plt.subplot2grid((1, 1), (0, 0))
ax1.set_color_cycle(Set1_9.mpl_colors)
lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)


for n, lmin_iter in enumerate(lmin):
    pass

    # ======
    fid = pickle.load(open(base_dir + 'data/{}/run{}/fid_values.p'.format(data_type, str(run_idx)), "rb"))
    values = pickle.load(open(base_dir + 'data/{}/run{}/grid_values.p'.format(data_type, str(run_idx)), "rb"))
    par_gaps = pickle.load(open(base_dir + 'data/{}/run{}/par_gaps.p'.format(data_type, str(run_idx)), "rb"))
    fisher_mat = np.loadtxt(
        base_dir + 'data/{}/run{}/output/fisher_mat_joint_lmin={}_lmax={}.txt'.format(data_type, str(run_idx), lmin_iter, lmax))
    par_gaps, values, fid, fisher_mat = utils.exclude_parameters_from_fisher(
        excluded_parameters, par_gaps, values, fid, fisher_mat)
    plot_param = list(set(fid.keys()) -set(excluded_parameters))

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
    font_size = 14
    plt.rcParams['font.size'] = font_size
    plt.rcParams['axes.labelsize'] = font_size / 1.3
    plt.rcParams['axes.titlesize'] = font_size / 1.1
    plt.rcParams['legend.fontsize'] = font_size / 1.9
    plt.rcParams['xtick.labelsize'] = font_size / 2.
    plt.rcParams['ytick.labelsize'] = font_size / 2.


    plt.rcParams['lines.linewidth'] = font_size / 20.


    # ============================================
    # ============================================
    # TICKS

    plt.rcParams['xtick.major.width'] = 0.13 / 1.1
    plt.rcParams['xtick.major.size'] = 5 / 1.1
    plt.rcParams['xtick.minor.width'] = 0.13 / 1.1
    plt.rcParams['xtick.minor.size'] = 2.8 / 1.1
    plt.rcParams['ytick.major.width'] = 0.13 / 1.1
    plt.rcParams['ytick.major.size'] = 5 / 1.1
    plt.rcParams['ytick.minor.width'] = 0.13 / 1.1
    plt.rcParams['ytick.minor.size'] = 2.8 / 1.1


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

    plt.rcParams['axes.linewidth'] = font_size / 18.5
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
    label['w'] = 'w'
    fisher_inv = np.linalg.inv(fisher_mat)
    fisher_inplace = fisher_mat.copy()


    # CYCLE ON PARAMETERS (KEYS HERE)


    for y, key_y in enumerate(plot_now):
        print key_y

        sigma_just_CMB_y = (np.sqrt(fisher_inv[fid.keys().index(key_y), fid.keys().index(key_y)]))

        for i, key in enumerate(priors_to_plot):

            if key == key_y:
                continue
                # DO NOT TEST PRIOR ON ONE PARMS IN ITSELF; TRIVIAL
            # ERROR ON Y IN THE CMB S4
            # RELATIVE ERROR ON X IN THE CMB S4 we need relative cause this is how prior are used
            sigma_just_CMB_x = (np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / np.abs(fid[key]))
            # the abs value abbove is taken to deal with the negative fiducial value of w
            prior_value = np.linspace(sigma_just_CMB_x / 10., sigma_just_CMB_x * 4.5, 900)
            # print sigma_just_CMB_x
            # sigma_just_CMB_y_percent =sigma_just_CMB_y_percent /fid[key_y]
            normalize_x = prior_value / sigma_just_CMB_x  # prior respect to actual error
            # compute the new sigma on y given the prior
            new_sigma = utils.return_simgax_y_prior(fid, fisher_mat, key_y, key, prior_value)
            normalize_y = new_sigma / sigma_just_CMB_y  # make the new sigma y relative.
            # plot
            line_plot = ax1.plot(normalize_x, new_sigma, label=r'${0}={1:.1f}\% , lmin = {2:.1f}$'.format(
                str(label[key]), np.abs(sigma_just_CMB_x * 100.), lmin_iter ), linestyle=next(linecycler))

        # new_sigma_all = utils.return_simgax_all_prior(fid, fisher_mat, key_y)

        # plt.plot(normalize_x, new_sigma_all, label=r'All',
        #          linestyle=next(linecycler), linewidth=font_size / 10., alpha=0.6)



ax1.legend(loc=0)

ax1.minorticks_on()
ax1.set_ylim((0.8 * np.amin(new_sigma), 1.1 * np.amax(new_sigma)))
ax1.set_xlim((0.1, 3.1))
# ax1.set_title(r'$\sigma({0})={1:.1f}\%$'.format(str(label[key_y]), np.abs(sigma_just_CMB_y / fid[key_y] * 100.)))
ax1.set_ylabel(r'$\sigma(' + label[key_y] + r')$')
ax1.set_xlabel(r'$\rm{prior}/\sigma(x)_{\rm old}$')
y1, y2 = ax1.get_ylim()
ax2 = ax1.twinx()
minor_loc = ax1.yaxis.get_minor_locator()
ax2.minorticks_on()
ax2.set_ylim(y1 / np.abs(fid[key_y]) * 100., y2 / np.abs(fid[key_y]) * 100.)
new_ticks = ax2.get_yticks().tolist()
for i,tick in enumerate(ax2.get_yticks().tolist()):
    new_ticks[i] = str(tick)+ r'$\%$'
ax2.set_yticklabels(new_ticks)
        # ============================================
        # FINALLY SAVE
        # ============================================
        # ============================================

        # ============================================

plt.savefig(base_dir + 'data/{}/run{}/output/prior_{}_snow_mass_lmax={}.png'.format(data_type, str(run_idx), str(key_y), lmax), dpi=400, papertype='Letter',
            format='png', bbox_inches='tight')

plt.close()
