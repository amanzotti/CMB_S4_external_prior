
#!/usr/bin/env python


'''

FAST VERSION

Here we produce the Planck Fisher matrix we will lately add to our fisher.


TODO:
low importance: ini file. Think about PCA to understand what is more important for N_eff

# DATA FROM CAMB

we use the _lenspotentialcls so

 l CTT CEE CBB CTE Cdd CdT CdE
 0  1  2    3   4   5   6  7

  CX are l(l+1)Cl/2pi and Cdd=[l(l+1)]^2 Clphi/2pi, CdT=[l(l+1)]^(3/2) ClphiT/2pi, CdE=[l(l+1)]^(3/2)ClphiE/2pi


The matrix is going to run an 4 parameters:

CONVENTIONS:
alphabetical in CAMB description. Use the ordered dict tools to explore it.
hubble,massless_neutrinos,omnuh2,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)


GOAL

1) reproduce Wu paper and beyond/ DONE. Cross checked with Zhen Pan UC Davis

'''
__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"

import numpy as np
import utils
import pickle
import sys
import collections


# Auxiliary functions.
# ===================

def years2sec(years):
    ''' years to sec '''
    return years * 365 * 24. * 60. * 60.


def fsky2arcmin(fsky):
    '''convert fsky in fraction of unity to arcmin^2'''
    # 41253 square degrees in all sky
    return 41253. * fsky * 60. * 60.
# ===================


def calc_deriv_vectorial(fisher_index, dats, pargaps, values, order=5):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE
    '''
    if order == 5 and np.shape(values[values.keys()[0]])[0] == 4:
        ci = (-calc_c_general(dats, fisher_index * 4 + 4) + 8. * calc_c_general(dats, fisher_index * 4 + 3) - 8. *
              calc_c_general(dats, fisher_index * 4 + 2) + calc_c_general(dats, fisher_index * 4 + 1)) / (12. * pargaps[values.keys()[fisher_index]])

    elif order == 7 and np.shape(values[values.keys()[0]])[0] == 6:
        ci = (calc_c_general(dats, fisher_index * 6 + 6) - 9. * calc_c_general(dats, fisher_index * 6 + 5) + 45. *
              calc_c_general(dats, fisher_index * 6 + 4) - 45. * calc_c_general(dats, fisher_index * 6 + 3) + 9. * calc_c_general(dats, fisher_index * 6 + 2) - calc_c_general(dats, fisher_index * 6 + 1)) / (60. * pargaps[values.keys()[fisher_index]])

    elif order == 9 and np.shape(values[values.keys()[0]])[0] == 8:
        ci = (-3. * calc_c_general(dats, fisher_index * 8 + 8) + 32. * calc_c_general(dats, fisher_index * 8 + 7) - 168. *
              calc_c_general(dats, fisher_index * 8 + 6) + 672. * calc_c_general(dats, fisher_index * 8 + 5) - 672. * calc_c_general(dats, fisher_index * 8 + 4) + 168. * calc_c_general(dats, fisher_index * 8 + 3) - 32. * calc_c_general(dats, fisher_index * 8 + 2) + 3. * calc_c_general(dats, fisher_index * 8 + 1)) / (840. * pargaps[values.keys()[fisher_index]])
    else:
        sys.exit()
    # if eight points
    # ci = (-calc_c_general(dats, fisher_index * 4 + 4) + 8. * calc_c_general(dats, fisher_index * 4 + 3) - 8. *
    # calc_c_general(dats, fisher_index * 4 + 2) + calc_c_general(dats,
    # fisher_index * 4 + 1)) / (12. * pargaps[values.keys()[fisher_index]])
    return ci


def calc_c_fiducial(data):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    THE LEVEL OF NOISE IN THE CHANNEL IS TAKEN FROM ALLISON ET AL (DUNKLEY)

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE


    '''

    # THIS is a planck like experiment different from the usual formalism
    # to combine different channel see eq 6 of 1403.5271v1

    global Y, sec_of_obs, arcmin_from_fsky, lmax_index
    ell = data[lmin_index:lmax_index, 0, 0]
    s = {}
    s_pol = {}
    t = {}
    N = {}
    N_pol = {}
    fac = (ell * (ell + 1.) / 2. / np.pi) / (7.4311 * 10 ** 12)
    N_tot = np.zeros_like(ell)
    N_pol_tot = np.zeros_like(ell)


    # 30,44,70,100,143,217,353
    s['30'] = 145  # this is supposed to be in muk arcmin to fit in the later used formulas
    s['44'] = 149
    s['70'] = 137
    s['100'] = 65
    s['143'] = 43
    s['217'] = 66
    s['353'] = 200

    t['30'] = 33. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['44'] = 23. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['70'] = 14. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['100'] = 10. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['143'] = 7. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['217'] = 5. / 60. / 180. * np.pi  # 2arcmin to rads beam
    t['353'] = 5. / 60. / 180. * np.pi  # 2arcmin to rads bea


# NOISE IN DIFFERENT CHANNELS IS ADDED IN QUADRATURE.

    for nu in ['30', '44', '70', '100', '143', '217', '353']:

        # Final CMB noise definition
        N[nu] = (s[nu] * np.pi / 180. / 60.) ** 2 * np.exp(ell * (ell + 1.) * t[nu] ** 2 / 8. / np.log(2))
        # print 'Nnu',N[nu]
        N_tot += 1 / N[nu]

    N_tot = 1. / N_tot

    # 30,44,70,100,143,217,353
    s_pol['30'] = 0.  # this is supposed to be in muk arcmin to fit in the later used formulas
    s_pol['44'] = 0.
    s_pol['70'] = 450.
    s_pol['100'] = 103.
    s_pol['143'] = 81.
    s_pol['217'] = 134.
    s_pol['353'] = 406.

#Use this for Planck below for planck pol
    for nu in ['70']:
    # for nu in ['70', '100', '143', '217', '353']:

        # Final CMB noise definition
        N_pol[nu] = (s_pol[nu] * np.pi / 180. / 60.) ** 2 * np.exp(ell * (ell + 1.) * t[nu] ** 2 / 8. / np.log(2))
        # print 'Nnu',N[nu]
        N_pol_tot += 1 / N_pol[nu]

    N_pol_tot = 1. / N_pol_tot

    # np.save('noise.npy',N_tot*fac)

    # print 'N_tot', N_tot
    # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    # C_inv = [[data[:lmax_index, 1, 0] + fac * N
    #     , data[:lmax_index, 4, 0], data[:lmax_index, 6, 0]],
    #                  [data[:lmax_index, 4, 0], data[:lmax_index, 2, 0] + fac * N * 2.,               0.],
    #                  [data[:lmax_index, 6, 0],          0.,         data[:lmax_index, 5, 0] + N_phi_l[:lmax_index, 1]]]

    return np.array([[data[lmin_index:lmax_index, 1, 0] + fac * N_tot, data[lmin_index:lmax_index, 4, 0]],

                     [data[lmin_index:lmax_index, 4, 0], data[lmin_index:lmax_index, 2, 0] +
                         fac * N_pol_tot]
                     ])


def calc_c_general(data, parabin):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE

    '''

    return np.array([[data[lmin_index:lmax_index, 1, parabin], data[lmin_index:lmax_index, 4, parabin]],
                     [data[lmin_index:lmax_index, 4, parabin], data[lmin_index:lmax_index,
                                                                    2, parabin]]
                     ])


def C(iell, ell, parbin, data):
    '''

    Here is used to compute derivatives so we do not need noise

    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE


    '''
    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    return np.array([[data[iell, 1, parbin], data[iell, 4, parbin], data[iell, 6, parbin]],
                     [data[iell, 4, parbin], data[iell, 2, parbin],               0.]
                     ]
                    )

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================

# \THIS IS  A PLANCK \TIPE EXPERIMENT


l_t_max = 2500  # this is the multipole you want to cut the temperature Cl at, to simulate the effect of foreground
lmax = 2500
lmin = 2
N_det = 'Planck'
# No lensing
# N_phi_l = np.loadtxt('data/noise/wu_cdd_noise_6.txt')
data_folder = 'varying_all/run7'
output_folder = ''
fsky = 0.44
lensed = False
# exclude = ['massless_neutrinos', 'w','omnuh2']  # None
# exclude = ['helium_fraction', 'scalar_nrun(1)', 'omk', 'wa','massless_neutrinos']  # None
exclude = None

# =============================
# DERIVED
arcmin_from_fsky = fsky2arcmin(fsky)
sec_of_obs = years2sec(5)
Y = 0.25  # 25% yeld
# ===================
header = 'Joint fisher CMB T E + phi lensing use Planck noise\n'
header += 'lmax={} \n lmin={} \n l_t_max={} \n fsky={} \n lensed={} \n data_folder={} \n N_det={} \n'.format(
    lmax, lmin, l_t_max, fsky, lensed, data_folder, N_det)


# READ PARAMS
# load fiducial data
# load fiducial parameters used
fid = pickle.load(open('data/{}/fid_values.p'.format(data_folder), "rb"))

print "fid ", fid

# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/{}/grid_values.p'.format(data_folder), "rb"))
par_gaps = pickle.load(open('data/{}/par_gaps.p'.format(data_folder), "rb"))
# order = 5

# select values to change gaps in derivative
new_value = {}

# ===============
# order = 5

# i = 0  # or 1 or 2 to select different gaps for 5 point formula
# for key in values.keys():
#     new_value[key] = [values[key][i], values[key][i + 1], values[key][-i - 2], values[key][-i - 1]]

#     if i != 2:
#         par_gaps[key] = par_gaps[key] * (4 - 2 * i)

# new_value = collections.OrderedDict(sorted(new_value.items(), key=lambda t: t[0]))
# values = new_value
# ===============

# ===============

# use different order formula same gap
# use different order formula same gap
order = 5
# step = np.array([-8,-7,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8])

step = np.array([-2, -1, 1, 2])

for key in values.keys():
    # if key == 'massless_neutrinos':
    #     new_value[key] = par_gaps[key] * step + fid[key]
    #     continue
    # if key=='omnuh2':
    #   new_value[key] = par_gaps[key] * np.array([-2,-1,1,2]) + fid[key]
    #   print new_value[key]
    #   continue

    new_value[key] = par_gaps[key] * step + fid[key]
for key in values.keys():
    par_gaps[key] = np.abs(new_value[key][0] - new_value[key][1])
new_value = collections.OrderedDict(sorted(new_value.items(), key=lambda t: t[0]))
values = new_value
print values, par_gaps

# # ===============

par_gaps, values, fid = utils.exclude_parameters(exclude, par_gaps, values, fid)

print 'loading files'
dats = utils.load_data(data_folder, values, lensed)
n_values = np.size(values.keys())
lmax_index = np.where(dats[:, 0, 0] == lmax)[0][0]
ltmax_index = np.where(dats[:, 0, 0] == l_t_max)[0][0]
lmin_index = np.where(dats[:, 0, 0] == lmin)[0][0]

# cut Cl^T at ells bigger than l_t_max
dats[ltmax_index:, 1, 1:] = 0.
# phi_T has oscillations in it.
# dats[900:, 6, 0:] = 0.

# creating the n_values by n_values matrix
fisher = np.zeros((n_values, n_values))
fisher_save = np.zeros((n_values, n_values, np.size(dats[lmin_index:lmax_index, 0, 0])))

fisher_inv = np.zeros((n_values, n_values))

no_marginalized_ell = np.zeros((np.size(dats[lmin_index:lmax_index, 0, 0]), n_values))
marginalized_ell = np.zeros((np.size(dats[lmin_index:lmax_index, 0, 0]), n_values))


print 'fisher_size', fisher.shape

# generate C for fiducial at all ell
C_inv_array = calc_c_fiducial(dats)

derivatives = np.ndarray((2, 2, np.size(dats[lmin_index:lmax_index, 0, 0]), n_values), dtype='float64')
for i in range(0, n_values):
    # computing derivatives.
    # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
                  # ---------------------------------
                              #   12h
    derivatives[:, :, :, i] = calc_deriv_vectorial(i, dats, par_gaps, values, order)[:, :, :]

# you may want not to compute the derivatives but to impose a known degeneracy.

# ratio = derivatives[:2, :2, 300:, fid.keys().index('w')]/derivatives[:2, :2, 300:, fid.keys().index('hubble')]
derivatives[:2, :2, 1000:, fid.keys().index('w')] = 1 / 3.41 * derivatives[:2, :2, 1000:, fid.keys().index('hubble')]


for iell, ell in enumerate(dats[lmin_index:lmax_index, 0, 0]):

    #  filling it the matrix l goes from l_min =2 to l_max =5000

    # c0 = np.zeros((3, 3))
    # c0 = C(iell, ell, 0, dats)  # 3x3 matrix in the fiducial cosmology
    # this is the covariance matrix of the data. So in this case we have C^T C^E C^phi
    c0 = C_inv_array[:, :, iell]

    cinv = np.linalg.inv(c0)

    for i in range(0, n_values):

        for j in range(0, n_values):
            # computing derivatives.
            # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
                  # ---------------------------------
                              #   12h

            ci = derivatives[:, :, iell, i]

            cj = derivatives[:, :, iell, j]

            # Eq 4.
            # tot = np.dot(np.dot(np.dot(cinv, ci),  cinv), cj)
            trace = np.sum(np.dot(np.dot(cinv, ci),  cinv) * cj.T)

            # See http://arxiv.org/pdf/1509.07471v1.pdf tabel IV THIS IS what to do if to be combined with S4 experiment.
            if ell < 50:
                fsky = 0.44
            else:
                fsky = 0.2
            fisher[i, j] += (2. * ell + 1.) / 2. * fsky * trace
            fisher_save[i, j, iell] = (2. * ell + 1.) / 2. * fsky * trace
            # print np.sum(fisher_save[:,:,:], axis =2)[0,0],fisher[0,0]

    no_marginalized_ell[iell, :] = 1. / np.sqrt(np.diag(fisher))
    fisher_inv = np.linalg.inv(fisher)
    marginalized_ell[iell, :] = np.sqrt(np.diag(fisher_inv))

print 'lmax =', ell


BAO_fisher = np.loadtxt('/home/manzotti/n_eff-dependence-on-prior/n_priors_code/data/fisher_mat_BAO.txt')
BAO_fisher_DESI = np.loadtxt('/home/manzotti/n_eff-dependence-on-prior/n_priors_code/data/fisher_mat_BAO_DESI.txt')

print 'ADDING BAO'
fisher += BAO_fisher


print np.shape(fisher)
np.savetxt('/home/manzotti/n_eff-dependence-on-prior/n_priors_code/data/fisher_mat_Planck.txt', fisher, header=header)

print 'fisher=', fisher

no_lcdm_parameters = ['massless_neutrinos', 'w', 'omnuh2']
plot_now = ['']
excluded_parameters = list(set(no_lcdm_parameters) - set(plot_now))

par_gaps, values, fid, fisher_single = utils.exclude_parameters_from_fisher(
    excluded_parameters, par_gaps, values, fid, fisher)

utils.print_resume_stat(fisher_single, fid)
