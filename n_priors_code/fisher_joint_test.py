
#!/usr/bin/env python


'''

FAST VERSION
Fisher code to test prior effect on N effective neutrino estimates

TODO:
HIGH: lensing noise (almost there: Komatsu and quicklens.)

low importance: ini file. Think about PCA to understand what is more important for N_eff


create classes or at least external utilities. Maybe a better way to interface with cosmosis

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
import math
import matplotlib.pyplot as plt
import utils
import pickle
import sys
import collections

# util.nl(noise_uK_arcmin, fwhm_arcmin, lmax)


def years2sec(years):
    ''' years to sec '''
    return years * 365 * 24. * 60. * 60.


def fsky2arcmin(fsky):
    '''convert fsky in fraction of unity to arcmin^2'''
    # 41253 square degrees in all sky
    return 41253. * fsky * 60. * 60.


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

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE


    '''
    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    # Y = 0.25  # 25% yeld
    # N_det = 10 ** 6  # 1 milion of detectors
    # these are taken from global
    global Y, sec_of_obs, arcmin_from_fsky, lmax_index
    ell = data[lmin_index:lmax_index, 0, 0]
    s = 350. * np.sqrt(arcmin_from_fsky) / np.sqrt(N_det * Y * sec_of_obs)  # half sky in arcmin^2
    t = 1. / 60. / 180. * np.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / np.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    N = (s * np.pi / 180. / 60.) ** 2 * np.exp(ell * (ell + 1.) * t ** 2 / 8. / np.log(2))
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

    return np.array([[data[lmin_index:lmax_index, 1, 0] + fac * N, data[lmin_index:lmax_index, 4, 0], data[lmin_index:lmax_index, 6, 0]],

                     [data[lmin_index:lmax_index, 4, 0], data[lmin_index:lmax_index, 2, 0] +
                         fac * N * 2.,  data[lmin_index:lmax_index, 2, 0] * 0.],

                     [data[lmin_index:lmax_index, 6, 0], data[lmin_index:lmax_index, 2, 0] * 0.,
                      data[lmin_index:lmax_index, 5, 0] + N_phi_l[lmin_index:lmax_index, 1]]

                     ])


def calc_c_general(data, parabin):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE

    '''
    return np.array([[data[lmin_index:lmax_index, 1, parabin], data[lmin_index:lmax_index, 4, parabin], data[lmin_index:lmax_index, 6, parabin]],
                     [data[lmin_index:lmax_index, 4, parabin], data[lmin_index:lmax_index,
                                                                    2, parabin],  data[lmin_index:lmax_index, 2, parabin] * 0.],
                     [data[lmin_index:lmax_index, 6, parabin], data[lmin_index:lmax_index, 2, parabin] * 0.,
                      data[lmin_index:lmax_index, 5, parabin]]
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
                     [data[iell, 4, parbin], data[iell, 2, parbin],               0.],
                     [data[iell, 6, parbin],          0.,         data[iell, 5, parbin]]]
                    )

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================
l_t_max = 3000  # this is the multipole you want to cut the temperature Cl at, to simulate the effect of foregrounds
lmax = 4499
lmin = 4
N_det = 10 ** 6
N_phi_l = np.loadtxt('data/noise/wu_cdd_noise_6.txt')
data_folder = 'varying_all/run4'
output_folder = 'varying_all/run4/output'
fsky = 0.75
lensed = False
# exclude = ['helium_fraction', 'scalar_nrun(1)', 'massless_neutrinos', 'omk', 'w','wa']  # None
# exclude = ['massless_neutrinos','w']
exclude = None
# =============================
# DERIVED
arcmin_from_fsky = fsky2arcmin(fsky)
sec_of_obs = years2sec(5)
Y = 0.25  # 25% yeld
# ===================
header = 'Joint fisher CMB T E + phi lensing used \n'
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
order = 5
# # step = np.array([-4,-3,-2,-1,1,2,3,4])
# step = np.array([-8,-6,-4,-2,2,4,6,8])

# for key in values.keys():
#   new_value[key] = par_gaps[key] * step + fid[key]
# for key in values.keys():
#   par_gaps[key]= np.abs(new_value[key][0]-new_value[key][1])
# new_value = collections.OrderedDict(sorted(new_value.items(), key=lambda t: t[0]))
# values = new_value
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

derivatives = np.ndarray((3, 3, np.size(dats[lmin_index:lmax_index, 0, 0]), n_values), dtype='float64')

for i in range(0, n_values):
    # computing derivatives.
    # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
                  # ---------------------------------
                              #   12h
    derivatives[:, :, :, i] = calc_deriv_vectorial(i, dats, par_gaps, values, order)[:, :, :]


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

            # assuming f Eq.4
            fisher[i, j] += (2. * ell + 1.) / 2. * fsky * trace
            fisher_save[i, j, iell] = (2. * ell + 1.) / 2. * fsky * trace
            # print np.sum(fisher_save[:,:,:], axis =2)[0,0],fisher[0,0]

    no_marginalized_ell[iell, :] = 1. / np.sqrt(np.diag(fisher))
    fisher_inv = np.linalg.inv(fisher)
    marginalized_ell[iell, :] = np.sqrt(np.diag(fisher_inv))


print 'ADDING PLANCK'
planck_fisher = np.loadtxt('/home/manzotti/n_eff-dependence-on-prior/n_priors_code/data/fisher_mat_joint_lmin=2_lmax=2500_ndet=Planck_fsky=0.2.txt')
# print planck_fisher
# print ''
# print fisher

# print ''
# print ''
# fisher +=planck_fisher

print 'lmax =', ell
# print fisher_inv

np.savetxt('data/{}/no_marginalized_ell_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder, lmin, lmax, N_det, fsky),
           np.column_stack((dats[lmin_index:lmax_index, 0, 0], no_marginalized_ell)), header=header)
np.savetxt('data/{}/marginalized_ell_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder, lmin, lmax, N_det, fsky),
           np.column_stack((dats[lmin_index:lmax_index, 0, 0], marginalized_ell)), header=header)
np.save('data/{}/full_fisher_mat_joint_lmin={}_lmax={}_ndet={}_fsky={}.npy'.format(output_folder, lmin, lmax, N_det, fsky),
        fisher_save)

np.savetxt('data/{}/ell_indeces_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder, lmin, lmax, N_det, fsky),
           dats[lmin_index:lmax_index, 0, 0], header=header)

# utils.study_prior_tau_on_N_eff(fid, fisher, 'data/' + output_folder, header)


print 'finally how much constraint on parameters without prior?'
print ''
fisher_single = fisher.copy()

fisher_inv = np.linalg.inv(fisher_single)

utils.save_cov_matrix(
    fisher_inv, 'data/{}/param_cov_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder, lmin, lmax, N_det, fsky))


np.savetxt('data/{}/invetered_sqrt_fisher_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder,
                                                                                            lmin, lmax, N_det, fsky), np.sqrt(fisher_inv), header=header)
np.savetxt('data/{}/fisher_mat_joint_lmin={}_lmax={}_ndet={}_fsky={}.txt'.format(output_folder,
                                                                                 lmin, lmax, N_det, fsky), fisher_single, header=header)

print 'fisher=', fisher

utils.print_resume_stat(fisher_single, fid)
