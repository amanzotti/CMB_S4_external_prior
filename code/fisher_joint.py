'''

SECOND VERSION TO SUBSTITUTE THE first one (AM)

Simple fisher code to test prior effect on N effective neutrino estimates

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
alphabetical in CAMB description
hubble,massless_neutrinos,omnuh2,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
PARAMETER ORDER = H0,Neff,Omega_nuh^2,tau,As,ns
                    0  1      2        3   4  5



GOAL

1) reproduce Wu paper and beyond

'''


import numpy as np
import math
import matplotlib.pyplot as plt
import utils
import pickle
import sys

# util.nl(noise_uK_arcmin, fwhm_arcmin, lmax)


def years2sec(years):
    ''' years to sec '''
    return years * 365 * 24. * 60. * 60.


def fsky2arcmin(fsky):
    '''convert fsky in fraction of unity to arcmin^2'''
    # 41253 square degrees in all sky
    return 41253. * fsky * 60. * 60.


def C(iell, ell, parbin, data):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE

    '''

    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    Y = 0.25  # 25% yeld
    N_det = 10 ** 6  # 1 milion of detectors
    s = 350. * np.sqrt(fsky2arcmin(0.75)) / np.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # s = 0.48 as in table from paper so it is ok.
    t = 1. / 60. / 180. * np.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / np.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    # overwrite s to check
    s=1.5 #1.5 muK
    N = (s * np.pi / 180. / 60.) ** 2 * np.exp(ell * (ell + 1.) * t ** 2 / 8. / np.log(2))
    # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    C = np.array([[data[iell, 1, parbin] +  fac*N, data[iell, 4, parbin], data[iell, 6, parbin]],
                  [data[iell, 4, parbin], data[iell, 2, parbin] + fac * N * 2.,               0.],
                  [data[iell, 6, parbin],          0.,         data[iell, 5, parbin] + N_phi_l[iell, 1]]]
                 )

    return C

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================
l_t_max = 3000  # this is the multipole you want to cut the temperature Cl at, to simulate the effect of foregrounds
lmax = 3000
lmin = 2
N_phi_l = np.loadtxt('data/noise/wu_cdd_noise_6.txt')
run_idx = 'test_fisher'
data_folder = 'test_fisher/run1'
fsky = 0.5
lensed= False
exclude = None
# =============================


# READ PARAMS
# load fiducial data
# load fiducial parameters used
fid = pickle.load(open('data/{}/fid_values.p'.format(data_folder), "rb"))
print "fid ", fid
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/{}/grid_values.p'.format(data_folder), "rb"))
par_gaps = pickle.load(open('data/{}/par_gaps.p'.format(data_folder), "rb"))

# exclude = ['omnuh2','w','ombh2','omch2']
# exclude = ['w']
# par_gaps,values,fid = utils.exclude_parameters(exclude,par_gaps,values,fid)

print 'loading files'
dats = utils.load_data(data_folder,values, lensed)
n_values = np.size(values.keys())

# # Load data for all parameters variations
# for key, value in values.iteritems():
#     for i in np.arange(0, 4):
#         print key, values[key][i]
#         filename = 'data/run{}/'.format(run_idx)
#         filename += key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
#         newdat = np.genfromtxt(filename)
#         dats = np.dstack((dats, newdat))




# cut Cl^T at ells bigger than l_t_max
dats[l_t_max:, 1, 1:] = 0.
# phi_T has oscillations in it.
dats[900:, 6, 0:] = 0.

# creating the n_values by n_values matrix
fisher = np.zeros((n_values, n_values))
fisher_inv = np.zeros((n_values, n_values))

no_marginalized_ell = np.zeros((np.size(range(lmin, lmax)), n_values))
marginalized_ell = np.zeros((np.size(range(lmin, lmax)), n_values))


print 'fisher_size', fisher.shape
pargaps = par_gaps

for iell, ell in enumerate(range(lmin, lmax)):
    #  filling it the matrix l goes from l_min =2 to l_max =5000

    ell_index = np.where(dats[:,0,0]==ell)[0][0]

    c0 = np.zeros((3, 3))
    c0 = C(ell_index, ell, 0, dats)  # 3x3 matrix in the fiducial cosmology
    # this is the covariance matrix of the data. So in this case we have C^T C^E C^phi

    # sys.exit()

    cinv = np.linalg.inv(c0)

    for i in range(0, n_values):

        for j in range(0, n_values):
            # computing derivatives.
            # ci = (C(iell, ell, i * 4 + 3,dats) - C(iell, ell, i * 4 + 2,dats)) /2./ pargaps[values.keys()[i]]
            # cj = (C(iell, ell, j * 4 + 3,dats) - C(iell, ell, j * 4 + 2,dats)) /2./ pargaps[values.keys()[j]]

            # print j, cj[1,1]

            # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
                  # ---------------------------------
                              #   12h

            ci = (-C(ell_index, ell, i * 4 + 4, dats) + 8. * C(ell_index, ell, i * 4 + 3, dats) - 8. *
                  C(ell_index, ell, i * 4 + 2, dats) + C(ell_index, ell, i * 4 + 1, dats)) / (12. * pargaps[values.keys()[i]])
            cj = (-C(ell_index, ell, j * 4 + 4, dats) + 8. * C(ell_index, ell, j * 4 + 3, dats) - 8. *
                  C(ell_index, ell, j * 4 + 2, dats) + C(ell_index, ell, j * 4 + 1, dats)) / (12. * pargaps[values.keys()[j]])
            # Eq 4.
            tot = np.dot(np.dot(np.dot(cinv, ci),  cinv), cj)
            # assuming f Eq.4
            fisher[i, j] += (2. * ell + 1.) / 2. * fsky * np.trace(tot)

    no_marginalized_ell[iell, :] = 1. / np.sqrt(np.diag(fisher))
    fisher_inv = np.linalg.inv(fisher)
    marginalized_ell[iell, :] = np.sqrt(np.diag(fisher_inv))

print ell

np.savetxt('output/no_marginalized_ell.txt', np.column_stack((np.arange(lmin, lmax), no_marginalized_ell)))
np.savetxt('output/marginalized_ell.txt', np.column_stack((np.arange(lmin, lmax), marginalized_ell)))


# # =======================================================
# #  N_eff with H0 priors
# # =======================================================
# d = []
# d2 = []
# d3 = []
# for i in np.arange(-3, -1, 0.1):

#     # '''alphabetical in CAMB description
#     # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)


#     fisher1 = fisher.copy()
#     # Cicle on H0 priors
#     fisher1[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (10 ** i * fid['hubble']) ** 2
#     # Invert and get Neff error with these priors

#     d.append(
#         math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher2 = fisher.copy()
#     # Cicle on H0 priors

#     fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (10 ** i * fid['hubble']) ** 2

#     # add 1% prior on ns
#     fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
#         (0.01 * fid['scalar_spectral_index(1)']) ** 2
#     # add 1% prior on As
#     fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
#         (0.01 * fid['scalar_amp(1)']) ** 2
#     fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
#         (0.01 * fid['re_optical_depth']) ** 2

#     # Invert and get Neff error with these priors
#     d2.append(
#         math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher3 = fisher.copy()[[fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')], :][
#         :, [fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')]]

#     fisher3[0, 0] += 1 / (10 ** i * fid['hubble']) ** 2

#     # Invert and get Neff error with these priors
#     d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

# np.savetxt('output/sigma_H0_1percent.txt', d2)
# np.savetxt('output/sigma_H0_noPrior.txt', d)
# np.savetxt('output/sigma_H0_perfect_prior.txt', d3)


# # =======================================================
# #  N_eff with tau priors
# # =======================================================

# d = []
# d2 = []
# d3 = []

# for i in np.arange(-3, -1, 0.1):

#     # '''alphabetical in CAMB description
#     # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
#     # PARAMETER ORDER = H0,Neff,tau,As,ns
#     #                     0  1   2  3   4'''
#     fisher1 = fisher.copy()
#     # Cicle on H0 priors
#     fisher1[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
#         (10 ** i * fid['re_optical_depth']) ** 2
#     # Invert and get Neff error with these priors

#     # print 'test = ', fisher1[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')], fisher[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] / (1 /
#     #                                                                                                                                                                                             (10 ** i * fid['re_optical_depth']) ** 2)

#     d.append(
#         math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher2 = fisher.copy()
#     # Cicle on H0 priors

#     fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
#         (10 ** i * fid['re_optical_depth']) ** 2

#     # add 1% prior on ns
#     fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
#         (0.01 * fid['scalar_spectral_index(1)']) ** 2
#     # add 1% prior on As
#     fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
#         (0.01 * fid['scalar_amp(1)']) ** 2
#     fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (0.01 * fid['hubble']) ** 2

#     # Invert and get Neff error with these priors
#     d2.append(
#         math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher3 = fisher.copy()[[fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')], :][
#         :, [fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')]]
#     # Cicle on H0 priors
#     # in the cut matrix tau is in the 0 place
#     fisher3[0, 0] += 1 / (100 ** i * fid['re_optical_depth']) ** 2

#     # Invert and get Neff error with these priors
#     d3.append(
#         math.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

# np.savetxt('output/sigma_tau_1percent.txt', d2)
# np.savetxt('output/sigma_tau_noPrior.txt', d)
# np.savetxt('output/sigma_tau_perfect_prior.txt', d3)


# # =======================================================
# #  N_eff with ns priors
# # =======================================================


# d = []
# d2 = []
# d3 = []

# for i in np.arange(-3, -1, 0.1):

#     # '''alphabetical in CAMB description
#     # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
#     # PARAMETER ORDER = H0,Neff,tau,As,ns
#     #                     0  1   2  3   4'''
#     fisher1 = fisher.copy()
#     # Cicle on H0 priors
#     fisher1[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
#         (10 ** i * fid['scalar_spectral_index(1)']) ** 2
#     # Invert and get Neff error with these priors

#     d.append(
#         math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher2 = fisher.copy()
#     # Cicle on H0 priors

#     fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
#         (10 ** i * fid['scalar_spectral_index(1)']) ** 2

#     # add 1% prior on ns
#     fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
#         (0.01 * fid['re_optical_depth']) ** 2
#     # add 1% prior on As
#     fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
#         (0.01 * fid['scalar_amp(1)']) ** 2
#     fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (0.01 * fid['hubble']) ** 2

#     # Invert and get Neff error with these priors
#     d2.append(
#         math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

#     fisher3 = fisher.copy()[[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('massless_neutrinos')], :][
#         :, [fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('massless_neutrinos')]]
#     # Cicle on H0 priors
#     fisher3[0, 0] += 1 / (10 ** i * fid['scalar_spectral_index(1)']) ** 2

#     # Invert and get Neff error with these priors
#     d3.append(
#         math.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

# np.savetxt('output/sigma_ns_1percent.txt', d2)
# np.savetxt('output/sigma_ns_noPrior.txt', d)
# np.savetxt('output/sigma_ns_perfect_prior.txt', d3)

print 'finally how much constraint on parameters without prior?'
print ''
fisher_single = fisher.copy()

fisher_inv = np.linalg.inv(fisher_single)

utils.save_cov_matrix(fisher_inv)


# param_cov = np.zeros((6, 6))
# for i in range(6):
#     for j in range(6):
#         if i != j:
#             param_cov[i, j] = fisher_inv[i, j] / np.sqrt(fisher_inv[i, i] * fisher_inv[j, j])
# # print param_cov
# np.savetxt('output/param_cov.txt', param_cov)

np.savetxt('output/invetered_sqrt_fisher.txt', np.sqrt(fisher_inv))
np.savetxt('output/fisher_mat.txt', fisher_single)

print fisher

utils.print_resume_stat(fisher,fid)

# for key, value in values.iteritems():

#     print 'sigma(',key,')', np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / fid[key], '%' ,"with no degeneracies", 1./np.sqrt(fisher[fid.keys().index(key), fid.keys().index(key)])
