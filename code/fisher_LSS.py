'''

Fisher matrix code. This is intended for LSS so it will load and take derivatives of the power spectrum.

Simple fisher code to test prior effect on N effective neutrino estimates

TODO:
HIGH:

low importance: ini file. Think about PCA to understand what is more important for N_eff


create classes or at least external utilities. Maybe a better way to interface with cosmosis

# DATA FROM CAMB
Look at eq (24) of http://arxiv.org/pdf/0805.4238.pdf
or somewhere else fot the formula

Fij = \int _kmin^kmax (4pik^2dk)/(2pi)^3 partial P(k)/parameter i partial P(k)/parameter j w(k)

w(k) = 1/2(ngPg(1+ngPg))Vsurvey = 1/2 Veff

BE VERY CAREFULL AS USUALLY AT the values you get from CAMB there are h factor everywhere

The matrix is going to run an 4 parameters:

CONVENTIONS:
alphabetical in CAMB description
hubble,massless_neutrinos,omnuh2,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)




GOAL

1) reproduce Wu paper and beyond

'''


import numpy as np
import math
import matplotlib.pyplot as plt
import util
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


# def C(iell, ell, parbin, data):
#     '''
#     Given CMB data dats it normalize them anf create a 3x3 matrix

#     ell is the multiple
#     iell is the index in the data ell corresponds to

#     remember the order from CAMB
#      l CTT CEE CBB CTE Cdd CdT CdE

#     '''

# noise definition from the number of observations and time
# eq 1 of W.hu et al snowmass paper 10^6 detectors
# Y = 0.25  # 25% yeld
# N_det = 10 ** 6  # 1 milion of detectors
# s = 350. * math.sqrt(fsky2arcmin(0.75)) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
# s = 0.48 as in table from paper so it is ok.
# t = 1. / 60. / 180. * math.pi  # 2arcmin to rads beam
#     fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
#     fac2 = (ell * (ell + 1.))
# Final CMB noise definition
#     N = (s * np.pi / 180. / 60.) ** 2 * math.exp(ell * (ell + 1.) * t ** 2 / 8. / math.log(2))
# this noise is in mu_K so check CMB accordingly
# N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
# is it a 3x3 matrix? with
# TT,TE,Tphi
# TE,EE,Ephi
# phiT,phiE,phiphi
#     C = np.array([[data[iell, 1, parbin] / fac + N, data[iell, 4, parbin], data[iell, 6, parbin]],
#                   [data[iell, 4, parbin], data[iell, 2, parbin] / fac + N * 2.,               0.],
#                   [data[iell, 6, parbin],          0.,         data[iell, 5, parbin] + N_phi_l[iell, 1]]]
#                  )

#     return C


# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
run_idx = 3

# READ PARAMS
# load fiducial data
dats = np.genfromtxt('data/run{}/fiducial_matterpower.dat'.format(run_idx))
lenght = dats.shape[0]
print dats.shape[0]
# load fiducial parameters used
fid = pickle.load(open('data/run{}/fid_values.p'.format(run_idx), "rb"))
print "fid ", fid
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/run{}/grid_values.p'.format(run_idx), "rb"))
par_gaps = pickle.load(open('data/run{}/par_gaps.p'.format(run_idx), "rb"))

# exclude = ['re_optical_depth', 'ombh2', 'w', 'omnuh2']
exclude = ['re_optical_depth']

for e in exclude:
    par_gaps.pop(e)
    values.pop(e)
    fid.pop(e)

# Load data for all parameters variations
for key, value in values.iteritems():
    for i in np.arange(0, 4):
        print dats.shape
        print key, values[key][i]
        filename = 'data/run{}/'.format(run_idx)
        filename += key + '_{:.13f}'.format(values[key][i]) + '_matterpower.dat'
        newdat = np.genfromtxt(filename)
        print newdat.shape[0]
        lenght_temp = newdat.shape[0]
        if lenght_temp < lenght:
            lenght = lenght_temp
            # print dats[:lenght_temp,:,:].shape[0],newdta.shape[0]
            try:
                dats = np.dstack((dats[:lenght_temp, :], newdat))
            except:
                dats = np.dstack((dats[:lenght_temp, :, :], newdat))

        elif lenght_temp > lenght:
            try:
                dats = np.dstack((dats, newdat[:lenght, :, :]))
            except:
                dats = np.dstack((dats, newdat[:lenght, :]))

        elif lenght_temp == lenght:
            dats = np.dstack((dats, newdat))
# =============================
# =============================
# =============================
# k_ov_h_max = dats[-1,0,0]
# k_ov_h_min = dats[0,0,0]
V_eff = 1e9  # in h^-3 Mpc^3 unit so this is 1 Gpc^h^-3
# =============================
n_values = np.size(values.keys())


print 'data loaded'

k_ov_h_bins =dats[10:-1, 0, 0]

# creating the n_values by n_values matrix
fisher = np.zeros((n_values, n_values))
fisher_inv = np.zeros((n_values, n_values))

no_marginalized_k = np.zeros((np.size(k_ov_h_bins), n_values))

marginalized_k = np.zeros((np.size(k_ov_h_bins), n_values))


print 'fisher_size', fisher.shape
pargaps = par_gaps
# this is k_ov_h
# sum over k you are mimicking an integral here

print 'start fisher with fid', fid
k=0
k_old = 0
for ik, k in enumerate(k_ov_h_bins):
    #  filling it the matrix k goes from k_min k_max

    delta_k = k-k_old
    k_old=k
    print ik, k, 'delta_k' ,delta_k

    print ik, "out of", np.size(dats[3:-100, 0, 0])
    covariance = (2. * dats[ik, 1, 0] ** 2) / V_eff  # cov = 2*P(k)^2/N_k with N_k number of mode

    cinv = 1 / covariance

    for i in range(0, n_values):
        # cicling on cosmological parameters

        for j in range(0, n_values):
            # computing derivatives.
            # ci = (C(iell, ell, i * 4 + 3,dats) - C(iell, ell, i * 4 + 2,dats)) /2./ pargaps[values.keys()[i]]
            # cj = (C(iell, ell, j * 4 + 3,dats) - C(iell, ell, j * 4 + 2,dats)) /2./ pargaps[values.keys()[j]]

            # print j, cj[1,1]

            # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
                  # ---------------------------------
                              #   12h

            ci = (-dats[ik, 1, i * 4 + 4] + 8. * dats[ik, 1, i * 4 + 3] - 8. *
                  dats[ik, 1, i * 4 + 2] + dats[ik, 1, i * 4 + 1]) / (12. * pargaps[values.keys()[i]])

            cj = (-dats[ik, 1, j * 4 + 4] + 8. * dats[ik, 1, j * 4 + 3] - 8. *
                  dats[ik, 1, j * 4 + 2] + dats[ik, 1, j * 4 + 1]) / (12. * pargaps[values.keys()[j]])
            # Eq 4.
            # print ci, cj
            tot = cinv * ci * cj
            # assuming f Eq.4
            fisher[i, j] += tot * delta_k * k ** 2 / (2. * np.pi ** 2)  # ! h^3 Mpc^-3
            # print fisher

    no_marginalized_k[ik, :] = 1. / np.sqrt(np.diag(fisher))
    fisher_inv = np.linalg.inv(fisher)
    marginalized_k[ik, :] = np.sqrt(np.diag(fisher_inv))

print k

np.savetxt('output/no_marginalized_k.txt', np.column_stack((k_ov_h_bins, no_marginalized_k)))
np.savetxt('output/marginalized_k.txt', np.column_stack((k_ov_h_bins, marginalized_k)))

# sys.exit()

# =======================================================
#  N_eff with H0 priors
# =======================================================
d = []
d2 = []
d3 = []
for i in np.arange(-3, -1, 0.1):

    # '''alphabetical in CAMB description
    # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)

    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (10 ** i * fid['hubble']) ** 2
    # Invert and get Neff error with these priors

    d.append(
        math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (10 ** i * fid['hubble']) ** 2

    # add 1% prior on ns
    fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
        (0.01 * fid['scalar_spectral_index(1)']) ** 2
    # add 1% prior on As
    fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
        (0.01 * fid['scalar_amp(1)']) ** 2
    # fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
    #     (0.01 * fid['re_optical_depth']) ** 2

    # Invert and get Neff error with these priors
    d2.append(
        math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher3 = fisher.copy()[[fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')], :][
        :, [fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')]]

    fisher3[0, 0] += 1 / (10 ** i * fid['hubble']) ** 2

    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

np.savetxt('output/sigma_H0_1percent_LSS.txt', d2)
np.savetxt('output/sigma_H0_noPrior_LSS.txt', d)
np.savetxt('output/sigma_H0_perfect_prior_LSS.txt', d3)


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

# np.savetxt('output/sigma_ns_1percent_LSS.txt', d2)
# np.savetxt('output/sigma_ns_noPrior_LSS.txt', d)
# np.savetxt('output/sigma_ns_perfect_prior_LSS.txt', d3)

print 'finally how much constraint on parameters without prior?'
print ''
fisher_single = fisher.copy()

fisher_inv = np.linalg.inv(fisher_single)

param_cov = np.zeros((n_values, n_values))
for i in range(n_values):
    for j in range(n_values):
        if i != j:
            param_cov[i, j] = fisher_inv[i, j] / np.sqrt(fisher_inv[i, i] * fisher_inv[j, j])
# print param_cov
np.savetxt('output/param_cov_LSS.txt', param_cov)
np.savetxt('output/invetered_sqrt_fisher_LSS.txt', np.sqrt(fisher_inv))


# print fisher_inv
print 'sigma(H0)', np.sqrt(fisher_inv[fid.keys().index('hubble'), fid.keys().index('hubble')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('hubble'), fid.keys().index('hubble')]) / fid['hubble'], '%'

print ''
print "sigma(Neff)", np.sqrt(fisher_inv[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]) / fid['massless_neutrinos'], '%'


print ''
# print "sigma(tau)", np.sqrt(fisher_inv[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')]) / fid['re_optical_depth'], '%'
# print ''
# print ''
# print "sigma(omnuh2)", np.sqrt(fisher_inv[fid.keys().index('omnuh2'), fid.keys().index('omnuh2')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('omnuh2'), fid.keys().index('omnuh2')]) / fid['omnuh2'], '%'

print ''
print "sigma(As)", np.sqrt(fisher_inv[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')]) / fid['scalar_amp(1)'], '%'

print ''
print "sigma(ns)", np.sqrt(fisher_inv[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')]) / fid['scalar_spectral_index(1)'], '%'
print ''
