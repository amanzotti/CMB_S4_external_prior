'''

SECOND VERSION TO SUBSTITUTE THE first one A. Manzotti

Simple fisher code to test prior effect on N effective neutrino estimates

TODO:
HIGH: lensing noise (almost there: Komatsu and quicklens. MV not yet)

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

1) reproduce Wu paper

'''


import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
import sys
import utils

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
    Given CMB data dats it normalize them and create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    remember the order from CAMB
     l CTT CEE CBB CTE Cdd CdT CdE

    '''

    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    Y = 0.25  # 25%yeld
    N_det = 10 ** 4  # 1 milion of detectors
    s = 350. * math.sqrt(fsky2arcmin(0.75)) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # s = 0.48 as in table from paper so it is ok.
    t = 1. / 60. / 180. * math.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    N = (s * np.pi / 180. / 60.) ** 2 * math.exp(ell * (ell + 1.) * t ** 2 / 8. / math.log(2))
    # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    # C = np.array([data[iell, 5, parbin] + N_phi_l[iell, 1]])

    C = np.array([[data[iell, 1, parbin] / fac + N, data[iell, 4, parbin]],
                  [data[iell, 4, parbin], data[iell, 2, parbin] / fac + N * 2.],
                  ]
                 )

    return C

def Cov_matrix(iell, ell, data):
    '''
    Given CMB data dats it normalize them and create the covariance matrix.
    This was create to test the different formulas used in the fisher formalism
    ell is the multiple
    iell is the index in the data ell corresponds to
    remember the order from CAMB
    l CTT CEE CBB CTE Cdd CdT CdE

    '''

    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    Y = 0.25  # 25%yeld
    N_det = 10 ** 4  # 1 milion of detectors
    s = 350. * math.sqrt(fsky2arcmin(0.75)) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # s = 0.48 as in table from paper so it is ok.
    t = 1. / 60. / 180. * math.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    N = (s * np.pi / 180. / 60.) ** 2 * math.exp(ell * (ell + 1.) * t ** 2 / 8. / math.log(2))
    # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    # C = np.array([data[iell, 5, parbin] + N_phi_l[iell, 1]])
    CTTTT = (data[iell, 1, 0] / fac + N) ** 2
    CEEEE = (data[iell, 2, 0] / fac + 2. * N) ** 2
    CTETE = 0.5 * (data[iell, 4, 0] ** 2 + (data[iell, 1, 0] / fac + N) * (data[iell, 2, 0] / fac + 2. * N))

    CTTTE = data[iell, 4, 0] * (data[iell, 1, 0] / fac + N)


    CTTEE = (data[iell, 4, 0]) ** 2

    CEETE = data[iell, 4, 0] * (data[iell, 2, 0] / fac + 2. * N)


    Cov = np.array([[CTTTT,CTTEE,CTTTE],
                  [CTTEE,CEEEE,CEETE],[CTTTE,CEETE,CTETE]
                  ]
                 )

    return Cov

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================
l_t_max = 3000  # this is the multipole you want to cut the temperature Cl at, to simulate the effect of foregrounds
lmax = 2500
lmin = 30
N_phi_l = np.loadtxt('data/noise/wu_cdd_noise_4.txt')
run_idx = 3
fsky = 0.75
lensed =False
# =============================



# READ PARAMS
# load fiducial data
fid = pickle.load(open('data/run{}/fid_values.p'.format(run_idx), "rb"))
print "fid ", fid
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/run{}/grid_values.p'.format(run_idx), "rb"))
par_gaps = pickle.load(open('data/run{}/par_gaps.p'.format(run_idx), "rb"))
dats = utils.load_data(run_idx,lensed,values)

# cmb = np.genfromtxt('data/run{}/fiducial_lensedcls.dat'.format(run_idx),usecols=(0,1,2,3,4)) if lensed else np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx),usecols=(0,1,2,3,4))
# lensing = np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx),usecols=(5,6,7))
# dats = np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx))
# dats = np.concatenate((cmb,lensing[:cmb.shape[0],:]),axis=1)

# # Load data for all parameters variations
# for key, value in values.iteritems():
#     for i in np.arange(0, 4):
#         print key, values[key][i]
#         filename = 'data/run{}/'.format(run_idx)
#         # filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat'
#         filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat' if lensed else filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'

#         filename_lensing = filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
#         cmb = np.genfromtxt(filename_cmb.format(run_idx),usecols=(0,1,2,3,4))
#         lensing = np.genfromtxt(filename_lensing.format(run_idx),usecols=(5,6,7))
#         newdat = np.concatenate((cmb,lensing[:dats.shape[0],:]),axis=1)
#         # newdat = np.genfromtxt(filename)
#         dats = np.dstack((dats, newdat))

print fid.keys().index('omnuh2')

# cut Cl^T at ells bigger than l_t_max
dats[l_t_max:, 1, 1:] = 0.
dats[900:, 6, 0:] = 0.

n_values = np.size(values.keys())
print n_values
# creating the n_values by n_values matrix
fisher = np.zeros((n_values, n_values))
fisher_inv = np.zeros((n_values, n_values))

no_marginalized_ell = np.zeros((np.size(range(lmin, lmax)), n_values))

marginalized_ell = np.zeros((np.size(range(lmin, lmax)), n_values))


print 'fisher_size', fisher.shape
# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = par_gaps  # h0, ns, As, Neff,tau

for iell, ell in enumerate(range(lmin, lmax)):
    #  filling it the matrix l goes from l_min =2 to l_max =5000

    ell_index = np.where(dats[:, 0, 0] == ell)[0][0]

    c0 = np.zeros((3, 3))
    Cov = np.zeros((3, 3))

    c0 = C(ell_index, ell, 0, dats)  # 3x3 matrix in the fiducial cosmology
    Cov = Cov_matrix(ell_index, ell, dats)
    # this is the covariance matrix of the data. So in this case we have C^T C^E C^phi

    # sys.exit()

    cinv = np.linalg.inv(c0)
    Cov_inv = np.linalg.inv(Cov)

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
            # Compute Total value for this ell, the fsky 2l+1/2 correction is done later.
            tot = (ci[0,0] * cj[0,0] * Cov_inv[0,0] + ci[0,0] * cj[1,1] * Cov_inv[0,1] + ci[0,0] * cj[0,1] * Cov_inv[0,2]
                 +ci[1,1] * cj[0,0] * Cov_inv[0,1] + ci[1,1] * cj[1,1] * Cov_inv[1,1] + ci[1,1] * cj[0,1] * Cov_inv[1,2]
                 +ci[0,1] * cj[0,0] * Cov_inv[0,2] + ci[0,1] * cj[1,1] * Cov_inv[1,2] + ci[0,1] * cj[0,1] * Cov_inv[2,2])

            # assuming f Eq.4
            fisher[i, j] += (2. * ell + 1.) / 2. * fsky * tot
    no_marginalized_ell[iell, :] = 1. / np.sqrt(np.diag(fisher))
    fisher_inv = np.linalg.inv(fisher)
    marginalized_ell[iell, :] = np.sqrt(np.diag(fisher_inv))


np.savetxt('output_cmb/no_marginalized_ell.txt', np.column_stack((np.arange(lmin, lmax), no_marginalized_ell)))
np.savetxt('output_cmb/marginalized_ell.txt', np.column_stack((np.arange(lmin, lmax), marginalized_ell)))

d = []
d2 = []
d3 = []
for i in np.arange(-3, -1, 0.1):

    # '''alphabetical in CAMB description
    # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
    # PARAMETER ORDER = H0,Neff,tau,As,ns
    #                     0  1   2  3   4'''

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
    fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
        (0.01 * fid['re_optical_depth']) ** 2

    # Invert and get Neff error with these priors
    d2.append(
        math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher3 = fisher.copy()[[fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')], :][
        :, [fid.keys().index('hubble'), fid.keys().index('massless_neutrinos')]]

    fisher3[0, 0] += 1 / (10 ** i * fid['hubble']) ** 2

    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))
print ell

np.savetxt('output_cmb/sigma_H0_1percent.txt', d2)
np.savetxt('output_cmb/sigma_H0_noPrior.txt', d)
np.savetxt('output_cmb/sigma_H0_perfect_prior.txt', d3)

# plt.clf()
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d) * 100., label='No Priors')
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d2) * 100., label=r'1$\%$ Priors')
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d3) * 100., label='Perfect Priors')
# plt.xscale('log')
# plt.xlabel(r'$\Delta H_0 / H_0$', fontsize=16)
# plt.ylabel(r'$10^{2} ~ \sigma(N_\mathrm{eff}) $', fontsize=16)
# plt.legend(loc=0)

# plt.savefig('../images/h0_fisher.pdf')

# DO the same for tau

d = []
d2 = []
d3 = []

for i in np.arange(-3, -1, 0.1):

    # '''alphabetical in CAMB description
    # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
    # PARAMETER ORDER = H0,Neff,tau,As,ns
    #                     0  1   2  3   4'''
    fid_tau = fid['re_optical_depth']
    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
        (10 ** i * fid['re_optical_depth']) ** 2
    # Invert and get Neff error with these priors

    # print 'test = ', fisher1[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')], fisher[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] / (1 / (10 ** i * fid['re_optical_depth']) ** 2)

    d.append(
        math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
        (10 ** i * fid['re_optical_depth']) ** 2

    # add 1% prior on ns
    fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
        (0.01 * fid['scalar_spectral_index(1)']) ** 2
    # add 1% prior on As
    fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
        (0.01 * fid['scalar_amp(1)']) ** 2
    fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (0.01 * fid['hubble']) ** 2

    # Invert and get Neff error with these priors
    d2.append(
        math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher3 = fisher.copy()[[fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')], :][
        :, [fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')]]
    # Cicle on H0 priors
    # in the cut matrix tau is in the 0 place
    fisher3[0, 0] += 1 / (10 ** i * fid['re_optical_depth']) ** 2

    # Invert and get Neff error with these priors
    d3.append(
        math.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

np.savetxt('output_cmb/sigma_tau_1percent.txt', d2)
np.savetxt('output_cmb/sigma_tau_noPrior.txt', d)
np.savetxt('output_cmb/sigma_tau_perfect_prior.txt', d3)


# ===========================
# DO the same for n_s
# ===========================


d = []
d2 = []
d3 = []

for i in np.arange(-3, -1, 0.1):

    # '''alphabetical in CAMB description
    # hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
    # PARAMETER ORDER = H0,Neff,tau,As,ns
    #                     0  1   2  3   4'''
    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
        (10 ** i * fid['scalar_spectral_index(1)']) ** 2
    # Invert and get Neff error with these priors

    d.append(
        math.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 / \
        (10 ** i * fid['scalar_spectral_index(1)']) ** 2

    # add 1% prior on ns
    fisher2[fid.keys().index('re_optical_depth'), fid.keys().index('re_optical_depth')] += 1 / \
        (0.01 * fid['re_optical_depth']) ** 2
    # add 1% prior on As
    fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 / \
        (0.01 * fid['scalar_amp(1)']) ** 2
    fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (0.01 * fid['hubble']) ** 2

    # Invert and get Neff error with these priors
    d2.append(
        math.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

    fisher3 = fisher.copy()[[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('massless_neutrinos')], :][
        :, [fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('massless_neutrinos')]]
    # Cicle on H0 priors
    fisher3[0, 0] += 1 / (10 ** i * fid['scalar_spectral_index(1)']) ** 2

    # Invert and get Neff error with these priors
    d3.append(
        math.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

np.savetxt('output_cmb/sigma_ns_1percent.txt', d2)
np.savetxt('output_cmb/sigma_ns_noPrior.txt', d)
np.savetxt('output_cmb/sigma_ns_perfect_prior.txt', d3)

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
np.savetxt('output_cmb/param_cov.txt', param_cov)
np.savetxt('output_cmb/invetered_sqrt_fisher.txt', np.sqrt(fisher_inv))


for key, value in values.iteritems():

    print 'sigma(', key, ')', np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / fid[key], '%', "with no degeneracies", 1. / np.sqrt(fisher[fid.keys().index(key), fid.keys().index(key)])
