'''

SECOND VERSION TO SUBSTITUTE THE first one

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
hubble,massless_neutrinos,re_optical_depth,scalar_amp(1),scalar_spectral_index(1)
PARAMETER ORDER = H0,Neff,tau,As,ns
                    0  1   2  3   4



GOAL

1) reproduce Wu paper

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
    Y = 0.25  # 25%yeld
    N_det = 10 ** 6  # 1 milion of detectors
    s = 350. * math.sqrt(fsky2arcmin(0.5)) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # s = 0.48 as in table from paper so it is ok.
    t = 2. / 60. / 180. * math.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))
    # Final CMB noise definition
    N = (s * np.pi / 180. / 60.) ** 2 * math.exp(ell * (ell + 1) * t ** 2 / 8. / math.log(2))
    # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2
    # is it a 3x3 matrix? with
    # TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    C = np.array([[data[iell, 1, parbin] / fac + N, data[iell, 4, parbin], data[iell, 6, parbin]],
                  [data[iell, 4, parbin], data[iell, 2, parbin] / fac + N * 2.,                0],
                  [data[iell, 6, parbin],          0,         data[iell, 5, parbin] + 1e-8]]
                 )
    return C

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================
l_t_max = 3000  # this is the multipole you want to cut the temperature Cl at, to simulate the effect of foregrounds
lmax = 5000
N_phi_l = np.loadtxt('multipole_noisebias.txt')
run_idx = 3
fsky = 0.5

# =============================


# READ PARAMS
# load fiducial data
dats = np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx))
fid = np.genfromtxt('data/run{}/fiducial_pars.txt'.format(run_idx))
# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/run{}/grid_values.p'.format(run_idx), "rb"))
par_gaps = pickle.load(open('data/run{}/par_gaps.p'.format(run_idx), "rb"))
n_values = np.size(values.keys())
# Load data for all parameters variations
for key, value in values.iteritems():
    for i in np.arange(0, 4):
        print key, values[key][i]
        filename = 'data/run{}/'.format(run_idx)
        filename += key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
        newdat = np.genfromtxt(filename)
        dats = np.dstack((dats, newdat))

# cut Cl^T at ells bigger than l_t_max
dats[l_t_max:, 1, 1:] = 0.
# creating the n_values by n_values matrix
fisher = np.zeros((n_values, n_values))
# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = par_gaps  # h0, ns, As, Neff,tau

for iell, ell in enumerate(range(2, lmax)):
    #  filling it the matrix l goes from l_min =2 to l_max =5000

    print ell
    c0 = np.zeros((3, 3))
    c0 = C(iell, ell, 0, dats)  # 3x3 matrix in the fiducial cosmology
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

            ci = (-C(iell, ell, i * 4 + 4, dats) + 8. * C(iell, ell, i * 4 + 3, dats) - 8. *
                  C(iell, ell, i * 4 + 2, dats) + C(iell, ell, i * 4 + 1, dats)) / (12. * pargaps[values.keys()[i]])
            cj = (-C(iell, ell, j * 4 + 4, dats) + 8. * C(iell, ell, j * 4 + 3, dats) - 8. *
                  C(iell, ell, j * 4 + 2, dats) + C(iell, ell, j * 4 + 1, dats)) / (12. * pargaps[values.keys()[j]])
            # Eq 4.
            tot = np.dot(np.dot(np.dot(cinv, ci),  cinv), cj)
            # assuming f Eq.4
            fisher[i, j] += (2. * ell + 1.) / 2. * fsky * np.trace(tot)

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
    fisher1[0, 0] += 1 / (10 ** i * 67.04346) ** 2
    # Invert and get Neff error with these priors

    d.append(math.sqrt(np.linalg.inv(fisher1)[1, 1]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[0, 0] += 1 / (10 ** i * 67.04346) ** 2

    # add 1% prior on ns
    fisher2[4, 4] += 1 / (0.01 * 0.96) ** 2
    # add 1% prior on As
    fisher2[3, 3] += 1 / (0.01 * 2.2e-9) ** 2
    fisher2[2, 2] += 1 / (0.01 * 0.0924518) ** 2

    # Invert and get Neff error with these priors
    d2.append(math.sqrt(np.linalg.inv(fisher2)[1, 1]))

    fisher3 = fisher.copy()[[0, 1], :][:, [0, 1]]
    # Cicle on H0 priors
    fisher3[0, 0] += 1 / (10 ** i * 67.04346) ** 2

    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

np.savetxt('sigma_H0_1percent.txt',d)
np.savetxt('sigma_H0_noPrior.txt',d2)
np.savetxt('sigma_H0_perfect_prior.txt',d3)

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
    fid_tau=fid[2]
    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[2, 2] += 1 / (10 ** i * 0.0924518) ** 2
    # Invert and get Neff error with these priors

    d.append(math.sqrt(np.linalg.inv(fisher1)[1, 1]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[2, 2] += 1 / (10 ** i * 0.0924518) ** 2

    # add 1% prior on ns
    fisher2[4, 4] += 1 / (0.01 * fid[4]) ** 2
    # add 1% prior on As
    fisher2[3, 3] += 1 / (0.01 * fid[3]) ** 2
    fisher2[0, 0] += 1 / (0.01 * fid[0]) ** 2

    # Invert and get Neff error with these priors
    d2.append(math.sqrt(np.linalg.inv(fisher2)[1, 1]))



    fisher3 = fisher.copy()[[2, 1], :][:, [2, 1]]
    # Cicle on H0 priors
    fisher3[0, 0] += 1 / (10 ** i * 0.0924518) ** 2

    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

np.savetxt('sigma_tau_1percent.txt',d)
np.savetxt('sigma_tau_noPrior.txt',d2)
np.savetxt('sigma_tau_perfect_prior.txt',d3)


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
    fid_ns=fid[4]
    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[4, 4] += 1 / (10 ** i * fid_ns) ** 2
    # Invert and get Neff error with these priors

    d.append(math.sqrt(np.linalg.inv(fisher1)[1, 1]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[4, 4] += 1 / (10 ** i * fid_ns) ** 2

    # add 1% prior on ns
    fisher2[2, 2] += 1 / (0.01 * fid[2]) ** 2
    # add 1% prior on As
    fisher2[3, 3] += 1 / (0.01 * fid[3]) ** 2
    fisher2[0, 0] += 1 / (0.01 * fid[0]) ** 2

    # Invert and get Neff error with these priors
    d2.append(math.sqrt(np.linalg.inv(fisher2)[1, 1]))



    fisher3 = fisher.copy()[[4, 1], :][:, [4, 1]]
    # Cicle on H0 priors
    fisher3[0, 0] += 1 / (10 ** i * fid_ns) ** 2

    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

np.savetxt('sigma_ns_1percent.txt',d)
np.savetxt('sigma_ns_noPrior.txt',d2)
np.savetxt('sigma_ns_perfect_prior.txt',d3)

print 'finally how much constraint on parameters without prior?'
print ''
fisher_single = fisher.copy()
fisher_inv = np.sqrt(np.linalg.inv(fisher_single))
print 'sigma(H0)', fisher_inv[0,0],'=',100.*fisher_inv[0,0]/fid[0],'%'
print ''
print "sigma(Neff)", fisher_inv[1,1],'=',100.*fisher_inv[1,1]/fid[1],'%'
print ''
print "sigma(tau)", fisher_inv[2,2],'=',100.*fisher_inv[2,2]/fid[2],'%'
print ''
print "sigma(As)", fisher_inv[3,3],'=',100.*fisher_inv[3,3]/fid[3],'%'
print ''
print "sigma(ns)", fisher_inv[4,4],'=',100.*fisher_inv[4,4]/fid[4],'%'
print ''


# plt.clf()
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d) * 100., label='No Priors')
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d2) * 100., label=r'1$\%$ Priors')
# plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d3) * 100., label='Perfect Priors')
# plt.xscale('log')
# plt.xlabel(r'$\Delta \tau / \tau$', fontsize=16)
# plt.ylabel(r'$10^{2} ~ \sigma(N_\mathrm{eff}) $', fontsize=16)
# plt.legend(loc=0)
# plt.savefig('../images/tau_fisher.pdf')

