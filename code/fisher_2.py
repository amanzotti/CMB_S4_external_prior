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

PARAMETER ORDER = Neff,H0,ns,As
                    0  1   2  3



GOAL

1) reproduce Wu paper

'''

__author__ = "Y. Park"


import numpy as np
import math
import matplotlib.pyplot as plt
import util
import pickle
import sys

# util.nl(noise_uK_arcmin, fwhm_arcmin, lmax)


def years2sec(years):
    return years * 365 * 24. * 60. * 60.


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
    s = 350 * math.sqrt(20626 * 60 * 60) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2

    # s = 0.48 as in table from paper so it is ok.

    t = 2 / 60 / 180 * math.pi  # 2arcmin to rads beam
    fac = (ell * (ell + 1.) / 2. / math.pi) / (7.4311 * 10 ** 12)
    fac2 = (ell * (ell + 1.))

    # Final CMB noise definition
    N = s ** 2 * math.exp(ell * (ell + 1) * t ** 2 / 8 / math.log(2))  # this noise is in mu_K so check CMB accordingly
    # N_phi = 0. * N_phi_l[iell, 1] * ell ** 2

    # is it a 3x3 matrix? with    TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi
    C = np.array([[data[iell, 1, parbin] / fac + N, data[iell, 4, parbin], data[iell, 6, parbin]],
                  [data[iell, 4, parbin], data[iell, 2, parbin] / fac + N * 2.,                0],
                  [data[iell, 6, parbin],          0,         data[iell, 5, parbin]]]
                 )
    return C

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively


# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================
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

# Load data for all parameters variations
for key, value in values.iteritems():
    for i in np.arange(0, 4):
        # print key, values[key][i]
        filename = 'data/run{}/'.format(run_idx)
        filename += key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
        newdat = np.genfromtxt(filename)
        dats = np.dstack((dats, newdat))


# creating the 4 by 4 matrix
fisher = np.zeros((4, 4))
# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = par_gaps  # h0, ns, As, Neff

for iell, ell in enumerate(range(2, lmax)):
    #  filling it the matrix l goes from l_min =2 to l_max =5000

    print ell
    c0 = np.zeros((3, 3))
    c0 = C(iell, ell, 0, dats)  # 3x3 matrix in the fiducial cosmology
    # this is the covariance matrix of the data. So in this case we have C^T C^E C^phi

    # sys.exit()

    cinv = np.linalg.inv(c0)

    for i in range(0, 4):

        for j in range(0, 4):
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
            fisher[i, j] += (2. * ell + 1.) / 2. * 0.5 * np.trace(tot)

d = []
d2 = []
d3 = []

for i in np.arange(-3, -1, 0.1):

    fisher1 = fisher.copy()
    # Cicle on H0 priors
    fisher1[1, 1] += 1 / (10 ** i * 67.04346) ** 2
    # Invert and get Neff error with these priors

    d.append(math.sqrt(np.linalg.inv(fisher1)[0, 0]))

    fisher2 = fisher.copy()
    # Cicle on H0 priors

    fisher2[1, 1] += 1 / (10 ** i * 67.04346) ** 2

    # add 1% prior on ns
    fisher2[2, 2] += 1 / (0.01 * 0.96) ** 2
    # add 1% prior on As
    fisher2[3, 3] += 1 / (0.01 * 2.2e-9) ** 2
    # Invert and get Neff error with these priors
    d2.append(math.sqrt(np.linalg.inv(fisher2)[0, 0]))

    fisher3 = fisher.copy()[[0, 1], :][:, [0, 1]]
    # Cicle on H0 priors

    fisher3[1, 1] += 1 / (10 ** i * 67.04346) ** 2
    # Invert and get Neff error with these priors
    d3.append(math.sqrt(np.linalg.inv(fisher3)[0, 0]))

plt.clf()
plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d) * 100., label='No Priors')
plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d2) * 100., label=r'1$\%$ Priors')
plt.plot(10 ** np.arange(-3, -1, 0.1), np.array(d3) * 100., label='Perfect Priors')
plt.xscale('log')
plt.xlabel(r'$\Delta H_0 / H_0$', fontsize=16)
plt.ylabel(r'$10^{2} ~ \sigma(N_\mathrm{eff}) $', fontsize=16)
plt.legend(loc=0)

plt.savefig('h0_fisher2.pdf')
plt.show()
