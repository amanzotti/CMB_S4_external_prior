'''

Simple fisher code to test prior effect on N effective neutrino estimates

TODO:
HIGH: lensing noise

low importance: ini file, automate the data creation process and derivative.

what are this data? documents better

create classes or at least external utilities. Maybe a better way to interface with cosmosis


The matrix is going to run an 4 parameters:

CONVENTIONS:

PARAMETER ORDER = h0, ns, As, Neff

'''

__author__ = "Y. Park"


import numpy as np
import math
import matplotlib.pyplot as plt
import util

# util.nl(noise_uK_arcmin, fwhm_arcmin, lmax)


def years2sec(years):
    return years * 365 * 24. * 60. * 60.

lmax = 3000

N_phi_l = np.loadtxt('multipole_noisebias.txt')


def C(iell, ell, parbin):
    '''
    Given CMB data dats it normalize them anf create a 3x3 matrix

    ell is the multiple
    iell is the index in the data ell corresponds to

    '''

    # noise definition from the number of observations and time
    # eq 1 of W.hu et al snowmass paper 10^6 detectors
    Y = 0.25  # 25%yeld
    N_det = 10 ** 6  # 1 milion of detectors
    s = 350 * math.sqrt(20626 * 60 * 60) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
    # print 'noise is ',s
    # s = 0.48 as in table from paper so it is ok.
    t = 2 / 60 / 180 * math.pi  # 2arcmin beam
    fac = ell * (ell + 1) / 2 / math.pi
    fac2 = ell ** 4 / ell / (ell + 1)

    # Final CMB noise definition
    N = s ** 2 * math.exp(ell * (ell + 1) * t ** 2 / 8 / math.log(2))
    N_phi = 0.* N_phi_l[iell, 1] * ell**2

    # Check again in particular cosmosis ouptut lensing
    # is it a 3x3 matrix? with    TT,TE,Tphi
    # TE,EE,Ephi
    # phiT,phiE,phiphi

    C = np.array([[dats[iell, 1, parbin] / fac + N, dats[iell, 2, parbin] / fac, 0],
                  [dats[iell, 2, parbin] / fac, dats[iell, 3, parbin] / fac + N * 2 ** 0.5, 0],
                  [0, 0, dats[iell, 4, parbin] + N_phi]])
    return C

# loading data. Each of this is a cmb Spectrum? probably cmb Tand E plus lensing
#  so the structure is data(:,:,i) is the i change in the parameters.
# rgw 2 :,: are the multiples and the column and the type of data respectively

dats = np.genfromtxt('data/dat0.txt')

for i in range(1, 9):
    newdat = np.genfromtxt('data/dat{}.txt'.format(i))
    dats = np.dstack((dats, newdat))

print dats[:, :, 0]

# creating the 4 by 4 matrix
fisher = np.zeros((4, 4))
# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = np.array([0.002, 0.02, 2e-11, 0.02])  # h0, ns, As, Neff

for iell, ell in enumerate(range(2, 4800)):
    #  filling it the matrix l goes from l_min =2 to l_max =5000

    print ell

    c0 = C(iell, ell, 0)  # 3x3 matrix in the fiducial cosmology
    # this is the covariance matrix of the data.

    cinv = np.linalg.inv(c0)

    for i in range(0, 4):
        for j in range(0, 4):
            # computing derivatives.
            ci = (C(iell, ell, i * 2 + 1) - C(iell, ell, i * 2 + 2)) / pargaps[i]
            cj = (C(iell, ell, j * 2 + 1) - C(iell, ell, j * 2 + 2)) / pargaps[j]

            # check and see here
            tot = np.dot(np.dot(np.dot(cinv, ci),  cinv), cj)

            # assuming f
            fisher[i, j] += (2 * ell + 1) / 2 * 0.5 * np.trace(tot)

d = []
d2 = []
d3 = []

for i in np.arange(-3, -1, 0.1):
    fisher1 = fisher.copy()
    fisher1[0, 0] += 1 / (10 ** i * 0.673) ** 2
    d.append(math.sqrt(np.linalg.inv(fisher1)[3, 3]))

    fisher2 = fisher.copy()
    fisher2[0, 0] += 1 / (10 ** i * 0.673) ** 2
    fisher2[1, 1] += 1 / (0.01 * 0.96) ** 2
    fisher2[2, 2] += 1 / (0.01 * 2.2e-9) ** 2

    d2.append(math.sqrt(np.linalg.inv(fisher2)[3, 3]))

    fisher3 = fisher.copy()[[0, 3], :][:, [0, 3]]
    fisher3[0, 0] += 1 / (10 ** i * 0.673) ** 2

    d3.append(math.sqrt(np.linalg.inv(fisher3)[1, 1]))

plt.plot(10 ** np.arange(-3, -1, 0.1), d, label='No Priors')
plt.plot(10 ** np.arange(-3, -1, 0.1), d2, label='1% Priors')
plt.plot(10 ** np.arange(-3, -1, 0.1), d3, label='Perfect Priors')
plt.xscale('log')
plt.xlabel(r'$\Delta H_0 / H_0$', fontsize=16)
plt.ylabel(r'$\sigma(N_\mathrm{eff})$', fontsize=16)
plt.legend(loc=0)

plt.savefig('h01.pdf')
plt.show()
