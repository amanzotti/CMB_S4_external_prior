#/usr/bin/env python
# --
# quicklens/examples/plot_lens_reconstruction_noise_levels.py
# --
# calculates  reconstruction noise levels
# for an idealized experiment with symmetric beam and white
# pixel noise calculation is done in both the
# full-sky and flat-sky limits for comparison.

# TODO load details from ini

import numpy as np
import pylab as pl
import configparser
import quicklens as ql
import sys
import math


def years2sec(years):
    return years * 365 * 24. * 60. * 60.


def fsky2arcmin(fsky):
    '''convert fsky in fraction of unity to arcmin^2'''
    # 41253 square degrees in all sky
    return 41253. * fsky * 60. * 60.

parser = configparser.ConfigParser()
configfile = './lensing_noise.ini'
parser.read(configfile)
N_det = parser.getfloat('lensing_noise', 'N_det')
beam_fwhm = parser.getfloat('lensing_noise', 'beam_fwhm')
lmin = parser.getint('lensing_noise', 'lmin')
lmax = parser.getint('lensing_noise', 'lmax')
fsky = 0.5
# calculation parameters.
# lmax = 3000  # maximum multipole for T, E, B and \phi.
nx = 512  # number of pixels for flat-sky calc.
dx = 1. / 60. / 180. * np.pi  # pixel width in radians.
Y = 0.25
print N_det, beam_fwhm, lmin, lmax, Y

# nlev_t = 5.   # temperature noise level, in uK.arcmin.
# nlev_p = 5.   # polarization noise level, in uK.arcmin.
# beam_fwhm = 1.   # Gaussian beam full-width-at-half-maximum.

# unlensed theory spectra.
cl_unl = ql.spec.get_camb_scalcl(
    fname='/home/manzotti/n_eff-dependence-on-prior/code/data/run3/fiducial_scalcls.dat', lmax=lmax)
# lensed theory spectra.
cl_len = ql.spec.get_camb_lensedcl(
    fname='/home/manzotti/n_eff-dependence-on-prior/code/data/run3/fiducial_lensedcls.dat', lmax=lmax)
fac = (7.4311 * 10 ** 12)
print cl_unl.cltt
bl = ql.spec.bl(beam_fwhm, lmax)  # transfer function.
pix = ql.maps.pix(nx, dx)


# noise definition from the number of observations and time
# eq 1 of W.hu et al snowmass paper 10^6 detectors
# Y = 0.25  # 25%yeld
# N_det = 10 ** 6  # 1 milion of detectors
# math.sqrt(20626. * 60. * 60.)  # half sky in arcmin^2

nlev_t = 350. * math.sqrt(fsky2arcmin(0.5)) / math.sqrt(N_det * Y * years2sec(5))  # half sky in arcmin^2
nlev_p = nlev_t * 2.
# print 'error' , nlev_t,nlev_p

# noise spectra
nltt = (np.pi / 180. / 60. * nlev_t) ** 2 / bl ** 2
nlee = nlbb = (np.pi / 180. / 60. * nlev_p) ** 2 / bl ** 2


# signal spectra
sltt = cl_len.cltt * (7.4311 * 10 ** 12)
slte = cl_len.clte * (7.4311 * 10 ** 12)
slee = cl_len.clee * (7.4311 * 10 ** 12)
slbb = cl_len.clbb * (7.4311 * 10 ** 12)
zero = np.zeros(lmax + 1)

# signal+noise spectra
cltt = sltt + nltt
clee = slee + nlee
clbb = slbb + nlbb

# filter functions
flt = np.zeros(lmax + 1)
flt[2:] = 1. / cltt[2:]
fle = np.zeros(lmax + 1)
fle[2:] = 1. / clee[2:]
flb = np.zeros(lmax + 1)
flb[2:] = 1. / clbb[2:]

# intialize quadratic estimators
qest_TT = ql.qest.lens.phi_TT(sltt)
qest_EE = ql.qest.lens.phi_EE(slee)
qest_TE = ql.qest.lens.phi_TE(slte)
qest_TB = ql.qest.lens.phi_TB(slte)
qest_EB = ql.qest.lens.phi_EB(slee)

# calculate the noise spectra.
watch = ql.util.stopwatch()


def calc_nlqq(qest, clXX, clXY, clYY, flX, flY):
    errs = np.geterr()
    np.seterr(divide='ignore', invalid='ignore')

# Commenting flat sky, uncomment if you need it
    # print "[%s]" % watch.elapsed(), "calculating flat-sky noise level for estimator of type", type(qest)
    # clqq_flatsky = qest.fill_clqq(ql.maps.cfft(nx, dx), clXX * flX * flX, clXY * flX * flY, clYY * flY * flY)
    # resp_flatsky = qest.fill_resp(qest, ql.maps.cfft(nx, dx), flX, flY)
    # nlqq_flatsky = clqq_flatsky / resp_flatsky ** 2

    print "[%s]" % watch.elapsed(), "calculating full-sky noise level for estimator of type", type(qest)
    clqq_fullsky = qest.fill_clqq(
        np.zeros(lmax + 1, dtype=np.complex), clXX * flX * flX, clXY * flX * flY, clYY * flY * flY)
    resp_fullsky = qest.fill_resp(qest, np.zeros(lmax + 1, dtype=np.complex), flX, flY)
    nlqq_fullsky = clqq_fullsky / resp_fullsky ** 2

    np.seterr(**errs)
    # return nlqq_flatsky, nlqq_fullsky
    return nlqq_fullsky


def compute_mv(nlpp_TT, nlpp_EE, nlpp_TE, nlpp_TB, nlpp_EB):
    '''
you need a 5x5 wit all the estimator and cross correlation between them

    '''

    nlpp_TT = np.nan_to_num(nlpp_TT.astype(float))
    nlpp_EE = np.nan_to_num(nlpp_EE.astype(float))
    nlpp_TE = np.nan_to_num(nlpp_TE.astype(float))
    nlpp_TB = np.nan_to_num(nlpp_TB.astype(float))
    nlpp_EB = np.nan_to_num(nlpp_EB.astype(float))

    ells = np.shape(nlpp_TT)[0]
    temp = np.zeros((3, 3, ells))
    nlpp_mv = np.zeros_like(nlpp_TT)
    temp[0, 0, :] = nlpp_TT[:]
    temp[0, 1, :] = nlpp_TE[:]
    temp[1, 0, :] = temp[0, 1, :]
    temp[0, 2, :] = nlpp_TB[:]
    temp[2, 0, :] = temp[0, 2, :]
    temp[1, 1, :] = nlpp_EE[:]
    temp[1, 2, :] = nlpp_EB[:]
    temp[2, 1, :] = temp[1, 2, :]

    # invert at a given l
    for ell in np.arange(1, ells):
        # invert at a given l
        nlpp_mv[ell] = 1. / np.sum(np.linalg.inv(temp[:, :, ell]).flatten())

    return nlpp_mv


nlpp_TT_fullsky = calc_nlqq(qest_TT, cltt, cltt, cltt, flt, flt)
# nlpp_EE_flatsky, nlpp_EE_fullsky = calc_nlqq(qest_EE, clee, clee, clee, fle, fle)
nlpp_EE_fullsky = calc_nlqq(qest_EE, clee, clee, clee, fle, fle)
nlpp_TE_fullsky = calc_nlqq(qest_TE, cltt, slte, clee, flt, fle)
nlpp_TB_fullsky = calc_nlqq(qest_TB, cltt, zero, clbb, flt, flb)
nlpp_EB_fullsky = calc_nlqq(qest_EB, clee, zero, clbb, fle, flb)


# convert nan

nlpp_TT_fullsky = np.nan_to_num(nlpp_TT_fullsky.astype(float))
nlpp_EE_fullsky = np.nan_to_num(nlpp_EE_fullsky.astype(float))
nlpp_TE_fullsky = np.nan_to_num(nlpp_TE_fullsky.astype(float))
nlpp_TB_fullsky = np.nan_to_num(nlpp_TB_fullsky.astype(float))
nlpp_EB_fullsky = np.nan_to_num(nlpp_EB_fullsky.astype(float))

# compute_mv()

# make plot
ls = np.arange(0, lmax + 1)
t = lambda l: (l * (l + 1.)) ** 2 / (2. * np.pi)  # scaling to apply to cl_phiphi when plotting.
lbins = np.linspace(2, lmax, 1)       # multipole bins.

np.savetxt('TT_lensing_noise.txt', np.vstack((ls, nlpp_TT_fullsky)).T)
np.savetxt('EE_lensing_noise.txt', np.vstack((ls, nlpp_EE_fullsky)).T)
np.savetxt('TE_lensing_noise.txt', np.vstack((ls, nlpp_TE_fullsky)).T)
np.savetxt('TB_lensing_noise.txt', np.vstack((ls, nlpp_TB_fullsky)).T)
np.savetxt('EB_lensing_noise.txt', np.vstack((ls, nlpp_EB_fullsky)).T)


# ql.spec.cl2cfft(cl_unl.clpp, ql.maps.cfft(nx, dx)).get_ml(lbins, t=t).plot(color='gray', ls='--')

# nlpp_TE_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')
