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


parser = configparser.ConfigParser()
configfile='./lensing_noise.ini'
parser.read(configfile)
nlev_t = parser.getfloat('lensing_noise', 'nlev_t')
nlev_p = parser.getfloat('lensing_noise', 'nlev_p')
beam_fwhm = parser.getfloat('lensing_noise', 'beam_fwhm')
lmin = parser.getint('lensing_noise', 'lmin')
lmax = parser.getint('lensing_noise', 'lmax')


# calculation parameters.
# lmax = 3000  # maximum multipole for T, E, B and \phi.
nx = 512  # number of pixels for flat-sky calc.
dx = 1. / 60. / 180. * np.pi  # pixel width in radians.

# nlev_t = 5.   # temperature noise level, in uK.arcmin.
# nlev_p = 5.   # polarization noise level, in uK.arcmin.
# beam_fwhm = 1.   # Gaussian beam full-width-at-half-maximum.

cl_unl = ql.spec.get_camb_scalcl(fname='/home/manzotti/n_eff-dependence-on-prior/code/fiducial_scalcls.dat',lmax=lmax)  # unlensed theory spectra.
cl_len = ql.spec.get_camb_lensedcl(fname='/home/manzotti/n_eff-dependence-on-prior/code/fiducial_lensedcls.dat',lmax=lmax)  # lensed theory spectra.

bl = ql.spec.bl(beam_fwhm, lmax)  # transfer function.
pix = ql.maps.pix(nx, dx)

# noise spectra
nltt = (np.pi / 180. / 60. * nlev_t) ** 2 / bl ** 2
nlee = nlbb = (np.pi / 180. / 60. * nlev_p) ** 2 / bl ** 2

# signal spectra
sltt = cl_len.cltt
slte = cl_len.clte
slee = cl_len.clee
slbb = cl_len.clbb
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


# make plot
ls = np.arange(0, lmax + 1)
t = lambda l: (l * (l + 1.)) ** 2 / (2. * np.pi)  # scaling to apply to cl_phiphi when plotting.
lbins = np.linspace(2, lmax, 1)       # multipole bins.

np.savetxt('TT_lensing_noise.txt', np.vstack((ls,nlpp_TT_fullsky)).T )
np.savetxt('EE_lensing_noise.txt', np.vstack((ls,nlpp_EE_fullsky)).T )
np.savetxt('TE_lensing_noise.txt', np.vstack((ls,nlpp_TE_fullsky)).T )
np.savetxt('TB_lensing_noise.txt', np.vstack((ls,nlpp_TB_fullsky)).T )
np.savetxt('EB_lensing_noise.txt', np.vstack((ls,nlpp_EB_fullsky)).T )




# ql.spec.cl2cfft(cl_unl.clpp, ql.maps.cfft(nx, dx)).get_ml(lbins, t=t).plot(color='gray', ls='--')

# nlpp_TE_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')
