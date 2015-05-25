'''

Simple fisher code to test prior effect on N effective neutrino estimates

TODO:

what are this data? documents better

create classes or at least external utilities. Maybe a better way to interface with cosmosis


Noise inspired by W. Hu et al 1402.4108 SNOWMASS paper

TODO:

noise calculatiions routines

'''

__author__ = "Y. Park"


import numpy as np
import math
import matplotlib.pyplot as plt

import os, copy, glob

# Simple noise generation

def bl(fwhm_arcmin, lmax):
    """ returns the map-level transfer function for a symmetric Gaussian beam.
         * fwhm_arcmin      - beam full-width-at-half-maximum (fwhm) in arcmin.
         * lmax             - maximum multipole.
    """
    ls = np.arange(0, lmax+1)
    return np.exp( -(fwhm_arcmin * np.pi/180./60.)**2 / (16.*np.log(2.)) * ls*(ls+1.) )

def nl(noise_uK_arcmin, fwhm_arcmin, lmax):
    """ returns the beam-deconvolved noise power spectrum in units of uK^2 for
          * noise_uK_arcmin - map noise level in uK.arcmin
          * fwhm_arcmin     - beam full-width-at-half-maximum (fwhm) in arcmin.
          * lmax            - maximum multipole.
    """
    return (noise_uK_arcmin * np.pi/180./60.)**2 / bl(fwhm_arcmin, lmax)**2

def noise_uK_arcmin(noise_uK_arcmin, fwhm_arcmin, lmax):
    """ returns noise_uK_arcmin given the number of detectors and the time of obsrvation
    """
    return (noise_uK_arcmin * np.pi/180./60.)**2 / bl(fwhm_arcmin, lmax)**2






#