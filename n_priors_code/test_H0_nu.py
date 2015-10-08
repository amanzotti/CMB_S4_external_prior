#!/usr/bin/env python

__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"

import configparser
import subprocess
import numpy as np
import sys
import pickle
import collections
import os
import warnings

import cosmolopy
h = 0.73
cosmo = {'omega_M_0': 0.142968/h**2, 'omega_lambda_0': 1-0.142968/h**2, 'omega_k_0': 0.0, 'h': h}
# from eq

# astro-ph/0006089v3
# g==D\a


def g_hamilton(om, oml):
    return 5 * om / 2. / (om ** (4 / 7) - oml + (1 + om / 2.)(1 + oml / 70.))

# this is actually the same used by Hu in POWER SPECTRA FOR CDM AND VARIANTS


def g(z, om, oml):
    return om * (1 + z) ** 3 + (1 - om - ol) * (1 + z) ** 2 + oml


def D_hu(z):
    def om(z):
        return om(1 + z) ** 3 * g(z, om, oml) ** -2

    def omz(z):
        oml * g(z, om, oml) ** -2

    return (1 + 1100) / (1 + z) * 5 * om(z) / 2. * (om(z) ** (4 / 7) - oml(z) + (1 + om(z) / 2.) * (1 + oml(z) / 70)) ** -1
