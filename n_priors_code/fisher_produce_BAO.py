
#!/usr/bin/env python


'''
produce Fisher matrix for BAO experiments
'''
__author__ = "A.Manzotti"
__license__ = "GPL"
__version__ = "2.0"
__maintainer__ = "A.Manzotti"
__email__ = "manzotti.alessandro@gmail.com"
__status__ = "Production"

import numpy as np
import math
import matplotlib.pyplot as plt
import utils
import pickle
import sys
import collections
from scipy.integrate import quad
from astropy.constants import hbar, k_B, G, c
import NuDensity as nd


# # ROUTINES TO COMPUTE WHAT BAO CONSTRAIN, rs/DV(z)

# # following EISENSTEIN & HU


# def nu_plus_r_density():
#     factor = (8.* np.pi**3 * k_B**4 * G ) / (45. * hbar**3 * c**5 *h*100)

#  from http://arxiv.org/pdf/1411.1074v3.pdf
# def integral_I(z,omnuh2,T_CMB):
#     T_nu = T_CMB *(4./11.)**(1/3.) * (3.046/3)**(1/4.)
#     m_nu =
#     r_variable= (m_nu  * c**2 /k_B/T_nu)
#     def integrand(x, r):
#         return 15./np.pi**4 * np.sqrt(x**2+r**2)*x**2/(np.exp(x)+1.)
#     return quad(integrand, 0, np.inf , args=(r_variable), epsrel=1.49e-11)[0]


def z_d(omch2, ombh2, omnuh2=0.0009):
    '''
    redshift drag from Eisenstein hu

    om is the total matter density omc+omb
    '''
    om0h2 = omch2 + ombh2 + omnuh2
    b1 = 0.313 * (om0h2)**(-0.419) * (1. + 0.607 * (om0h2)**(0.674))
    b2 = 0.238 * (om0h2)**(0.223)
    return 1291 * (om0h2**(0.251)) / (1. + 0.659 * om0h2**0.828) * (1. + b1 * ombh2**b2)


def k_eq(ombh2, omch2, Theta_cmb=2.725 / 2.7):
    '''
    From http://arxiv.org/pdf/astro-ph/9709112v1.pdf
    in Mpc-1 (no h in there)
    '''
    om0h2 = ombh2 + omch2
    return 7.461e-2 * om0h2 * Theta_cmb ** -2


def z_eq(ombh2, omch2, Theta_cmb=2.725 / 2.7):
    '''
    From http://arxiv.org/pdf/astro-ph/9709112v1.pdf
    '''
    om0h2 = ombh2 + omch2
    return 2.50 * 10**4 * om0h2 * Theta_cmb ** -4


def R(z, ombh2, T_cmb=2.725 / 2.7):
    return 31.5 * ombh2 * T_cmb**-4 * (z / 10**3)**-1


def r_s_approx(ombh2, omch2):
    '''
    from eq 26 of http://arxiv.org/pdf/astro-ph/9709112v1.pdf Eisenstein hu
    result is in Mpc and om0h2 is omc+omb

    possibly correct by 1.0395 empirically from CAMB
    '''
    om0h2 = omch2 + ombh2
    return 44.5 * np.log(9.83 / om0h2) / np.sqrt(1 + 10 * (ombh2)**(3. / 4.)) / 1.0395


def rs_over_Dv(z, ombh2, omch2, h0, w=-1, omnuh2=0):
    '''
    '''
    c = 1.
    return r_s_approx(ombh2, omch2) / Dv(z, ombh2, omch2, h0, w, omnuh2)


def H_z_r_s(z, ombh2, omch2, h0, w=-1, omnuh2=0):
    '''
    '''
    c = 1.
    return r_s_approx(ombh2, omch2) * H_over_h0(z, ombh2, omch2, h0, w, omnuh2) * 100.* h0

def Dv(z, ombh2, omch2, h0, w=-1, omnuh2=0):
    '''
    from eq 2 of
    1304.6984v2
    '''
    c = 1.
    return ((1. + z)**2 * Da(z, ombh2, omch2, h0, w, omnuh2)**2 * (3000 / h0) * z / (H_over_h0(z, ombh2, omch2, h0, w, omnuh2)))**(1. / 3.)


def cs(z, ombh2, T_cmb=2.725 / 2.7):
    c = 1.
    return c / np.sqrt(3. * (1. + R(ombh2, z, T_cmb=2.725 / 2.7)))


def r_d_approx_N_eff(ombh2, omch2, omnuh2=0., N_eff=3.046):
    '''
    from eq 17 of http://arxiv.org/pdf/1411.1074v3.pdf
    result is in Mpc and omcbh2 is omc+omb

    which is accurate to 0.119% if we restrict to neutrino
    masses in the range 0 <
    mv < 0.6 eV and 3 < Neff <5.
     '''
    omcbh2 = omch2 + ombh2
    return 56.067 * np.exp(-49.7 * (omnuh2 + 0.002)**2) / omcbh2**0.2436 / ombh2**0.128876 / (1 + (N_eff - 3.046) / 30.60)


def H_over_h0(z, ombh2, omch2, h0, w=-1, omnuh2=0.):
    om = (omch2 + ombh2 + omnuh2) / h0**2
    om_rad = 5.46 * 10**-5
    return np.sqrt((om) * (1 + z)**3 + om_rad ** (1 + z)**4 + (1 - om) * (1. + z) ** (3 * (1 + w)))


def Da(z, ombh2, omch2, h0, w=-1, omnuh2=0.):
    '''comoving distance'''
    from scipy.integrate import quad

    def integrand(z, ombh2, omch2, h0, w, omnuh2=0.):
        return (3000 / h0) / (H_over_h0(z, ombh2, omch2, h0, w, omnuh2))
    return quad(integrand, 0, z, args=(ombh2, omch2, h0, w, omnuh2), epsrel=1.49e-11)[0] / (1. + z)


def partial_derivative(func, var, gaps, point):
    from scipy.misc import derivative
    args = point[:]

    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx=point[var] * 0.1, order=7)


# BAO experiment
# errors on rs/Dv
sigma_BAO = [0.00071, 0.0015, 0.00034, 0.00032]
# redshift of rs/Dv measurametns
z_BAO = [0.44, 0.106, 0.60, 0.73]
# parameters

# TODO LOAD EVERYTHING FROM INI
# =============================

data_folder = 'varying_all/run4'
output_folder = 'varying_all/run4/output'


# ===================
header = 'BAO fisher from fisher_produce_BAO.py \n'
# header += 'lmax={} \n lmin={} \n l_t_max={} \n fsky={} \n lensed={} \n data_folder={} \n N_det={} \n'.format(
#     lmax, lmin, l_t_max, fsky, lensed, data_folder, N_det)


# READ PARAMS
# load fiducial data
# load fiducial parameters used
fid = pickle.load(open('data/{}/fid_values.p'.format(data_folder), "rb"))

print "fid ", fid

# load parameter grid dictionary. The format is a pickle
values = pickle.load(open('data/{}/grid_values.p'.format(data_folder), "rb"))
par_gaps = pickle.load(open('data/{}/par_gaps.p'.format(data_folder), "rb"))


# you just need ombh2,omch2,h0
exclude = ['helium_fraction', 'massless_neutrinos', 'omk', 're_optical_depth',
           'scalar_amp(1)', 'scalar_nrun(1)', 'scalar_spectral_index(1)', 'wa','omnuh2','w']
par_gaps, values, fid = utils.exclude_parameters(exclude, par_gaps, values, fid)
n_values = np.size(values.keys())


fisher = np.zeros((n_values, n_values))
fisher_inv = np.zeros((n_values, n_values))

print 'fisher_size', fisher.shape
var = {}

# WIGGLE Z + SDSS + 6dFGS covariance from
# http://arxiv.org/pdf/1212.5226.pdf
# http://arxiv.org/pdf/1304.6984v2.pdf

# errors on rs/Dv
# ==================================================================
#========================================================================


# redshift of rs/Dv measurametns
z_BAO = np.array([0.1, 0.35, 0.57, 0.44, 0.60, 0.73])  # table 1 of 1212.5226 6dFGS + SDSS DR7,DR9 + 3 WiggleZ point
# parameters
cov_BAO_inv = np.array([[4444.4, 0, 0, 0, 0, 0], [0, 250000, 0, 0, 0, 0], [0, 0, 1000000.0, 0, 0, 0], [
                       0, 0, 0, 24532.1, -25137.7, 12099.1], [0, 0, 0, -25137.7, 134598.4, -64783.9], [0, 0, 0, 12099.1, -64783.9, 128837.6]])


var['hubble'] = 3
var['omch2'] = 2
var['ombh2'] = 1
var['w'] = 4
var['omnuh2'] = 5


for i in range(0, n_values):

    for j in range(0, n_values):
        # computing derivatives.
        # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
              # ---------------------------------
                          #   12h
        if fid.keys()[i] != 'hubble' and fid.keys()[i] != 'ombh2' and fid.keys()[i] != 'omch2' :#and fid.keys()[i] != 'w' :#and fid.keys()[i] != 'omnuh2':
            fisher[i, j] += 0
            continue

        if fid.keys()[j] != 'hubble' and fid.keys()[j] != 'ombh2' and fid.keys()[j] != 'omch2' :#and fid.keys()[j] != 'w' :#and fid.keys()[i] != 'omnuh2':
            fisher[i, j] += 0
            continue
        else:

            # print fid.keys()[j],var[fid.keys()[j]],par_gaps[fid.keys()[j]]
            # print fid.keys()[i],var[fid.keys()[i]],par_gaps[fid.keys()[i]]

            c_i = np.array([partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
                i]], gaps=par_gaps[fid.keys()[i]], point=[z, fid['ombh2'], fid['omch2'], fid['hubble']]) for z in z_BAO])
            print c_i
            c_j = np.array([partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
                j]], gaps=par_gaps[fid.keys()[j]], point=[z, fid['ombh2'], fid['omch2'], fid['hubble']]) for z in z_BAO])

            fisher[i, j] += np.sum(np.dot(c_i,  cov_BAO_inv) * c_j.T)


# # BOSS CMASS BE carefull

# BOSS is technically part of the fisher above




# #redshift of rs/Dv measurametns
# z_boss = 0.57
# # #parameters

# args = [z_boss,fid['ombh2'], fid['omch2'], fid['hubble']]
# sigma_DA= 45/1408.* Da(args[0],args[1],args[2],args[3])
# sigma_H= 7.8/98. * H_over_h0(args[0],args[1],args[2],args[3])
# cov_DA_H = 0.55 * sigma_DA*sigma_H # correlation taken from footnotes in Table 1 of http://arxiv.org/pdf/1304.6984v2.pdf
# cov_BOSS = np.array([[sigma_DA**2,cov_DA_H],[cov_DA_H,sigma_H**2]]) # order is DA,H
# cov_BOSS_inv = np.linalg.inv(cov_BOSS)
# var['hubble'] = 3
# var['omch2'] = 2
# var['ombh2'] = 1
# for z,sigma in zip(z_BAO,sigma_BAO):
#     print z,sigma
#     args = [z,fid['ombh2'], fid['omch2'], fid['hubble']]

#     for i in range(0, n_values):

#         for j in range(0, n_values):
#             # computing derivatives.
#             # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
#                   # ---------------------------------
#                               #   12h
#             if fid.keys()[i] != 'hubble' and fid.keys()[i] != 'ombh2' and fid.keys()[i] != 'omch2':
#                 fisher[i, j] += 0
#                 continue

#             if fid.keys()[j] != 'hubble' and fid.keys()[j] != 'ombh2' and fid.keys()[j] != 'omch2':
#                 fisher[i, j] += 0
#                 continue
#             else:

#                 # print fid.keys()[j],var[fid.keys()[j]],par_gaps[fid.keys()[j]]
#                 # print fid.keys()[i],var[fid.keys()[i]],par_gaps[fid.keys()[i]]

#                 c_i = np.array( [ partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
#                                            i]], gaps=par_gaps[fid.keys()[i]], point=[z_boss,fid['ombh2'], fid['omch2'], fid['hubble']]),
#                 partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
# i]], gaps=par_gaps[fid.keys()[i]], point=[z_boss,fid['ombh2'],
# fid['omch2'], fid['hubble']])])

#                 c_j = np.array( [ partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
#                                            j]], gaps=par_gaps[fid.keys()[j]], point=[z_boss,fid['ombh2'], fid['omch2'], fid['hubble']]),
#                 partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
# j]], gaps=par_gaps[fid.keys()[j]], point=[z_boss,fid['ombh2'],
# fid['omch2'], fid['hubble']])])

#                 fisher[i, j] += np.sum(np.dot(np.dot(cov_BOSS_inv, c_i),  cov_BOSS_inv) * c_j.T)


# var['hubble'] = 3
# var['omch2'] = 2
# var['ombh2'] = 1
# var['w'] = 4
# var['omnuh2'] = 5

# for z, sigma in zip(z_BAO, sigma_BAO):
#     print z, sigma
#     args = [z, fid['ombh2'], fid['omch2'], fid['hubble']]

#     for i in range(0, n_values):

#         for j in range(0, n_values):
#             # computing derivatives.
#             # f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
#                   # ---------------------------------
#                               #   12h
#             if fid.keys()[i] != 'hubble' and fid.keys()[i] != 'ombh2' and fid.keys()[i] != 'omch2' :#and fid.keys()[i] != 'w' :#and fid.keys()[i] != 'omnuh2':
#                 fisher[i, j] += 0
#                 continue

#             if fid.keys()[j] != 'hubble' and fid.keys()[j] != 'ombh2' and fid.keys()[j] != 'omch2' :#and fid.keys()[j] != 'w' :#and fid.keys()[i] != 'omnuh2':
#                 fisher[i, j] += 0
#                 continue
#             else:

#                 # print fid.keys()[j],var[fid.keys()[j]],par_gaps[fid.keys()[j]]
#                 # print fid.keys()[i],var[fid.keys()[i]],par_gaps[fid.keys()[i]]

#                 c_i = np.array([partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
#                     i]], gaps=par_gaps[fid.keys()[i]], point=[z, fid['ombh2'], fid['omch2'], fid['hubble']]) for z in z_BAO])
#                 c_j = np.array([partial_derivative(func=rs_over_Dv, var=var[fid.keys()[
#                     j]], gaps=par_gaps[fid.keys()[j]], point=[z, fid['ombh2'], fid['omch2'], fid['hubble']]) for z in z_BAO])

#                 fisher[i, j] += np.sum(np.dot(c_i,  cov_BAO_inv) * c_j.T)



print fisher

# print 'finally how much constraint on parameters without prior?'
# print ''
fisher_single = fisher.copy()
fisher_single[fid.keys().index('ombh2'), fid.keys().index('ombh2')] += 1 / \
    (0.04 * fid['ombh2']) ** 2
# fisher_inv = np.linalg.inv(fisher_single)

# print 'fisher=', fisher

utils.print_resume_stat(fisher_single, fid)
