import numpy as np
import math
import matplotlib.pyplot as plt
import util
import pickle
import sys



def load_data(run_idx,lensed,values):

    '''
    build the dats matrix of data used in the fisher code. If lensed it is composed by lesned CMB cls + CMB lensing cls like cldd clde cldt

    the order is
    l CTT CEE CBB CTE Cdd CdT CdE
    0  1  2    3   4   5   6  7
    '''

    cmb = np.genfromtxt('data/run{}/fiducial_lensedcls.dat'.format(run_idx),usecols=(0,1,2,3,4)) if lensed else np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx),usecols=(0,1,2,3,4))
    lensing = np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx),usecols=(5,6,7))
    dats = np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx))
    dats = np.concatenate((cmb,lensing[:cmb.shape[0],:]),axis=1)
    # Load data for all parameters variations
    for key, value in values.iteritems():
        for i in np.arange(0, 4):
            print key, values[key][i]
            filename = 'data/run{}/'.format(run_idx)
            # filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat'
            filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat' if lensed else filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'

            filename_lensing = filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
            cmb = np.genfromtxt(filename_cmb.format(run_idx),usecols=(0,1,2,3,4))
            lensing = np.genfromtxt(filename_lensing.format(run_idx),usecols=(5,6,7))
            newdat = np.concatenate((cmb,lensing[:dats.shape[0],:]),axis=1)
            # newdat = np.genfromtxt(filename)
            dats = np.dstack((dats, newdat))

    return dats

def study_prior_H0_on_N_eff():
    pass

def study_prior_tau_on_N_eff():
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

def study_prior_ns_on_N_eff():
    pass

def save_cov_matrix(filename='output_cmb/param_cov.txt'):

    param_cov = np.zeros((n_values, n_values))
    for i in range(n_values):
        for j in range(n_values):
            if i != j:
                param_cov[i, j] = fisher_inv[i, j] / np.sqrt(fisher_inv[i, i] * fisher_inv[j, j])
    # print param_cov
    np.savetxt(filename, param_cov)

def print_resume_stat(values):
    for key, value in values.iteritems():

        print 'sigma(', key, ')', np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / fid[key], '%', "with no degeneracies", 1. / np.sqrt(fisher[fid.keys().index(key), fid.keys().index(key)])


