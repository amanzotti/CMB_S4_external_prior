import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
import sys


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

def load_data(run_idx,  values,lensed = True):
    '''
    build the dats matrix of data used in the fisher code. If lensed it is composed by lesned CMB cls + CMB lensing cls like cldd clde cldt

    the order is
    l CTT CEE CBB CTE Cdd CdT CdE
    0  1  2    3   4   5   6  7
    '''

    cmb = np.genfromtxt('../data/run{}/fiducial_lensedcls.dat'.format(run_idx), usecols=(0, 1, 2, 3, 4)
                        ) if lensed else np.genfromtxt('data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx), usecols=(0, 1, 2, 3, 4))
    lensing = np.genfromtxt('../data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx), usecols=(5, 6, 7))
    dats = np.genfromtxt('../data/run{}/fiducial_lenspotentialcls.dat'.format(run_idx))
    dats = np.concatenate((cmb, lensing[:cmb.shape[0], :]), axis=1)
    # Load data for all parameters variations
    for key, value in values.iteritems():
        for i in np.arange(0, 4):
            print key, values[key][i]
            filename = '../data/run{}/'.format(run_idx)
            # filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat'
            filename_cmb = filename + key + \
                '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat' if lensed else filename + \
                key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'

            filename_lensing = filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
            cmb = np.genfromtxt(filename_cmb.format(run_idx), usecols=(0, 1, 2, 3, 4))
            lensing = np.genfromtxt(filename_lensing.format(run_idx), usecols=(5, 6, 7))
            newdat = np.concatenate((cmb, lensing[:dats.shape[0], :]), axis=1)
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

        # print 'test = ', fisher1[fid.keys().index('re_optical_depth'),
        # fid.keys().index('re_optical_depth')],
        # fisher[fid.keys().index('re_optical_depth'),
        # fid.keys().index('re_optical_depth')] / (1 / (10 ** i *
        # fid['re_optical_depth']) ** 2)

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


def print_resume_stat(fisher,fid):
    fisher_inv =np.linalg.inv(fisher)
    for key, value in fid.iteritems():

        print 'sigma(', key, ')', np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]), '=', 100. * np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / fid[key], '%', "with no degeneracies", 1. / np.sqrt(fisher[fid.keys().index(key), fid.keys().index(key)])


def fisher_marginalize(fisher, marginalize_parameters_list, fid):
    '''
    This function is used to marginalize the fisher matrix following http://wfirst.gsfc.nasa.gov/science/fomswg/fomswg_technical.pdf eq (14) DETF.
    It is used to avoid ill-conditioned matrix to be inverted.

    Input:

    fisher= fisher matrix to marginalize
    list: parameters to be marginalized
    fid : describe the parameters


    Note:
    Fqq containes the parameters we want to keep.

    '''
    import scipy.linalg as linalg
    parameters_to_keep = list(set(fid.keys()) - set(marginalize_parameters_list))
    indeces_to_keep = []
    indeces_to_marg = []
    for key in parameters_to_keep:
        indeces_to_keep.append(fid.keys().index(key))
    #  build the 2 matrices
    indeces_to_marg = list(set(np.arange(0, np.size(fid.keys()))) - set(indeces_to_keep))

    Fqq = fisher[indeces_to_keep, :][:, indeces_to_keep]

    Frr = fisher[indeces_to_marg, :][:, indeces_to_marg]
    # U orthogonal matrix diagonalising V its transpose
    # s array of eigenvalues
    Fqr = fisher[indeces_to_keep, :][:, indeces_to_marg]
    U, s, V = linalg.svd(Frr, full_matrices=True)
    print 'marginalize matrix eig', s
    s_inv = np.diag(1/s)
    G = Fqq - np.dot( (np.dot(Fqr,U)) , np.dot(s_inv,((np.dot(Fqr,U)).T)))
    return G,parameters_to_keep

# ROUTINES TO STUDY INSTABILITY
def pert_element(fisher,pos,pert_ampl):
    b= np.zeros_like(fisher)
    b[pos] = np.random.rand()*pert_ampl
    return b

def cond_removing_el(fisher,fid):
    temp =np.linalg.cond(fisher)
    for i in np.arange(0,np.shape(fisher)[0]):
        survived = list(set(np.arange(0,np.shape(fisher)[0]))-set([i]))
        print 'element removed',i,fid.keys()[i], 'condition number',np.linalg.cond(fisher[survived,:][:,survived]),'improvement', np.linalg.cond(fisher[survived,:][:,survived])/temp

class Fisher_matrix(object):

    """docstring for Fisher_matrix"""

    def __init__(self, arg):
        super(Fisher_matrix, self).__init__()
        self.arg = arg
