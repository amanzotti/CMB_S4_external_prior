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
    ls = np.arange(0, lmax + 1)
    return np.exp(-(fwhm_arcmin * np.pi / 180. / 60.) ** 2 / (16. * np.log(2.)) * ls * (ls + 1.))


def nl(noise_uK_arcmin, fwhm_arcmin, lmax):
    """ returns the beam-deconvolved noise power spectrum in units of uK^2 for
          * noise_uK_arcmin - map noise level in uK.arcmin
          * fwhm_arcmin     - beam full-width-at-half-maximum (fwhm) in arcmin.
          * lmax            - maximum multipole.
    """
    return (noise_uK_arcmin * np.pi / 180. / 60.) ** 2 / bl(fwhm_arcmin, lmax) ** 2


def noise_uK_arcmin(noise_uK_arcmin, fwhm_arcmin, lmax):
    """ returns noise_uK_arcmin given the number of detectors and the time of obsrvation
    """
    return (noise_uK_arcmin * np.pi / 180. / 60.) ** 2 / bl(fwhm_arcmin, lmax) ** 2


def load_data(data_folder,  values, lensed=False, verbose = True):
    '''
    build the dats matrix of data used in the fisher code. If lensed it is composed by lesned CMB cls + CMB lensing cls like cldd clde cldt

    the order is
    l CTT CEE CBB CTE Cdd CdT CdE
    0  1  2    3   4   5   6  7
    '''

    cmb = np.genfromtxt('data/{}/fiducial_lensedcls.dat'.format(data_folder), usecols=(0, 1, 2, 3, 4)
                        ) if lensed else np.genfromtxt('data/{}/fiducial_lenspotentialcls.dat'.format(data_folder), usecols=(0, 1, 2, 3, 4))
    lensing = np.genfromtxt('data/{}/fiducial_lenspotentialcls.dat'.format(data_folder), usecols=(5, 6, 7))
    dats = np.loadtxt('data/{}/fiducial_lenspotentialcls.dat'.format(data_folder))
    dats = np.concatenate((cmb, lensing[:cmb.shape[0], :]), axis=1)
    # Load data for all parameters variations
    for key, value in values.iteritems():
        if np.shape(values[values.keys()[0]])[0]!=4:
            if verbose: print "you are not using a 5 point formula are you sure?"
        for i in np.arange(0, np.shape(values[values.keys()[0]])[0]):
            if verbose: print key, values[key][i],i
            filename = 'data/{}/'.format(data_folder)
            # filename_cmb = filename + key + '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat'
            filename_cmb = filename + key + \
                '_{:.13f}'.format(values[key][i]) + '_lensedcls.dat' if lensed else filename + \
                key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'

            filename_lensing = filename + key + '_{:.13f}'.format(values[key][i]) + '_lenspotentialcls.dat'
            cmb = np.genfromtxt(filename_cmb.format(data_folder), usecols=(0, 1, 2, 3, 4))
            lensing = np.genfromtxt(filename_lensing.format(data_folder), usecols=(5, 6, 7))
            newdat = np.concatenate((cmb, lensing[:dats.shape[0], :]), axis=1)
            # newdat = np.genfromtxt(filename)
            dats = np.dstack((dats, newdat))
    np.save('data/{}/fisher_data.npy'.format(data_folder),dats)
    return dats


def study_prior_H0_on_N_eff():
    '''
    TO be finished we want to study the effect of priors on parameters manipulating the fisher matrix
    '''

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
            np.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

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
            np.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

        fisher3 = fisher.copy()[[fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')], :][
            :, [fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')]]
        # Cicle on H0 priors
        # in the cut matrix tau is in the 0 place
        fisher3[0, 0] += 1 / (10 ** i * fid['re_optical_depth']) ** 2

        # Invert and get Neff error with these priors
        d3.append(
            np.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

        np.savetxt('output_cmb/sigma_tau_1percent.txt', d2)
        np.savetxt('output_cmb/sigma_tau_noPrior.txt', d)
        np.savetxt('output_cmb/sigma_tau_perfect_prior.txt', d3)


def study_prior_tau_on_N_eff(fid, fisher, output_folder, header):
    '''
    TO be finished we want to study the effect of priors on parameters manipulating the fisher matrix
    '''
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
            np.sqrt(np.linalg.inv(fisher1)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

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
            np.sqrt(np.linalg.inv(fisher2)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

        fisher3 = fisher.copy()[[fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')], :][
            :, [fid.keys().index('re_optical_depth'), fid.keys().index('massless_neutrinos')]]
        # Cicle on H0 priors
        # in the cut matrix tau is in the 0 place
        fisher3[0, 0] += 1 / (10 ** i * fid['re_optical_depth']) ** 2

        # Invert and get Neff error with these priors
        d3.append(
            np.sqrt(np.linalg.inv(fisher3)[fid.keys().index('massless_neutrinos'), fid.keys().index('massless_neutrinos')]))

        np.savetxt('{}/sigma_tau_1percent.txt'.format(output_folder), d2, header=header)
        np.savetxt('{}/sigma_tau_noPrior.txt'.format(output_folder), d, header=header)
        np.savetxt('{}/sigma_tau_perfect_prior.txt'.format(output_folder), d3, header=header)


def return_simgax_y_prior(fid, fisher, x, y, prior_val):
    '''
    Concept return sigma x given a possible prior on y

    As an input it gets the fisher matrix but what it spits out is not a relative sigma.
    just sigma.


    '''

    print "CHECK"
    if isinstance(prior_val, (int, long, float)):
        assert (fid.has_key(x))
        assert (fid.has_key(y))

        fisher1 = fisher.copy()
        # Cicle on H0 priors
        fisher1[fid.keys().index(y), fid.keys().index(y)] += 1 / \
            (prior_val * fid[y]) ** 2
        # print np.diag(fisher1)

        return np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])

    elif (isinstance(prior_val, (np.ndarray))) or (isinstance(prior_val, list)):

        sigma_x_prior = np.zeros(np.size(prior_val))
        for i, prior in enumerate(prior_val):
            # print 'prior in single',prior
            # print 'fisher', fisher[0,0]
            fisher1 = fisher.copy()
            # Cicle on H0 priors

            fisher1[fid.keys().index(y), fid.keys().index(y)] += 1 / \
                (prior * fid[y]) ** 2

            # print prior,y
            print np.linalg.inv(fisher1).diagonal()/np.linalg.inv(fisher).diagonal()
            # print ''
            sigma_x_prior[i] = np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])
            # print 'fisher', fisher[0,0]
            # print ''
#
        return sigma_x_prior
    else:
        sys.exit('type not recognized')


# def return_simgax_y_prior(fid, fisher, x, y, prior_val):
#     '''
#     Concept return sigma x given a possible prior on y

#     As an input it gets the fisher matrix but what it spits out is not a relative sigma.
#     just sigma.


#     '''

#     if isinstance(prior_val, (int, long, float)):
#         assert (fid.has_key(x))
#         assert (fid.has_key(y))

#         fisher1 = fisher.copy()
#         # Cicle on H0 priors
#         fisher1[fid.keys().index(y), fid.keys().index(y)] += 1 / \
#             (prior_val * fid[y]) ** 2
#         print np.diag(fisher1)

#         return np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])

#     elif (isinstance(prior_val, (np.ndarray))) or (isinstance(prior_val, list)):

#         sigma_x_prior = np.zeros(np.size(prior_val))
#         for i, prior in enumerate(prior_val):
#             # print 'prior in single',prior
#             fisher1 = fisher.copy()
#             # Cicle on H0 priors

#             fisher1[fid.keys().index(y), fid.keys().index(y)] += 1 / \
#                 (prior * fid[y]) ** 2

#             sigma_x_prior[i] = np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])
#         return sigma_x_prior
#     else:
#         sys.exit('type not recognized')


def apply_simgax_y_prior(fid, fisher, x, y, prior_val):
    '''
    Concept return sigma x given a possible prior on y

    As an input it gets the fisher matrix but what it spits out is not a relative sigma.
    just sigma.

    IMPORTANT Note: here the prior is applied in place so the fisher matrix itsels id modified

    '''

    fisher_mat_multi = np.array([fisher for i in range(np.size(prior_val))])

    fisher_mat_multi[:, fid.keys().index(y), fid.keys().index(y)] += 1 / \
        (prior_val * fid[y]) ** 2

    return np.sqrt(np.linalg.inv(fisher_mat_multi)[:, fid.keys().index(x), fid.keys().index(x)])

    # if isinstance(prior_val, (int, long, float)):
    #     assert (fid.has_key(x))
    #     assert (fid.has_key(y))

    # fisher1 = fisher.copy()
    # Cicle on H0 priors
    #     fisher[fid.keys().index(y), fid.keys().index(y)] += 1 / \
    #         (prior_val * fid[y] ) ** 2
    # print np.diag(fisher)

    #     return np.sqrt(np.linalg.inv(fisher)[fid.keys().index(x), fid.keys().index(x)])

    # elif (isinstance(prior_val, (np.ndarray))) or (isinstance(prior_val, list)):

    #     sigma_x_prior = np.zeros(np.size(prior_val))
    #     for i, prior in enumerate(prior_val):
    # print 'prior in single',prior
    # fisher1 = fisher.copy()
    # Cicle on H0 priors

    #         fisher[fid.keys().index(y), fid.keys().index(y)] += 1 / \
    #             (prior * fid[y] ) ** 2

    #         sigma_x_prior[i] = np.sqrt(np.linalg.inv(fisher)[fid.keys().index(x), fid.keys().index(x)])
    #     return sigma_x_prior
    # else:
    #     sys.exit('type not recognized')


def return_simgax_all_prior(fid, fisher_mat, key_y):
    '''
    Concept return sigma x given a possible prior on *ALL* the other parameters

    As an input it gets the fisher matrix but what it spits out is not a relative sigma.
    just sigma.

    NOTE : *IMPORTANT*  here the applieded prior is the same for everyone


    '''
    fisher_inv = np.linalg.inv(fisher_mat)
    prior_size = 900
    fisher_mat_multi = np.array([fisher_mat for i in range(prior_size)])
    # this is  s [prior_size,9,9] 3d matrix
    for i, key in enumerate(fid.keys()):

        if key == key_y:
            continue
        # print key

        sigma_just_CMB_x = (np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / np.abs(fid[key]))
        prior_value = np.linspace(sigma_just_CMB_x / 10., sigma_just_CMB_x * 4.5, prior_size)
        # if you have an array of prior to apply, cycle on them.
        fisher_mat_multi[:, fid.keys().index(key), fid.keys().index(key)] += 1 / \
            (prior_value * fid[key]) ** 2

    return np.sqrt(np.linalg.inv(fisher_mat_multi)[:, fid.keys().index(key_y), fid.keys().index(key_y)])


def return_simgax_list_prior(fid, fisher_mat, key_prior_list,key_y):
    '''
    Concept return sigma x given a possible prior on list of the other parameters

    As an input it gets the fisher matrix but what it spits out is not a relative sigma.
    just sigma.

    NOTE : *IMPORTANT*  here the applieded prior is the same for everyone


    '''
    fisher_inv = np.linalg.inv(fisher_mat)
    prior_size = 900
    fisher_mat_multi = np.array([fisher_mat for i in range(prior_size)])
    # this is  s [prior_size,9,9] 3d matrix
    for i, key in enumerate(key_prior_list):

        if key == key_y:
            continue
        # print key

        sigma_just_CMB_x = (np.sqrt(fisher_inv[fid.keys().index(key), fid.keys().index(key)]) / np.abs(fid[key]))
        prior_value = np.linspace(sigma_just_CMB_x / 10., sigma_just_CMB_x * 4.5, prior_size)
        # if you have an array of prior to apply, cycle on them.
        fisher_mat_multi[:, fid.keys().index(key), fid.keys().index(key)] += 1 / \
            (prior_value * fid[key]) ** 2

    return np.sqrt(np.linalg.inv(fisher_mat_multi)[:, fid.keys().index(key_y), fid.keys().index(key_y)])




def return_simgax_y_prior2D(fid, fisher, x, y, prior_val1, prior_val2):
    '''
    Concept return sigma x given a possible prior on y
    '''

    par1 = y[0]
    par2 = y[1]
    if isinstance(prior_val1, (int, long, float)):
        assert (fid.has_key(x))
        assert (fid.has_key(y))
        fisher1 = fisher.copy()
        fisher1[fid.keys().index(par1), fid.keys().index(par1)] += 1 / \
            (prior_val1 * fid[par1]) ** 2
        fisher1[fid.keys().index(par2), fid.keys().index(par2)] += 1 / \
            (prior_val2 * fid[par1]) ** 2

        return np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])

    elif (isinstance(prior_val1, (np.ndarray))) or (isinstance(prior_val1, list)):

        sigma_x_prior = np.zeros((np.size(prior_val1), np.size(prior_val2)))
        for i, prior_1 in enumerate(prior_val1):
            for j, prior_2 in enumerate(prior_val2):

                fisher1 = fisher.copy()
                # Cicle on H0 priors
                fisher1[fid.keys().index(par1), fid.keys().index(par1)] += 1 / \
                    (prior_1 * fid[par1]) ** 2
                fisher1[fid.keys().index(par2), fid.keys().index(par2)] += 1 / \
                    (prior_2 * fid[par1]) ** 2

                sigma_x_prior[i, j] = np.sqrt(np.linalg.inv(fisher1)[fid.keys().index(x), fid.keys().index(x)])
        return sigma_x_prior


def study_prior_ns_on_N_eff():
    '''
    TO be finished we want to study the effect of priors on parameters manipulating the fisher matrix
    '''
    pass


# def study_prior_x_on_y(x, y, fid, fisher, output_folder):
#     '''
#     TO be finished we want to study the effect of priors on parameters manipulating the fisher matrix
#     '''

#     d = []
#     d2 = []
#     d3 = []

# cycling on x prior
#     for i in np.arange(-3, -1, 0.1):

#         fisher1 = fisher.copy()
# Cicle on x priors
#         fisher1[fid.keys().index(x)), fid.keys().index(x)] += 1 /
#             (10 ** i * fid[x]) ** 2
# Invert and get y error with these priors

#         d.append(
#             math.sqrt(np.linalg.inv(fisher1)[fid.keys().index(y), fid.keys().index(y)]))

#         fisher2=fisher.copy()
# Cicle on x priors

#         fisher2[fid.keys().index(x), fid.keys().index(x)] += 1 /
#             (10 ** i * fid[x]) ** 2

# add 1% prior on ns
#         fisher2[fid.keys().index('scalar_spectral_index(1)'), fid.keys().index('scalar_spectral_index(1)')] += 1 /
#             (0.01 * fid['scalar_spectral_index(1)']) ** 2
# add 1% prior on As
#         fisher2[fid.keys().index('scalar_amp(1)'), fid.keys().index('scalar_amp(1)')] += 1 /
#             (0.01 * fid['scalar_amp(1)']) ** 2
#         fisher2[fid.keys().index('hubble'), fid.keys().index('hubble')] += 1 / (0.01 * fid['hubble']) ** 2

# Invert and get y error with these priors
#         d2.append(
#             math.sqrt(np.linalg.inv(fisher2)[fid.keys().index(y), fid.keys().index(y)]))

#         fisher3=fisher.copy()[[fid.keys().index(x), fid.keys().index(y)], :][
#             :, [fid.keys().index(x), fid.keys().index(y)]]
# Cicle on x priors
# in the cut matrix tau is in the 0 place
#         fisher3[0, 0] += 1 / (10 ** i * fid[x]) ** 2

# Invert matrix and get error  on y with perfect priors
#         d3.append(
#             math.sqrt(np.linalg.inv(fisher3)[fid.keys().index(y), fid.keys().index(y)]))

#         np.savetxt('{}/sigma_{}_1percent.txt'.format(output_folder, y), d2, header = 'test_header')
#         np.savetxt('{}/sigma_{}_noPrior.txt'.format(output_folder, y), d, header = 'test_header')
#         np.savetxt('{}/sigma_{}_perfect_prior.txt'.format(output_folder, y), d3, header = 'test_header')


def save_cov_matrix(fisher_inv, filename='output_cmb/param_cov.txt'):
    n_values = np.shape(fisher_inv)[0]
    param_cov = np.zeros((n_values, n_values))
    for i in range(n_values):
        for j in range(n_values):
            if i != j:
                param_cov[i, j] = fisher_inv[i, j] / np.sqrt(fisher_inv[i, i] * fisher_inv[j, j])
    np.savetxt(filename, param_cov)


def exclude_parameters(excluded_parameters, par_gaps, values, fid):
    '''
    easily exclude parameters for the analysis
    '''
    if excluded_parameters == None:
        return par_gaps, values, fid
    else:
        for e in excluded_parameters:
            par_gaps.pop(e)
            values.pop(e)
            fid.pop(e)

        return par_gaps, values, fid


def exclude_parameters_from_fisher(excluded_parameters, par_gaps, values, fid, fisher_mat):
    '''
    easily exclude parameters for the analysis
    '''
    if excluded_parameters == None:
        return par_gaps, values, fid, fisher_mat
    else:
        for e in excluded_parameters:
            fisher_mat = np.delete(np.delete(fisher_mat, fid.keys().index(e), 0), fid.keys().index(e), 1)
            par_gaps.pop(e)
            values.pop(e)
            fid.pop(e)
        return par_gaps, values, fid, fisher_mat


def print_resume_stat(fisher, fid):
    fisher_inv = np.linalg.inv(fisher)
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
    s_inv = np.diag(1 / s)
    G = Fqq - np.dot((np.dot(Fqr, U)), np.dot(s_inv, ((np.dot(Fqr, U)).T)))
    red_fid = {}
    for key in parameters_to_keep:
        red_fid[key] = fid[key]

    return G, red_fid

# ROUTINES TO STUDY INSTABILITY


def pert_element(fisher, pos, pert_ampl):
    b = np.zeros_like(fisher)
    b[pos] = np.random.rand() * pert_ampl
    return b


def cond_removing_el(fisher, fid):
    temp = np.linalg.cond(fisher)
    for i in np.arange(0, np.shape(fisher)[0]):
        survived = list(set(np.arange(0, np.shape(fisher)[0])) - set([i]))
        print 'element removed', i, fid.keys()[i], 'condition number', np.linalg.cond(fisher[survived, :][:, survived]), 'improvement', temp / np.linalg.cond(fisher[survived, :][:, survived])


class Fisher_matrix(object):

    """docstring for Fisher_matrix"""

    def __init__(self, arg):
        super(Fisher_matrix, self).__init__()
        self.arg = arg
