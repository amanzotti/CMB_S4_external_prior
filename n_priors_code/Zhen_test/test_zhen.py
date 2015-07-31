'''
Script written to test Zhen Pan Fisher matrix and compare with mine.
'''

import numpy as np
import math
import matplotlib.pyplot as plt
import utils
import pickle
import sys
import n_priors_code.utils


# my fid
base_dir ='/home/manzotti/n_eff-dependence-on-prior/n_priors_code/'

fid = pickle.load(open(base_dir+'data/{}/run5/fid_values.p'.format('test_fisher'), "rb"))
values = pickle.load(open(base_dir+'data/{}/run5/grid_values.p'.format('test_fisher'), "rb"))
par_gaps = pickle.load(open(base_dir+'data/{}/run5/par_gaps.p'.format('test_fisher'), "rb"))
index = {}

index['hubble'] = 3
index['scalar_spectral_index(1)'] = 1
index['scalar_amp(1)'] = 0
index['re_optical_depth'] = 2
index['ombh2'] = 4
index['ombch2'] = 5
index['mnu'] = 6

new_order = [index[fid.keys()[0]], index[fid.keys()[1]], index[fid.keys()[2]], index[fid.keys()[3]],
             index[fid.keys()[4]], index[fid.keys()[5]], index[fid.keys()[6]]]
# import zhen fisher
fisher_zhen = np.genfromtxt(base_dir+'data/test_fisher/F_CMBS4.dat')

reordered_fisher = np.zeros_like(fisher_zhen)
reordered_fisher = fisher_zhen[:, new_order][new_order]
fisher_zhen = reordered_fisher

my_fisher_1 = np.genfromtxt(base_dir+'output/fisher_mat_test_run1.txt')
my_fisher_2 = np.genfromtxt(base_dir+'output/fisher_mat_test_run2.txt')
my_fisher_3 = np.genfromtxt(base_dir+'output/fisher_mat_test_run3.txt')
my_fisher_4 = np.genfromtxt(base_dir+'output/fisher_mat_test_run4.txt')
my_fisher_5 = np.genfromtxt(base_dir+'output/fisher_mat_test_run5.txt')


no_marg_1 = np.genfromtxt(base_dir+'output/no_marginalized_ell_test_run1.txt')
no_marg_2 = np.genfromtxt(base_dir+'output/no_marginalized_ell_test_run2.txt')
no_marg_4 = np.genfromtxt(base_dir+'output/no_marginalized_ell_test_run4.txt')
no_marg_3 = np.genfromtxt(base_dir+'output/no_marginalized_ell_test_run3.txt')
no_marg_5 = np.genfromtxt(base_dir+'output/no_marginalized_ell_test_run5.txt')

marg_1 = np.genfromtxt(base_dir+'output/marginalized_ell_test_run1.txt')
marg_2 = np.genfromtxt(base_dir+'output/marginalized_ell_test_run2.txt')
marg_4 = np.genfromtxt(base_dir+'output/marginalized_ell_test_run4.txt')
marg_3 = np.genfromtxt(base_dir+'output/marginalized_ell_test_run3.txt')
marg_5 = np.genfromtxt(base_dir+'output/marginalized_ell_test_run5.txt')

my_fisher_inv_1 = np.genfromtxt(base_dir+'output/invetered_fisher_test_run1.txt')
my_fisher_inv_2 = np.genfromtxt(base_dir+'output/invetered_fisher_test_run2.txt')
my_fisher_inv_3 = np.genfromtxt(base_dir+'output/invetered_fisher_test_run3.txt')
my_fisher_inv_4 = np.genfromtxt(base_dir+'output/invetered_fisher_test_run4.txt')
my_fisher_inv_5 = np.genfromtxt(base_dir+'output/invetered_fisher_test_run5.txt')

my_fisher_inv_1 = (np.linalg.inv(my_fisher_1))
my_fisher_inv_2 = (np.linalg.inv(my_fisher_2))
my_fisher_inv_3 = (np.linalg.inv(my_fisher_3))
my_fisher_inv_4 = (np.linalg.inv(my_fisher_4))
my_fisher_inv_5 = (np.linalg.inv(my_fisher_5))

fisher_zhen_inv = (np.linalg.inv(reordered_fisher))
sys.exit()
plt.figure()
plt.imshow(np.abs((my_fisher_1 - my_fisher_2) / my_fisher_1) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=3)
plt.title(r'Difference Fisher mine 1-2 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_12.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((my_fisher_1 - my_fisher_3) / my_fisher_1) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=10)
plt.title(r'Difference Fisher mine 1-3 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_13.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((my_fisher_1 - my_fisher_4) / my_fisher_1) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=10)
plt.title(r'Difference  Fisher mine 1-4 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_14.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((my_fisher_2 - my_fisher_3) / my_fisher_2) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=10)
plt.title(r'Difference Fisher mine 2-3 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_23.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((my_fisher_3 - my_fisher_4) / my_fisher_3) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=10)
plt.title(r'Difference  Fisher mine 3-4 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_34.pdf')
plt.close()

# Zhen comparison

plt.figure()
plt.imshow(np.abs((fisher_zhen - my_fisher_1) / fisher_zhen) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=100)
plt.title(r'Difference Fisher Zhen mine 1 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_zhen_1.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((fisher_zhen - my_fisher_2) / fisher_zhen) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=100)
plt.title(r'Difference Fisher Zhen mine 2 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_zhen_2.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((fisher_zhen - my_fisher_3) / fisher_zhen) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=100)
plt.title(r'Difference Fisher Zhen mine 3 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_zhen_3.pdf')
plt.close()

plt.figure()
plt.imshow(np.abs((fisher_zhen - my_fisher_4) / fisher_zhen) * 100.,
           cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=100)
plt.title(r'Difference Fisher Zhen mine 4 in $\%$')
plt.colorbar()
plt.savefig('Diff_mine_zhen_4.pdf')
plt.close()

# Now invert

plt.figure()
plt.imshow(np.sqrt(my_fisher_inv_1 / my_fisher_inv_2), cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
plt.title(r'Difference sqrt(inv(Fisher)) mine 1-2 ratio)')
plt.colorbar()
plt.savefig('Diff_ratio_inv_mine_12.pdf')
plt.close()

# plt.figure()
# plt.imshow(my_fisher_inv_1 / my_fisher_inv_3, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) mine 1-3 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_mine_13.pdf')
# plt.close()

# plt.figure()
# plt.imshow(my_fisher_inv_1 / my_fisher_inv_4, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) mine 1-4 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_mine_14.pdf')
# plt.close()

# plt.figure()
# plt.imshow(my_fisher_inv_2 / my_fisher_inv_3, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) mine 3-2 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_mine_23.pdf')
# plt.close()

# plt.figure()
# plt.imshow(my_fisher_inv_3 / my_fisher_inv_4, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) mine 3-4 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_mine_34.pdf')
# plt.close()

plt.figure()
plt.imshow(np.sqrt(fisher_zhen_inv / my_fisher_inv_1), cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
plt.title(r'Difference sqrt(inv(Fisher)) Zhen-1 ratio)')
plt.colorbar()
plt.savefig('Diff_ratio_inv_zhen_1.pdf')
plt.close()


plt.figure()
plt.imshow(np.sqrt(fisher_zhen_inv / my_fisher_inv_2), cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
plt.title(r'Difference sqrt(inv(Fisher)) Zhen-2 ratio)')
plt.colorbar()
plt.savefig('Diff_ratio_inv_zhen_2.pdf')
plt.close()

# plt.figure()
# plt.imshow(fisher_zhen_inv / my_fisher_inv_3, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) Zhen-3 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_zhen_3.pdf')
# plt.close()

# plt.figure()
# plt.imshow(fisher_zhen_inv / my_fisher_inv_4, cmap=plt.get_cmap('bwr'), interpolation='nearest', vmin=0, vmax=2)
# plt.title(r'Difference sqrt(inv(Fisher)) Zhen-4 ratio)')
# plt.colorbar()
# plt.savefig('Diff_ratio_inv_zhen_4.pdf')
# plt.close()