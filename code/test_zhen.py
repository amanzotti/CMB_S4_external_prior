'''
Script written to test Zhen Pan Fisher matrix and compare with mine.
'''

import numpy as np
import math
import matplotlib.pyplot as plt
import utils
import pickle
import sys

# my fid

fid = pickle.load(open('data/{}/run1/fid_values.p'.format('test_fisher'), "rb"))

index = {}

index['hubble'] = 3
index['scalar_spectral_index(1)'] = 1
index['scalar_amp(1)'] = 0
index['re_optical_depth'] = 2
index['ombh2'] = 4
index['omch2'] = 5
index['mnu'] = 6

new_order = [index[fid.keys()[0]], index[fid.keys()[1]], index[fid.keys()[2]], index[fid.keys()[3]],
             index[fid.keys()[4]], index[fid.keys()[5]], index[fid.keys()[6]]]
# import zhen fisher
fisher_zhen = np.genfromtxt('data/test_fisher/F_CMBS4.dat')

reordered_fisher = np.zeros_like(fisher_zhen)
reordered_fisher = fisher_zhen[:, new_order][new_order]


my_fisher_1 = np.genfromtxt('output/fisher_mat.txt')
my_fisher = np.genfromtxt('output/fisher_mat_run2.txt')
my_fisher = np.genfromtxt('output/fisher_mat_run3.txt')


