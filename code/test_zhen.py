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

#         As      ns      tau       h    omgb      omgbc      mnu
#par  = [2.215, 0.9619, 0.0925, 0.6704, 0.022068, 0.142968, 0.085]
# fid = {}
# fid['hubble'] =  0.6704
# fid['scalar_spectral_index(1)'] =  0.9619
# fid['scalar_amp(1)'] = np.log(10**10*config.getfloat('camb', 'scalar_amp(1)'))
# fid['re_optical_depth'] = 0.0925
# fid['ombh2'] = 0.022068
# fid['omch2'] = 0.142968
# fid['omnuh2'] = 0.009104624166
# fid = collections.OrderedDict(sorted(fid.items(), key=lambda t: t[0]))

index = {}

index['hubble'] =  3
index['scalar_spectral_index(1)'] =  1
index['scalar_amp(1)'] = 0
index['re_optical_depth'] = 2
index['ombh2'] = 4
index['omch2'] = 5
index['omnuh2'] = 6

new_order = [index[fid.keys()[0]],index[fid.keys()[1]],index[fid.keys()[2]],index[fid.keys()[3]],index[fid.keys()[4]],index[fid.keys()[5]],index[fid.keys()[6]]]
# import zhen fisher
fisher_zhen = np.genfromtxt('data/test_fisher/F_CMBS4.dat')

reordered_fisher = np.zeros_like(fisher)
reordered_fisher= fisher[new_order,new_order]
