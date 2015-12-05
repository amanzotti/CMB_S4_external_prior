
#!/usr/bin/env python

'''
Script to produced the different Cls with different parameters needed to compute derivative for the Fisher matrices.

It changes the ini file of CAMB at each step with configparser and it then runs camb.

For better accuracy the derivative are computed with a 5 point stencil algorithms

f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
      ---------------------------------
                    12h

CONVENTIONS:

PARAMETER ORDER = we are using ordered dictionaries. So the order is the alphabetical order of the name of the variable sin the camb ini.

THIS IA DIFFERENT FROM THE GENERATE DATA BECAUSE IN THIS CASE WE ARE CHANGING LAMBDA TO GET THE UNIVERSE FLAT AND NOT OMEGA M


'''
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

output_folder = 'varying_all'
output_folder_2 = 'dunkley'
camb_location = '/home/manzotti/local/camb2013/camb'

config = configparser.ConfigParser()
configfile = './fiducial_dunkley.ini'
config.read(configfile)
config.set('camb', 'output_root',
           './data/{}/{}/'.format(output_folder, output_folder_2) + 'fiducial')

if not os.path.exists('./data/{}/{}/'.format(output_folder, output_folder_2)):
    os.makedirs('./data/{}/{}/'.format(output_folder, output_folder_2))


# Change output in fiducial to save it in the right folder
with open(configfile, 'w') as confile:
    config.write(confile)

# run fiducial
subprocess.call([camb_location, configfile])

# Load it for later
omch2 = config.getfloat('camb', 'omch2')

# ================================================
# get fiducial values to figure out where to compute next

# ================================================
#  FIDUCIALS
# ================================================

fid = {}
fid['hubble'] = config.getfloat('camb', 'hubble') / 100.
fid['scalar_spectral_index(1)'] = config.getfloat('camb', 'scalar_spectral_index(1)')
# fid['helium_fraction'] = config.getfloat('camb', 'helium_fraction')
fid['scalar_amp(1)'] = 10 ** 9 * config.getfloat('camb', 'scalar_amp(1)')
fid['massless_neutrinos'] = config.getfloat('camb', 'massless_neutrinos') + 1.
fid['re_optical_depth'] = config.getfloat('camb', 're_optical_depth')
fid['w'] = config.getfloat('camb', 'w')  # DE W parameters
# fid['wa'] = config.getfloat('camb', 'wa')  # DE W parameters
fid['ombh2'] = config.getfloat('camb', 'ombh2')
fid['omch2'] = config.getfloat('camb', 'omch2')
# fid['omk'] = config.getfloat('camb', 'omk')
# fid['scalar_nrun(1)'] = config.getfloat('camb', 'scalar_nrun(1)')
fid['omnuh2'] = config.getfloat('camb', 'omnuh2')
fid = collections.OrderedDict(sorted(fid.items(), key=lambda t: t[0]))


print fid

print "./data/{}/{}/fid_values.p".format(output_folder, output_folder_2)
with open("./data/{}/{}/fid_values.p".format(output_folder, output_folder_2), "wb") as output_file:
    pickle.dump(fid, output_file)

# ================================================
#  Gaps dx of x + dx etc, usually percentage of parameters
# ================================================

pargaps_dict = {}
pargaps_dict['hubble'] = fid['hubble'] * 0.05
# pargaps_dict['helium_fraction'] = fid['helium_fraction'] * 0.02
# for the value with fiducial = 0 we can not take a percentage. So fraction of 1509.07471 Table III
# pargaps_dict['omk'] = 0.01
# pargaps_dict['scalar_nrun(1)'] = 5e-2  # for the value with fiducial = 0 we can not take a percentage
pargaps_dict['scalar_spectral_index(1)'] = 0.0075
pargaps_dict['scalar_amp(1)'] = 0.05
pargaps_dict['re_optical_depth'] = 0.008
pargaps_dict['omnuh2'] = 0.00015
pargaps_dict['ombh2'] = 0.0012
pargaps_dict['omch2'] = 0.003
pargaps_dict['massless_neutrinos'] = 0.16
pargaps_dict['w'] = 0.15
# for the value with fiducial = 0 we can not take a percentage. So fraction of 1509.07471 Table III
# pargaps_dict['wa'] = 0.3


pargaps_dict = collections.OrderedDict(sorted(pargaps_dict.items(), key=lambda t: t[0]))
# save datagaps
print pargaps_dict
with open("./data/{}/{}/par_gaps.p".format(output_folder, output_folder_2), "wb") as output_file:
    pickle.dump(pargaps_dict, output_file)


# ================================================
#  values x + dx  x-dx etc we will run CAMB on these
# ================================================
# generate values to compute Cls
# step = np.array([-8,-7,-6,-5])
# step = np.array([-4,-3,-2,-1])
step = np.array([-4, -3, -2, -1.5, -1., -0.75, -0.5, 0.5, 0.75, 1, 1.5, 2, 3, 4])
# step = np.array([5,6,7,8])
# step = np.array([-12,-10,10,12])
# step = np.array([-1.5,-0.75,-0.5,0.5,0.75,1.5])

# step = np.array([-8, -7, -6, -5, -4, -3, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.4, -0.3, -0.25, -
#                  0.2, -0.10, 0.10, 0.2, 0.25, 0.3, 0.4, 0.5,0.75 ,1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6, 7, 8])


values = {}
values['hubble'] = pargaps_dict['hubble'] * step + fid['hubble']
# values['helium_fraction'] = pargaps_dict['helium_fraction'] * step + fid['helium_fraction']
values['scalar_spectral_index(1)'] = pargaps_dict['scalar_spectral_index(1)'] * \
    step + fid['scalar_spectral_index(1)']
values['scalar_amp(1)'] = pargaps_dict['scalar_amp(1)'] * step + fid['scalar_amp(1)']
values['massless_neutrinos'] = pargaps_dict['massless_neutrinos'] * step + fid['massless_neutrinos']
values['re_optical_depth'] = pargaps_dict['re_optical_depth'] * step + fid['re_optical_depth']
values['omnuh2'] = pargaps_dict['omnuh2'] * step + fid['omnuh2']
values['w'] = pargaps_dict['w'] * step + fid['w']
# values['wa'] = pargaps_dict['wa'] * step + fid['wa']
values['ombh2'] = pargaps_dict['ombh2'] * step + fid['ombh2']
values['omch2'] = pargaps_dict['omch2'] * step + fid['omch2']
# values['omk'] = pargaps_dict['omk'] * step + fid['omk']
# values['scalar_nrun(1)'] = pargaps_dict['scalar_nrun(1)'] * step + fid['scalar_nrun(1)']


values = collections.OrderedDict(sorted(values.items(), key=lambda t: t[0]))
# save a pickle of data values to be re-used
with open("./data/{}/{}/grid_values.p".format(output_folder, output_folder_2), "wb") as output_file:
    pickle.dump(values, output_file)

print fid, values
# ================================================
# start loop on them and generate.

for key, value in values.iteritems():
    for i in np.arange(0, np.size(step)):
        print ''

        print 'modifying', key, 'values=', values[key][i]
        print ''
        config.read('./fiducial.ini')
        omch2 = config.getfloat('camb', 'omch2')

        # SPECIAL CONDITIONS FOR SOME VALUES FLATNESS IS ALWAYS ENFORCED. BUT THAT
        # IS IT EVERYTHING ELSE NEED TO BE INSERTED BY HAND

        # if key == 'scalar_spectral_index(1)':
        #     continue

        # if key == 're_optical_depth':
        #     continue

        # if key != 'massless_neutrinos':
        #     continue

        if (key == 'hubble'):
            config.set('camb', key, str(100. * values[key][i]))

        elif (key == 'massless_neutrinos'):
            # from eq (14) of http://arxiv.org/pdf/1402.4108v1.pdf

            config.set('camb', key, str(values[key][i] - 1.))
            S = np.sqrt(1. + 7. * (values[key][i] - 3.046) / 43.)
            eta10 = 273.9 * fid['ombh2']
            Y_p = 0.2486 + 0.0016 * ((eta10 - 6.) + 100. * (S - 1.))
            config.set('camb', 'helium_fraction', str(Y_p))

        elif (key == 'scalar_amp(1)'):
            config.set('camb', key, str(values[key][i] * 10 ** (-9)))

        else:

            config.set('camb', key, str(values[key][i]))

        # print type(config.getfloat('camb', key))
        config.set('camb', 'output_root', './data/{}/{}/'.format(output_folder, output_folder_2) +
                   key + '_{:.13f}'.format(values[key][i]))
        configfile_temp = './data/{}/{}/'.format(output_folder, output_folder_2) + \
            key + '_{:.13f}'.format(values[key][i]) + '.ini'
        # if os.path.isfile(path): sys.exit('this inifile already exist, delete if you want to ovwrwrite')
        if os.path.isfile(configfile_temp):
            print 'folder', configfile_temp
            print('THE DATA ALREADY EXIST IN {} I WILL SKIP IT. IF YOU WANT TO REGENERATE IT, DELETE.')
            continue

        print ''
        print config.getfloat('camb', 'hubble')
        print config.getfloat('camb', 'scalar_amp(1)')
        print config.getfloat('camb', 'scalar_spectral_index(1)')
        print config.getfloat('camb', 'massless_neutrinos')
        print config.getfloat('camb', 're_optical_depth')
        print config.getfloat('camb', 'omnuh2') * 93.359

        with open(configfile_temp, 'w') as confile:
            config.write(confile)

        subprocess.call([camb_location, configfile_temp])


sys.exit()
