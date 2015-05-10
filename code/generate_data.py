'''
Script to produced the different Cls with different parameters needed to compute derivative for the Fisher matrices.

It chenge the ini file of CAMB at each step with configparser
and it then run camb.

f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
      ---------------------------------
                    12h
For better accuracy the derivative  should be computed with a 5 point stencil algorithms

CONVENTIONS:

PARAMETER ORDER = h0, ns, As, Neff


'''


import configparser
import subprocess
import numpy as np
import sys
import pickle

# deltax to be used with different parameters

# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = np.array([1., 0.05, 5e-11, 0.08])  # h0, ns, As, Neff


config = configparser.ConfigParser()
configfile = './fiducial.ini'
config.read(configfile)

# run fiducial
subprocess.call(['/Users/alessandromanzotti/Work/Software/camb2013/camb',configfile])

# get fiducial values to figure out where to compute next

h = config.getfloat('camb', 'hubble')
As = config.getfloat('camb', 'scalar_amp(1)')
ns = config.getfloat('camb', 'scalar_spectral_index(1)')
N_eff = config.getfloat('camb', 'massless_neutrinos')  # is it true? what do we want to keep fixed?

np.savetxt("./data/run2/fiducial_pars.txt", np.array([h, ns, As, N_eff]))

# generate values to compute Cls
values = {}
values['hubble'] = pargaps[0] * np.array([-2, -1, 1, 2]) + h
values['scalar_spectral_index(1)'] = pargaps[1] * np.array([-2, -1, 1, 2]) + ns
values['scalar_amp(1)'] = pargaps[2] * np.array([-2, -1, 1, 2]) + As
values['massless_neutrinos'] = pargaps[3] * np.array([-2, -1, 1, 2]) + N_eff

pargaps_dict = {}
pargaps_dict['hubble'] = pargaps[0]
pargaps_dict['scalar_spectral_index(1)'] = pargaps[1]
pargaps_dict['scalar_amp(1)'] = pargaps[2]
pargaps_dict['massless_neutrinos'] = pargaps[3]
# save a pickle of data values to be re-used
with open("./data/run2/grid_values.p", "wb") as output_file:
    pickle.dump(values, output_file)

# save datagaps

with open("./data/run2/par_gaps.p", "wb") as output_file:
    pickle.dump(pargaps_dict, output_file)



# start loop on them and generate.

for key, value in values.iteritems():
    for i in np.arange(0, 4):
        print key, values[key][i]
        print ''
        config.read('./fiducial.ini')
        # print config.getfloat('camb', 'massless_neutrinos')
        # print key=='massless_neutrinos',
        config.set('camb', key, str(values[key][i]))
        # print type(config.getfloat('camb', key))
        config.set('camb', 'output_root', './data/run3/' + key + '_{:.13f}'.format(values[key][i]))
        configfile_temp = key + '_{:.13f}'.format(values[key][i]) + '.ini'
        print config.getfloat('camb', 'hubble')
        print config.getfloat('camb', 'scalar_amp(1)')
        print config.getfloat('camb', 'scalar_spectral_index(1)')
        print config.getfloat('camb', 'massless_neutrinos')

        with open(configfile_temp, 'w') as confile:
            config.write(confile)

        subprocess.call(['/Users/alessandromanzotti/Work/Software/camb2013/camb', configfile_temp])


sys.exit()
