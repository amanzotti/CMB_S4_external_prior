'''
Script to produced the different Cls with different parameters needed to compute derivative for the Fisher matrices.

It changes the ini file of CAMB at each step with configparser and it then runs camb.

For better accuracy the derivative are computed with a 5 point stencil algorithms

f' = -f(x+2h) + 8f(x+h) -8f(x-h)+f(x-2h)
      ---------------------------------
                    12h

CONVENTIONS:

PARAMETER ORDER = we are using ordered dictionaries. So the order is the alphabetical order of the name of the variable sin the camb ini.


'''


import configparser
import subprocess
import numpy as np
import sys
import pickle
import collections


output_folder = 'run4'

# load fiducial camb ini file
config = configparser.ConfigParser()
configfile = './fiducial.ini'
config.read(configfile)

# run fiducial
subprocess.call(['/home/manzotti/local/camb2013/camb', configfile])


# Load it for later
omch2 = config.getfloat('camb', 'omch2')

# ================================================
# get fiducial values to figure out where to compute next


# ================================================
#  FIDUCIALS
# ================================================

fid = {}
fid['hubble'] =  config.getfloat('camb', 'hubble')
fid['scalar_spectral_index(1)'] =  config.getfloat('camb', 'scalar_spectral_index(1)')
fid['scalar_amp(1)'] = config.getfloat('camb', 'scalar_amp(1)')
fid['massless_neutrinos'] =  config.getfloat('camb', 'massless_neutrinos')
fid['re_optical_depth'] = config.getfloat('camb', 're_optical_depth')
fid['w'] = config.getfloat('camb', 'w') #DE W parameters
fid['ombh2'] = config.getfloat('camb', 'ombh2')
fid['omch2'] = config.getfloat('camb', 'omch2')
fid['omnuh2'] = config.getfloat('camb', 'omnuh2')
fid = collections.OrderedDict(sorted(fid.items(), key=lambda t: t[0]))


print "./data/{}/fid_values.p".format(output_folder)
with open("./data/{}/fid_values.p".format(output_folder), "wb") as output_file:
    pickle.dump(fid, output_file)

# ================================================
#  Gaps dx of x + dx etc, usually percentage of parameters
# ================================================

pargaps_dict = {}
pargaps_dict['hubble'] = fid['hubble']*0.06
pargaps_dict['scalar_spectral_index(1)'] = fid['scalar_spectral_index(1)']*0.06
pargaps_dict['scalar_amp(1)'] = fid['scalar_amp(1)']*0.06
pargaps_dict['massless_neutrinos'] = fid['massless_neutrinos']*0.1
pargaps_dict['re_optical_depth'] = fid['re_optical_depth']*0.06
pargaps_dict['omnuh2'] = fid['omnuh2']*0.3
pargaps_dict['w'] = fid['w']*0.06
pargaps_dict['ombh2'] = fid['ombh2']*0.06
pargaps_dict['omch2'] = fid['omch2']*0.06


pargaps_dict = collections.OrderedDict(sorted(pargaps_dict.items(), key=lambda t: t[0]))
# save datagaps

with open("./data/{}/par_gaps.p".format(output_folder), "wb") as output_file:
    pickle.dump(pargaps_dict, output_file)



# ================================================
#  values x + dx  x-dx etc we will run CAMB on these
# ================================================
# generate values to compute Cls
values = {}
values['hubble'] = pargaps_dict['hubble'] * np.array([-2, -1, 1, 2]) + fid['hubble']
values['scalar_spectral_index(1)'] = pargaps_dict['scalar_spectral_index(1)'] * np.array([-2, -1, 1, 2]) + fid['scalar_spectral_index(1)']
values['scalar_amp(1)'] = pargaps_dict['scalar_amp(1)'] * np.array([-2, -1, 1, 2]) + fid['scalar_amp(1)']
values['massless_neutrinos'] = pargaps_dict['massless_neutrinos'] * np.array([-2, -1, 1, 2]) + fid['massless_neutrinos']
values['re_optical_depth'] = pargaps_dict['re_optical_depth'] * np.array([-2, -1, 1, 2]) + fid['re_optical_depth']
values['omnuh2'] = pargaps_dict['omnuh2'] * np.array([-2, -1, 1, 2]) + fid['omnuh2']
values['w'] = pargaps_dict['w'] * np.array([-2, -1, 1, 2]) + fid['w']
values['ombh2'] = pargaps_dict['ombh2'] * np.array([-2, -1, 1, 2]) + fid['ombh2']
values['omch2'] = pargaps_dict['omch2'] * np.array([-2, -1, 1, 2]) + fid['omch2']

values = collections.OrderedDict(sorted(values.items(), key=lambda t: t[0]))
# save a pickle of data values to be re-used
with open("./data/{}/grid_values.p".format(output_folder), "wb") as output_file:
    pickle.dump(values, output_file)


# ================================================


# start loop on them and generate.

for key, value in values.iteritems():
    for i in np.arange(0, 4):
        print ''

        print 'modifying', key, 'values=', values[key][i]
        print ''
        config.read('./fiducial.ini')
        omch2 = config.getfloat('camb', 'omch2')

        # SPECIAL CONDITIONS FOR SOME VALUES FLATNESS IS ALWAYS ENFORCED. BUT THAT IS IT EVERYTHING ELSE NEED TO BE INSERTED BY HAND


        if key=='hubble':
            # CHANGE hubble-> chenge all the h^2 quantity to keep lambda fixed
            # set omega
            config.set('camb', key, str( (fid['hubble']/values[key][i])) )
            # reduce omega_m
            print 'test_hubble',values[key][i]/fid['hubble'], fid['omch2']/fid['hubble']**2, fid['omch2'] * (fid['hubble']/values[key][i])**2
            config.set('camb', 'omch2', str(fid['omch2'] * (values[key][i]/fid['hubble'])**2 ) )
            config.set('camb', 'ombh2', str(fid['ombh2'] * (values[key][i]/fid['hubble'])**2) )
            config.set('camb', 'omnuh2', str(fid['omnuh2']* (values[key][i]/fid['hubble'])**2) )



        if key=='omnuh2':
            # CHANGE OMEGA NU but keeping lambda fixed
            # set omega
            config.set('camb', key, str(values[key][i]))
            # reduce omega_m
            print 'test_omeganu',omch2, (values[key][i]-fid[key]), fid['omnuh2'], values[key][i]

            config.set('camb', 'omch2', str(omch2-(values[key][i]-fid[key])) )

        if key=='ombh2':
            # set omega
            config.set('camb', key, str(values[key][i]))
            # reduce omega_m
            print 'test_ombh2 ',omch2, (values[key][i]-fid[key]), fid['ombh2'], values[key][i]
            # reduce omega_c to keep lambda constant
            config.set('camb', 'omch2', str(omch2-(values[key][i]-fid[key])) )


        # print config.getfloat('camb', 'massless_neutrinos')
        # print key=='massless_neutrinos',
        config.set('camb', key, str(values[key][i]))

        # print type(config.getfloat('camb', key))
        config.set('camb', 'output_root', './data/{}/'.format(output_folder) + key + '_{:.13f}'.format(values[key][i]))
        configfile_temp = './data/{}/'.format(output_folder) + key + '_{:.13f}'.format(values[key][i]) + '.ini'
        print ''
        print config.getfloat('camb', 'hubble')
        print config.getfloat('camb', 'scalar_amp(1)')
        print config.getfloat('camb', 'scalar_spectral_index(1)')
        print config.getfloat('camb', 'massless_neutrinos')
        print config.getfloat('camb', 're_optical_depth')
        print config.getfloat('camb', 'omnuh2')

        with open(configfile_temp, 'w') as confile:
            config.write(confile)

        subprocess.call(['/home/manzotti/local/camb2013/camb', configfile_temp])


sys.exit()
