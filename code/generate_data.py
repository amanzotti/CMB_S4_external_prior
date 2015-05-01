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

# deltax to be used with different parameters

# gaps beween  x1 x_-1 these three are used to get the value of the derivative in the middle
pargaps = np.array([0.002, 0.02, 2e-11, 0.02])  # h0, ns, As, Neff



config = configparser.ConfigParser()
configfile='./fiducial.ini'

# run fiducial
subprocess.call(['/Users/alessandromanzotti/Work/Software/camb2013/camb',configfile])




# config.read(configfile)
# # config.set('camb','massless_neutrinos',2.046)
# config.set('camb','output_root','N_eff_{:f}'.format(2.046))

# with open(configfile, 'w') as confile:
#     config.write(confile)


# subprocess.call(['/Users/alessandromanzotti/Work/Software/camb2013/camb',configfile])