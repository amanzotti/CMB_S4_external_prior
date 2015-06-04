import matplotlib.pyplot as plt
import numpy as np

dat = np.genfromtxt('../output/param_cov_LSS.txt')

plt.imshow(dat,cmap=plt.get_cmap('bwr'),interpolation='nearest',vmin=-1,vmax=1)
plt.xticks([0,1,2,3,4,5],[r'$H_0$',r'$N_\mathrm{eff}$',r'$\Omega_\nu h^2 $',r'$\tau $',r'$A_s $',r'$n_s $'])
plt.yticks([0,1,2,3,4,5],[r'$H_0$',r'$N_\mathrm{eff}$',r'$\Omega_\nu h^2 $',r'$\tau $',r'$A_s $',r'$n_s $'])
plt.gca().xaxis.set_ticks_position('top')

plt.colorbar()

plt.savefig('../../images/parcov.pdf')
plt.show()
