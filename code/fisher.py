import numpy as np
import math
import matplotlib.pyplot as plt

def C(iell,ell,parbin):

	s = 350*math.sqrt(20626*60*60)/math.sqrt(10**6*0.25*5*365*24*60*60)
	t = 2/60/180*math.pi

	fac = ell*(ell+1)/2/math.pi
	fac2 = ell**4/ell/(ell+1)

	N = s**2 * math.exp( ell * (ell+1) * t**2 / 8 / math.log(2) )

	C = np.array([[dats[iell,1,parbin]/fac+N,dats[iell,2,parbin]/fac,0],\
		[dats[iell,2,parbin]/fac,dats[iell,3,parbin]/fac+N*2**0.5,0],\
		[0,0,dats[iell,4,parbin]]])
	return C


dats = np.genfromtxt('data/dat0.txt')

for i in range(1,9):
	newdat = np.genfromtxt('data/dat{}.txt'.format(i))
	dats = np.dstack((dats,newdat))

print dats[:,:,0]

fisher = np.zeros((4,4))
pargaps = np.array([0.002,0.02,2e-11,0.02]) #h0, ns, As, Neff

for iell, ell in enumerate(range(2,5000)):

	print ell

	c0 = C(iell,ell,0)
	cinv = np.linalg.inv(c0)

	for i in range(0,4):
		for j in range(0,4):
			ci = (C(iell,ell,i*2+1)-C(iell,ell,i*2+2))/pargaps[i]
			cj = (C(iell,ell,j*2+1)-C(iell,ell,j*2+2))/pargaps[j]

			tot = np.dot(np.dot(np.dot(cinv,ci),cinv),cj)

			fisher[i,j] += (2*ell+1)/2*0.5*np.trace(tot)

d = []
d2 = []
d3= []

for i in np.arange(-3,-1,0.1):
	fisher1 = fisher.copy()
	fisher1[0,0] += 1/(10**i*0.673)**2
	d.append(math.sqrt(np.linalg.inv(fisher1)[3,3]))

	fisher2 = fisher.copy()
	fisher2[0,0] += 1/(10**i*0.673)**2
	fisher2[1,1] += 1/(0.01*0.96)**2
	fisher2[2,2] += 1/(0.01*2.2e-9)**2

	d2.append(math.sqrt(np.linalg.inv(fisher2)[3,3]))

	fisher3 = fisher.copy()[[0,3],:][:,[0,3]]
	fisher3[0,0] += 1/(10**i*0.673)**2

	d3.append(math.sqrt(np.linalg.inv(fisher3)[1,1]))

plt.plot(10**np.arange(-3,-1,0.1),d,label='No Priors')
plt.plot(10**np.arange(-3,-1,0.1),d2,label='1% Priors')
plt.plot(10**np.arange(-3,-1,0.1),d3,label='Perfect Priors')
plt.xscale('log')
plt.xlabel(r'$\Delta H_0 / H_0$',fontsize=16)
plt.ylabel(r'$\sigma(N_\mathrm{eff})$',fontsize=16)
plt.legend(loc='middle right')

plt.savefig('h0.pdf')
plt.show()