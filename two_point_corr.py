from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from scipy.special import legendre
import time
start = time.time()

ang_sep = np.arange(2*np.pi/360, (181/360)*2*np.pi, 2*np.pi/360)
no_seps = len(ang_sep)
l_max = 1000
no_sims = 5

c_thetas = np.zeros((no_sims, no_seps))

for j in range(no_sims):
	cmb_map_og = hp.fitsfunc.read_map("../sims/dx12_v3_smica_cmb_mc_000"+str(j).zfill(2)+"_raw.fits")
	NSIDE = hp.npix2nside(len(cmb_map_og))

	alms = hp.sphtfunc.map2alm(cmb_map_og, lmax = l_max)
	alms_sq = (np.abs(alms))**2

	alm_sep = np.zeros((l_max+1,l_max+1))
	leg_poly = np.zeros((no_seps,l_max+1))

	for i in range(no_seps):
		for l in range(l_max+1):
			leg_poly[i,l] = legendre(l)(np.cos(ang_sep[i]))
			for m in range(l+1):
				alm_sep[l,m] = alms_sq[hp.sphtfunc.Alm.getidx(l_max,l,m)]
		print("i="+str(i), time.time() - start)

	alm_sep = alm_sep*2
	alm_sep[:,0]/2
	alm = np.sum(alm_sep, axis = 1)
	leg_alm = leg_poly*alm
	c_theta = ((1/4)*np.pi)*np.sum(leg_alm, axis = 1)
	c_theta = c_theta*(10**12)
	c_thetas[j] = c_theta
	print("j="+str(j), time.time()-start)

np.savetxt('../data/two_point_corr.csv', c_thetas, delimiter = ',')

ave_c_theta = np.mean(c_thetas, axis = 0)

ang_sep = ang_sep*(360/(2*np.pi))

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

fig, ax = plt.subplots()
ax.plot(ang_sep, ave_c_theta, color = 'k')
ax.set_xlabel(r'$\theta/^\circ$')
ax.set_ylabel(r'$C(\theta)/\mu$K$^2$')
ax.axhline(0, color = 'black')
plt.xlim(0,180)
plt.xticks(np.arange(0, 181, 30))
plt.tight_layout()
fig.savefig('../figures/c_theta.png', overwrite = True)
plt.show()
