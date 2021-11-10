from __future__ import division
import numpy as np
import time
from strip import strip_finder
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

cmb_map = cmb_map_og*mask_ring

CMB_DIST = 14000
CELL_SIZE = 320

bins = 360
m_max = 720

ang_rad = np.linspace((2*np.pi)/360, np.pi/2, 90)[:-1]
no_pix = np.zeros(len(ang_rad))

for i in range(len(ang_rad)):
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	hdul = fits.open('/opt/local/l4astro/rbbg94/cmb_maps/strip_a')
	circ_dat = hdul[1].data
	hdul.close()
	T = circ_dat['T']
	lon = circ_dat['long']
	lat = circ_dat['lat']
	index = circ_dat['index']
		
	circle = np.zeros((len(T),4))
	circle[:,3] = index
	circle[:,0] = T
	circle[:,1] = lon
	circle[:,2] = lat
	
	sort_ind = circle[:,1].argsort()

	circle = circle[sort_ind]

	array = np.linspace(0,bins-1,bins)
	lon_bin = np.array(np.cumsum(pd.cut(circle[:,1], bins = array, right = False, include_lowest = True).value_counts()))
	T_bin = np.array(np.split(circle[:,0], lon_bin))
	no_pix_circle = np.zeros(len(T_bin))
	for j in range(len(T_bin)):
		no_pix_circle[j] = len(T_bin[j])
	no_pix[i] = np.mean(no_pix_circle)

print no_pix[0]
print no_pix[-1]

ang_rad = ang_rad*(360/(2*np.pi))

fig, ax = plt.subplots()
ax.plot(ang_rad, no_pix, color = 'k')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('Average number of pixels per bin')
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/ave_pix_bin_ngp.png', overwrite = True)

plt.show()
