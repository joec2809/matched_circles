from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle_r import match_circle_r
from match_circle_s import match_circle_s
from match_circle_s_fast import match_circle_s_fast
from rotate import rotate_to_top
import time

start_time = time.time()

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

lon = 207.8
lat = -56.3

#cmb_map = np.multiply(cmb_map_og,mask_ring)
cmb_map = rotate_to_top(np.multiply(cmb_map_og,mask_ring), lon, lat)

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/360)
corr = np.zeros(len(ang_rad), dtype=complex)

for i in range(len(ang_rad)):
	rad_lag = 0
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	corr[i] = match_circle_s_fast(circle_a, circle_b, rad_lag)

print time.time()-start_time
	

ang_rad = ang_rad*(360/(2*np.pi))

fig, ax = plt.subplots()

ax.plot(ang_rad, corr)
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.set_title('Fast S: CS')
ax.annotate('Lag=0$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_cs_fast.png', overwrite = True)

plt.show()
