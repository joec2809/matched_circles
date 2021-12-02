from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import mc_functions

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000

lon = 297.8
lat = 33.7

bins = 360

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/360)
x_corr = np.zeros(len(ang_rad), dtype=complex)
	
cmb_map = mc_functions.rotate_to_top(np.multiply(cmb_map_og,mask_ring), lon, lat)

for i in range(len(ang_rad)):
	mc_functions.strip_finder(cmb_map, ang_rad[i], NSIDE)

	circle_a = mc_functions.load_file('strip_a', bins)
	circle_b = mc_functions.load_file('strip_b', bins)
	x_corr[i] = mc_functions.match_circle_s(circle_a, circle_b, 0, bins, 720)

np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_corr_lon_297_lat_33.csv', x_corr, delimiter = ',')

ang_rad_deg = ang_rad*(360/(2*np.pi))	

fig, ax = plt.subplots()

ax.plot(ang_rad_deg, x_corr, color = 'k')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel('$S$')
ax.annotate('Phase=0$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
plt.tight_layout()
