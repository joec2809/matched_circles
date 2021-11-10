from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle_s import match_circle_s
from T_m_setup import T_m_setup

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

bins = 360
m_max = 720

T_m = T_m_setup(m_max, bins)

cmb_map = cmb_map_og

x_corr_1 = np.zeros((bins), dtype = complex)

rad_lag = 0

strip_finder(cmb_map, 80*(2*np.pi/360), NSIDE)

circle_a = load_file('strip_a', bins)
circle_b = load_file('strip_b', bins)

for i in range(bins):
	x_corr_1[i] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	rad_lag += 2*np.pi/360

rand_ang = np.pi/np.random.randint(1,20)
rand_vec = np.array((np.random.random(), np.random.random(), np.random.random()))

cmb_map = hpt.rotate_map(cmb_map, rand_vec, rand_ang)
cmb_map = hpt.rotate_map(cmb_map, rand_vec, -rand_ang)

x_corr_2 = np.zeros((bins), dtype = complex)

rad_lag = 0

strip_finder(cmb_map, 80*(2*np.pi/360), NSIDE)

circle_a = load_file('strip_a', bins)
circle_b = load_file('strip_b', bins)

for i in range(bins):
	x_corr_2[i] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	rad_lag += 2*np.pi/360

diff = (x_corr_1-x_corr_2)

phase = np.arange(0,360,1)

fig, ax = plt.subplots()
ax.plot(phase, diff, color = 'k')
plt.xlim(0,360)
plt.axhline(y=0, color = 'k')
ax.set_xlabel('Phase/$^\circ$')
ax.set_ylabel('$\Delta{}S$')
plt.tight_layout()

#fig.savefig("/opt/local/l4astro/rbbg94/figures/rotate_diff.png", overwrite = True)

plt.show()
