from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import time
import mc_functions

cmb_map_og = hp.fitsfunc.read_map("../cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("../cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

lon = 207.8
lat = -56.3

bins = 360
m_max = 720

T_m = mc_functions.T_m_setup(m_max, bins)

cmb_map = cmb_map_og*mask_ring
#cmb_map = rotate_to_top(cmb_map_og*mask_ring, lon, lat)


# Corr
ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 89)
x_corr = np.zeros(len(ang_rad), dtype=complex)


start = time.time()
for i in range(len(ang_rad)):
	rad_lag = 0
	mc_functions.strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = mc_functions.load_file('strip_a', bins)
	circle_b = mc_functions.load_file('strip_b', bins)
	x_corr[i] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	print(i, time.time()-start)

np.savetxt('../data/ngp_corr_0_thin.csv', x_corr, delimiter = ',')
