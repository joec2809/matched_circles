from __future__ import division
import numpy as np
import healpy as hp
import pandas as pd
from astropy.io import fits
import mc_functions
import time

start = time.time()

cmb_map_og = hp.fitsfunc.read_map("../cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("../cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

lon = 144.0
lat = 85.0

mode = 'p3'

bins = 360
m_max = 720

pred = pd.read_csv(f"../data/pred_angs/omega_{mode}_reflection.txt", delim_whitespace=True, header = 1, index_col=0)

ang_rads = np.array(pred['alpha'])[~np.isnan(pred['alpha'])]

x_corrs = np.genfromtxt(f'../data/grid_search/mode_{mode}/{mode}_lon_{lon}_lat_{lat}.csv', dtype=complex, delimiter = ',')


ang_rad = ang_rads[np.argmax(np.real(x_corrs))]
x_corr = np.zeros(bins, dtype=complex)

cmb_map = mc_functions.rotate_to_top(cmb_map_og*mask_ring, lon, lat)

T_m = mc_functions.T_m_setup(m_max, bins)

rad_lag = 0
mc_functions.strip_finder(cmb_map, ang_rad, NSIDE)
circle_a = mc_functions.load_file("strip_a", bins)
circle_b = mc_functions.load_file("strip_b", bins)
for j in range(bins):
	x_corr[j] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, np.radians(rad_lag))
	rad_lag += 360/bins				

np.savetxt(f'../data/p3_lon_{lon}_lat_{lat}_lags.csv', x_corr, delimiter = ',')