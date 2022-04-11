from __future__ import division
import numpy as np
import healpy as hp
from astropy.io import fits
from astrotools import healpytools as hpt
from tqdm import tqdm
import mc_functions
import time

cmb_map_og = hp.fitsfunc.read_map("../cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("../cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320
bins = 360
m_max = 720

og_lon = 238.2
og_lat = -28.8
og_vec = hp.pixelfunc.ang2vec(og_lon, og_lat, lonlat=True)
vec_lon = 238.2
vec_lat = og_lat + 90
vec = hp.pixelfunc.ang2vec(vec_lon, vec_lat, lonlat=True)

no_dir = 1
rot_ang = np.pi/no_dir

ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros((no_dir,no_circles), dtype = complex)

T_m = mc_functions.T_m_setup(m_max, bins)

cmb_map = mc_functions.rotate_to_top(cmb_map_og*mask_ring, og_lon, og_lat)

for i in tqdm(range(no_dir)):
	for j in range(len(ang_rad)):
		rad_lag = 0
		mc_functions.strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = mc_functions.load_file('strip_a', bins)
		circle_b = mc_functions.load_file('strip_b', bins)
		x_corr[i,j] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	np.savetxt(f'../data/dipole_corr.csv', x_corr[i], delimiter = ',')
	cmb_map = hpt.rotate_map(cmb_map, og_vec, rot_ang)