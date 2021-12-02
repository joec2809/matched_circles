from __future__ import division
import numpy as np
import healpy as hp
from astropy.io import fits
from astrotools import healpytools as hpt
import mc_functions
import time

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320
bins = 360

cs_lon = 207.8
cs_lat = -56.3
cs_vec = hp.pixelfunc.ang2vec(cs_lon, cs_lat, lonlat=True)
vec_lon = 207.8
vec_lat = cs_lat + 90
vec = hp.pixelfunc.ang2vec(vec_lon, vec_lat, lonlat=True)

no_dir = 12
rot_ang = np.pi/no_dir

ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros((no_dir,no_circles), dtype = complex)

start = time.time()

cmb_map = mc_functions.rotate_to_top(cmb_map_og*mask_ring, vec_lon, vec_lat)

for i in range(no_dir):
	for j in range(len(ang_rad)):
		rad_lag = 0
		mc_functions.strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = mc_functions.load_file('strip_a', bins)
		circle_b = mc_functions.load_file('strip_b', bins)
		x_corr[i,j] = mc_functions.match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
	np.savetxt('/opt/local/l4astro/rbbg94/data/cs_dir_'+str(i)+'.csv', x_corr[i], delimiter = ',')
	print(time.time()-start)
	cmb_map = hpt.rotate_map(cmb_map, cs_vec, rot_ang)
