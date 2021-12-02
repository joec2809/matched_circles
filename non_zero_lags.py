from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import mc_functions
import time

start = time.time()

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

bins = 360

cmb_map = cmb_map_og*mask_ring

ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros((no_circles,bins), dtype=complex)

for i in range(no_circles):
	rad_lag = 0
	mc_functions.strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = mc_functions.load_file("strip_a", bins)
	circle_b = mc_functions.load_file("strip_b", bins)
	for j in range(bins):
		x_corr[i,j] = mc_functions.match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
		rad_lag += 2*np.pi/360				
	print(i, time.time()-start)
	
	
ang_rad = ang_rad*(360/(2*np.pi))

x_corr = np.swapaxes(x_corr,0,1)
for i in range(bins):
	np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_corr_'+str(i)+'.csv', x_corr[i], delimiter = ',')
