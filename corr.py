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
from rotate import rotate_to_top
from T_m_setup import T_m_setup
import time

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
m_max = 720

T_m = T_m_setup(m_max, bins)

cmb_map = cmb_map_og*mask_ring
#cmb_map = rotate_to_top(cmb_map_og*mask_ring, lon, lat)

# Predict Circles
"""ang_rad = np.zeros(100)
x_corr = np.zeros(100, dtype = complex)
a = -1
rad_lag = 0


for i in range(len(ang_rad)):
	xi = (i+1)*1
	yi = (i+1)*0
	if np.sqrt(xi**2+yi**2) <= (2*CMB_DIST/CELL_SIZE):
		strip_finder(cmb_map, circle_finder(CELL_SIZE, xi, yi), NSIDE)

		circle_a = load_file('strip_a', bins)
		circle_b = load_file('strip_b', bins)
	
		ang_rad[i] = (360/(2*np.pi))*circle_finder(CELL_SIZE, xi, yi)
		x_corr[i] = match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
		a += 1
	else:
		break

ang_rad = ang_rad[:a]
x_corr = x_corr[:a]

np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_ang_rad_pred.csv', ang_rad, delimiter = ',')
np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_corr_pred.csv', x_corr, delimiter = ',')"""

# vs lag
"""strip_finder(cmb_map, 87.9*(2*np.pi/360), NSIDE)

circle_a = load_file('strip_a', bins)
circle_b = load_file('strip_b', bins)

rad_lag = 0
	
x_corr = np.zeros((bins), dtype = complex)
for i in range(bins):
	x_corr[i] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	rad_lag += 2*np.pi/360

np.savetxt('/opt/local/l4astro/rbbg94/data/87.9_s_corr.csv', x_corr, delimiter = ',')"""

# Errors
"""ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros((no_circles,bins), dtype=complex)

for i in range(no_circles):
	rad_lag = 0
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = load_file("strip_a", bins)
	circle_b = load_file("strip_b", bins)
	for j in range(bins):
		x_corr[i,j] = match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
		rad_lag += 2*np.pi/360				
	print i, time.time()-start
	
	
ang_rad = ang_rad*(360/(2*np.pi))

for i in range(bins):
	err_corr = np.delete(x_corr, i, axis = 1)
	std_dev = np.std(err_corr, axis=1)
	np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_lag_err_'+str(i)+'.csv', std_dev, delimiter = ',')"""

# Use errors
"""ang_rad = linspace((1/360)*2*np.pi, np.pi/2, 360)
x_corr = np.zeros(len(ang_rad), dtype=complex)
#std_dev = np.genfromtxt('/opt/local/l4astro/rbbg94/data/cs_sim_err_100.csv', delimiter = ',')
std_dev = 0

for i in range(len(ang_rad)):
	rad_lag = 0
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = load_file('strip_a', bins)
	circle_b = load_file('strip_b', bins)
	x_corr[i] = match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
	
np.savetxt('/opt/local/l4astro/rbbg94/data/cs_90_corr.csv', x_corr, delimiter = ',')"""

# Corr
ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 89)
x_corr = np.zeros(len(ang_rad), dtype=complex)
T_m = T_m_setup(m_max, bins)


start = time.time()
for i in range(len(ang_rad)):
	rad_lag = 0
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	circle_a = load_file('strip_a', bins)
	circle_b = load_file('strip_b', bins)
	x_corr[i] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	print(i, time.time()-start)

np.savetxt('/opt/local/l4astro/rbbg94/data/ngp_corr_0_thin.csv', x_corr, delimiter = ',')
