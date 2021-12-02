from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
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

no_sims = 100
bins = 360
m_max = 720

lon = 0
lat = 90

ang_rad = np.linspace((2*np.pi)/360, np.pi/2, 90)[:-1]
no_circles = len(ang_rad)
T_m = mc_functions.T_m_setup(m_max, bins)

start = time.time()

for i in range(1):
	"""x_corr_0 = np.zeros(no_circles, dtype=complex)
	cmb_map = rotate_to_top(cmb_map_og*mask_ring, lon, lat)
	for j in range(no_circles):
		rad_lag = 0
		strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = load_file('strip_a', bins)
		circle_b = load_file('strip_b', bins)
		x_corr_0[j] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	np.savetxt('/opt/local/l4astro/rbbg94/data/plane_corrs_and_errs/lon_'+str(lon)+'_lat_'+str(lat)+'_corr_0.csv', x_corr_0, delimiter = ',')
	lat += 45
	x_corr_45 = np.zeros(no_circles, dtype=complex)
	cmb_map = rotate_to_top(cmb_map_og*mask_ring, lon, lat)
	for j in range(no_circles):
		rad_lag = 0
		strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = load_file('strip_a', bins)
		circle_b = load_file('strip_b', bins)
		x_corr_45[j] = match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	np.savetxt('/opt/local/l4astro/rbbg94/data/plane_corrs_and_errs/lon_'+str(lon)+'_lat_'+str(lat)+'_corr_0.csv', x_corr_45, delimiter = ',')
	#lat -= 45"""

	x_corr_sims_0 = np.zeros((no_sims,no_circles), dtype=complex)
	x_corr_sims_180 = np.zeros((no_sims,no_circles), dtype=complex)
	for j in range(no_sims):
		cmb_map_sim = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/sims/dx12_v3_smica_cmb_mc_000"+str(j).zfill(2)+"_raw.fits")
		NSIDE = hp.npix2nside(len(cmb_map_sim))
		
		
		for k in range(no_circles):
			rad_lag = 0
			mc_functions.strip_finder(cmb_map_sim, ang_rad[k], NSIDE)
			circle_a = mc_functions.load_file('strip_a', bins)
			circle_b = mc_functions.load_file('strip_b', bins)
			x_corr_sims_0[j,k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
		
		
		cmb_map = mc_functions.rotate_to_top(cmb_map_sim, 0, 0)
		for k in range(no_circles):
			rad_lag = 0
			mc_functions.strip_finder(cmb_map_sim, ang_rad[k], NSIDE)
			circle_a = mc_functions.load_file('strip_a', bins)
			circle_b = mc_functions.load_file('strip_b', bins)
			x_corr_sims_180[j,k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
		#lat -= 45
		print("j="+str(j), time.time() - start)
	np.savetxt("/opt/local/l4astro/rbbg94/data/x_corr_sims_ngp.csv", x_corr_sims_0, delimiter = ',')
	np.savetxt("/opt/local/l4astro/rbbg94/data/x_corr_sims_0.csv", x_corr_sims_180, delimiter = ',')
	"""std_dev_0 = np.std(x_corr_sims_0, axis = 0)
	np.savetxt("/opt/local/l4astro/rbbg94/data/plane_corrs_and_errs/lon_"+str(lon)+"_lat_0_sim_err_"+str(no_sims)+"_sims.csv", std_dev_0, delimiter = ',')
	std_dev_45 = np.std(x_corr_sims_45, axis = 0)
	np.savetxt("/opt/local/l4astro/rbbg94/data/plane_corrs_and_errs/lon_"+str(lon)+"_lat_45_sim_err_"+str(no_sims)+"_sims.csv", std_dev_45, delimiter = ',')
	lon += 30"""
	print("i="+str(i), time.time() - start)
