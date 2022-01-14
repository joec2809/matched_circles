from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import mc_functions
import time

CMB_DIST = 14000
CELL_SIZE = 320

no_sims = 100
bins = 360
m_max = 720
rad_lag = 0

lon = 0
lat = 88
no_dir = 10
rot_ang = 2*np.pi/no_dir

ang_rad = np.linspace((2*np.pi)/360, np.pi/2, 90)[:-1]
no_circles = len(ang_rad)
T_m = mc_functions.T_m_setup(m_max, bins)

start = time.time()

for i in range(no_dir):
	x_corr_sims = np.zeros((no_sims,no_circles), dtype=complex)
	for j in range(no_sims):
		cmb_map_sim = hp.fitsfunc.read_map("../sims/dx12_v3_smica_cmb_mc_000"+str(j).zfill(2)+"_raw.fits")
		NSIDE = hp.npix2nside(len(cmb_map_sim))
		
		for k in range(no_circles):
			mc_functions.strip_finder(cmb_map_sim, ang_rad[k], NSIDE)
			circle_a = mc_functions.load_file('strip_a', bins)
			circle_b = mc_functions.load_file('strip_b', bins)
			x_corr_sims[j,k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
		
		print("j="+str(j), time.time() - start)

	std_dev = np.std(x_corr_sims, axis = 0)
	np.savetxt(f"../data/near_ngp/near_ngp_lon_{180*lon/np.pi}_{no_sims}_sims.csv", std_dev, delimiter = ',')
	
	lon += rot_ang
	print("i="+str(i), time.time() - start)
