from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import mc_functions
import time
from tqdm import tqdm
from scipy.spatial.transform import Rotation as R

CMB_DIST = 14000
CELL_SIZE = 320

no_sims = 100
bins = 360
m_max = 720
rad_lag = 0

lon = 238.2
lat = -28.8
og_vec = hp.pixelfunc.ang2vec(lon, lat, lonlat=True)
vec_lon = 238.2
vec_lat = lat + 90
vec = hp.pixelfunc.ang2vec(vec_lon, vec_lat, lonlat=True)

no_dir = 12
rot_ang = 2*np.pi/no_dir

ang_rad = np.linspace((2*np.pi)/360, np.pi/2, 90)[:-1]
no_circles = len(ang_rad)
T_m = mc_functions.T_m_setup(m_max, bins)

start = time.time()

for i in tqdm(range(no_dir)):
	x_corr_sims = np.zeros((no_sims,no_circles), dtype=complex)
	"""for j in range(no_sims):
		cmb_map_og = hp.fitsfunc.read_map("../sims/dx12_v3_smica_cmb_mc_000"+str(j).zfill(2)+"_raw.fits")
		cmb_map_sim = mc_functions.rotate_to_top(cmb_map_og, vec_lon, vec_lat)
		NSIDE = hp.npix2nside(len(cmb_map_sim))
		
		for k in range(no_circles):
			mc_functions.strip_finder(cmb_map_sim, ang_rad[k], NSIDE)
			circle_a = mc_functions.load_file('strip_a', bins)
			circle_b = mc_functions.load_file('strip_b', bins)
			x_corr_sims[j,k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
		
		print("j="+str(j), time.time() - start)

	std_dev = np.std(x_corr_sims, axis = 0)
	np.savetxt(f"../data/dipole_dir_{i}_{no_sims}_sims.csv", std_dev, delimiter = ',')"""

	rot_vec = og_vec * rot_ang*(i+1)
	rotation = R.from_rotvec(rot_vec)
	new_vec = rotation.apply(vec)
	vec_lat, vec_lon = hp.pixelfunc.vec2ang(new_vec, lonlat = True)
	print(vec_lat, vec_lon)