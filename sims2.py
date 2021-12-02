from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
import time
import mc_functions

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

lon = 207.8
lat = 33.7

cs_lon = 207.8
cs_lat = -56.3
cs_vec = hp.pixelfunc.ang2vec(cs_lon, cs_lat, lonlat=True)

cs_hat = (1/np.sqrt((cs_vec[0])**2+(cs_vec[1])**2+(cs_vec[2])**2))*cs_vec
rot_ang = np.pi/6
n1 = cs_hat[0]
n2 = cs_hat[1]
n3 = cs_hat[2]
rot_mat = np.zeros((3,3))

vec_lon = 207.8
vec_lat = cs_lat + 90
vec = hp.pixelfunc.ang2vec(vec_lon, vec_lat, lonlat=True)

cmb_map_start = mc_functions.rotate_to_top(cmb_map_og*mask_ring, vec_lon, vec_lat)

rot_vec = hp.pixelfunc.ang2vec(207.8, 0, lonlat=True)

ang_rad = np.linspace((2*np.pi)/360, np.pi/2, 90)[:-1]
no_circles = len(ang_rad)
T_m = mc_functions.T_m_setup(m_max, bins)

start = time.time()

a=4

for i in range(6):
	c = np.cos(np.pi/12+i*rot_ang)
	s = np.sin(np.pi/12+i*rot_ang)
	rot_mat[0,0] = c+n1**2*(1-c)
	rot_mat[0,1] = n1*n2*(1-c)-n3*s
	rot_mat[0,2] = n1*n3*(1-c)+n2*s
	rot_mat[1,0] = n1*n2*(1-c)+n3*s
	rot_mat[1,1] = c+n2**2*(1-c)
	rot_mat[1,2] = n2*n3*(1-c)-n1*s
	rot_mat[2,0] = n1*n3*(1-c)-n2*s
	rot_mat[2,1] = n2*n3*(1-c)+n1*s
	rot_mat[2,2] = c+n3**2*(1-c)
	new_vec = np.dot(rot_mat, vec)
	lon_lat = hp.pixelfunc.vec2ang(new_vec, lonlat = True)
	lon = lon_lat[0]
	lat = lon_lat[1]
	cmb_map = hpt.rotate_map(cmb_map_start, rot_vec, np.pi/12+i*rot_ang)
	x_corr_0 = np.zeros(no_circles, dtype=complex)
	for j in range(no_circles):
		rad_lag = 0
		mc_functions.strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = mc_functions.load_file('strip_a', bins)
		circle_b = mc_functions.load_file('strip_b', bins)
		x_corr_0[j] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
	np.savetxt('/opt/local/l4astro/rbbg94/data/cs_corrs_and_errs/lon_'+str(lon)+'_lat_'+str(lat)+'_corr_0_'+str(a)+'.csv', x_corr_0, delimiter = ',')

	x_corr_sims_0 = np.zeros((no_sims,no_circles), dtype=complex)
	for j in range(no_sims):
		cmb_map_sim = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/sims/dx12_v3_smica_cmb_mc_000"+str(j).zfill(2)+"_raw.fits")
		NSIDE = hp.npix2nside(len(cmb_map_sim))
		
		cmb_map_sim = mc_functions.rotate_to_top(cmb_map_sim, vec_lon, vec_lat)
		cmb_map_sim = hpt.rotate_map(cmb_map_sim, rot_vec, np.pi/12+i*rot_ang)
		for k in range(no_circles):
			rad_lag = 0
			mc_functions.strip_finder(cmb_map_sim, ang_rad[k], NSIDE)
			circle_a = mc_functions.load_file('strip_a', bins)
			circle_b = mc_functions.load_file('strip_b', bins)
			x_corr_sims_0[j,k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
		
		print("j="+str(j), time.time() - start)
	std_dev_0 = np.std(x_corr_sims_0, axis = 0)
	np.savetxt("/opt/local/l4astro/rbbg94/data/cs_corrs_and_errs/lon_"+str(lon)+"_lat_"+str(lat)+"_sim_err_"+str(no_sims)+"_sims_"+str(a)+".csv", std_dev_0, delimiter = ',')
	a+=1
	if a == 7:
		a=4
	print("i="+str(i), time.time() - start)
