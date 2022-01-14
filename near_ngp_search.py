from __future__ import division
import numpy as np
import healpy as hp
from astropy.io import fits
from astrotools import healpytools as hpt
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
rad_lag = 0
m_max = 720

lon = 0
lat = 88
no_dir = 10
rot_ang = 2*np.pi/no_dir



ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros(no_circles, dtype = complex)

start = time.time()

T_m = mc_functions.T_m_setup(m_max, bins)

for i in range(no_dir):
    cmb_map = mc_functions.rotate_to_top(cmb_map_og*mask_ring, lon, lat)
    for j in range(len(ang_rad)):
        mc_functions.strip_finder(cmb_map, ang_rad[j], NSIDE)
        circle_a = mc_functions.load_file('strip_a', bins)
        circle_b = mc_functions.load_file('strip_b', bins)
        x_corr[j] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, bins)
    np.savetxt(f'../data/near_ngp/near_ngp_lon_{180*lon/np.pi}.csv', x_corr, delimiter = ',')
    lon += rot_ang
    print(time.time()-start)