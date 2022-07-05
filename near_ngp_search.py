from __future__ import division
import numpy as np
import healpy as hp
from astropy.io import fits
from astrotools import healpytools as hpt
import awkward as ak
import pandas as pd
import mc_functions
import time

mode = "1"

if mode == "p3":
    pred = pd.read_csv("../data/pred_angs/omega_p3_reflection.txt", delim_whitespace=True, header = 1, index_col=0)
elif mode == "1":
    pred = pd.read_csv("../data/pred_angs/omega_1_reflection.txt", delim_whitespace=True, header = 1, index_col=0)

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

lat = np.radians(np.arange(75,90,2))
init_no_dir = 360/2
no_dir = (init_no_dir*np.cos(lat)).astype(int)
lons = ak.ArrayBuilder()
for i in range(len(no_dir)):
    lon = np.radians(np.linspace(0,360,no_dir[i]+1))
    lons.append(np.delete(lon,-1))



ang_rad = np.radians(np.array(pred['alpha'])[~np.isnan(pred['alpha'])])
x_corr = np.empty(len(ang_rad), dtype = complex)

start = time.time()

T_m = mc_functions.T_m_setup(m_max, bins)

for i in range(len(lat)):
    for j in range(len(lons[i])):
        cmb_map = mc_functions.rotate_to_top(cmb_map_og*mask_ring, lons[i][j], lat[i])
        for k in range(len(ang_rad)):
            mc_functions.strip_finder(cmb_map, ang_rad[k], NSIDE)
            circle_a = mc_functions.load_file('strip_a', bins)
            circle_b = mc_functions.load_file('strip_b', bins)
            x_corr[k] = mc_functions.match_circle_s(circle_a, circle_b, T_m, m_max, rad_lag)
        np.savetxt(f'../data/grid_search/mode_{mode}/{mode}_lon_{round(np.degrees(lons[i][j]), 1)}_lat_{round(np.degrees(lat[i]), 1)}.csv', x_corr, delimiter = ',')
        print(i, j, time.time()-start)