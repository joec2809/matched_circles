from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt

def rotate_to_top(cmb_map, lon, lat):
	og_vec = hp.pixelfunc.ang2vec(lon, lat, lonlat=True)
	z_vec = np.array([0,0,1])

	rot_vec = np.cross(z_vec, og_vec)

	rot_ang = np.arccos(np.dot(og_vec, z_vec)/(np.linalg.norm(og_vec)*np.linalg.norm(z_vec)))

	cmb_map_2 = hpt.rotate_map(cmb_map, rot_vec, rot_ang)

	return cmb_map_2
