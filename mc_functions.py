from __future__ import division
import numpy as np
import healpy as hp
import pandas as pd
from astropy.io import fits
from astrotools import healpytools as hpt

# Finds angular radii of circles based on cell size
def circle_finder(cell_size, cmb_dist, x, y):
	d = cell_size*(np.sqrt(x**2+y**2))
	ang_rad = np.arccos(d/(2*cmb_dist))
	return ang_rad

# Selects strips of CMB data with certain angular radii
def strip_finder(data, ang_rad, nside):

	ipix_strip1 = hp.query_strip(nside, ang_rad-0.1*(np.pi/360), ang_rad+0.1*(np.pi/360))
	ipix_strip2 = hp.query_strip(nside, np.pi-ang_rad-0.1*(np.pi/360), np.pi-ang_rad+0.1*(np.pi/360))

	strip1_data = np.zeros((len(ipix_strip1), 3))
	strip2_data = np.zeros((len(ipix_strip2), 3))

	lon1, lat1 = hp.pixelfunc.pix2ang(2048, ipix_strip1, lonlat=True)
	lon2, lat2 = hp.pixelfunc.pix2ang(2048, ipix_strip2, lonlat=True)

	strip1_data[:,0] = data[ipix_strip1]
	strip1_data[:,1] = lon1
	strip1_data[:,2] = lat1
	strip2_data[:,0] = data[ipix_strip2]
	strip2_data[:,1] = lon2
	strip2_data[:,2] = lat2

	fname1 = 'strip_a'
	fname2 = 'strip_b'

	col11 = fits.Column(name='index', array = ipix_strip1,format='D')
	col12 = fits.Column(name='T', array = strip1_data[:,0],format='D')
	col13 = fits.Column(name='long', array=lon1, format='D')
	col14 = fits.Column(name='lat', array=lat1, format='D')
	u=fits.BinTableHDU.from_columns([col11,col12,col13,col14])
	u.writeto('../cmb_maps/'+fname1, overwrite=True)

	col21 = fits.Column(name='index', array = ipix_strip2,format='D')
	col22 = fits.Column(name='T', array = strip2_data[:,0],format='D')
	col23 = fits.Column(name='long', array=lon2, format='D')
	col24 = fits.Column(name='lat', array=lat2, format='D')
	v=fits.BinTableHDU.from_columns([col21,col22,col23,col24])
	v.writeto('../cmb_maps/'+fname2, overwrite=True)

# Provides strips to use with auto correlation
def strip_auto(data, ang_rad, nside):

	ipix_strip1 = hp.query_strip(nside, ang_rad-(np.pi/360), ang_rad+(np.pi/360))
	ipix_strip2 = hp.query_strip(nside, ang_rad-(np.pi/360), ang_rad+(np.pi/360))

	strip1_data = np.zeros((len(ipix_strip1), 3))
	strip2_data = np.zeros((len(ipix_strip2), 3))

	lon1, lat1 = hp.pixelfunc.pix2ang(2048, ipix_strip1, lonlat=True)
	lon2, lat2 = hp.pixelfunc.pix2ang(2048, ipix_strip2, lonlat=True)

	strip1_data[:,0] = data[ipix_strip1]
	strip1_data[:,1] = lon1
	strip1_data[:,2] = lat1
	strip2_data[:,0] = data[ipix_strip2]
	strip2_data[:,1] = lon2
	strip2_data[:,2] = lat2

	fname1 = 'strip_a'
	fname2 = 'strip_b'

	col11 = fits.Column(name='index', array = ipix_strip1,format='D')
	col12 = fits.Column(name='T', array = strip1_data[:,0],format='D')
	col13 = fits.Column(name='long', array=lon1, format='D')
	col14 = fits.Column(name='lat', array=lat1, format='D')
	u=fits.BinTableHDU.from_columns([col11,col12,col13,col14])
	u.writeto('/opt/local/l4astro/rbbg94/cmb_maps/'+fname1, overwrite=True)

	col21 = fits.Column(name='index', array = ipix_strip2,format='D')
	col22 = fits.Column(name='T', array = strip2_data[:,0],format='D')
	col23 = fits.Column(name='long', array=lon2, format='D')
	col24 = fits.Column(name='lat', array=lat2, format='D')
	v=fits.BinTableHDU.from_columns([col21,col22,col23,col24])
	v.writeto('/opt/local/l4astro/rbbg94/cmb_maps/'+fname2, overwrite=True)


# Loads a saved fits file ands puts the temperature in bins defined by galactic longitude
def load_file(fname, bins):
	hdul = fits.open('../cmb_maps/'+fname)
	circ_dat = hdul[1].data
	hdul.close()
	T = circ_dat['T']
	lon = circ_dat['long']
	lat = circ_dat['lat']
	index = circ_dat['index']
		
	circle = np.zeros((len(T),4))
	circle[:,3] = index
	circle[:,0] = T
	circle[:,1] = lon
	circle[:,2] = lat
	
	sort_ind = circle[:,1].argsort()

	circle = circle[sort_ind]

	array = np.linspace(0,bins-1,bins)
	lon_bin = np.array(np.cumsum(pd.cut(circle[:,1], bins = array, right = False, include_lowest = True).value_counts()))
	T_bin = np.array(np.split(circle[:,0], lon_bin))
	T_binn = [np.mean(arr) for arr in T_bin]

	return T_binn


# Calculate the S correlation between two strips 
def match_circle_s(data1, data2, T_m, m_max, phase):
	m = np.arange(0,m_max,1)
	T_m_1 = (1/(2*np.pi))*np.sum(data1*T_m, axis = 1)
	T_m_2 = (1/(2*np.pi))*np.sum(data2*T_m, axis = 1)
	T_m_2_star = np.conj(T_m_2)
	phase_factor = (np.full(m_max, np.exp(-1j*phase), dtype = complex))**m
	s_num = 2*np.sum(m*T_m_1*T_m_2_star*phase_factor)
	s_den = np.sum(m*((np.abs(T_m_1))**2+(np.abs(T_m_2))**2))

	if s_den != 0:
		s = s_num/s_den
	else:
		s=0
	
	return s

# Rotate CMB map so vector defined by longitude and 
def rotate_to_top(cmb_map, lon, lat):
	og_vec = hp.pixelfunc.ang2vec(lon, lat, lonlat=True)
	z_vec = np.array([0,0,1])

	rot_vec = np.cross(z_vec, og_vec)

	rot_ang = np.arccos(np.dot(og_vec, z_vec)/(np.linalg.norm(og_vec)*np.linalg.norm(z_vec)))

	cmb_map_2 = hpt.rotate_map(cmb_map, rot_vec, rot_ang)

	return cmb_map_2

# Setup of arrays used for S statistic calculation
def T_m_setup(m_max, bins):
	T_m = np.full((m_max,bins), np.exp(-1j), dtype = complex)
	m = np.arange(0,m_max,1)
	phi = np.arange(0,2*np.pi,2*np.pi/bins)
	T_m = T_m**phi
	T_m = np.swapaxes(T_m, 0, 1)
	T_m = T_m**m
	T_m = np.swapaxes(T_m, 0, 1)
	return T_m

def stack(*args):
	n = len(args)
	array = np.zeros(len(args[0]), dtype=complex)
	for arg in args:
		array += arg

	return array/n

def stack_errors(*args):
	n = len(args)
	array = np.zeros(len(args[0]), dtype=complex)
	for arg in args:
		array += arg**2

	return np.sqrt(array)/np.sqrt(n)