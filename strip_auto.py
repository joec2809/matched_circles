from __future__ import division
import numpy as np
import healpy as hp
from astropy.io import fits



NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

def strip_auto(data, ang_rad, nside):

	ipix_strip1 = hp.query_strip(NSIDE, ang_rad-(np.pi/360), ang_rad+(np.pi/360))
	ipix_strip2 = hp.query_strip(NSIDE, ang_rad-(np.pi/360), ang_rad+(np.pi/360))

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
