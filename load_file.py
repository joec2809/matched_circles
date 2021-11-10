from __future__ import division
from astropy.io import fits
import numpy as np
import pandas as pd


def load_file(fname, bins):
	hdul = fits.open('/opt/local/l4astro/rbbg94/cmb_maps/'+fname)
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
