from __future__ import division
import numpy as np

def T_m_setup(m_max, bins):
	T_m = np.full((m_max,bins), np.exp(-1j), dtype = complex)
	m = np.arange(0,m_max,1)
	phi = np.arange(0,2*np.pi,2*np.pi/bins)
	T_m = T_m**phi
	T_m = np.swapaxes(T_m, 0, 1)
	T_m = T_m**m
	T_m = np.swapaxes(T_m, 0, 1)
	return T_m
