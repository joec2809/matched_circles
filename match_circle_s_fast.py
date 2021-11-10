from __future__ import division
import numpy as np

def match_circle_s_fast(data1, data2, phase):	
	s_den = 0
	m = np.arange(0,360,1)
	T_m_a = np.fft.fft((data1[0]))
	T_m_b = np.fft.fft((data2[1]))
	T_m_a_sum = np.abs(np.sum(T_m_a))
	T_m_b_sum = np.abs(np.sum(T_m_b))
	lag = np.full(360, np.exp(-1j*phase))
	m_lag = np.multiply(m,lag)
	T_m_a = np.multiply(T_m_a, m)
	for m in range(360):
		s_den += m*((T_m_a_sum)**2+(T_m_b_sum)**2)

	if s_den != 0:
		s_m = (2*np.multiply(T_m_a,np.conj(T_m_b)))/s_den
		s = np.sum(np.multiply(s_m, m_lag))
	else:
		s=0
	
	return s
