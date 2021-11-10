from __future__ import division
import numpy as np

def match_circle_r(data1, data2):
	
	T1_bin = [[] for _ in range(360)]
	T2_bin = [[] for _ in range(360)]

	for i in range(len(data1[:,1])):
		iy = int(data1[i,1])
		if data1[i,0] != 0:
			T1_bin[iy].append(data1[i,0])
	T1_bin = np.array([np.array(T1i) for T1i in T1_bin])
	
			
	for i in range(len(data2[:,2])):
		iy = int(data2[i,1])
		if data2[i,0] != 0:      	  
			T2_bin[iy].append(data2[i,0])
	T2_bin = np.array([np.array(T2i) for T2i in T2_bin])
	
	T1_binn = np.zeros((2,len(T1_bin)))
	for i in range(len(T1_bin)):
		T1_binn[1,i] = np.sum(T1_bin[i])
		T1_binn[0,i] = T1_bin[i].size

	T2_binn = np.zeros((2,len(T2_bin)))
	for i in range(len(T2_bin)):
		T2_binn[1,i] = np.sum(T2_bin[i])
		T2_binn[0,i] = T2_bin[i].size

	T = np.zeros((2,len(T1_binn[0,:])))
	

	for i in range(len(T1_binn[1,:])):
		if T1_binn[0,i] != 0:
			T[0,i] = T1_binn[1,i]/T1_binn[0,i]
		else:
			T[0,i] = 0
		if T2_binn[0,i] != 0:
			T[1,i] = T2_binn[1,i]/T2_binn[0,i]
		else:
			T[1,i] = 0


	sum_val = np.zeros(360)
	ind_val = np.zeros((360, len(T[0])))
	circle1_tot = 0
	circle2_tot = 0

	for i in range(360):
		circle1_tot += T[0,i]**2
		circle2_tot += T[1,i]**2
		for j in range(len(T[0])):
		    if i + j >= 360:
		        x = (i + j - 360)
		    else:
		        x = (i + j)
			ind_val[i][j] = T[0,x]*T[1,j]
		    sum_val[i] = sum_val[i] + ind_val[i][j]
	
	if circle1_tot != 0 or circle2_tot != 0:
		norm_all = sum_val/(np.sqrt(circle1_tot)*np.sqrt(circle2_tot))    
	else:
		norm_all = np.arange(0,360,1)

	return norm_all
    
