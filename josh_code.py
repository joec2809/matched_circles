#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 10:06:36 2019

@author: qgcs22
"""

import numpy as np
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as pyplot

if __name__ == '__main__':

	'''Healpy'''
	data = hp.fitsfunc.read_map('/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits', field = None, nest = True)

def circle_selector(data, theta_i, vector, i):
    rad1 = 2*np.pi*(theta_i+0.35)/360
    rad2 = 2*np.pi*(theta_i-0.35)/360
    vectorb = -vector
    
    '''Healpy'''
    circle1a = hp.query_disc(2048, vector, rad1, nest = True)
    circle1b = hp.query_disc(2048, vector, rad2, nest = True)
    circle2a = hp.query_disc(2048, vectorb, rad1, nest = True)
    circle2b = hp.query_disc(2048, vectorb, rad2, nest = True)
    
    strip1 = np.setdiff1d(circle1a, circle1b)
    strip2 = np.setdiff1d(circle2a, circle2b)
    strip1_data = np.zeros((len(strip1),3))
    strip2_data = np.zeros((len(strip2),3))
    
    '''Healpy'''
    lon1, lat1 = hp.pixelfunc.pix2ang(2048, strip1, nest=True, lonlat=True)
    lon2, lat2 = hp.pixelfunc.pix2ang(2048, strip2, nest=True, lonlat=True)
    
    strip1_data[:,0] = data[strip1]
    strip1_data[:,1] = lon1
    strip1_data[:,2] = lat1
    strip2_data[:,0] = data[strip2]
    strip2_data[:,1] = lon2
    strip2_data[:,2] = lat2
    
    x = str(i)
    fname1 = 'circle'+x+'a'
    fname2 = 'circle'+x+'b'
    
    col11 = fits.Column(name='index', array = strip1,format='D')
    col12 = fits.Column(name='T', array = strip1_data[:,0],format='D')
    col13 = fits.Column(name='long', array=lon1, format='D')
    col14 = fits.Column(name='lat', array=lat1, format='D')
    u=fits.BinTableHDU.from_columns([col11,col12,col13,col14])
    u.writeto('/opt/local/l4astro/rbbg94/cmb_maps/'+fname1, overwrite=True)
    
    col21 = fits.Column(name='index', array = strip2,format='D')
    col22 = fits.Column(name='T', array = strip2_data[:,0],format='D')
    col23 = fits.Column(name='long', array=lon2, format='D')
    col24 = fits.Column(name='lat', array=lat2, format='D')
    v=fits.BinTableHDU.from_columns([col21,col22,col23,col24])
    v.writeto('/opt/local/l4astro/rbbg94/cmb_maps/'+fname2, overwrite=True)


def load_file(fname1,fname2):
    hdul = fits.open(fname1)
    circ_dat = hdul[1].data
    hdul.close()
    T1 = circ_dat['T']
    lon1 = circ_dat['long']
    lat1 = circ_dat['lat']
    index1 = circ_dat['index']
    
    hdul = fits.open(fname2)
    circ_dat2 = hdul[1].data
    hdul.close()
    T2 = circ_dat2['T']
    lon2 = circ_dat2['long']
    lat2 = circ_dat2['lat']
    index2 = circ_dat2['index']
    
    circlea = np.zeros((len(T1),4))
    circlea[:,3] = index1
    circlea[:,0] = T1
    circlea[:,1] = lon1
    circlea[:,2] = lat1
    circleb = np.zeros((len(T2),4))
    circleb[:,3] = index2
    circleb[:,0] = T2
    circleb[:,1] = lon2
    circleb[:,2] = lat2
    
    return circlea, circleb


def match_circle(data1, data2):

    T1_bin = np.zeros((360, 2))
    T2_bin = np.zeros((360, 2))
    for i in range(len(data1[:,1])):
        iy = int(data1[i,1])
        
        T1_bin[iy,1] += data1[i,0]
        T1_bin[iy,0] += 1
        
    for i in range(len(data2[:,2])):
        iy = int(data2[i,1])
                
        T2_bin[iy,1] += data2[i,0]
        T2_bin[iy,0] += 1
    
    T1 = np.zeros(len(T1_bin[:,0]))
    T2 = np.zeros(len(T2_bin[:,0]))
    for i in range(len(T1_bin[:,1])):
        T1[i] = T1_bin[i,1]/T1_bin[i,0]

    for i in range(len(T2_bin[:,1])):
        T2[i] = T2_bin[i,1]/T2_bin[i,0]

    sum_val = np.zeros(360)
    ind_val = np.zeros((360, len(T1)))
    circle1_tot = 0
    circle2_tot = 0
    

    for i in range(360):
        circle1_tot += T1[i]**2
        circle2_tot += T2[i]**2
        for j in range(len(T1)):
            if i + j >= 360:
                x = (i + j - 360)
            else:
                x = (i + j)
            ind_val[i][j] = T1[x]*T2[j]
            sum_val[i] = sum_val[i] + ind_val[i][j]
            
    norm_all = sum_val/(np.sqrt(circle1_tot)*np.sqrt(circle2_tot))    
    m = np.amax(norm_all)
    lag = np.argwhere(norm_all == m)
    max_norm = norm_all[lag]    
    
    
    count1 = 0    
    count2 = 0
    smooth_factor = 10
    smooth_dat1 = np.zeros(len(T1)/smooth_factor)
    smooth_dat2 = np.zeros(len(T2)/smooth_factor)
    
    for i in range(len(smooth_dat1)):
        smooth_dat1[i] = np.sum(T1[count1:smooth_factor*(i+1)])/smooth_factor
        count1 += smooth_factor
    
    for i in range(len(smooth_dat2)):
        smooth_dat2[i] = np.sum(T2[count2:smooth_factor*(i+1)])/smooth_factor
        count2 += smooth_factor    
    
    #num = np.arange(0, len(smooth_dat1), 1)

    #for i in range(360/smooth_factor):
    #    pyplot.scatter(i, smooth_dat1[i], s=10, c='g')
    #pyplot.plot(num, smooth_dat1, c = 'g')
    #
    #for i in range(360/smooth_factor):
    #    pyplot.scatter(i, smooth_dat2[i], s=10, c='r')
    #pyplot.plot(num, smooth_dat2, c = 'r')

    #pyplot.ylim(-0.0002, 0.0002)
    #pyplot.title('Smoothed T')
    #pyplot.xlabel('Longitude/smoothing factor')
    #pyplot.ylabel('T')
    #pyplot.show()
    
    return max_norm, norm_all

"""domains = np.arange(1, 100, 1)
theta = y = np.zeros(len(domains))
for i in range(len(domains)):
    y[i] = domains[i]*160
    if y[i]<13000:
        theta[i] = 360*np.arccos(y[i]/13000)/(2*np.pi)   
    else:
        theta = theta[:i]
        
for i in range(len(theta)):
    vector = np.zeros(3)
    vector[0] = 0
    vector[1] = 0
    vector[2] = 1
    circle_selector(theta[i],vector,i)
    
circle_norms = np.zeros(len(theta))
norm_all = np.zeros((len(theta),360))
for i in range(len(theta)):
    if 90-theta[i] < 5:
        circle_norms[i] = 0
    else:
        x = str(i)
        circle_norms[i] = match_circle(load_file('circle'+x+'a','circle'+x+'b')[0],load_file('circle'+x+'a','circle'+x+'b')[1])[0]
    #print [circle_norms[i], theta[i], i]
#print np.amax(circle_norms)
#np.savetxt('N-S_axis_xcor_vals', circle_norms)
    norm_all[i] = match_circle(load_file('circle'+x+'a','circle'+x+'b')[0],load_file('circle'+x+'a','circle'+x+'b')[1])[1]    
    j = np.arange(0,360,1)
    pyplot.plot(j, norm_all[i])
    pyplot.xlabel('Lag')
    pyplot.ylabel('Xcor')   
    pyplot.title('Circle pair at latitude +-'+str(90-round(theta[i],2)))
    pyplot.savefig(str(round(theta[i],2))+'.png')    
    pyplot.show()"""

