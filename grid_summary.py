from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm, truncnorm
import pandas as pd
import awkward as ak


SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

mode = "1"

pred = pd.read_csv(f"../data/pred_angs/omega_{mode}_reflection.txt", delim_whitespace=True, header = 1, index_col=0)

ang_rad = np.array(pred['alpha'])[~np.isnan(pred['alpha'])]

lat = np.radians(np.arange(75,90,2))
init_no_dir = 360/2
no_dir = (init_no_dir*np.cos(lat)).astype(int)
lons = ak.ArrayBuilder()
for i in range(len(no_dir)):
    lon = np.radians(np.linspace(0,360,no_dir[i]+1))
    lons.append(np.delete(lon,-1))


probs = []
no_pos_points = []
no_neg_points = []

table_lat = []

table_lon = []

signal_detectors = []


for i in range(len(lat)):
    for j in range(len(lons[i])):
        table_lat.append(lat[i])
        table_lon.append(lons[i][j])
        x_corr = np.genfromtxt(f'../data/grid_search/mode_{mode}/{mode}_lon_{round(np.degrees(lons[i][j]),1)}_lat_{round(np.degrees(lat[i]),1)}.csv', dtype=complex, delimiter = ',')
        x_corr_real = np.zeros(len(x_corr))

        pos_points = 0
        neg_points = 0

        for k, point in enumerate(x_corr):
            x_corr_real[k] = float(np.real(point))
            
            if np.real(point) >= 0:
                pos_points += 1
            elif np.real(point) < 0:
                neg_points += 1

        signal_detectors.append(np.sum(x_corr_real))
        mean = np.mean(x_corr_real)
        var = np.std(x_corr_real)**2

        print(var)

        dist_mean = norm.stats(moments = 'm')


        z_score = (mean-dist_mean)/np.sqrt(var/len(x_corr_real))

        prob = 2*norm.cdf(-np.abs(z_score))

        probs.append(prob)
        no_pos_points.append(pos_points)
        no_neg_points.append(neg_points)


signal_detectors = np.array(signal_detectors)
probs = np.array(probs)
no_pos_points = np.array(no_pos_points)
no_neg_points = np.array(no_neg_points)

table_lats = np.degrees(table_lat)
table_lons = np.round(np.degrees(table_lon), 1)

data = {'Lat.': table_lats, 'Lon.': table_lons, 'Signal Detector': signal_detectors, 'p values': probs, 'no. +ve points': no_pos_points, 'no. -ve points': no_neg_points}

df = pd.DataFrame(data)

print(df[['Signal Detector', 'p values', 'no. +ve points', 'no. -ve points']])

df.to_csv(f'../data/grid_search/{mode}_grid_summary.csv',index = False)

print(signal_detectors[np.argwhere(no_pos_points-no_neg_points >= 30)])

print(df.loc[df['Signal Detector'] >= 2])

fig, ax = plt.subplots(constrained_layout = True)

ax.scatter(no_pos_points-no_neg_points, signal_detectors)

plt.show()