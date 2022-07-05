from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import truncnorm
import pandas as pd
import awkward as ak
import mc_functions


SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

mode = "p3"

lon = 144.0
lat = 85.0

pred = pd.read_csv(f"../data/pred_angs/omega_{mode}_reflection.txt", delim_whitespace=True, header = 1, index_col=0)

ang_rad = np.array(pred['alpha'])[~np.isnan(pred['alpha'])]

x_corr = np.genfromtxt(f'../data/grid_search/mode_{mode}/{mode}_lon_{lon}_lat_{lat}.csv', dtype=complex, delimiter = ',')
errors = np.genfromtxt(f'../data/{mode}_lon_{int(lon)}_lat_{int(lat)}_100_sims.csv', dtype=complex, delimiter = ',')
#errors_real = np.zeros(len(errors))
#errors = 0

print(ang_rad[np.argmax(x_corr)])

#x_corr, errors = mc_functions.array_size_match(x_corr, errors)


fig, ax1 = plt.subplots(figsize = (14,8), constrained_layout = True)
#ax1.errorbar(ang_rad, x_corr, yerr = errors, color = 'k', ecolor = 'k', elinewidth = 0.8, capsize = 3)
ax1.plot(ang_rad, x_corr, color = 'k')
ax1.axhline(y=0, color = 'k')
ax1.set_xticks(np.arange(0, 91, 10))
ax1.set_xlim(0,90)
ax1.set_xlabel(r'$\alpha/^\circ$')
ax1.set_ylabel('$S$')
ax1.annotate('Phase=0$^\circ$', xy = (0.8,0.9), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax1.fill_between(ang_rad, 2*errors, -2*errors, color = 'lightgrey')
ax1.fill_between(ang_rad, errors, -errors, color = 'darkgrey')

plt.show()
fig.savefig(f"../figures/{mode}_lon_{lon}_lat_{lat}.png")