from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from circle_finder import circle_finder

CMB_DIST = 14000
CELL_SIZE = 320

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/(360))
data = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_corr.csv', dtype = complex, delimiter = ',')
peaks_ind= np.argwhere(np.logical_and(ang_rad>=((40/360)*2*np.pi), ang_rad<=((65/360)*2*np.pi)))
peaks_rad_pred = circle_finder(CELL_SIZE, peaks_ind, 0)
"""ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/(360*3))
peaks_ind= np.argwhere(np.logical_and(ang_rad>=((40/360)*2*np.pi), ang_rad<=((65/360)*2*np.pi)))"""
peaks_rad = ang_rad[peaks_ind]
peaks_data = data[peaks_ind]
#std_dev = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_lag_err.csv', delimiter = ',')
#peaks_err = std_dev[peaks_ind]
peaks_err = 0

peaks_rad = peaks_rad*(360/(2*np.pi))
peaks_rad_pred = peaks_rad_pred*(360/(2*np.pi))

fig, ax = plt.subplots(figsize=(18, 10))

ax.errorbar(peaks_rad, peaks_data, yerr = peaks_err, ecolor = 'red')
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.annotate('Phase=0$^\circ$', xy = (0.02,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax.axhline(0, color = 'black')
#plt.xticks(np.arange(0, 91, 10))
plt.xlim(40,65)
#for i in range(len(peaks_rad_pred)):
#	ax.axvline(peaks_rad_pred[i], color = 'red', ls = '--')
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_ngp_no_err_peaks.png', overwrite = True)

plt.show()
