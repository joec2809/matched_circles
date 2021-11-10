from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

ang_rad = np.linspace(1,89,89)
x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_corr_0_thin.csv', dtype=complex, delimiter = ',')
errors = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_sim_err_100_sims.csv', dtype=complex, delimiter = ',')
errors = 0

fig, ax = plt.subplots(figsize = (14,8))
#ax.errorbar(ang_rad, x_corr, yerr = errors, color = 'k', ecolor = 'k', elinewidth = 0.8, capsize = 3)
ax.plot(ang_rad, x_corr, color = 'k')
plt.axhline(y=0, color = 'k')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.annotate('Phase=0$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax.fill_between(ang_rad, 2*errors, -2*errors, color = 'lightgrey')
ax.fill_between(ang_rad, errors, -errors, color = 'darkgrey')
plt.tight_layout()

#fig.savefig("/opt/local/l4astro/rbbg94/figures/corr_0_ngp_sim_err_signif.png", overwrite = True)

plt.show()
