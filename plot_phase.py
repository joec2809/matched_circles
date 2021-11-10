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

phase = np.linspace(0,359,360)
x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/87.9_s_corr.csv', dtype=complex, delimiter = ',')

fig, ax = plt.subplots(figsize = (14,8))
ax.plot(phase, x_corr, color = 'k')
plt.axhline(y=0, color = 'k')
plt.xlim(0,360)
plt.ylim(-0.3, 0.3)
ax.set_xlabel('Phase/$^\circ$')
ax.set_ylabel('$S$')
plt.tight_layout()

fig.savefig("/opt/local/l4astro/rbbg94/figures/87.9_s_corr.png", overwrite = True)

plt.show()
