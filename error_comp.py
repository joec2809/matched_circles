from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


sim = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_sim_err_100_sims.csv', dtype = complex, delimiter = ',')
lag = np.genfromtxt('/opt/local/l4astro/rbbg94/data/all_ngp_lag_errs/ngp_lag_err.csv', dtype = complex, delimiter = ',')
ang_rad = np.linspace(1,89,89)

ratio = lag/sim

mean = np.mean(ratio[:75])
print(mean)

fig, ax = plt.subplots()
plt.plot(ang_rad, ratio, color = 'k')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('Phase error/Simulation error')
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/error_ratio.png', overwrite = True)

plt.show()
