from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import random
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip_auto import strip_auto
from load_file import load_file
from match_circle_s import match_circle_s


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/sims/dx12_v3_smica_cmb_mc_00000_raw.fits")

NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320
bins = 360

strip_finder(cmb_map, 88.6*(2*np.pi/360), NSIDE)

circle_a = load_file('strip_a', bins)
circle_b = load_file('strip_b', bins)

x_corr = np.zeros((bins), dtype = complex)
for i in range(bins):
	x_corr[i] = match_circle_s(circle_a, circle_b, i*(2*np.pi/bins), bins, 720)

fig, ax = plt.subplots()

ax.plot(np.arange(0,360,1), x_corr)
ax.set_title('Auto-Correlation: Simulation')
ax.set_xlabel('Phase / $\circ$')
ax.set_ylabel('X-Correlation')
ax.annotate(r'$\alpha$ = 88.6$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
fig.savefig("/opt/local/l4astro/rbbg94/figures/auto_corr_sim.png", overwrite = True)

plt.show()
