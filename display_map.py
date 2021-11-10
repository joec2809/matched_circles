from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

SMICA_MAP = hp.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048

ang_rad = np.radians(40)

ipix_strip1 = hp.query_strip(NSIDE, ang_rad-(np.pi/360), ang_rad+(np.pi/360))
ipix_strip2 = hp.query_strip(NSIDE, np.pi-ang_rad-(np.pi/360), np.pi-ang_rad+(np.pi/360))

SMICA_MAP[ipix_strip1] = SMICA_MAP.max()
SMICA_MAP[ipix_strip2] = SMICA_MAP.max()

hp.mollview(SMICA_MAP, title = '', cbar = False)

plt.savefig("/opt/local/l4astro/rbbg94/figures/cmb_map.png")

plt.show()
