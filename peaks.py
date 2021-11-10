from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

CMB_DIST = 14000 #Mpc

ang_rad = np.linspace(1,89,89)
x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/all_ngp_corrs/ngp_corr.csv', dtype=complex, delimiter = ',')[:-1]
errors = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_sim_err_100_sims.csv', dtype=complex, delimiter = ',')

CMB_DIST = 14000 #Mpc

peaks_ind = find_peaks(x_corr, height=0.12)[0][-3:]

peaks_ang_rad = ang_rad[peaks_ind]*(2*np.pi)/360
print(peaks_ang_rad*(360/(2*np.pi)))

domain_seps = 2*CMB_DIST*np.cos(peaks_ang_rad)

print(domain_seps)

guess = 2850 #Mpc
gap = 50 #Mpc

for i in range(5):

	peaks_dist_low = np.array((7*guess, 6*guess, 5*guess))

	low_per_diff = np.mean((np.abs(peaks_dist_low-domain_seps)/domain_seps)*100)

	peaks_dist_high = np.array((7*(guess+gap), 6*(guess+gap), 5*(guess+gap)))

	high_per_diff = np.mean((np.abs(peaks_dist_high-domain_seps)/domain_seps)*100)

	if high_per_diff < low_per_diff:
		guess = (guess+gap)
		peaks_dist = peaks_dist_high
	elif high_per_diff > low_per_diff:
		peaks_dist = peaks_dist_low
		gap = gap/2


print(guess, gap, low_per_diff, high_per_diff)


data = np.genfromtxt('/opt/local/l4astro/rbbg94/data/all_ngp_corrs/ngp_corr.csv', dtype = complex, delimiter = ',')

print(len(ang_rad), len(data))
peaks_ind= np.argwhere(np.logical_and(ang_rad>=40, ang_rad<=65))
peaks_rad = ang_rad[peaks_ind]
peaks_data = data[peaks_ind]
peaks_err = 0

print(peaks_rad, peaks_data)

fig1, ax = plt.subplots(figsize=(14, 8))

ax.plot(peaks_rad, peaks_data, color = 'k')
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(40,65)
for i in range(90):
	ax.axvline((np.arccos((i+1)*320/(2*CMB_DIST)))*(360/(2*np.pi)), color = 'red', ls = '--')
ax.annotate('Phase=0$^\circ$', xy = (0.01,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black'))
ax.annotate('Cell Size = '+str(int(320))+' Mpc', xy = (0.75,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black'))
plt.tight_layout()

fig1.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_ngp_no_err_peaks_320.png', overwrite = True)

fig2, ax = plt.subplots(figsize=(14, 8))

ax.plot(ang_rad, data, color = 'k')
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
for i in range(int(2*CMB_DIST/guess)):
	ax.axvline((np.arccos((i+1)*guess/(2*CMB_DIST)))*(360/(2*np.pi)), color = 'red', ls = '--')
ax.annotate('Phase=0$^\circ$', xy = (0.03,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black'))
ax.annotate('Cell Size = '+str(int(guess))+' Mpc', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black'))
plt.tight_layout()

fig2.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_ngp_no_err_new_size.png', overwrite = True)

plt.show()
