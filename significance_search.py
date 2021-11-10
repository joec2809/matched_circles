from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

ang_rad = np.linspace(1,89,89)

sig_sets_corr = np.zeros((360,89), dtype = complex)
sig_sets_err = np.zeros((360,89), dtype = complex)
a = 0

for i in range(360):
	x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/all_ngp_corrs/ngp_corr_'+str(i)+'.csv', dtype=complex, delimiter = ',')[:-1]
	#errors = np.genfromtxt('/opt/local/l4astro/rbbg94/data/all_ngp_lag_errs/ngp_lag_err_'+str(i)+'.csv', dtype=complex, delimiter = ',')[:-1]
	errors = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_sim_err_100.csv', dtype=complex, delimiter = ',')
	for j in range(len(x_corr)):
		if x_corr[j] > 2*errors[j]:
			sig_sets_corr[i,:] = x_corr
			sig_sets_err[i,:] = errors
			a += 1
			break

phases = np.where(np.any(sig_sets_corr>0,axis=1))[0]
print(phases, len(phases))
sig_sets_corr = sig_sets_corr[~np.all(sig_sets_corr == 0, axis=1)]
sig_sets_err = sig_sets_err[~np.all(sig_sets_err == 0, axis=1)]

for i in range(len(sig_sets_corr)):
	fig, ax = plt.subplots(figsize = (14,10))
	ax.plot(ang_rad, sig_sets_corr[i], color = 'k')
	ax.fill_between(ang_rad, 2*sig_sets_err[i], -2*sig_sets_err[i], color = 'lightgrey')
	ax.fill_between(ang_rad, sig_sets_err[i], -sig_sets_err[i], color = 'darkgrey')
	ax.set_xlabel(r'$\alpha/^\circ$')
	ax.set_ylabel('$S$')
	ax.annotate('Phase='+str(phases[i])+'$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
	ax.axhline(0, color = 'black')
	plt.xticks(np.arange(0, 91, 10))
	plt.xlim(0,90)
	plt.show()
