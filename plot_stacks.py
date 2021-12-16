from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
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

ang_rad = np.linspace(1,89,89)

lon = 0
lat = 0

no_dir = 6

for i in range(2*no_dir):
    x_corr_1 = np.genfromtxt('/home/joe/documents/rbbg94/data/all_ngp_corrs/ngp_corr_0.csv', dtype=complex, delimiter = ',')[:-1]
    x_corr_2 = np.genfromtxt(f'/home/joe/documents/rbbg94/data/plane_corrs_and_errs/lon_{lon}_lat_{lat}_corr_0.csv', dtype=complex, delimiter = ',')

    errors_1 = np.genfromtxt('/home/joe/documents/rbbg94/data/ngp_sim_err_100_sims.csv', dtype=complex, delimiter = ',')
    errors_2 = np.genfromtxt(f'/home/joe/documents/rbbg94/data/plane_corrs_and_errs/lon_{lon}_lat_{lat}_sim_err_100_sims.csv', dtype=complex, delimiter = ',')

    x_corr = mc_functions.stack(x_corr_1, x_corr_2)
    errors = mc_functions.stack_errors(errors_1, errors_2)

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

    fig.savefig(f'/home/joe/documents/rbbg94/figures/ngp_lon_{lon}_lat_{lat}_stack.png')

    if (-1)**i == 1:
        lat += 45
    elif (-1)**i == -1:
        lat -= 45
        lon += 30