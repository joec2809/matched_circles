from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
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

for i in range(1):
    x_corr_1 = np.genfromtxt('../data/lon_297_lat_0_corr_0.csv', dtype=complex, delimiter = ',')
    x_corr_2 = np.genfromtxt(f'../data/ngp_corrs/ngp_corr.csv', dtype=complex, delimiter = ',')

    errors_1 = np.genfromtxt('../data/lon_297_lat_0_sim_err_100_sims.csv', dtype=complex, delimiter = ',')
    errors_2 = np.genfromtxt(f'../data/ngp_corrs/ngp_sim_err_100_sims.csv', dtype=complex, delimiter = ',')

    x_corr_1, x_corr_2, errors_1, errors_2 = mc_functions.array_size_match(x_corr_1, x_corr_2, errors_1, errors_2)

    x_corr = mc_functions.stack(x_corr_1, x_corr_2)
    errors = mc_functions.stack_errors(errors_1, errors_2)

    pos_points = 0
    neg_points = 0

    for i, point in enumerate(x_corr):
        if point >= 0:
            pos_points += 1
        elif point < 0:
            neg_points += 1

    std_dev = np.std(x_corr)
    mean = np.real(np.mean(x_corr))
    
    z_score = (mean - 0)/(std_dev)

    norm_dist = norm(0, 1)

    prob = norm_dist.cdf(-z_score)

    print(prob)



    fig, ax = plt.subplots(figsize = (14,8))
    #ax.errorbar(ang_rad, x_corr, yerr = errors, color = 'k', ecolor = 'k', elinewidth = 0.8, capsize = 3)
    ax.plot(ang_rad, x_corr, color = 'k')
    plt.axhline(y=0, color = 'k')
    plt.xticks(np.arange(0, 91, 10))
    plt.xlim(0,90)
    ax.set_xlabel(r'$\alpha/^\circ$')
    ax.set_ylabel('$S$')
    ax.annotate('Phase=0$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
    ax.annotate(f"Number of positive points = {pos_points} \nNumber of negative points = {neg_points}", xy = (0.6,0.2), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
    ax.fill_between(ang_rad, 2*errors, -2*errors, color = 'lightgrey')
    ax.fill_between(ang_rad, errors, -errors, color = 'darkgrey')

    plt.tight_layout()

    plt.show()

    fig.savefig(f'../figures/ngp_lon_297_lat_0_stack.png')