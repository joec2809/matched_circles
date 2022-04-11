from __future__ import division
from audioop import avg
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm, truncnorm
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
diff = 36

for i in range(1):
    x_corr = np.genfromtxt(f'../data/lon_297_lat_0_corr_0.csv', dtype=complex, delimiter = ',')
    x_corr_real = np.zeros(len(x_corr))
    errors = np.genfromtxt(f'../data/lon_297_lat_0_sim_err_100_sims.csv', dtype=complex, delimiter = ',')[:-1]
    errors_real = np.zeros(len(errors))
    #errors = 0

    print(len(x_corr), len(errors))

    x_corr, errors = mc_functions.array_size_match(x_corr, errors)

    pos_points = 0
    neg_points = 0

    for i, point in enumerate(x_corr):
        x_corr_real[i] = float(np.real(point))
        errors_real[i] = float(np.real(errors[i]))
        if np.real(point) >= 0:
            pos_points += 1
        elif np.real(point) < 0:
            neg_points += 1

    mean = np.mean(x_corr_real)
    
    """
    z_score = (mean - 0)/(std_dev/np.sqrt(len(x_corr)))

    norm_dist = norm(0, 1)

    prob = norm_dist.cdf(-z_score)"""

    a, b = -1,1

    truncnorm_dist = truncnorm(a,b)

    dist_mean, var = truncnorm_dist.stats(moments = 'mv')

    z_score = (mean-dist_mean)/np.sqrt(var/len(x_corr))

    prob = truncnorm_dist.cdf(-z_score)

    print(prob)

    #print(prob)

    fig, ax1 = plt.subplots(figsize = (14,8))
    ax1.errorbar(ang_rad, x_corr_real, yerr = errors_real, color = 'k', ecolor = 'k', elinewidth = 0.8, capsize = 3)
    ax1.plot(ang_rad, x_corr_real, color = 'k')
    ax1.axhline(y=0, color = 'k')
    ax1.set_xticks(np.arange(0, 91, 10))
    ax1.set_xlim(0,90)
    ax1.set_xlabel(r'$\alpha/^\circ$')
    ax1.set_ylabel('$S$')
    ax1.annotate('Phase=0$^\circ$', xy = (0.8,0.9), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
    ax1.annotate(f"Number of positive points = {pos_points} \nNumber of negative points = {neg_points}", xy = (0.6,0.2), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
    """ax1.fill_between(ang_rad, 2*errors, -2*errors, color = 'lightgrey')
    ax1.fill_between(ang_rad, errors, -errors, color = 'darkgrey')"""
    plt.tight_layout()

    #plt.show()
    #fig.savefig(f"../figures/lon_297_lat_0.png")
    lon += diff




"""

non_zero_x_corr = np.delete(x_corr, np.where(x_corr == 0))
no_points = len(non_zero_x_corr)

avg_value = np.mean(np.real(non_zero_x_corr))
expected_avg = 0

std_dev = np.std(np.real(non_zero_x_corr))

z_score = (avg_value-expected_avg)/std_dev

p_value = 1 - norm.cdf(abs(z_score))

print(avg_value)

print(p_value)"""