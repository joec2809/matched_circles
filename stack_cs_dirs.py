from __future__ import division
import numpy as np
from stack import stack

CMB_DIST = 14000
CELL_SIZE = 320

no_plots = 6

for i in range(no_plots):
	x_corr_1 = np.genfromtxt('/opt/local/l4astro/rbbg94/data/cs_dir_'+str(i)+'.csv', dtype = complex, delimiter = ',')
	x_corr_2 = np.genfromtxt('/opt/local/l4astro/rbbg94/data/cs_dir_'+str(i+no_plots)+'.csv', dtype = complex, delimiter = ',')
	x_corr_1 = x_corr_1[:-1]
	x_corr_2 = x_corr_2[:-1]
	errors= 0

	x_corr = stack(x_corr_1, x_corr_2)

	np.savetxt('/opt/local/l4astro/rbbg94/data/cs_dir_stack_'+str(i)+'.csv', x_corr, delimiter = ',')
