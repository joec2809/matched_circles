import numpy as np
import mc_functions

array_1 = np.arange(1,10,1)
array_2 = np.arange(0.1,1,0.1)

print(mc_functions.stack_errors(array_1, array_2))