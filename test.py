from click import argument
import numpy as np
import mc_functions

def func(*args):
    arguments = np.array(args, dtype = object)
    return arguments

arrays = func(np.arange(1,10,1), np.arange(5,10,1))

lens = []

for i, array in enumerate(arrays):
    lens.append(len(array))

min_len = min(lens)

out = []

for j, array in enumerate(arrays):
    out.append(array[:min_len])

final_out = np.array(out)

print(arrays, final_out)