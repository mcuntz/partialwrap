# File: logli.py
import sys
import numpy as np
from partialwrap import standard_parameter_reader

# log-likelihood
def log_prob(theta):
    return -0.5 * np.sum(theta**2)

# get pid
if len(sys.argv) > 1:
    pid = int(sys.argv[1])
else:
    pid = None

# read parameters
x = standard_parameter_reader('params.txt', pid=pid)

# calc function
y = log_prob(x)

# write output file
if pid is None:
    fname = 'out.txt'
else:
    fname = f'out.txt.{pid}'
with open(fname, 'w') as ff:
    print(y, file=ff)
