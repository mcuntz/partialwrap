# File: logli1.py
import numpy as np

# get pid
import sys
if len(sys.argv) > 1:
    pid = int(sys.argv[1])
else:
    pid = None

# log-likelihood
def log_prob(theta):
    return -0.5 * np.sum(theta**2)

# read parameters
from partialwrap import standard_parameter_reader
x = standard_parameter_reader('params.txt', pid=pid)

# calc function
y = log_prob(x)

# write output file
if pid:
    fname = 'out.txt'+'.'+str(pid)
else:
    fname = 'out.txt'
with open(fname, 'w') as ff:
    print(y, file=ff)
