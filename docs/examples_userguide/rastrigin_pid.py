import sys
import numpy as np
from partialwrap import standard_parameter_reader


# Rastrigin function a=10, b=2*pi
def rastrigin(x):
    return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))


# get pid
if len(sys.argv) > 1:
    pid = int(sys.argv[1])
else:
    pid = None

# read parameters
x = standard_parameter_reader('params.txt', pid=pid)

# calc function
y = rastrigin(x)

# write output file
if pid is None:
    fname = 'out.txt'
else:
    fname = f'out.txt.{pid}'
with open(fname, 'w') as ff:
    print(y, file=ff)
