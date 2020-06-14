# Rastrigin function a=10, b=2*pi
import numpy as np
def rastrigin1(x):
    return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))

import sys
if len(sys.argv) > 1:
    pid = sys.argv[1]
else:
    pid = None

# read parameters
from partialwrap import standard_parameter_reader
pfile = 'params.txt'
if pid:
    pfile = pfile + '.' + pid
x = standard_parameter_reader(pfile)

# calc function
y = rastrigin1(x)

# write output file
ofile = 'out.txt'
if pid:
    ofile = ofile + '.' + pid
with open(ofile, 'w') as ff:
    print(y, file=ff)
