# Rastrigin function a=10, b=2*pi
import numpy as np
def rastrigin1(x):
    return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))

# read parameters
from partialwrap import standard_parameter_reader
x = standard_parameter_reader('params.txt')

# calc function
y = rastrigin1(x)

# write output file
with open('out.txt', 'w') as ff:
    print(y, file=ff)
