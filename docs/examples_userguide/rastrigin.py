import numpy as np
from partialwrap import standard_parameter_reader


# Rastrigin function a=10, b=2*pi
def rastrigin(x):
    return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))


# read parameters
x = standard_parameter_reader('params.txt')

# calc function
y = rastrigin(x)

# write output file
with open('out.txt', 'w') as ff:
    print(y, file=ff)
