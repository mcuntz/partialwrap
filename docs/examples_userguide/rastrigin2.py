# Rastrigin function a=10, b=2*pi
import numpy as np
def rastrigin1(x):
    return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))


# read parameters
with open('params.txt', 'r') as fi:
    pdict = {}
    for line in fi:
        ll = line.split()
        if (len(ll) == 0) or ll[0].startswith('#'):
            continue
        pdict[ll[0]] = float(ll[2])
x = np.array([ pdict[kk] for kk in sorted(pdict.keys()) ])

# calc function
y = rastrigin1(x)

# write output file
with open('out.txt', 'w') as ff:
    print(y, file=ff)
