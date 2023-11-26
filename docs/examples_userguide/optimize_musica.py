import sys
from functools import partial
import numpy as np
import scipy.optimize as opt
import netCDF4 as nc
from partialwrap import exe_wrapper, sub_params_ja

#
# Functions
#


# Read NEE from MuSICA's netCDF output
def read_model_nee(ofile):
    with nc.Dataset(ofile, 'r') as fi:
        nee = fi.variables['nee'][:].squeeze()
    return nee


# The objective: RMSE(obs, mod)
def rmse(obs, mod):
    return np.sqrt(((obs - mod)**2).mean())


# RMSE given the output file and the observations
def calc_rmse(ofile, obs):
    mod = read_model_nee(ofile)
    return rmse(obs, mod)


# Read NEE observations
# Assume a csv file such as submitted to europe-fluxdata.org
#   TIMESTAMP_END,H_1_1_1,LE_1_1_1,FC_1_1_1,FC_SSITC_TEST_1_1_1,TA_1_1_1,...
#   201601010030,-20.4583,-1.8627,1.9019,0,5.5533,...
#   ...
# Assume that timestamps are the same as MuSICA output file
def read_obs_nee(ofile):
    with open(ofile, 'r') as fi:
        head = fi.readline().strip().split(',')
    iivar = head.index('FC_1_1_1')
    iiflag = head.index('FC_SSITC_TEST_1_1_1')
    dat = np.loadtxt(ofile, delimiter=',', skiprows=1)
    nee = np.ma.array(dat[:, iivar], mask=(dat[:, iiflag] > 0))
    return nee


# RMSE is around 1-10 (umol m-2 s-1).
# Use a large random number in case of model error because of
# odd parameter combinations.
def err(x):
    return (1. + np.random.random()) * 1000.


#
# Setup
#

# namelist files
nfiles = ['musica.nml', 'musica_soil.nml', 'fagus_sylvatica.nml']
# parameter names (not used, only for info)
names = ['GS_SLOPE', 'GS_HX_HALF', 'GS_HX_SHAPE']
# lower bounds
lb = [1.0, -4.5, 2.0]
# upper bounds
ub = [13.0, -1.0, 20.]

# observations
obsfile = 'FR-Hes_europe-fluxdata_2017.txt'
obs = read_obs_nee(obsfile)

# Use a Python wrapper for MuSICA that deals with pid to run musica.exe
exe             = ['./run_musica.py']
parameterfile   = nfiles
parameterwriter = sub_params_ja
outputfile   = 'musica_out.nc'
outputreader = calc_rmse
oargs        = [obs]  # outputreader arguments in addition to outputfile
# use exe_wrapper for run_musica.py
wrap = partial(exe_wrapper, exe,
               parameterfile, parameterwriter,
               outputfile, outputreader,
               {'oargs': oargs,
                'pid': True, 'error': err})

#
# Optimize
#

print('Start optimization')
ncpu = 4
bounds = list(zip(lb, ub))
res = opt.differential_evolution(wrap, bounds, workers=ncpu)
print('Best parameters:', res.x, ' with objective:', res.fun)

# write parameter files with optimized parameters and suffix .opt
print('Write parameter files with optimized parameters')
sub_params_ja('fagus_sylvatica.nml', res.x, pid='opt')
