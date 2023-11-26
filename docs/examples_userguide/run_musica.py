import os
import shutil
import subprocess
import sys

# get pid
if len(sys.argv) > 1:
    pid = sys.argv[1]
else:
    pid = None

exe    = './musica.exe'
pfiles = ['musica.nml', 'musica_soil.nml', 'fagus_sylvatica.nml']
ofile  = 'musica_out.nc'

# make individual run directory
if pid is None:
    rundir = 'tmp'
else:
    rundir = f'tmp.{pid}'
os.mkdir(rundir)

# copy individual parameter files
for pfile in pfiles:
    if pid is None:
        os.rename(f'{pfile}', f'{rundir}/{pfile}')
    else:
        os.rename(f'{pfile}.{pid}', f'{rundir}/{pfile}')

# run in individual directory
shutil.copyfile(exe, f'{rundir}/{exe}')
os.chdir(rundir)
err = subprocess.check_output([exe], stderr=subprocess.STDOUT)

# make output available to exe_wrapper
if pid is None:
    os.rename(ofile, f'../{ofile}')
else:
    os.rename(ofile, f'../{ofile}.{pid}')

# clean up
os.chdir('..')
shutil.rmtree(rundir)
