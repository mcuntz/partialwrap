# File: run_rastrigin.py
import os
import shutil
import subprocess
import sys

# get pid
if len(sys.argv) > 1:
    pid = sys.argv[1]
else:
    raise IOError('This scripts needs a process identifier (pid) as'
                  ' command line argument.')

exe   = 'rastrigin.py'
pfile = 'params.txt'
ofile = 'out.txt'

# make individual run directory
rundir = f'tmp.{pid}'
os.mkdir(rundir)

# copy individual parameter file
os.rename(f'{pfile}.{pid}', f'{rundir}/{pfile}')

# run in individual directory
shutil.copyfile(exe, f'{rundir}/{exe}')
os.chdir(rundir)
err = subprocess.check_output(['python3', exe],
                              stderr=subprocess.STDOUT)

# make output available to exe_wrapper
os.rename(ofile, f'../{ofile}.{pid}')

# clean up
os.chdir('..')
shutil.rmtree(rundir)
