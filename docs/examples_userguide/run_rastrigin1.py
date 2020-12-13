# get pid
import sys
if len(sys.argv) > 1:
    pid = sys.argv[1]
else:
    estr  = 'This scripts needs a process identifier (pid) as command line'
    estr += ' argument.'
    raise IOError(estr)

import os
import shutil
import subprocess

exe   = 'rastrigin1.py'
pfile = 'params.txt'
ofile = 'out.txt'

# make individual run directory
rundir = 'tmp.'+pid
os.mkdir(rundir)

# copy individual parameter file
os.rename(pfile+'.'+pid, rundir+'/'+pfile)

# run in individual directory
shutil.copyfile(exe, rundir+'/'+exe)
os.chdir(rundir)
err = subprocess.check_output(['python3', exe], stderr=subprocess.STDOUT)

# make output available to exe_wrapper
os.rename(ofile, '../'+ofile+'.'+pid)

# clean up
os.chdir('..')
shutil.rmtree(rundir)
