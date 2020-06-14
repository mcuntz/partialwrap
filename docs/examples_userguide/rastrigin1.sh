#!/bin/bash

set -e

# get pid
pid=${1}

exe=rastrigin1.py
pfile=params.txt
ofile=out.txt

# make individual run directory
rundir=tmp.${pid}
mkdir ${rundir}

# copy individual parameter file
mv ${pfile}.${pid} ${rundir}/${pfile}

# run in individual directory
cd ${rundir}
ln -s ../${exe}
python3 ${exe}

# individualize output file
mv ${ofile} ../${ofile}.${pid}

# clean up
cd ..
rm -r ${rundir}
