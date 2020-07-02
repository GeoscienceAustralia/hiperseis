#!/bin/bash
#PBS -P vy72
#PBS -N cwb2mseed
#PBS -q hugemem
#PBS -l walltime=10:00:00,mem=2950GB,ncpus=96,jobfs=1024GB
#PBS -l storage=gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae

module purge
module load python3/3.7.4
module load openmpi/2.1.6-mt
module load hdf5/1.10.5p
#export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
#export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PYTHONPATH=/home/547/fxz547/github/hiperseis/:$PYTHONPATH
export PATH=$HOME/.local/bin:$PATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

cd /home/547/fxz547/github/hiperseis/seismic/ASDFdatabase/cwb2asdf
mpirun -np 48 --map-by ppr:24:node  python3 demultiplex.py /g/data/ha3/GASeisDataArchive/2019 /g/data/ha3/GASeisDataArchive/2019_demultiplex

