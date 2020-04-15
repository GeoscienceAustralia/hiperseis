#!/bin/bash

#PBS -P vy72
#PBS -N cwb2asdf
#PBS -q normal
#PBS -l walltime=20:00:00,mem=384GB,ncpus=96,jobfs=250GB
#PBS -l storage=scratch/fxz547+gdata/ha3

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

mpirun -np 96 python3 cwb2asdf.py /g/data/ha3/GASeisDataArchive/2019_demultiplex /g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml /g/data/ha3/GASeisDataArchive/2019.h5 > /g/data/ha3/GASeisDataArchive/2019_cwb2asdf_run.out

#mpirun -np 96 python3 cwb2asdf.py /g/data/ha3/GASeisDataArchive/2020_demultiplexed  /g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml /g/data/ha3/GASeisDataArchive/2020.h5 > /g/data/ha3/GASeisDataArchive/2020_cwb2asdf_run.out

