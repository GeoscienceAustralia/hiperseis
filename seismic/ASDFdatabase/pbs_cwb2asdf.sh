#!/bin/bash

#PBS -P vy72
#PBS -N cwb2asdf
#PBS -q normal
#PBS -l walltime=40:00:00,mem=32GB,ncpus=2,jobfs=250GB
#PBS -l storage=scratch/fxz547+gdata/ha3

#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae

# this is tested working in NCI gadi system.
# Alternative is to use anaconda python3 env "conda activate hiperseispy37" (vdi)
module purge
module load pbs
module load python3-as-python
module load openmpi/3.1.4/
module load hdf5/1.10.5p

#export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
#export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PYTHONPATH=/home/547/fxz547/github/hiperseis/:$PYTHONPATH
export PATH=$HOME/.local/bin:$PATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


cd /home/547/fxz547/github/hiperseis/seismic/ASDFdatabase/cwb2asdf

# In order to use more than 1-processor h5py pyasdf must be parallel-IO support.
mpirun -np 2 python3 cwb2asdf.py /g/data/ha3/GASeisDataArchive/2019_demultiplex /g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml /g/data/ha3/GASeisDataArchive/2019.h5 > /g/data/ha3/GASeisDataArchive/2019_cwb2asdf_run.out

