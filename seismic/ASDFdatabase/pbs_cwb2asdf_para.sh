#!/bin/bash

# Fei Zhang 2020-05-03
# collate many mseed files into a single ASDF file by mpi parallel writing

### PBS directives
#PBS -P vy72
#PBS -N cwb2asdf_para
#PBS -q normal
#PBS -l walltime=10:00:00,mem=96GB,ncpus=16,jobfs=250GB
#PBS -l storage=scratch/fxz547+gdata/ha3

#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae


# this is tested working in NCI gadi system.
module purge
module load python3/3.7.4
module loat hdf5/1.10.5p 
module load openmpi/3.1.4
#export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
#export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PYTHONPATH=/home/547/fxz547/github/hiperseis/:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


# the virtual env  with parallel h5py together with many other python libs:
# pip install ordered_set tqdm  click  ujson  psutil
# pip install pykml  pyyaml  joblib pywavelets
source /g/data/ha3/Passive/Software/VENV/para_h5py/bin/activate

# the cmdline 

cd /home/547/fxz547/github/hiperseis/seismic/ASDFdatabase/cwb2asdf

# In order to use more than 1-processor h5py pyasdf must be parallel-IO support.
# with parallel pyasdf installed in venv, test run shown that the writing is not scaled-up. Less than 2 CPU is used even though 16 are requested!!
mpirun -np  $PBS_NCPUS python3 cwb2asdf.py /g/data/ha3/GASeisDataArchive/2019_demultiplex /g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml /g/data/ha3/GASeisDataArchive/2019.h5 > /g/data/ha3/GASeisDataArchive/2019_cwb2asdf_run.out

