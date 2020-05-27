#!/bin/bash

# Convert an input ASDF file into multiple mseed waveform data file (Preparation for Local Earthquake picking workflow)
# Fei Zhang 2020-05-21

#PBS -P vy72
#PBS -N gadi_asdf2mseed
#PBS -q normal
#PBS -l walltime=10:00:00,mem=384GB,ncpus=96,jobfs=100GB
#PBS -l storage=scratch/fxz547+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae


# Environment tested working in NCI gadi system.
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
cd /home/547/fxz547/github/hiperseis/seismic/ASDFdatabase

# In preparation for Local Earthquake Picking, convert ASDF into mseed files
mpirun -np $PBS_NCPUS python asdf2mseed.py /g/data/ha3/Passive/STRIPED_DATA/TEMP/OA_AUSARRAY_Yr2_S1.h5 /g/data/ha3/Passive/SHARED_DATA/Scratch/local_picking_workflow/Miniseed/ --start-date 2014-09-01T00:00:00 --end-date 2021-12-01T00:00:00 


## mpirun -np  $PBS_NCPUS python3 FederatedASDFDataSet.py /g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt
