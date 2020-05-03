#!/bin/bash

# Create sqlite database index for a list of federeated ASDF files
# Fei Zhang 2020-05

#PBS -P vy72
#PBS -N para_make_fed_asdf_index
#PBS -q normal
#PBS -l walltime=20:00:00,mem=160GB,ncpus=16,jobfs=250GB
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
cd /home/547/fxz547/github/hiperseis/seismic/ASDFdatabase
mpirun -np  $PBS_NCPUS python3 FederatedASDFDataSet.py /g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt
# FZ2020-05 after installation of parallel h5py, the making of fed-asdf can work with mpi run. 
# However, in this particular case when only 2 asdf files in the list. It appears only reduced walltime by half.
# The expected spead-up was not achieved by $PBS_NCPUS=16.

# Without support for parallel h5py I/O, we got Error:  
# File "/home/547/fxz547/.local/lib/python3.7/site-packages/pyasdf/asdf_data_set.py", line 459, in mpi
#    raise RuntimeError(msg)
# RuntimeError: Running under MPI requires HDF5/h5py to be complied with support for parallel I/O
# mpirun -np $PBS_NCPUS python3 FederatedASDFDataSet.py /g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt

