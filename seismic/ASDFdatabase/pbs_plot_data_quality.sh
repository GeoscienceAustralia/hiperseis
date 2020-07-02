#!/bin/bash

# Check waveform data in a Federated ASDF Database 
# Fei Zhang 2020-05-28

#PBS -P vy72
#PBS -N gadi_plot_data_q
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

# The pointer to Fed ASDF
FED_ASDF_DB_TXT=/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt

mpirun -np $PBS_NCPUS python plot_data_quality.py $FED_ASDF_DB_TXT  2018:01:01 2019:01:01 'AU' '*' '*' test2018_plot_data_quality_gadi

