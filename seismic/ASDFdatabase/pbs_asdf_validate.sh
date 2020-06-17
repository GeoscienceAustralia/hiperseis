#!/bin/bash

# Fei Zhang 2020-06-05  PBS run asdf-validate

### PBS directives
#PBS -P vy72
#PBS -N run_asdf_validator 
#PBS -q normal
#PBS -l walltime=48:00:00,mem=32GB,ncpus=1,jobfs=350GB
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

echo "asdf-validate /g/data/ha3/Passive/STRIPED_DATA/TEMP/OA_AUSARRAY_Yr2_S1.h5 " > validation_OA_AUSARRAY_Yr2_S1.h5.log
asdf-validate /g/data/ha3/Passive/STRIPED_DATA/TEMP/OA_AUSARRAY_Yr2_S1.h5 >> validation_OA_AUSARRAY_Yr2_S1.h5.log

#####################################################################################################
# cat /g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt

# permanent
#/g/data/ha3/Passive/STRIPED_DATA/GA_PERM/1995-1999.h5
#/g/data/ha3/Passive/STRIPED_DATA/GA_PERM/2018-2019.h5
# temporary
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/7G(2013-2015).h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/7J_CAPRAL2006-2007.h5
# OA
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/OA_AUSARRAY1_rev1.h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/OA_AUSARRAY_Yr2_S1.h5
# AQT
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/AQT_2015-2018.h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/7Q_SEAL1.h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/7T_SEAL2_2007.h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/7U_SEAL3.h5
#/g/data/ha3/Passive/STRIPED_DATA/TEMP/W1_2012-2017.h5

