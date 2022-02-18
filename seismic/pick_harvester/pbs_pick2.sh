#!/bin/bash

# Fei Zhang 2020-05-01

# PBS directives
#PBS -P vy72
#PBS -N pick2_pbs_job
#PBS -q normal
#PBS -l walltime=48:00:00,mem=950GB,ncpus=240,jobfs=160GB
#PBS -l storage=scratch/fxz547+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae


# this is tested working in NCI gadi system.
module purge
module load pbs
module load python3-as-python
module load openmpi/3.1.4/
module load hdf5/1.10.5p

#export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
#export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PYTHONPATH=/home/547/fxz547/github/hiperseis/:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


# the virtual env  with parallel h5py together with many other python libs:
# pip install ordered_set tqdm  click  ujson  psutil
# pip install pykml  pyyaml  joblib pywavelets
source /g/data/ha3/Passive/Software/VENV/para_h5py/bin/activate


# files paths 
ASDF_FILES="/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"    # "/g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt"
PICK_OUTPUT="/g/data/ha3/GASeisDataArchive/DevSpace/pick_workflow2"
CATALOGD="/g/data/ha3/Passive/Events/Unified/"  # $PICK_OUTPUT/catalog/

# the cmdline 
mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/pick.py $ASDF_FILES $CATALOGD $PICK_OUTPUT > $PICK_OUTPUT/run_output_0501.txt
# mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/pick.py $ASDF_FILES $CATALOGD $PICK_OUTPUT > $PICK_OUTPUT/run_output.txt 2>&1

