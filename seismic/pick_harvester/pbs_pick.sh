#!/bin/bash

## Fei Zhang 2020-04-20

#PBS -P vy72
#PBS -N pick_pbs_job
#PBS -q normal
#PBS -l walltime=24:00:00,mem=190GB,ncpus=48,jobfs=160GB
#PBS -l storage=scratch/fxz547+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae

# this is tested working in NCI gadi system.
# Alternative is to use anaconda python3 env "conda activate hiperseispy37" (vdi)
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

ASDF_FILES="/g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt"
CATALOGD="/g/data/ha3/Passive/Events/Unified/"
PICK_OUTPUT="/g/data/ha3/GASeisDataArchive/DevSpace/pick_workflow"

cd $PICK_OUTPUT

mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/pick.py $ASDF_FILES $CATALOGD $PICK_OUTPUT > $PICK_OUTPUT/run_output.txt 2>&1

