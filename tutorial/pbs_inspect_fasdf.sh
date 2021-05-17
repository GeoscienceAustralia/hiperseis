#!/bin/bash

# Fei Zhang 2021-05-17

# PBS directives
#    #PBS -q express 
#PBS -q normal 
#PBS -P vy72
#PBS -N pbs_job
#PBS -l walltime=30:10:00,mem=50GB,ncpus=2,jobfs=160GB
#PBS -l storage=scratch/fxz547+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae

export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

export PYTHONPATH=/g/data/ha3/fxz547/Githubz/hiperseis:$PYTHONPATH

HIPERSEIS_HOME=/g/data/ha3/fxz547/Githubz/hiperseis
# the virtual env  with parallel h5py together with many other python libs:
source $HIPERSEIS_HOME/gadi_env.sh 


# Run the cmdline

cd $HIPERSEIS_HOME/tutorial
python3 inspect_fasdf.py

###################################################################

# mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/pick.py $ASDF_FILES $CATALOGD $PICK_OUTPUT > $PICK_OUTPUT/run_output_0501.txt
# mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/pick.py $ASDF_FILES $CATALOGD $PICK_OUTPUT > $PICK_OUTPUT/run_output.txt 2>&1

