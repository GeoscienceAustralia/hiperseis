#!/bin/bash
#PBS -P vy72
#PBS -N oa_mseed
#PBS -q normal
#PBS -l walltime=24:00:00,mem=570GB,ncpus=144,jobfs=160GB
#PBS -l storage=scratch/rxh562+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M rakib.hassan@ga.gov.au
#PBS -m bae

module purge
module load pbs
module load python3-as-python
module load openmpi/2.1.6-mt
module load hdf5/1.10.5p

export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=/g/data/ha3/rakib/seismic/pst/hiperseis/:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

mpirun -np $PBS_NCPUS python /g/data/ha3/rakib/seismic/pst/hiperseis/seismic/ASDFdatabase/asdf2event_mseed.py /g/data/ha3/rakib/seismic/data/db/asdf_files.txt /g/data/ha3/Passive/tmp/OA_event_mseeds/events.xml OA /g/data/ha3/Passive/tmp/OA_event_mseeds/mseeds

