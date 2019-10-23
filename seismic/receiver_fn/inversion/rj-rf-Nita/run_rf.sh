#!/bin/bash
#PBS -P vy72
#PBS -N rfinv
#PBS -q normalbw
#PBS -l walltime=10:00:00,mem=160GB,ncpus=896,jobfs=1GB
#PBS -l other=hyperthread
#PBS -l wd
#PBS -j oe
#PBS -M andrew.medlin@ga.gov.au
#PBS -m bae

module load openmpi/1.10.2-mt
mpirun ./run

