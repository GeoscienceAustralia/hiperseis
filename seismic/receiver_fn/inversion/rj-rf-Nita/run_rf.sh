#!/bin/bash
# To run:
#   qsub -M my_email_address -N jobname -v INFILE=input_filename,OUT=output_folder ./run_rf.sh
#PBS -P vy72
#PBS -q normalbw
#PBS -l walltime=24:00:00,mem=160GB,ncpus=56,jobfs=500MB
#PBS -l other=hyperthread
#PBS -l wd
#PBS -j oe
#PBS -m bae

module load openmpi/1.10.2-mt
echo $INFILE
echo $OUT
mpirun ./run $INFILE $OUT

