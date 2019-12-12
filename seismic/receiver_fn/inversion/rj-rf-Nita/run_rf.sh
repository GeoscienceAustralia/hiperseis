#!/bin/bash
# To run:
#   qsub -M my_email_address -N jobname -v INFILE=input_filename,OUT=output_folder ./run_rf.sh
#PBS -P vy72
#PBS -q normalbw
#PBS -l walltime=48:00:00,mem=16GB,ncpus=28,jobfs=200MB
#PBS -l other=hyperthread
#PBS -l wd
#PBS -j oe
#PBS -m bae

module purge
module load openmpi/1.10.2
echo $INFILE
echo $OUT
# Map by node so that if using OpenMP threads, processes will be assigned
# to node in round-robin fashion.
mpirun --report-bindings --map-by node ./run $INFILE $OUT
