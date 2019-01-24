#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -q express
#PBS -l walltime=3:00:00,mem=128GB
#PBS -l ncpus=32
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

python generate_rf_parallel.py > generate_rf_parallel.log
