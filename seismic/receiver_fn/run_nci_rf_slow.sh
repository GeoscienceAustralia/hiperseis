#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -q express
#PBS -l walltime=2:00:00,mem=128GB
#PBS -l ncpus=16
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

python rf_slow_stack.py 
