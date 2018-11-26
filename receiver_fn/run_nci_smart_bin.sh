#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -l walltime=30:00:00,mem=257GB
#PBS -l ncpus=7
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

python rf_smart_bin.py > rf_smart_bin.log
