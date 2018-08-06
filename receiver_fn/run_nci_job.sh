#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -l walltime=15:00:00,mem=257GB
#PBS -l ncpus=7
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

python generate.py > output.txt

