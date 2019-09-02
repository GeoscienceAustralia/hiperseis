#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -l walltime=30:00:00,mem=257GB
#PBS -l ncpus=7
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

python extract_event_traces.py > extract_event_traces.log
