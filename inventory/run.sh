#!/bin/bash
## To run, enter 'qsub run.sh'
#PBS -P ha3
#PBS -q normal
#PBS -l mem=8gb
#PBS -l jobfs=2gb
#PBS -l ncpus=1
#PBS -l software=python
## Tell it gdata1a filesystem needs to be available
#PBS -lother=gdata1a
## The job will be executed from current working directory instead of home.
#PBS -l wd
#PBS -N StationMetadataCleanup
python ./engd2stxml.py
