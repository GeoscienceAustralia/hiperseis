#!/bin/bash
#PBS -l ngpus=2
#PBS -l ncpus=6
#PBS -l mem=8GB
#PBS -l walltime=12:00:00
#PBS -l wd
#PBS -lother=gdata1
#PBS -q gpu
#PBS -P vy72

module purge
module load hdf5/1.8.14p
module load tensorflow/1.8-cudnn7.1-python2.7
module load python/2.7.13

python main.py
