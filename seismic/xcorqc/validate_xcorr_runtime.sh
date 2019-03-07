#!/bin/bash
module purge
module load python/2.7.13
module load openmpi/2.1.1
module load hdf5/1.10.2p
module load mpi4py/3.0.0-py2

export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=../../:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

mkdir -p validation_result
mpirun -np 2 python correlator.py test/test_data_ARMA.h5 test/test_data_CMSA.h5 validation_result 3600 60 --nearest-neighbours 1 --start-time 2011-03-11T00:00:00Z --end-time 2011-03-11T01:00:00Z --read-buffer-size 1 --ds1-dec-factor 1 --ds2-dec-factor 1 --fmin .01 --fmax 10 --clip-to-2std True --one-bit-normalize True --ignore-ds1-json-db --ignore-ds2-json-db
