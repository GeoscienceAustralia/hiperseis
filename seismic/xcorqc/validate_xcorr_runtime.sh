#!/bin/bash
module purge
module load python/2.7.13
module load openmpi/2.1.1
module load hdf5/1.10.2p
module load mpi4py/3.0.0-py2
module list

export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=../../:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

mkdir -p validation_result
rm -f validation_result/*.nc
rm -f validation_result/*.log
mpirun -np 2 python correlator.py ../../tests/seismic/xcorqc/data/asdf_file_list1.txt ../../tests/seismic/xcorqc/data/asdf_file_list2.txt validation_result 3600 60 --nearest-neighbours 1 --start-time 2011-03-11T00:00:00Z --end-time 2011-03-12T01:00:00Z --read-buffer-size 1  --fmin .01 --fmax 10 --clip-to-2std --one-bit-normalize

