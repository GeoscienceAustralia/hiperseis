#!/bin/bash

module purge
module load pbs
module load python3-as-python
module load openmpi/3.1.4/
module load hdf5/1.10.5p
module list

export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=../../:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8

mkdir -p validation_result
rm -f validation_result/*.nc
rm -f validation_result/*.log
rm -f validation_result/*.txt

echo "../../tests/test_seismic/xcorqc/data/test_data_ARMA.h5" > validation_result/asdf_file_list1.txt
echo "../../tests/test_seismic/xcorqc/data/test_data_CMSA.h5" > validation_result/asdf_file_list2.txt

mpirun -np 2 python correlator.py validation_result/asdf_file_list1.txt validation_result/asdf_file_list2.txt validation_result 3600 60 0.1 --nearest-neighbours 1 --start-time 2011-03-11T00:00:00Z --end-time 2011-03-12T01:00:00Z --read-buffer-size 1  --fmin .01 --fmax 10 --one-bit-normalize
