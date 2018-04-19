#!/usr/bin/env bash

# example PBS script provided in handover knowledge transfer

#PBS -P vy72
#PBS -q express
#PBS -l walltime=2:00:00,mem=64GB,ncpus=32,jobfs=100GB
#PBS -l wd
#PBS -j oe

module load python3/3.6.2
module load hdf5/1.8.14p openmpi/1.8
module load geos/3.5.0 netcdf/4.4.1.1

# setup environment
export PATH=/g/data/ha3/sudipta/venvs/seismic/bin:$HOME/.local/bin:$PATH
export PYTHONPATH=$HOME/.local/lib/python3.6/site-packages:$PYTHONPATH
export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.6.2/bin/python3
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8
export WORKON_HOME=/g/data/ha3/sudipta/venvs/
source $HOME/.local/bin/virtualenvwrapper.sh

# start the virtualenv
workon seismic

# ELLIPCORR env variable should point to `passive-seismic/ellip-corr` dir
export ELLIPCORR=$HOME/passive-seismic/ellip-corr


# gather
mpirun --mca mpi_warn_on_fork 0 cluster gather /g/data/ha3/events_xmls_sc3ml


# sort
cluster sort outfile_P.csv 5. -s sorted_P.csv
cluster sort outfile_S.csv 10. -s sorted_S.csv


# match
cluster match sorted_P.csv sorted_S.csv -p matched_P.csv -s matched_S.csv


# zones
cluster zone sorted_S.csv -z '0 -50.0 100 190' -r sorted_region_S.csv -g sorted_global_S.csv
cluster zone sorted_P.csv -z '0 -50.0 100 190' -r sorted_region_P.csv -g sorted_global_P.csv
cluster zone matched_P.csv -z '0 -50.0 100 190' -r region_P.csv -g global_P.csv
cluster zone matched_S.csv -z '0 -50.0 100 190' -r region_S.csv -g global_S.csv
