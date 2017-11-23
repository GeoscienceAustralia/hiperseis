#!/usr/bin/env bash
#PBS -P vy72
#PBS -q normalbw
#PBS -l walltime=2:00:00,mem=256GB,ncpus=56,jobfs=100GB
#PBS -l wd
#PBS -j oe

module load python3/3.4.3 python3/3.4.3-matplotlib
module load hdf5/1.8.14p openmpi/1.8
module load geos/3.5.0 netcdf/4.4.1.1

# setup environment
export PATH=/g/data/ha3/sudipta/venvs/seismic/bin:$HOME/.local/bin:$PATH
export PYTHONPATH=$HOME/.local/lib/python3.4/site-packages:$PYTHONPATH
export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.4.3/bin/python3
export LANG=en_AU.UTF-8
export WORKON_HOME=/g/data/ha3/sudipta/venvs/
source $HOME/.local/bin/virtualenvwrapper.sh

# start the virtualenv
workon seismic

# gather
mpirun --mca mpi_warn_on_fork 0 cluster gather /g/data/ha3/sudipta/event_xmls


# sort
cluster sort outfile_P.csv -s sorted_P.csv
cluster sort outfile_S.csv -s sorted_S.csv


# match
cluster match sorted_P.csv sorted_S.csv
