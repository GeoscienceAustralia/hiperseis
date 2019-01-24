#!/usr/bin/env bash

# run the unittest suite within a virtual environment
#

module load python3/3.6.2
module load hdf5/1.8.14p openmpi/1.8
module load geos/3.5.0 netcdf/4.4.1.1

# setup environment
#export PYTHONPATH=$HOME/.local/lib/python3.6/site-packages:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


# User please modify the below INSTALL_DIR, where pstvenv is installed (see install_venv.sh)
INSTALL_DIR=/g/data/ha3/PST2

# ELLIPCORR env variable should point to `passive-seismic/ellip-corr` dir
export ELLIPCORR=$INSTALL_DIR/passive-seismic/ellip-corr # check this dir is correct

# start the virtualenv which was installed through install_venv.sh
source $INSTALL_DIR/pstvenv/bin/activate

# cd into the source code dir
cd $INSTALL_DIR/passive-seismic/
pytest -v tests/