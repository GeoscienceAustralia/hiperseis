#!/usr/bin/env bash

# Install a python3 virtual enviroment for passive-seismic  (travel time, cluster to run)
# These instructions currently only work with gcc and not Intel compilers. 
# The script is tested working in raijin, but had issues in VDI.


# User can modify the below INSTALL_DIR=/g/data/ha3/PST2, where everything will be installed.
INSTALL_DIR=/g/data/ha3/PST2
VENV_NAME=pstvenv

[ -d $INSTALL_DIR ] ||   mkdir $INSTALL_DIR

# Clone the passive-seismic repository into install_dir directory of your choice:
cd $INSTALL_DIR
git clone https://github.com/GeoscienceAustralia/passive-seismic.git

module load git/2.9.5
cd passive-seismic
/apps/git/2.9.5/bin/git submodule update --init --recursive --remote


#Unload the icc compiler and default openmpi from the terminal:
module unload intel-cc intel-fc openmpi
#Load the modules required for installation and running:

module load python3/3.6.2
# module load python3  # default
module load hdf5/1.8.14p openmpi/1.8
module load geos/3.5.0 netcdf/4.4.1.1
#(Alternatively, you may wish to add the above lines to your ~/.profile file)

#Optional Install virtualenv and virtualenvwrapper
#pip3 install  --user virtualenv virtualenvwrapper
#which python3 /apps/python3/3.6.2/bin/python3

#Create a new virtualenv for passive-seismic:
mkdir $INSTALL_DIR/$VENV_NAME
python3 -m venv $INSTALL_DIR/$VENV_NAME

# Make sure the virtualenv is activated: workon seismic

source $INSTALL_DIR/$VENV_NAME/bin/activate

# Install mpi4py as required by h5py in the next step.

pip3 install mpi4py==3.0.0 --no-binary :all:

# Clone h5py from https://github.com/basaks/h5py.git:

cd $INSTALL_DIR
git clone https://github.com/basaks/h5py.git
cd ./h5py
export CC=mpicc
python setup.py configure --mpi --hdf5=/apps/hdf5/1.8.14p
python setup.py install


# Install passive-seismic:

cd $INSTALL_DIR/passive-seismic
export GEOS_DIR=$GEOS_BASE
pip3 install --process-dependency-links -e .[dev] --no-binary :all:

# alternative install commands tested in ubnuntu python2
#sudo pip2 install --process-dependency-links -e .[dev] --no-binary :all:

# additional stuffs:
pip3 install flake8

# Once installation has completed, you can run the tests to verify everything has gone correctly:

# pip3 install pytest

# pytest tests/test_cluster.py
# pytest tests/
