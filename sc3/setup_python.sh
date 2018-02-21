#!/usr/bin/env bash
set -e
set -o pipefail

# python stuff
sudo yum install -y geos-devel python-pip openjpeg2-devel
sudo pip install -U pip setuptools virtualenv virtualenvwrapper numpy

# install mpi4py
sudo env MPICC=/usr/lib64/openmpi/bin/mpicc pip install mpi4py==3.0.0


# "Installing parallel hdf5...."
# Build parallel HDF5
cd $HOME && \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && \
    tar -xzf hdf5-1.8.14.tar.gz && \
    cd hdf5-1.8.14 && \
    CC=/usr/lib64/openmpi/bin/mpicc ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    sudo make install && \
    cd $HOME && rm -rf hdf5-1.8.14 hdf5-1.8.14.tar.gz

# build parallel h5py
# /usr/include/openmpi-x86_64/mpi.h
cd $HOME && \
    git clone https://github.com/h5py/h5py.git && \
    cd h5py && git checkout tags/2.7.0  && \
    sed -i "52i \ \ \ \ COMPILER_SETTINGS['include_dirs'].extend(['/usr/include/openmpi-x86_64'])" setup_build.py && \
    export CC=/usr/lib64/openmpi/bin/mpicc && \
    python setup.py configure --mpi --hdf5=/usr/local/hdf5/ && \
    sudo python setup.py install && \
    cd $HOME && sudo rm -rf h5py

