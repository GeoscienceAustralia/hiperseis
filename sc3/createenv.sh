#!/usr/bin/env bash

sudo yum install -y development \
                    wget \
                    tar \
                    bzip2-devel \
                    m4 \
                    openmpi openmpi-devel \
                    git \
                    vim \
                    libpng \
                    libpng-devel \
                    epel-release \
                    blas blas-devel \
                    lapack lapack-devel \
                    libtiff-devel \
                    libjpeg-devel \
                    zlib-devel \
                    freetype-devel \
                    lcms2-devel \
                    libwebp-devel \
                    openjpeg2-devel \
                    tkinter \
                    tcl-devel \
                    tk-devel \
                    libffi-devel \
                    libxml2 \
                    libxml2-devel \
                    libxslt \
                    libxslt-devel

sudo yum install -y python-pip
sudo pip install -U pip virtualenv virtualenvwrapper numpy

# required due to h5py install
sudo pip install mpi4py==2.0.0

#Build parallel HDF5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz
tar -xzvf hdf5-1.8.14.tar.gz
cd hdf5-1.8.14 && \
    ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    sudo make install && \
    cd .. && sudo rm -rf /hdf5-1.8.14 /hdf5-1.8.14.tar.gz2

# build parallel h5py
# /usr/include/openmpi-x86_64/mpi.h
git clone https://github.com/h5py/h5py.git && \
    cd h5py && git checkout tags/2.7.0  && \
    sed -i "52i \ \ \ \ COMPILER_SETTINGS['include_dirs'].extend(['/usr/include/openmpi-x86_64'])" setup_build.py && \
    export CC=/usr/lib64/openmpi/bin/mpicc && \
    python setup.py configure --mpi --hdf5=/usr/local/hdf5/ && \
    sudo python setup.py install && \
    cd .. && sudo rm -rf h5py

# Setup virtualenv
echo "source /usr/bin/virtualenvwrapper.sh" >> $HOME/.bashrc
source $HOME/.bashrc

mkvirtualenv --system-site-packages seismic

workon seismic && \
    git clone https://github.com/GeoscienceAustralia/passive-seismic && \
    cd passive-seismic && \
    python setup.py develop
