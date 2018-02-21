#!/usr/bin/env bash
set -e
set -o pipefail
sudo yum update -y
sudo yum groupinstall 'Development Tools' -y
sudo yum install -y wget \
                   tar \
                   bzip2-devel \
                   m4 \
                   openmpi openmpi-devel \
                   git \
                   vim \
                   screen \
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
                   tkinter \
                   tcl-devel \
                   tk-devel \
                   libffi-devel \
                   libxml2 \
                   libxml2-devel \
                   libxslt \
                   libxslt-devel \
                   python-devel \
                   netcdf-devel.x86_64 \
                   cpan \
                   evince \
                   evince-nautilus \
                   java-1.7.0-openjdk


# manage envs
cat $HOME/passive-seismic/sc3/perl_envs.txt >> $HOME/.bashrc;
source $HOME/.bashrc

# install perl modules
perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit' &&\
    cpan App::cpanminus && \
    cpanm version && \
    cpanm PDL && \
    cpanm DBI && \
    rm -rf $HOME/.cpanm

# Install GMP with mapproject from source, yum install does not bring
# mapproject, required by Istvan's scripts
cd $HOME && \
    svn checkout svn://gmtserver.soest.hawaii.edu/gmt/trunk gmt-dev && \
    wget ftp://ftp.star.nesdis.noaa.gov/pub/sod/lsa/gmt/gshhg-gmt-2.3.7.tar.gz && \
    wget ftp://ftp.star.nesdis.noaa.gov/pub/sod/lsa/gmt/dcw-gmt-1.1.3.tar.gz && \
    tar -xzf gshhg-gmt-2.3.7.tar.gz && \
    tar -xzf dcw-gmt-1.1.3.tar.gz && \
    rm -rf $HOME/dcw-gmt-1.1.3.tar.gz $HOME/gshhg-gmt-2.3.7.tar.gz

sudo yum -y install cmake \
                   libcurl-devel \
                   netcdf-devel \
                   gdal-devel \
                   pcre-devel \
                   fftw3-devel

cd $HOME/gmt-dev && \
    cp cmake/ConfigUserTemplate.cmake cmake/ConfigUser.cmake && \
    echo "set (CMAKE_INSTALL_PREFIX /usr/local/)" >> cmake/ConfigUser.cmake && \
    echo "set (GSHHG_ROOT $HOME/gshhg-gmt-2.3.7)" >> cmake/ConfigUser.cmake && \
    echo "set (DCW_ROOT $HOME/dcw-gmt-1.1.3)" >> cmake/ConfigUser.cmake && \
    mkdir build && cd build && \
    cmake .. && sudo make install && \
    rm -rf $HOME/gmt-dev


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
