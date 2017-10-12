#!/usr/bin/env bash

CWD=`pwd`
cd $HOME
git clone https://github.com/GeoscienceAustralia/passive-seismic.git
cd $CWD

# manage env variables
if [ -z ${ILOCROOT+x} ]; then
    cat $HOME/passive-seismic/sc3/sc3_envs.txt >> $HOME/.bashrc;
    source $HOME/.bashrc
    env
fi

sudo yum update -y
sudo yum groupinstall 'Development Tools' -y

sudo yum install -y wget \
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
                    tkinter \
                    tcl-devel \
                    tk-devel \
                    libffi-devel \
                    libxml2 \
                    libxml2-devel \
                    libxslt \
                    libxslt-devel \
                    python-devel \
                    java-1.7.0-openjdk  # required due to cwbquery

# clone passive-seismic in home directory after git install
if [ ! -d "$HOME/passive-seismic" ]; then
	CWD=`pwd`
	cd $HOME
	git clone https://github.com/GeoscienceAustralia/passive-seismic.git
	cd $CWD
fi

sudo yum install -y geos-devel python-pip openjpeg2-devel
sudo pip install -U pip setuptools virtualenv virtualenvwrapper numpy

# install mpi4py
sudo env MPICC=/usr/lib64/openmpi/bin/mpicc pip install mpi4py==2.0.0


echo "Installing parallel hdf5...."
#Build parallel HDF5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz
tar -xzf hdf5-1.8.14.tar.gz
cd hdf5-1.8.14 && \
    CC=/usr/lib64/openmpi/bin/mpicc ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    sudo make install && \
    cd .. && sudo rm -rf hdf5-1.8.14 hdf5-1.8.14.tar.gz

# build parallel h5py
# /usr/include/openmpi-x86_64/mpi.h
echo "Installing parallel h5py...."
git clone https://github.com/h5py/h5py.git && \
    cd h5py && git checkout tags/2.7.0  && \
    sed -i "52i \ \ \ \ COMPILER_SETTINGS['include_dirs'].extend(['/usr/include/openmpi-x86_64'])" setup_build.py && \
    export CC=/usr/lib64/openmpi/bin/mpicc && \
    python setup.py configure --mpi --hdf5=/usr/local/hdf5/ && \
    sudo python setup.py install && \
    cd .. && sudo rm -rf h5py

# Setup virtualenv
echo '' >> $HOME/.bashrc
echo "source /usr/bin/virtualenvwrapper.sh" >> $HOME/.bashrc
source $HOME/.bashrc

echo "Installing passive seismic software....."


if [ -z ${CIRCLECI+x} ];
    then
    cd $HOME/passive-seismic &&
    mkvirtualenv --system-site-packages seismic && \
    # install python packages
    workon seismic && \
    pip install git+https://github.com/obspy/obspy.git && \
    pip install --ignore-installed ipython && \
    export PATH=$HOME/.virtualenvs/seismic/bin/:$PATH && \
    env GEOS_DIR=/usr/include/ pip install --process-dependency-links .[dev] && \
    $HOME/passive-seismic/iloc_rstt/install_iloc_rstt.sh;  # install iloc and rstt
fi
