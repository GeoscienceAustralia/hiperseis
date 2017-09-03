FROM ubuntu:16.04

MAINTAINER Sudipta Basak <basaks@gmail.com>

ENV HOME /root

RUN apt-get update
RUN apt-get -y install build-essential \
                       wget \
                       tar \
                       m4 \
                       libopenmpi-dev \
                       git \
                       vim \
                       openssh-server \
                       openssl \
                       libssl-dev \
                       libblas-dev \
                       liblapack-dev \
                       libpng-dev \
                       python-setuptools \
                       libffi-dev \
                       libxml2-dev \
                       libxslt1-dev \
                       libtiff4-dev \
                       libjpeg8-dev \
                       zlib1g-dev \
                       libfreetype6-dev \
                       liblcms2-dev \
                       libwebp-dev \
                       tcl8.5-dev \
                       tk8.5-dev \
                       python-tk \
                       python-dev \
                       python3-dev \
                       python3-pip

RUN pip install -U pip

# required due to h5py install
RUN pip install mpi4py==2.0.0

# Build parallel HDF5
# choose version 1.8.14 as it matches with NCI version
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && \
    tar -xzvf hdf5-1.8.14.tar.gz && \
    cd hdf5-1.8.14 && \
    ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    make install && \
    cd .. && \
    rm -rf /hdf5-1.8.14 /hdf5-1.8.14.tar.gz

RUN git clone https://github.com/h5py/h5py.git && \
    export CC=mpicc && \
    cd h5py && python setup.py configure --mpi --hdf5=/usr/local/hdf5/ && \
    python setup.py install && \
    cd .. && rm -rf h5py

RUN git clone https://github.com/GeoscienceAustralia/passive-seismic && \
    cd passive-seismic && \
    python setup.py develop


EXPOSE 22

# Start sshd
CMD ["/usr/sbin/sshd", "-D"]
