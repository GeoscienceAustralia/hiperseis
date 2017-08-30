FROM ubuntu:16.04

MAINTAINER Sudipta Basak <basaks@gmail.com>

ENV HOME /root

RUN apt-get update
RUN apt-get -y install build-essential \
                       wget \
                       bzip2 \
                       tar \
                       m4 \
                       libopenmpi-dev

#Build parallel HDF5
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && \
    tar -xzvf hdf5-1.8.14.tar.gz && \
    cd hdf5-1.8.14 && \
    ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    make install && \
    cd .. && \
    rm -rf /hdf5-1.8.14 /hdf5-1.8.14.tar.gz

