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
