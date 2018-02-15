FROM centos:centos7
MAINTAINER Sudipta Basak <basaks@gmail.com>
ENV HOME /root
RUN yum update -y
RUN yum groupinstall 'Development Tools' -y
RUN yum install -y wget \
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
                   java-1.7.0-openjdk \
                   cmake \
                   libcurl-devel \
                   netcdf-devel \
                   gdal-devel \
                   pcre-devel \
                   fftw3-devel


# manage envs
ENV PERL_LOCAL_LIB_ROOT="$PERL_LOCAL_LIB_ROOT:/root/perl5"
ENV PERL_MB_OPT="--install_base /root/perl5"
ENV PERL_MM_OPT="INSTALL_BASE=/root/perl5"
ENV PERL5LIB="/root/perl5/lib/perl5:$PERL5LIB"
ENV PATH="/root/perl5/bin:$PATH"

# install perl modules
RUN perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit'
RUN cpan App::cpanminus
RUN cpanm version
RUN cpanm PDL


# Install GMP with mapproject from source, yum install does not bring
# mapproject, required by Istvan's scripts
RUN cd $HOME && \
    svn checkout svn://gmtserver.soest.hawaii.edu/gmt/trunk gmt-dev && \
    wget ftp://ftp.star.nesdis.noaa.gov/pub/sod/lsa/gmt/gshhg-gmt-2.3.7.tar.gz && \
    wget ftp://ftp.star.nesdis.noaa.gov/pub/sod/lsa/gmt/dcw-gmt-1.1.3.tar.gz && \
    tar -xzf gshhg-gmt-2.3.7.tar.gz && \
    tar -xzf dcw-gmt-1.1.3.tar.gz && \
    rm -rf $HOME/dcw-gmt-1.1.3.tar.gz $HOME/gshhg-gmt-2.3.7.tar.gz

RUN yum -y install \
        cmake \
        libcurl-devel \
        netcdf-devel \
        gdal-devel \
        pcre-devel \
        fftw3-devel

RUN cd $HOME/gmt-dev && \
    cp cmake/ConfigUserTemplate.cmake cmake/ConfigUser.cmake && \
    echo "set (CMAKE_INSTALL_PREFIX /usr/local/)" >> cmake/ConfigUser.cmake && \
    echo "set (GSHHG_ROOT $HOME/gshhg-gmt-2.3.7)" >> cmake/ConfigUser.cmake && \
    echo "set (DCW_ROOT $HOME/dcw-gmt-1.1.3)" >> cmake/ConfigUser.cmake && \
    mkdir build && cd build && \
    cmake .. && make install

# python stuff
RUN yum install -y geos-devel python-pip openjpeg2-devel
RUN pip install -U pip setuptools virtualenv virtualenvwrapper numpy

# install mpi4py
RUN env MPICC=/usr/lib64/openmpi/bin/mpicc pip install mpi4py==3.0.0


# "Installing parallel hdf5...."
# Build parallel HDF5
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz && \
    tar -xzf hdf5-1.8.14.tar.gz && \
    cd hdf5-1.8.14 && \
    CC=/usr/lib64/openmpi/bin/mpicc ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5 && \
    make && \
    make install && \
    cd .. && rm -rf hdf5-1.8.14 hdf5-1.8.14.tar.gz

# build parallel h5py
# /usr/include/openmpi-x86_64/mpi.h
RUN git clone https://github.com/h5py/h5py.git && \
    cd h5py && git checkout tags/2.7.0  && \
    sed -i "52i \ \ \ \ COMPILER_SETTINGS['include_dirs'].extend(['/usr/include/openmpi-x86_64'])" setup_build.py && \
    export CC=/usr/lib64/openmpi/bin/mpicc && \
    python setup.py configure --mpi --hdf5=/usr/local/hdf5/ && \
    python setup.py install && \
    cd .. && rm -rf h5py
