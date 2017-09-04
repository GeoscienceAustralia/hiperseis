#!/usr/bin/env bash

LIB_DIR=$HOME/lib
BIN_DIR=$HOME/bin
ILOC=iLocRelease1.60.tar.gz
ILOC_DIR=$HOME/iLocRelease1.60
SLBM=$ILOC_DIR/SLBM_Root.3.0.5.Linux.tar.gz
SLBM_DIR=$HOME/SLBM_Root.3.0.5.Linux

sudo yum update -y

# install build essential for centos
sudo yum groupinstall 'Development Tools' -y

# Install the lapack libraries
sudo yum install -y blas \
                    lapack \
                    wget \
                    tar

# get iLoc and RSTT
wget http://www.seismology.hu/data/iLoc/$ILOC && \
    tar -xzf $ILOC --directory $HOME && \
    rm $ILOC

mkdir $LIB_DIR $BIN_DIR

# if iloc env vars are still not available
if [ -z ${ILOCROOT+x} ]; then
    cat $HOME/passive-seismic/sc3/sc3_envs.txt >> $HOME/.bashrc;
    source $HOME/.bashrc
fi

echo 'Installing RSTT ...'
# Install RSTT
cp $SLBM $HOME/ && \
    tar -xzf $SLBM --directory $HOME && \
    rm $HOME/SLBM_Root.3.0.5.Linux.tar.gz && \
    cd $SLBM_DIR \
        && make clean_objs \
        && make geotess \
        && make cc \
        && make c

# install mysql, since siescomp3 already comes with mysql, this may not be
# required.
sudo yum install mysql-devel -y

# links
sudo ln -s /usr/lib64/mysql/libmysqlclient.so /usr/lib/libmysqlclient.so
sudo ln -s /usr/lib64/liblapack.so.3.4 /usr/lib/liblapack.so
sudo ln -s /usr/lib64/libblas.so.3.4 /usr/lib/librefblas.so

echo 'Installing iLoc ...'
# install iLoc
cd $ILOC_DIR/src/ \
    && sed -i '120s/\#//g' Makefile \
    && sed -i '121s/\#//g' Makefile \
    && sed -i '122s/\#//g' Makefile \
    && sed -i '123s/\#//g' Makefile \
    && sed -i '124s/\#//g' Makefile \
    && make sc3db


# create a symlink to the executable
sudo ln -s $ILOC_DIR/src/iloc_sc3db /usr/bin/iloc


echo 'Copy mysql conf into home dir'
# copy the mysql conf
if [ -z ${CIRCLECI+x} ];
    then cp $HOME/passive-seismic/iloc_rstt/.my.cnf $HOME/;
    else cp /root/passive-seismic/iloc_rstt/.my.cnf $HOME/
fi

# use iloc with seiscomp3 db:
# echo "event_id update_db=1 do_gridsearch=0 depth=5" | iloc sc3db
