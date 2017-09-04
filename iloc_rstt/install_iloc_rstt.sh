#!/usr/bin/env bash
# get iLoc and RSTT
wget http://www.seismology.hu/data/iLoc/iLocRelease1.60.tar.gz
tar -xzf iLocRelease1.60.tar.gz
mkdir lib bin

# cat passive-seismic/iloc_rstt/iloc_envs.sh >> $HOME/.bashrc
# source $HOME/.bashrc

# install build essential for centos
# sudo yum groupinstall 'Development Tools' -y

# Install the lapack libraries
# sudo yum install blas lapack -y


# Install RSTT
cp $HOME/iLocRelease1.60/SLBM_Root.3.0.5.Linux.tar.gz $HOME/
tar -xzf $HOME/SLBM_Root.3.0.5.Linux.tar.gz
cd SLBM_Root.3.0.5.Linux/ \
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


# install iLoc
cd $HOME/iLocRelease1.60/src/ \
    && sed -i '120s/\#//g' Makefile \
    && sed -i '121s/\#//g' Makefile \
    && sed -i '122s/\#//g' Makefile \
    && sed -i '123s/\#//g' Makefile \
    && sed -i '124s/\#//g' Makefile \
    && make sc3db


# create a symlink to the executable
sudo ln -s $HOME/iLocRelease1.60/src/iloc_sc3db /usr/bin/iloc


# copy the mysql conf
cp $HOME/passive-seismic/iloc_rstt/.my.cnf $HOME/

# use iloc with seiscomp3 db:
# echo "event_id update_db=1 do_gridsearch=0 depth=5" | iloc sc3db
