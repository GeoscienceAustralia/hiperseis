# get iLoc and RSTT
wget http://www.seismology.hu/data/iLoc/iLocRelease1.60.tar.gz
tar -xzf iLocRelease1.60.tar.gz
mkdir lib bin

# export env variables
export ILOCROOT=~/iLocRelease1.60
export SLBMROOT=~/SLBM_Root.3.0.5.Linux
export LD_LIBRARY_PATH=~/lib:$SLBMROOT/lib:/opt/seiscomp3/lib:
export PATH=$PATH:~/bin:$SLBMROOT/bin

# install build essential
sudo yum groupinstall 'Development Tools' -y

# Install the lapack libraries
sudo yum install blas lapack -y


# Install RSTT
cp iLocRelease1.60/SLBM_Root.3.0.5.Linux.tar.gz .
tar -xzf SLBM_Root.3.0.5.Linux.tar.gz
cd SLBM_Root.3.0.5.Linux/
make clean_objs
make geotess
make cc
make c
cd


# install mysql, since siescomp3 already comes with mysql, this may not be
# required.
sudo yum install mysql-devel -y

#
sudo ln -s /usr/lib64/mysql/libmysqlclient.so /usr/lib/libmysqlclient.so
sudo ln -s /usr/lib64/liblapack.so.3.4 /usr/lib/liblapack.so
sudo ln -s /usr/lib64/libblas.so.3.4 /usr/lib/librefblas.so


# install iLoc
cd iLocRelease1.60/src/

# edit the Makefile in iLocRelease1.60/src/
# uncomment variables under RedHat:

# RedHat:
       IPGSQL = -I/usr/include/postgresql
       IMYSQL = -I/usr/include/mysql
       PGSQL = -lpq
       MYSQL = -lmysqlclient
       LAPACK = -llapack -lrefblas

make sc3db

cd

# create a symlink to the executable
sudo ln -s ~/iLocRelease1.60/src/iloc_sc3db /usr/bin/iloc

# use iloc with seiscomp3 db:
echo "event_id update_db=1 do_gridsearch=0 depth=5" | iloc sc3db


# May need to add this as a ~/.my.cnf as mysql config

#Sample MySQL config file for medium systems.
#
# This is for a system with little memory (32M - 64M) where MySQL plays
# an important part, or systems up to 128M where MySQL is used together with
# other programs (such as a web server)
#
# You can copy this file to
# /etc/my.cnf to set global options,
# mysql-data-dir/my.cnf to set server-specific options (in this
# installation this directory is /var/lib/mysql) or
# ~/.my.cnf to set user-specific options.
#
# In this file, you can use all long options that a program supports.
# If you want to know which options a program supports, run the program
# with the "--help" option.

# The following options will be passed to all MySQL clients
[client]
host=localhost
database=seiscomp3
user=sysop
password=sysop


