#! /bin/bash

sudo yum update
sudo yum install git

sudo yum install xorg-x11-xauth.x86_64 xorg-x11-server-utils.x86_64 dbus-x11.x86
sudo yum install xeyes

# exit and 
# ssh -X login again to test xeyes


sudo yum install mysql56*
sudo yum install  MySQL-python27.x86_64

sudo chkconfig mysqld on

sudo service mysqld start

# set mysql root passwd
mysqladmin -u root password 'CWBr00tpwd'
# Alternatively run   /usr/libexec/mysql56/mysql_secure_installation

#verify root pass setup, check version
mysqladmin -u root -pCWBr00tpwd version


# create users and set permissions
sudo adduser vdl
# make vdl home group accessible
sudo chmod -R 750 /home/vdl

sudo adduser davidk
sudo usermod -g vdl davidk
sudo adduser sudipta
sudo usermod -g vdl sudipta

sudo usermod -g vdl ec2-user


# become user vdl and get the software from USGS ftp site
sudo su - vdl
mkdir TEMP

# cd into the TEMP dir and do the following
cd TEMP
wget ftp://hazards.cr.usgs.gov/CWBQuery/EdgeCWBRelease.tar.gz

tar -zxvf EdgeCWBRelease.tar.gz
tar -xvf bin_release.tar 
tar -xvf scripts_release.tar 

mv bin ~vdl
cp -p Jars/* ~vdl/bin/
chmod -R 755 ~vdl/bin

# make sure java class lib has jars
ln -s  /home/vdl/bin bin/lib   

cp -r scripts ~vdl/

cp .bash* ~vdl/

cd ~vdl

# the scripts/installCWBRelease.bash has issues. as a reference only
dirs="bin log LOG EDGEMOM config SEEDLINK EW log/MSLOG log/CWBLOG log/q330";
# create any missing directories
for dir in $dirs; do
if [ ! -e ${dir} ]; then
        echo "Create ${dir} directory"
        mkdir ${dir} 
fi
done

tar -xvf TEMP/dbconn.tar 
#.dbconn/dbconn.conf
#.dbconn/dbconn_mysql_localhost_3306.conf
#.dbconn/dbconn_mysql_localhost.conf

# creat databases - require type in mysql user info: root, passwd, y
cd scripts/INSTALL/
bash ./makedb.bash   # if run success, this will create a few db

bash createketchum.bash
# type the mysqlroot passpwd

bash createvdlro.bash
# type the mysqlroot passpwd

# restart the server

sudo shutdown -r now

# after restarting, the mysql server should be up and running.
# ssh -X to enable X-windows

###################################################################
# More Configurations:

# In the /etc/my.cnf file always add to the [mysqld] section
# default-time-zone='+00:00'

# database ops users setup by X-win GUI:
#this will create 2 .dbconnn/*.conf  encryptedfile 
 ~vdl/scripts/INSTALL/dbsetup

# created proper ~vdl/*.prop files
#ip-172-31-2-68::vdl:~>ls *.prop
#edgecon.prop*  edge_ip-172-31-2-68.prop*  edge.prop*	metagui.prop  msread.prop*  query.prop*  queryserver.prop*
ip

# create a role file
ip-172-31-2-68::vdl:~>ls roles_ip-172-31-2-68*
roles_ip-172-31-2-68

# Finally install the crontab.
# install crontab under the user vdl


# Testings:

# query -ls
# query -lsc -b 2005,1
# query -h localhost -lsc -b 2005,1

# query -lsc -b 2016,1

