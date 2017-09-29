#! /bin/bash

sudo yum update -y
sudo yum install git -y
sudo yum install wget -y


sudo yum install xorg-x11-xauth.x86_64 xorg-x11-server-utils.x86_64 dbus-x11.x86 -y
sudo yum install xeyes -y

sudo yum install expect expectk -y

# exit and 
# ssh -X login again to test xeyes


sudo yum install mysql56* -y
sudo yum install  MySQL-python27.x86_64 -y

sudo chkconfig mysqld on

sudo service mysqld start

# set mysql root passwd
mysqladmin -u root password \"CWBr00tpwd\"
# Alternatively run   /usr/libexec/mysql56/mysql_secure_installation

#verify root pass setup, check version
mysqladmin -u root -pCWBr00tpwd version


# create users and set permissions
sudo adduser vdl
sudo passwd -d vdl

# make vdl home group accessible
sudo chmod -R 750 /home/vdl

sudo adduser ketchum
sudo usermod -g vdl ketchum
sudo adduser sudipta
sudo usermod -g vdl sudipta
sudo adduser niket
sudo usermod -g vdl niket
sudo adduser fei
sudo usermod -g vdl fei

sudo usermod -g vdl ec2-user

sudo cp /etc/sudoers /tmp/sudoers.bak

sudo chmod 755 /tmp/sudoers.bak
sudo chown ec2-user /tmp/sudoers.bak
sudo chgrp ec2-user /tmp/sudoers.bak

sudo echo "niket	ALL=(vdl)	ALL" >> /tmp/sudoers.bak
sudo echo "ketchum	ALL=(vdl)	ALL" >> /etc/sudoers.bak
sudo echo "sudipta	ALL=(vdl)	ALL" >> /etc/sudoers.bak
sudo echo "fei	ALL=(vdl)	ALL" >> /etc/sudoers.bak
sudo echo "vdl	ALL=(vdl)	NOPASSWD: ALL" >> /etc/sudoers.bak

sudo visudo -cf /tmp/sudoers.bak
if [ $? -eq 0 ]; then
	sudo cp /tmp/sudoers.bak /etc/sudoers
else
	echo "Could not modify /etc/sudoers file. Please do this manually."
fi

# become user vdl and get the software from USGS ftp site
sudo -u vdl bash install-vdl-bash

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

# created proper ~vdl/*.prop files0
#ip-172-31-2-68::vdl:~>ls *.prop
#edgecon.prop*  edge_ip-172-31-2-68.prop*  edge.prop*	metagui.prop  msread.prop*  query.prop*  queryserver.prop*

# create a role file
#ip-172-31-2-68::vdl:~>ls roles_ip-172-31-2-68*
#roles_ip-172-31-2-68

# Finally install the crontab.
# install crontab under the user vdl


# Testings:

# query -ls
# query -lsc -b 2005,1
# query -h localhost -lsc -b 2005,1

# query -lsc -b 2016,1

