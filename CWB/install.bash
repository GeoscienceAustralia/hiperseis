#! /bin/bash

sudo yum update
sudo yum install git

sudo yum install xorg-x11-xauth.x86_64 xorg-x11-server-utils.x86_64 dbus-x11.x86
sudo yum install xeyes

# exit and 
# ssh -X login again to test xeyes


sudo yum install mysql56*
sudo yum install  MySQL-python27.x86_64
