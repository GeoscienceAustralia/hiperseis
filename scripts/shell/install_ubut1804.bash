#!/usr/bin/env/ bash
# install passive-seispic into Fresh Ubuntu 18.04
# Fei Zhang
# 20180812

sudo apt-get update
sudo apt install python3

sudo apt install gfortran

sudo apt-get install python3-pip

sudo pip3 install pyasdf
sudo pip3 install obspy
sudo pip3 install pandas
sudo pip3 install matplotlib

git clone https://github.com/GeoscienceAustralia/passive-seismic
cd passive-seismic/
git submodule update --init --recursive --remote

cd seismic/cluster/
python3 cluster.py gather
sudo pip3 install basemap
export ELLIPCORR=$PWD/ellip-corr
echo $ELLIPCORR

mkdir pstven
#python3 -m venv  pstven/
sudo  apt-get install python3-venv
python3 -m venv  pstven/
# ls pstven/

source  pstven/bin/activate
sudo apt-get install mpich
# mpirun

pip3 install mpi4py==3.0.0 --no-binary :all:


sudo apt-get install libgeos-3.3.3 libgeos-c1 libgeos-dev
sudo apt-get install  libgeos-dev
echo $GEOS_DIR
export GEOS_DIR=$GEOS_BASE

sudo apt-get install libfreetype6-dev pkg-config libpng12-dev
sudo apt-get install libfreetype6-dev pkg-config
sudo apt-get install libpng*


pip3 install matplotlib
python3 seismic/cluster/cluster.py

sudo apt-get install libhdf5-serial-dev netcdf-bin libnetcdf-dev


sudo apt-get install libxml2*
sudo apt-get install libxml
sudo apt-get install libxml2-dev
STATIC_DEPS=true pip3 install lxml
sudo apt-get install libblas3 liblapack3 liblapack-dev libblas-dev

pip3 install --process-dependency-links -e .[dev] --no-binary :all:
echo $?


###################### Successful installation will have output look like
#
#(pstven) feizhang01@ubunt1804inst1:~/passive-seismic$ pip3 install --process-dependency-links -e .[dev] --no-binary :all:
#Obtaining file:///home/feizhang01/passive-seismic
#  DEPRECATION: Dependency Links processing has been deprecated and will be removed in a future release.
#Requirement already satisfied: Click>=6.0 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: Cython>=0.22.1 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: PyYAML>=3.11 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: basemap==1.1.0 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: chardet==3.0.4 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: h5py>=2.6.0 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: joblib in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: lxml>=3.3.5 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: matplotlib>=1.4.3 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: mpi4py==3.0.0 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: netCDF4>=1.3.0 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Requirement already satisfied: numpy>=1.9.2 in /home/feizhang01/pstven/lib/python3.6/site-packages (from Passive-Seismic==0.0.1)
#Collecting obspy>=1.1.0 (from Passive-Seismic==0.0.1)
#  Using cached https://files.pythonhosted.org/packages/c5/f5/9a6150bc707bfb82e1d0533f83d3a7ef00ae3bdb2b8c8f85f5ac2d3b2bd9/obspy-1.1.0.zip
#Collecting pandas (from Passive-Seismic==0.0.1)

# ....

#  Running setup.py install for sphinx ... done
#  Running setup.py install for sphinxcontrib-programoutput ... done
#  Running setup.py install for virtualenv ... done
#  Running setup.py install for tox ... done
#  Running setup.py develop for Passive-Seismic
#Successfully installed Jinja2-2.10 MarkupSafe-1.0 Passive-Seismic Pygments-2.2.0 alabaster-0.7.11 atomicwrites-1.1.5 attrs-18.1.0 babel-2.6.0 codecov-2.0.9 colorama-0.3.9 coverage-4.4.1 dill-0.2.8.2 docutils-0.14 flake8-3.5.0 flake8-docstrings-1.3.0 flake8-polyfill-1.0.2 ghp-import-0.5.5 imagesize-1.0.0 isodate-0.6.0 mccabe-0.6.1 more-itertools-4.3.0 networkx-2.1 obspy-1.1.0 packaging-17.1 pandas-0.23.4 pluggy-0.7.1 prov-1.5.2 py-1.5.4 pyasdf-0.4.0 pycodestyle-2.3.1 pydocstyle-2.1.1 pyflakes-1.6.0 pyqtgraph-0.10.0 pytest-3.7.1 pytest-cov-2.5.1 pytest-flake8-1.0.2 pytest-lazy-fixture-0.4.1 pytest-mock-1.10.0 pytest-regtest-1.3.0 rdflib-4.2.2 scipy-1.1.0 snowballstemmer-1.2.1 sphinx-1.7.6 sphinxcontrib-programoutput-0.11 sphinxcontrib-websupport-1.1.0 sqlalchemy-1.2.10 tox-3.2.1 virtualenv-16.0.0
#(pstven) feizhang01@ubunt1804inst1:~/passive-seismic$ echo $?
#0
#
# How to run????
#
# ssh feizhang01@35.189.33.175
#
#cd ~/testruns
#
#source  ~/pstven/bin/activate
#
#export ELLIPCORR=/home/feizhang01/passive-seismic/ellip-corr
#
#nohup cluster gather /home/feizhang01/testdata/  &
#
##################################################
#source  ~/pstven/bin/activate
#
#export ELLIPCORR=/home/feizhang01/passive-seismic/ellip-corr
#
## nohup cluster gather /home/feizhang01/testdata/  &
## mpi  using 2 processor
#mpirun -n 2 cluster gather /home/feizhang01/testdata/Passive_ANU_7W_2008-2011_extracted_events_isc/
#
