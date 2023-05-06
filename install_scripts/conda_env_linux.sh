#!/bin/bash

#set -e
#set -o pipefail

help()
{
    echo "-----------------------------------------"
    echo "Usage: source conda_env_linux.sh env_name"
    echo "Assumes Anaconda3-2021.11 is in PATH"
    echo "-----------------------------------------"
}
VALID_ARGUMENTS=$#
if [ "$VALID_ARGUMENTS" -eq 0 ]; then
    help
    return
fi

echo "==== Create a conda environment with requisite packages ===="
conda create -n $1 -c conda-forge python=3.6.8 gfortran_linux-64==7.5.0 gcc_linux-64==7.5.0 gxx_linux-64==7.5.0 proj4 geos

echo "==== Activating conda environment ===="
conda activate $1

echo "==== Installing gdal ===="
conda install -n $1 gdal

echo "==== Installing packages not available through conda ===="
gfortran=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gfortran;
gcc=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
pip3 install setuptools==45.0.0
pip3 install pip==21.1.2
pip3 install numpy==1.18.5
pip3 install mpi4py==3.1.3
pip3 install cython==0.29.22
pip3 install h5py==3.1.0
pip3 install pyfftw==0.12.0

pip3 install scipy==1.4.1
pip3 install joblib==0.14.1
pip3 install scikit-learn==0.22.2.post1
pip3 install tqdm==4.43.0
pip3 install sortedcontainers==2.3.0
pip3 install obspy==1.2.0
pip3 install click==7.1.2
pip3 install netCDF4==1.4.0
pip3 install pyasdf==0.5.1
pip3 install ordered_set ujson psutil
pip3 install obspyh5==0.5.0
pip3 install matplotlib==3.3.4
pip3 install PyPDF2==1.26.0
pip3 install shapely==1.8.1.post1 --no-binary shapely
pip3 install cartopy==0.19.0.post1 --no-binary cartopy
pip3 install PyWavelets==1.1.1
pip3 install rf==0.8.0
pip3 install affine==2.3.0
pip3 install future==0.18.2
pip3 install pandas==1.1.5
pip3 install pyproj==3.0.1
pip3 install rtree==0.9.7
pip3 install mtspec==0.3.2
pip3 install stockwell==1.0.7
pip3 install ipython==7.10.1
pip3 install opencv-python==4.5.3.56
pip3 install pillow==8.4.0
pip3 install basemap==1.3.0
pip3 install descartes==1.1.0
pip3 install PyYAML
pip3 install rasterio==1.2.10
pip3 install notebook==6.4.10
pip3 install ipython==7.10.0
pip3 install jedi==0.17

echo "#######################################################################"
echo "######################### VARIABLES TO EXPORT #########################"
echo "#######################################################################"
echo "Add the following to your .bashrc"
echo ""
echo "conda activate " $1
echo ""
echo "export PYTHONPATH=___P_A_T_H___T_O___/hiperseis/:\$PYTHONPATH"
echo ""
