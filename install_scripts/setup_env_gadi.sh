#!/bin/bash

set -e 
set -o pipefail

help()
{   
    echo "----------------------------"
    echo "Usage: setup_env.sh env_name"
    echo "----------------------------"
    exit 2
}

VALID_ARGUMENTS=$#
if [ "$VALID_ARGUMENTS" -eq 0 ]; then
    help
fi

ENV_DIR=`pwd`/$1

echo "================ Loading NCI Modules ================"
MODULES=$(cat << EOF
module purge
module load pbs
module load python3-as-python
module load openmpi/3.1.4
module load hdf5/1.10.5p
module load geos
module load proj/6.2.1
module load gdal/3.0.2/
module load fftw3/3.3.8
EOF
)
echo "$MODULES"
eval "$MODULES"

echo "================ Creating virtual-environment in $ENV_DIR ================"
python -m venv $ENV_DIR

echo "================ Activating virtual-environment ================"
source $ENV_DIR/bin/activate

echo "================ Installing pip 21.1.2 ================"
pip3.6 install pip==21.1.2

echo "================ Installing numpy 1.18.5 ================"
pip3.6 install numpy==1.18.5

echo "================ Building mpi4py against module openmpi-3.1.4 ================"
MPICC=/apps/openmpi/3.1.4/bin/mpicc pip3.6 install --no-binary=mpi4py mpi4py==3.1.3

echo "================ Installing cython 0.29.22 ================"
pip3.6 install cython==0.29.22

if [ 1 -eq 1 ]
then
    echo "================ Building parallel h5py against module hdf5/1.10.5p ================"
    echo "=== Fetching a version of h5py, tweaked for Gadi ==="
    # create src folder, remove 'h5py' if it exists and cd into src
    mkdir -p $ENV_DIR/src
    rm -rf $ENV_DIR/src/h5py
    cd $ENV_DIR/src

    git clone --single-branch --branch 3.1.0-gadi-tweaks https://github.com/rh-downunder/h5py.git
    cd h5py
    echo "=== Compile and install h5py ==="
    CC=mpicc HDF5_MPI="ON" HDF5_DIR=/apps/hdf5/1.10.5p/ python setup.py install
    cd ../../
fi

if [ 1 -eq 1 ]
then
    echo "================ Installing standard packages ================"
STD_PACKAGES=$(cat << EOF
pip3.6 install scipy==1.4.1 
pip3.6 install joblib==0.14.1 
pip3.6 install scikit-learn==0.22.2.post1 
pip3.6 install tqdm==4.43.0 
pip3.6 install sortedcontainers==2.3.0
pip3.6 install obspy==1.1.0 
pip3.6 install click==7.1.2 
pip3.6 install netCDF4==1.4.0 
pip3.6 install pyasdf==0.5.1 
pip3.6 install ordered_set ujson psutil 
pip3.6 install obspyh5==0.5.0 
pip3.6 install matplotlib==3.3.4 
pip3.6 install PyPDF2==1.26.0 
pip3.6 install shapely==1.8.1.post1 --no-binary shapely 
pip3.6 install cartopy==0.19.0.post1 --no-binary cartopy 
pip3.6 install PyWavelets==1.1.1 
pip3.6 install rf==0.8.0 
pip3.6 install affine==2.3.0 
pip3.6 install future==0.18.2 
pip3.6 install pandas==1.1.5 
pip3.6 install pyproj==3.0.1 
pip3.6 install rtree==0.9.7
pip3.6 install mtspec==0.3.2
pip3.6 install stockwell==1.0.7
pip3.6 install ipython==7.10.1
pip3.6 install opencv-python==4.5.3.56
pip3.6 install pillow==8.4.0
EOF
)
    echo "$STD_PACKAGES"
    eval "$STD_PACKAGES"
fi

echo "================ Donwloading cartopy coastlines ================"
python -c "import cartopy.crs as ccrs; import matplotlib.pyplot as plt; crs = ccrs.PlateCarree(); fig = plt.figure(); ax = plt.subplot(projection=crs); ax.coastlines('50m'); plt.savefig('/tmp/cartopy_$USER.pdf'); print('SUCCESS..\n');"

echo "Completed installing all dependencies within virtual environment in: $ENV_DIR .."
echo ""

echo "#######################################################################"
echo "######################### VARIABLES TO EXPORT #########################"
echo "#######################################################################"
echo "Add the following to your .bashrc and/or PBS job-script:"
echo ""
echo "$MODULES"
echo ""
echo "source $ENV_DIR/bin/activate"
echo ""
echo "export PYTHONPATH=___P_A_T_H___T_O___/hiperseis/:\$PYTHONPATH"
echo ""

