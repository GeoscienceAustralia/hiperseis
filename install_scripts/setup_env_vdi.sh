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
module load openmpi
module load mpi4py
module load hdf5/1.10.5p
module load geos
module load proj/6.2.1
module load gdal/3.0.2
module load fftw3/3.3.8
module load netcdf
EOF
)
echo "$MODULES"
eval "$MODULES"

echo "================ Creating virtual-environment in $ENV_DIR ================"
python3.7 -m venv $ENV_DIR

echo "================ Activating virtual-environment ================"
source $ENV_DIR/bin/activate

echo "================ Installing pip 21.1.2 ================"
python3.7 -m pip install pip==21.1.2

echo "================ Installing numpy 1.18.5 ================"
python3.7 -m pip install numpy==1.18.5

echo "================ Installing cython 0.29.22 ================"
python3.7 -m pip install cython==0.29.22

if [ 1 -eq 1 ]
then
    echo "================ Building pyfftw against fftw3/3.3.8 ================"
    mkdir -p $ENV_DIR/src
    rm -rf $ENV_DIR/src/pyFFTW-0.12.0
    cd $ENV_DIR/src
    wget https://github.com/pyFFTW/pyFFTW/archive/v0.12.0.tar.gz
    tar -zxvf v0.12.0.tar.gz
    cd pyFFTW-0.12.0/
    python3.7 setup.py build_ext --inplace
    python3.7 setup.py install
    cd ../../
fi

if [ 1 -eq 1 ]
then
    echo "================ Installing standard packages ================"
STD_PACKAGES=$(cat << EOF
python3.7 -m pip install scipy==1.4.1 
python3.7 -m pip install joblib==0.14.1 
python3.7 -m pip install scikit-learn==0.22.2.post1 
python3.7 -m pip install tqdm==4.43.0 
python3.7 -m pip install sortedcontainers==2.3.0
python3.7 -m pip install obspy==1.1.0 
python3.7 -m pip install click==7.1.2 
python3.7 -m pip install netCDF4==1.5.0 --only-binary netCDF4
python3.7 -m pip install pyasdf==0.5.1 
python3.7 -m pip install ordered_set ujson psutil 
python3.7 -m pip install obspyh5==0.5.0 
python3.7 -m pip install matplotlib==3.3.4 
python3.7 -m pip install PyPDF2==1.26.0 
python3.7 -m pip install shapely==1.8.1.post1 --no-binary shapely 
python3.7 -m pip install cartopy==0.19.0.post1 --no-binary cartopy 
python3.7 -m pip install PyWavelets==1.1.1 
python3.7 -m pip install rf==0.8.0 
python3.7 -m pip install affine==2.3.0 
python3.7 -m pip install future==0.18.2 
python3.7 -m pip install pandas==1.1.5 
python3.7 -m pip install pyproj==3.0.1 
python3.7 -m pip install rtree==0.9.7
python3.7 -m pip install mtspec==0.3.2
python3.7 -m pip install stockwell==1.0.7
python3.7 -m pip install ipython==7.10.1
python3.7 -m pip install opencv-python==4.5.3.56
python3.7 -m pip install pillow==8.4.0
python3.7 -m pip install basemap==1.3.2
python3.7 -m pip install descartes==1.1.0
EOF
)
    echo "$STD_PACKAGES"
    eval "$STD_PACKAGES"
fi

echo "================ Donwloading cartopy coastlines ================"
python3.7 -c "import cartopy.crs as ccrs; import matplotlib.pyplot as plt; crs = ccrs.PlateCarree(); fig = plt.figure(); ax = plt.subplot(projection=crs); ax.coastlines('50m'); plt.savefig('/tmp/cartopy_$USER.pdf'); print('SUCCESS..\n');"

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

