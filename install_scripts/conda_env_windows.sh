#!/bin/bash

#set -e
#set -o pipefail

help()
{
    echo "-------------------------------------------"
    echo "Usage: source conda_env_windows.sh env_name"
    echo "Assumes Anaconda3-2021.11 is in PATH"
    echo "-------------------------------------------"
    exit 0
}
VALID_ARGUMENTS=$#
if [ "$VALID_ARGUMENTS" -eq 0 ]; then
    help
fi

echo "==== Create a conda environment with requisite packages ===="
conda create -n $1 -c conda-forge python=3.6.8 obspy=1.1.0 cartopy shapely tqdm geographiclib cython h5py mpi4py==3.0.3 netCDF4==1.4.0 pyproj==3.0.1 scipy scikit-learn==0.22.2.post1 rtree==0.9.7 PyWavelets==1.1.1 matplotlib==3.3.4 pandas==1.1.5 opencv==4.5.3

echo "==== Activating conda environment ===="
conda activate $1

echo "==== Installing packages not available through conda ===="
pip3 install pip==21.1.2
pip3 install click==7.1.2 pyasdf==0.5.1 obspyh5==0.5.0 ordered_set ujson psutil PyPDF2==1.26.0 sortedcontainers stockwell==1.0.7 pillow==8.4.0 ipython==7.10.1 basemap==1.3.2 descartes==1.1.0 PyYAML rasterio==1.2.10

echo "==== Installing a whittled down version of the rf package ===="

pushd /tmp/
git clone --depth 1 --branch v0.8.0 https://github.com/trichter/rf.git
cd rf
# remove toeplitz as a dependency, since it requires a compiler to be available
sed -i "s/'cartopy', 'geographiclib', 'shapely', 'toeplitz', 'tqdm']/'cartopy', 'geographiclib', 'shapely', 'tqdm']/" setup.py
python setup.py install
cd ..
rm -rf rf
popd

echo "#######################################################################"
echo "######################### VARIABLES TO EXPORT #########################"
echo "#######################################################################"
echo "Add the following to your .bashrc"
echo ""
echo "conda activate " $1
echo ""
echo "export PYTHONPATH=___P_A_T_H___T_O___/hiperseis/:\$PYTHONPATH"
echo ""
