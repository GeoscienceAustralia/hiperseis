#!/bin/bash

# Fei Zhang 2020-05-06

#PBS -P vy72
#PBS -N create_ensemble
#PBS -q hugemem
#PBS -l walltime=2:00:00,mem=1000GB,ncpus=48,jobfs=200GB
#PBS -l storage=scratch/fxz547+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M fei.zhang@ga.gov.au
#PBS -m bae


# Load modules in NCI gadi system.
module purge
module load pbs
module load python3/3.7.4
module loat hdf5/1.10.5p 
module load openmpi/3.1.4
#export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
#export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PYTHONPATH=/home/547/fxz547/github/hiperseis/:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


# activate the virtual env  with parallel h5py together with many other python libs:
source /g/data/ha3/Passive/Software/VENV/para_h5py/bin/activate

# For new environment, need to install dependency libs:
# pip install ordered_set tqdm  click  ujson  psutil  pykml  pyyaml  joblib pywavelets


# files paths 
CATALOGD="/g/data/ha3/Passive/Events/Unified/"  # a csv file 
INVENTORY_XML="/g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml"  
ISC_STATION_COORD="/g/data/ha3/Passive/SHARED_DATA/Inventory/station_coords/stations.kml"
OUTPUT1="/g/data/ha3/GASeisDataArchive/DevSpace/pick_workflow2/step1"
OUTPUT2="/g/data/ha3/GASeisDataArchive/DevSpace/pick_workflow2/step2B"


# the cmdline 
mpirun -np $PBS_NCPUS python3 /home/547/fxz547/github/hiperseis/seismic/pick_harvester/createEnsembleXML.py $CATALOGD $INVENTORY_XML $ISC_STATION_COORD $OUTPUT2 --p-arrivals $OUTPUT1/p_combined.txt --s-arrivals $OUTPUT1/s_combined.txt

# example cmd by RH
# mpirun -np 96 python /g/data/ha3/rakib/seismic/pst/passive-seismic/seismic/pick_harvester/createEnsembleXML.py /g/data/ha3/Passive/Events/Unified/ /g/data/ha3/Passive/SHARED_DATA/Inventory/networks_fdsnstationxml/inventory.xml /g/data/ha3/Passive/SHARED_DATA/Inventory/station_coords/stations.kml /g/data/ha3/Passive/SHARED_DATA/Scratch/picking_workflow/step_2 --p-arrivals /g/data/ha3/Passive/SHARED_DATA/Scratch/picking_workflow/step_1/p_combined.txt --s-arrivals /g/data/ha3/Passive/SHARED_DATA/Scratch/picking_workflow/step_1/s_combined.txt

