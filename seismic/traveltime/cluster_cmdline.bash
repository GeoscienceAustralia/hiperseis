#!/usr/bin/env bash
# This can be run in raijin

WORKDIR=$1  # /g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/run3

if [ ! -d "$WORKDIR" ] 
then
    echo "Directory $WORKDIR DOES NOT exists." 
    exit 1 # die with error code 1 
fi

cd $WORKDIR

module load python3/3.6.2
module load hdf5/1.8.14p openmpi/1.8
module load geos/3.5.0 netcdf/4.4.1.1

# setup environment
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


# User please modify the below INSTALL_DIR, where pstvenv is installed (see install_venv.sh) 
INSTALL_DIR=/g/data/ha3/PST2

# ELLIPCORR env variable should point to `passive-seismic/ellip-corr` dir
export ELLIPCORR=/g/data/ha3/fxz547/Githubz/passive-seismic/ellip-corr

# start the virtualenv workon seismic
source $INSTALL_DIR/pstvenv/bin/activate

######################################################

# gather
# mpirun --mca mpi_warn_on_fork 0 cluster gather /g/data/ha3/events_xmls_sc3ml -w "P S"
#mpirun --mca mpi_warn_on_fork 0 cluster gather /g/data1a/ha3/fxz547/travel_time_tomography/run4_events/pbs_events_paths.txt  -w "P S"

# re-hash into a csv file with the right columns for sorting
# cat outfile_P_header.csv | cut -d , -f 1-10,14-17 >  P_out.csv
# cat outfile_S_header.csv | cut -d , -f 1-10,14-17 >  S_out.csv

# sort them
cluster sort2 S_out.csv 10. -s sorted_S.csv
cluster sort2 P_out.csv 5. -s sorted_P.csv

# zone local global 
cluster zone sorted_S.csv -z '0 -50.0 100 190' -r sorted_region_S.csv -g sorted_global_S.csv
cluster zone sorted_P.csv -z '0 -50.0 100 190' -r sorted_region_P.csv -g sorted_global_P.csv

# plots
cluster plot sorted_S.csv '0 -50.0 100 190'
cluster plot sorted_P.csv '0 -50.0 100 190'

