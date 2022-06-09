# Dependencies
The workflow requires MPI (mpi4py) and parallel HDF5 (h5py) capabilities on the NCI.
Installation instructions for NCI (Gadi ) are as follows:

### Load system modules:
  1. `module purge`
  2. `module load pbs` 
  3. `module load python3-as-python`
  4. `module load openmpi/3.1.4`
  5. `module load hdf5/1.10.5p`
  6. `module load geos`
  7. `module load proj/6.2.1`

### Remove old packages
1. `rm -rf ~/.local/lib/python3.6/site-packages/h5py*`
2. `rm -rf ~/.local/lib/python3.6/site-packages/mpi4py*`
3. `rm -rf ~/.local/lib/python3.6/site-packages/cartopy`
4. `rm -rf ~/.local/lib/python3.6/site-packages/Cartopy*`
5. `rm -rf ~/.local/lib/python3.6/site-packages/shapely`
6. `rm -rf ~/.local/lib/python3.6/site-packages/Shapely*`
7. `rm -rf ~/.local/lib/python3.6/site-packages/*geos*`

### Upgrade pip

 1. `pip3.6 install pip==21.1.2 --user`

### numpy

 1. `pip3.6 install numpy==1.18.5 --user`

### Install mpi4py that uses the correct OpenMPI libs

  1. `MPICC=/apps/openmpi/3.1.4/bin/mpicc pip3.6 install --no-binary=mpi4py mpi4py==3.1.3 --user` Note that we use `pip3.6`, the system-provided pip for python 3.6

### Build Parallel H5PY

1. `pip3.6 install cython==0.29.22 --user`
2. `git clone --single-branch --branch 3.1.0-gadi-tweaks https://github.com/rh-downunder/h5py.git` Pull a branch (based on version 3.1.0) from a fork of h5py, adapted for Gadi.
3. `cd h5py`
4. `CC=mpicc HDF5_MPI="ON" HDF5_DIR=/apps/hdf5/1.10.5p/ python setup.py install --user` Configure, build and install

### Setup standard packages

  1. `pip3.6 install obspy==1.1.0 --user`
  2. `pip3.6 install click==7.1.2 --user `
  3. `pip3.6 install netCDF4==1.4.0 --user`
  4. `pip3.6 install pyasdf==0.5.1 --user`
  5. `pip3.6 install ordered_set ujson psutil --user`
  6. `pip3.6 install obspyh5==0.5.0 --user`  
  7. `pip3.6 install matplotlib==3.3.4 --user`
  8. `pip3.6 install sortedcontainers --user`
  9. `pip3.6 install PyPDF2==1.26.0 --user`
  10. `pip3.6 install shapely==1.8.1.post1 --no-binary shapely --user`
  11. `pip3.6 install cartopy==0.19.0.post1 --no-binary cartopy --user`
  12. `python -c "import cartopy.crs as ccrs; import matplotlib.pyplot as plt; crs = ccrs.PlateCarree(); fig = plt.figure(); ax = plt.subplot(projection=crs); ax.coastlines('50m'); plt.savefig('"/tmp/cartopy_$USER.pdf"'); print('\nSUCCESS');"`  
  13. `pip3.6 install PyWavelets==1.1.1 --user`
  14. `pip3.6 install rf==0.8.0 --user`

Step 12 ensures coastline shapefiles used by Cartopy are downloaded and available for use before 
jobs are launched on NCI compute nodes that do not allow internet access.

### Initialize PhasePapy
Additionally, the workflow requires the PhasePapy package (a submodule of this repository) to be
initialized, if not already.
```
cd hiperseis
git submodule update --init --recursive
```

# Waveform Extraction and Orientation Analysis Tools

Scripts for extracting waveforms and generating station-orientation reports are 
described in the following sections:

## Event Trace Extraction

The first step in data preparation is to extract relevant waveform data for a set of 
registered teleseismic earthquakes obtained from the ISC web service.
The main program to extract raw waveforms from an FDSN web-service or the 
FederatedASDFDataSet infrastructure is `extract_event_traces.py`. It can be run in 
MPI-parallel mode with the FederatedASDFDataset infrastructure; otherwise it must be 
run in serial mode on an NCI login node or on the VDI, allowing access to an 
FDSN web-service.

The following parameters must be supplied to the script:

 - Waveform source (URL to an FDSN web-service or path to a text file to initialize a \
   FederatedASDFDataSet)
 - Event catalog to create or re-use
 - Network or list of networks to extract waveforms for
 - Arrival-type(s) to extract waveforms for: around P, S and Surface-waves
 - Output file name

All waveforms (except Surface-waves) are processed through an automatic (AIC-based) 
picker and the following metadata are added to all traces where a pick is found:

```
'arrival_time' # actual arrival time of a P/S-wave
'slope_ratio' # quality metric for arrival -- typically, values>5 indicate \
                robust arrivals
```
The P/S-arrival quality metric (`slope_ratio`) can be particularly useful for 
filtering out low-quality waveforms in downstream processing steps.

For usage details, see output of `python extract_event_traces.py --help`.

### Typical Usage:

Download events for OA stations (must be run on an NCI login node or on the VDI):
```
python extract_event_traces.py --waveform-database asdf_files.txt 
--event-catalog-file OA_events.xml --event-trace-datafile OA.h5 --catalog-only 
--network-list "OA"
```

Generate OA waveforms:
```
mpirun -np 48 python extract_event_traces.py --waveform-database asdf_files.txt 
--event-catalog-file OA_events.xml --event-trace-datafile OA_waveforms.h5 
--network-list "OA" --p-data --s-data --sw-data > out.txt 2>&1
```
## Orientation Analysis

Station-orientation reports are generated by an MPI-parallel driving script: 
`bulk_station_orientations.py`. The script generates two estimates of station-orientation 
based on (i) Receiver Functions and (ii) Surface-wave Polarization.

### Typical usage:
`mpirun -np 48 python bulk_station_orientations.py 
OA_waveforms_X-Y.h5 'OA' --output-basename OA_orientations`
