# Cross-Correlation Workflow

The cross-correlation functionality works off of two text-files, each containing a list of ASDF files. Users specify a subset of stations (or * for all stations) available in the ASDF files in each of the text-files to be cross-correlated. Depending on the parameterization, a cartesian product of two lists of station-names defines station-pairs that are to be cross-correlated -- see `correlator.py` for details.

Cross-correlation results are written out as NetCDF-4 files for each station-pair. Interrogating cross-correlation results requires interactive visualization capabilities. [Panoply], a freely available cross-platform tool, which is also available on NCI VDIs, can be used to visualize the results interactively. A quick intro to [Panoply] is available [here].

## Fetching Permanent Stations Data

For experimentation, `client_data.py` can be used to fetch permanent stations data, for a given time-range, at given location (lat, lon). Permanent station locations for AU can be found on the [FDSN website].

## Setting up NCI Gadi environment for running cross-correlations

These instructions are for users of the National Computational Infrastructure (NCI) facility at the
Australian National University. The following is the recommended Python environment setup process
within an individual user's account that has been validated and is a known working configuration.

### Order of operations

To set up the Python environment, the following high level order of operations must be observed:
* Load modules
* Build and install customized 3rd party Python libraries for MPI
* Install standard 3rd party Python libraries

### Library version summary

* Open MPI v2.1.6-mt
* HDF5 v1.10.5p (parallel build)
* mpi4py v3.0.3 (custom MPI build)
* osbpy 1.1.0
* click 7.0
* netCDF4 1.4.0
* h5py 2.10.0 (custom MPI build)
* pyasdf 0.5.1

### Limitations

In general virtual environments are preferred, but because of the limited number of python versions
available on Gadi and the advantages (natively compiled, etc.) they provide over precompiled versions
available through virtual environments, this setup uses the system-provided python3.6 and installs 
dependencies in user space (`--user` option of `pip`). 

### Setup process

#### Load requisite modules
  1. `module purge` is highly recommended if you have modules loaded
  2. `module load pbs` 
  3. `module load python3-as-python`
  4. `module load openmpi/2.1.6-mt`
  5. `module load hdf5/1.10.5p`

#### Setup custom packages

##### H5PY

  1. `git clone --single-branch --branch 2.10.0.gadi_tweaks https://github.com/rh-downunder/h5py.git` Pull a branch (based on version 2.10.0) from h5py repository from github fork of h5py, adapted for Gadi, for purpose of custom build
  2. `cd h5py`
  3. `CC=mpicc python setup.py configure --mpi --hdf5=/apps/hdf5/1.10.5p/` Configure with mpi enabled  
  4. `python setup.py build` Build h5py
  5. `python setup.py install --user` Install in user space

##### mpi4py
  1. `MPICC=/apps/openmpi/2.1.6-mt/bin/mpicc pip3.6 install mpi4py --user` Note that we use `pip3.6`, the system-provided pip for python 3.6

#### Setup standard packages
  1. `pip3.6 install obspy==1.1.0 --user`
  2. `pip3.6 install click --user `
  3. `pip3.6 install netCDF4==1.4.0 --user`
  4. `pip3.6 install pyasdf --user`

Date last validated: 10 December 2019

### Setup validation

The Python setup can be tested for running the cross-correlation on Raijin using scripts `validate_xcorr_setup.py`
(Python) and `validate_xcorr_runtime.sh` (shell script) in folder `hiperseis/seismic/xcorqc`. These scripts should
be run directly from that folder, not from another folder.

Firstly, run `python validate_xcorr_setup.py` from the command line. Various output will appear explaining the item
being tested and the result, plus output from Python libraries. If the test succeeds, the last line of output
should read `SUCCESS!`.

Next, run at the command line run `./validate_xcorr_runtime.sh`. You should see various output, mostly warnings. To
confirm successful completion of the test, you should see a file `validation_result/ARMA.CMSA.nc` with a current
file time stamp. If this file is not present or not with a current time stamp, then the test was not successful.

### Version dump of known good configuration

Known good `pip freeze` output:
```
atomicwrites==1.3.0
attrs==19.1.0
backports.shutil-get-terminal-size==1.0.0
certifi==2018.11.29
cftime==1.0.3.4
chardet==3.0.4
Click==7.0
colorama==0.4.1
configparser==3.7.3
cycler==0.10.0
Cython==0.26
decorator==4.1.2
dill==0.2.9
entrypoints==0.3
enum34==1.1.6
flake8==3.7.7
funcsigs==1.0.2
functools32==3.2.3.post2
future==0.17.1
h5py==2.7.0.post0
idna==2.8
ipython==5.4.1
ipython-genutils==0.2.0
isodate==0.6.0
lxml==4.3.2
matplotlib==2.0.2
mccabe==0.6.1
more-itertools==5.0.0
mpi4py==3.0.0
netCDF4==1.4.3
networkx==2.2
nose==1.3.7
numpy==1.14.2
obspy==1.1.0
pathlib2==2.3.0
pexpect==4.2.1
pickleshare==0.7.4
pluggy==0.9.0
prompt-toolkit==1.0.15
prov==1.5.3
ptyprocess==0.5.2
py==1.8.0
pyasdf==0.4.0
pycodestyle==2.5.0
pyflakes==2.1.1
Pygments==2.2.0
pyparsing==2.2.0
pytest==4.3.0
python-dateutil==2.6.1
pytz==2017.2
rdflib==4.2.2
requests==2.21.0
scandir==1.5
scipy==0.19.1
simplegeneric==0.8.1
six==1.10.0
SQLAlchemy==1.3.0
subprocess32==3.2.7
traitlets==4.3.2
typing==3.6.6
urllib3==1.24.1
wcwidth==0.1.7
```

# Visualizing Cross-Correlation Results

As there are many technical aspects to interpreting the time series cross-correlation results
generated by the above Cross-Correlation Workflow, functions are provided to convert `.nc` files
to a standard graphical visualization of the cross-correlation time series.

The following three functions from module `xcorr_station_clock_analysis` are used as entry points
for this purpose:
|Function | Purpose|
|---------|--------|
|`plot_xcorr_file_clock_analysis`| For a single `.nc` file. Does not overlay runtime options. |
|`batch_process_xcorr`| For an iterable set of `.nc` files, plot each. Adds traceability information (runtime configuration parameters).|
|`batch_process_folder`| Run `batch_process_xcorr` on all `.nc` files in a specified folder.|


[Panoply]:https://www.giss.nasa.gov/tools/panoply/
[here]:http://www.meteor.iastate.edu/classes/mt452/EdGCM/Documentation/EdGCM_Panoply.pdf
[FDSN website]:http://www.fdsn.org/networks/detail/AU/
