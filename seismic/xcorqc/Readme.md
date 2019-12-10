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
asn1crypto==0.24.0
attrs==19.3.0
backcall==0.1.0
bleach==3.1.0
certifi==2019.11.28
cffi==1.11.5
cftime==1.0.4.2
chardet==3.0.4
Click==7.0
colorama==0.4.3
configobj==5.0.6
cryptography==2.3
cycler==0.10.0
decorator==4.2.1
defusedxml==0.6.0
dill==0.3.1.1
entrypoints==0.3
flake8==3.7.9
future==0.18.2
gpg==1.10.0
h5py==2.10.0
idna==2.5
importlib-metadata==1.2.0
iniparse==0.4
iotop==0.6
ipykernel==5.1.3
ipython==7.10.1
ipython-genutils==0.2.0
ipywidgets==7.5.1
isc==2.0
isodate==0.6.0
jedi==0.15.1
Jinja2==2.10.3
jsonschema==3.2.0
jupyter==1.0.0
jupyter-client==5.3.4
jupyter-console==6.0.0
jupyter-core==4.6.1
kiwisolver==1.1.0
lxml==4.4.2
MarkupSafe==1.1.1
matplotlib==3.1.2
mccabe==0.6.1
mistune==0.8.4
more-itertools==8.0.2
mpi4py==3.0.3
nbconvert==5.6.1
nbformat==4.4.0
netCDF4==1.4.0
netifaces==0.10.6
networkx==2.4
notebook==6.0.2
numpy==1.16.4
obspy==1.1.0
ofed-le-utils==1.0.3
packaging==19.2
pandocfilters==1.4.2
parso==0.5.1
pciutils==2.3.6
perf==0.1
pexpect==4.7.0
pickleshare==0.7.5
pluggy==0.13.1
ply==3.9
prometheus-client==0.7.1
prompt-toolkit==2.0.10
prov==1.5.3
psutil==5.6.7
ptyprocess==0.6.0
py==1.8.0
pyasdf==0.5.1
pycodestyle==2.5.0
pycparser==2.14
pyflakes==2.1.1
Pygments==2.5.2
pygobject==3.28.3
pyOpenSSL==18.0.0
pyparsing==2.1.10
pyrsistent==0.15.6
pytest==5.3.1
python-dateutil==2.6.1
python-dmidecode==3.12.2
python-linux-procfs==0.6
pyudev==0.21.0
pyzmq==18.1.1
qtconsole==4.6.0
rdflib==4.2.2
requests==2.22.0
rhnlib==2.8.6
rpm==4.14.2
schedutils==0.6
scipy==1.3.3
Send2Trash==1.5.0
six==1.11.0
slip==0.6.4
slip.dbus==0.6.4
SQLAlchemy==1.3.11
SSSDConfig==2.0.0
syspurpose==1.23.8
terminado==0.8.3
testpath==0.4.4
tornado==6.0.3
traitlets==4.3.3
ujson==1.35
urllib3==1.25.7
wcwidth==0.1.7
webencodings==0.5.1
widgetsnbextension==3.5.1
zipp==0.6.0
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
