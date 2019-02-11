# Station Inventory Generation Tool

The purpose of this package is to provide tools and workflow for generating
a database of cleaned up station inventory files in Seiscomp3 format, based on
collections of station records (many of which are not in IRIS) in bespoke
`.STN` format.


## Requirements

The Python library requirements for using the tools and workflow of this package
are listed in `requirements.txt`. The package is written for Python version 2.7
or >=3.5. To install the library requirements, run

&nbsp;&nbsp;&nbsp;&nbsp;`pip install -r requirements.txt`

### Performance note

Version 1.1.0 of the obspy library has a severe performance problem scaling up
reading of large sc3ml inventory files at the size of data we typically process.
Geoscience Australia has produced a fix for this issue (TODO: add link), enabling
users to use the patched version. Until the release of future version of obspy
incorporating this patch, produce a custom build of obspy to allow runtime in the
order of minutes rather than many hours for loading sc3ml IRIS inventory.


## Workflow

1. Produce reference data file IRIS station inventory in sc3ml format:
   `python update_iris_inventory.py -o IRIS-ALL.xml`
2. Run station cleanup script, passing in IRIS inventory as input:
    `python engd2stxml.py -i IRIS-ALL.xml`

Note: in the call to `update_iris_inventory.py`, additional options are available to
narrow the scope of query to IRIS online database if desired.


## Input file formats

Station input files are a collection of bespoke .STN files generated and curated by
Geoscience Australia (see [future development roadmap][#future]). There are two file
formats (despite having common .STN filename extension):
- Networkless STN format: files `stations/BMG.STN` and `stations/ISC.STN`
- Networked STN format: files `stations/ehb.stn` and `stations/iscehb.stn`

The IRIS inventory file passed into `engd2stxml.py` can be any station inventory format
supported by Obspy. Current the `update_iris_inventory.py` script generates this in sc3ml
format.

### Networkless STN format

This file format consists of fixed width metadata fields having following format:
```
 SMPI                                                       -1.9810  138.7100      0.0   2005001  2286324  I
 SMRI  GEOFON Station SemJava, Indonesia                    -7.0492  110.4407      0.0   2005001  2286324  I
 SNAA  GEOFON/AWI StationAntarctica                        -71.6707   -2.8379      0.0   2005001  2286324  I
 SPSI                                                       -3.9646  119.7691      0.0   2005001  2286324  I
 SRBI                                                       -8.0848  115.2126      0.0   2005001  2286324  I
 STKA  Stephens Creek    NSW                               -31.8769  141.5952      0.0   2005001  2286324  I
```
Each row is a new station record with station code, coordinates and elevation. Other fields are
ignored, notably the dates which are not considered reliable. Since there are no network codes,
no dates and no channel data, a default network code of `GE` is applied, date fields are left empty
and a default channel code of `BHZ` is applied to each record. See source code of function
`engd2stxml.read_eng()` for more details.

### Networked STN format

This format consists of interleaved station header/location rows, followed optionally by an
arbitrary number of channel metadata rows. Each header row contains a station code, location,
elevation and station start/end dates. Each channel row adds alternate station code, a network
code, a channel code and start/end dates for that channel.
```
LMG      -8.9077  148.1500  1200.0 1971-01-10 07:17:04
LMGC     20.0637  -77.0050   307.0 1999-10-16 09:46:48 2000-01-20 06:13:03
LMGC     20.0637  -77.0050     0.0 2000-01-20 06:13:03
             FDSN LMGC   CW -- HHZ 1998-01-01 00:00:00 2599-12-31 23:59:59
             FDSN LMGC   CW 00 HHZ 1998-01-01 00:00:00 2010-09-01 00:00:00
             FDSN LMGC   CW 00 HHZ 2010-09-01 00:00:00 2599-12-31 23:59:59
LMHM     41.5780 -121.6570  2228.0 1985-04-10 20:37:37
LMI      54.2197   -3.3070   140.0 1992-09-26 05:45:52
LMK      53.4563   -0.3270   130.0 1992-05-21 05:00:00
             FDSN LMK    GB -- BHZ 2009-07-19 00:00:00 2599-12-31 23:59:59    50
             FDSN LMK    GB -- HHZ 2009-07-19 00:00:00 2599-12-31 23:59:59   100
LML      -6.7347  147.0090   100.0 1979-09-04 13:36:51
```
The channel rows (only) specify a network code, and repeat the station code which may or may not
match the station code stated in the preceding header row. Since standalone header rows do not
state a network code, all header rows are treated as distince records and assigned default network
code `IR` and channel code `BHZ`. If channel rows are present, they are treated as distinct channel
records with their own station and networks codes, but add channel dates and take their coordinates
from the preceding header row. For more details consult source code function
`engd2stxml.read_isc()`.


## Output artifacts

The station inventory cleanup script `engd2stxml.py` produces the following output artifacts:

| FILE NAME(S) | DESCRIPTION |
|--------------|-------------|
|LOG_LOCATION_DUPES_YYYYMMDDThhmmss.txt<BR>LOG_CODE_DUPES_YYYYMMDDThhmmss.txt | Log files of duplicate records removed|
|output/network_*.xml          | Station inventory files in FDSN station xml format|
|networks_sc3ml/network_*.xml  | Station inventory files in sc3ml format|
|INVENTORY_YYYYMMDDThhmmss.csv | Re-loadable inventory database in tabular CSV|
|INVENTORY_YYYYMMDDThhmmss.h5  | Re-loadable inventory database in tabular hdf5 under key "inventory"|
|INVENTORY_YYYYMMDDThhmmss.txt | Human readable inventory database in fixed with tabular text format|
|plots/networks/*.png          | Map graphic showing station locations for a given network|
|plots/networks/*.txt          | Per network station data in FDSN stationtxt format|
|station.txt                   | Entire inventory database in FDSN stationtxt format|
|stations/*.pkl                | As a side effect, cached `.pkl` files alongside `.STN` files used as input|

*Beware:* whenever a `.STN` file is changed, any corresponding `.pkl` file should be deleted manually
to force a regeneration of the `.pkl` file.


## Future development roadmap <a name="future"></a>

Agenda for future development:
* Make IRIS inventory update script `update_iris_inventory.py` use FDSN station XML format by default, rather
  than sc3ml, as it is more efficient to load, more widely supported and doesn't need
  seiscomp3 to be on host system.
* Since cleanup script `engd2stxml.py` no longer filters out records which are found in IRIS, the
  only purpose for loading IRIS-ALL inventory is to get a tree of valid network and
  station codes from which we can get instrument responses. For this purpose, we don't
  need sc3ml format, and we don't need channel level metadata stored in the file, so
  we could just store IRIS-ALL as FDSN xml only down to station level, which will be a
  much smaller file and could be queried in real-time directly from IRIS rather than
  cached in a file. This could remove the need for `update_iris_inventory.py` in the
  inventory cleanup workflow.
* Factor our the STN file readers into a separate module so that any reader can be
  substituted, as long as it generates Pandas DataFrame output, to support arbitrary
  input file formats.
* Performance improvements reading files in function `engd2stxml.read_isc()`. E.g. change
  header rows to include a number at the end indicating the number of subsequent channel
  rows for that station location. This will allow channel rows to be loaded in bulk instead
  of one at a time.


## Last known good configuration

Python:
```
Python 3.6.7 (default, Dec  5 2018, 15:02:05)
[GCC 4.8.5 20150623 (Red Hat 4.8.5-36)]
```

Obspy (standard release):
```
'1.1.0'
```

Obspy (custom):
```
'1.1.1rc1.post0+552.g8478dae0aa.obspy.master'
```

Python libraries (via `pip freeze`):
```
arrow==0.13.0
astroid==2.1.0
backcall==0.1.0
basemap==1.2.0
bleach==3.1.0
certifi==2018.11.29
chardet==3.0.4
ciso8601==2.1.1
Click==7.0
cycler==0.10.0
dask==1.1.0
decorator==4.3.2
defusedxml==0.5.0
dill==0.2.9
entrypoints==0.3
future==0.17.1
graphviz==0.10.1
h5py==2.9.0
idna==2.8
ipykernel==5.1.0
ipython==7.2.0
ipython-genutils==0.2.0
ipywidgets==7.4.2
isort==4.3.4
jedi==0.13.2
Jinja2==2.10
jsonschema==2.6.0
jupyter==1.0.0
jupyter-client==5.2.4
jupyter-console==6.0.0
jupyter-core==4.4.0
kiwisolver==1.0.1
lazy-object-proxy==1.3.1
line-profiler==2.1.2
llvmlite==0.27.0
lxml==4.3.0
MarkupSafe==1.1.0
matplotlib==3.0.2
mccabe==0.6.1
mistune==0.8.4
nbconvert==5.4.0
nbformat==4.4.0
notebook==5.7.4
numba==0.42.0
numexpr==2.6.9
numpy==1.16.1
-e git+https://github.com/obspy/obspy.git@8478dae0aac6d647171ab3092c48384ac666d5be#egg=obspy
pandas==0.24.1
pandoc==1.0.2
pandocfilters==1.4.2
parso==0.3.2
pexpect==4.6.0
pickleshare==0.7.5
ply==3.11
proj==0.1.0
prometheus-client==0.5.0
prompt-toolkit==2.0.8
ptyprocess==0.6.0
Pygments==2.3.1
pyimgur==0.6.0
pylint==2.2.2
pyparsing==2.3.1
pyproj==1.9.6
pyshp==2.0.1
python-dateutil==2.7.5
pytz==2018.9
pyzmq==17.1.2
qtconsole==4.4.3
requests==2.21.0
rope==0.11.0
Rtree==0.8.3
scikit-learn==0.20.2
scipy==1.2.0
seaborn==0.9.0
Send2Trash==1.5.0
six==1.12.0
SQLAlchemy==1.2.17
tables==3.4.4
terminado==0.8.1
testpath==0.4.2
tornado==5.1.1
tqdm==4.30.0
traitlets==4.3.2
typed-ast==1.2.0
urllib3==1.24.1
versioneer==0.18
wcwidth==0.1.7
webencodings==0.5.1
widgetsnbextension==3.4.2
wrapt==1.11.1
xlrd==1.2.0
```
