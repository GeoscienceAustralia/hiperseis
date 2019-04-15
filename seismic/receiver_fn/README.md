# Receiver Function (RF) tools

## Workflow

The RF workflow is based on the following main steps:

 - Data preparation. Extract waveforms associated with moderate events.
 - Calculate RF performing deconvolution to ZRT or QLT coordinate system
 - Discard bad results
 - Analyse RF to identify multiples
 - Perform 1D Earth structure inversion 

There are a number of manuals and literature about RF. You can refer to examples provided by
[Charles Ammon](http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/rftn01.html), manuals
distributed with [python RF libraries](https://rf.readthedocs.io/en/latest/) or [chapter 4.1 of
Tom Richter's dissertation](http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf).

## Data preparation

The first step in data preparation is to extract segments of waveforms stored in the H5 waveform
file containing all deployment data.
It is based on a set of registered teleseismic earthquakes extracted from the IRIS web service.
The main program to extract waveforms from H5 file is `prepare_rf_data.py`. The
following parameters must be set up in the program body:

 - File name of the H5 file that contains all waveform data
 - Start and end dates of the deployment
 - Station inventory XML file
 - Output file name

First, the program will create a catalogue of the events from 15 to 90 degrees of distance relative
to the centre of Australia, and then quit. The program should then be run a second time, either
directly or using `qsub` command on NCI. The latter is preferable because data extraction is a very
time consuming process. If there is catalogue with specified dates it will continue to extraction without quitting.

## Receiver function calculation

Receiver function calculation is a very simple strightforward process using the `generate_rf.py` program.
It reads H5 file with waveform segments and generates LQT receiver functions without moveout.
The parameters to specify are:

 - H5 input file name
 - H5 output file name
 - filtering parameters
 - stations to exclude if any
 - sampling rate to interpolate

There is a work-around to preserve H5 header information due to a obspy library bug [!REF?].
It was reported and will be fixed soon.

The parallelized version of the same program is `generate_rf_parallel.py`.


## Filtering and removing bad RF

There are number of ways to select good results and remove bad ones. 
Currently there are three methods implemented:
  - grouping by similarity - grouping method based on calculation of Euclidean distances and clustering by similarity (aca machine learning approach)
  - coherence - finding the coherent signals (in frequency domain) relative to median. Consequently, moveout should be applied to use this technique.
  - knive - analysing the change of RMS relative to median. Noisy stations will give higher input. Moveout should be applied to use this technique.

In the last two approaches, RFs are compared to the median and, consequently, the moveout should be applied to built one.
Currently it works with LQT convention that can be changed within the program. Existing libraries allow one to tailor
the process and more options can be easily added.

The example of how methods work can be seen running `rf_smart_tune_up.py` code. It uses `rf_test.dat` file to illustrate the process.

The actual program that works with RF is `rf_smart_bin.py`.

Following parameters should be specified within program:

- H5 Input file name
- H5 Output file name

Although there are number of methods coded, the current version outputs only the data grouped by similarity.
Please feel free to change and modify this code. Output RF are cut -5 to 60 seconds relative to arrival time.

## Analysing RF using vespagrams

Vespagrams are images of receiver functions ordered by distance or slowness. The character (shapes) of conversion distribution
could tell if they are direct conversions such as Pms or its multiples such as PpPms or PpSms.
The good reading about application of vespagrams is * Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885 *

The tool to plot vespagrams is `rf_slow_stack.py` and it takes set of RF stored in H5 file (python rf package). 
The file name must be specified within program body as well as set of parameters. Currently the plotting is set from -5 to 20 seconds and 5 to 9 s/degrees.
The filtering is bandpass from 0.03 to 0.5 Hz.

This program does not modify input data and produces only figures as a single PDF file. The parameters of the data processing are saved within PDF file in document properties.

The naming of the output PDF file is network-rf-vespagrams.pdf . The date stamp will be added into the middle of the PDF file name in case if such name already exists.
The date can be translated into human friendly format using Linux command line command `date -d @number` such as:
 
```
$ date -d @1542926631
Fri Nov 23 09:43:51 AEDT 2018
```

## Plotting station map and the ray piercing points

The `plot_map.py` draws a station map and piercing points of the receiver functions. It is highly adaptable plotting routine with use of high-resolution topography and automatic adjustment of the station name positions. It requires installation of `https://github.com/Phlya/adjustText` for better results.

## Visualization and extraction of RF for inversion

The program to visualize and extract RF for inversion is `extract_rf.py`.
The RF input file must be specified within the body of the program.
It shows different stacking options and allows to extract specific station as ASCII file for further inversion.


## Configuration and used modules on NCI

### Currently Loaded Modulefiles:
```
  1) pbs                    4) dot                    7) python/2.7.13         10) geos/3.5.0
  2) gcc/4.9.0              5) gmt/5.1.0              8) szip/2.1              11) proj/4.9.3
  3) openmpi/1.8.4-debug    6) intel-mkl/17.0.1.132   9) hdf5/1.8.14       
```
### PIP Installed modules (pip freeze --user)
```
adjustText==0.7.3
arrow==0.12.1
atomicwrites==1.2.1
attrs==18.2.0
backports.functools-lru-cache==1.5
basemap==1.1.0
Cartopy==0.16.0
certifi==2018.10.15
cftime==1.0.3.4
chardet==3.0.4
Click==7.0
colorama==0.4.0
configparser==3.5.0
dill==0.2.8.2
enum34==1.1.6
fastdtw==0.3.2
flake8==3.5.0
Flask==1.0.2
future==0.16.0
geographiclib==1.49
gmt-python==0.1a3
h5py==2.8.0
idna==2.7
intervaltree==2.1.0
isodate==0.6.0
itsdangerous==1.1.0
Jinja2==2.10
joblib==0.13.0
jobspy==0.26.1
lxml==4.2.5
MarkupSafe==1.1.0
mccabe==0.6.1
more-itertools==4.3.0
netCDF4==1.4.2
networkx==2.2
obspy==1.1.0
obspyh5==0.3.2
pathlib2==2.3.2
pluggy==0.8.0
prov==1.5.2
psutil==5.4.8
py==1.7.0
pyasdf==0.4.0
pycodestyle==2.3.1
pyflakes==1.6.0
pyproj==1.9.5.1
pyshp==1.2.12
rdflib==4.2.2
requests==2.20.0
rf==0.6.2
Rtree==0.8.3
scandir==1.9.0
scikit-learn==0.20.1
Shapely==1.6.4.post2
sklearn==0.0
sortedcontainers==2.1.0
SQLAlchemy==1.2.12
toeplitz==0.1.3
tqdm==4.27.0
ujson==1.35
urllib3==1.24
Werkzeug==0.14.1
```
