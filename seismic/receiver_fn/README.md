# Receiver Function (RF) tools

## Workflow

The RF workflow is based on the following main steps:

 - Data preparation: extract raw waveforms in time windows associated with selected seismic
   events from ASDF database for the deployment
 - Calculate Receiver Functions: performing transformation and deconvolution to ZRT or QLT coordinate system
 - Quality filtering on RFs: discard bad results by filtering
 - Analyse RF to identify multiples
 - Perform analysis such as 1D Earth structure inversion, H-k stacking, etc...

There are a number of manuals and literature about RF. You can refer to examples provided by
[Charles Ammon](http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/rftn01.html), manuals
distributed with [python RF libraries](https://rf.readthedocs.io/en/latest/) or [chapter 4.1 of
Tom Richter's dissertation](http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf).

## Data preparation

The first step in data preparation is to extract segments of waveforms stored in the H5 waveform
file containing all deployment data.
It is based on a set of registered teleseismic earthquakes extracted from the ISC web service.
The main program to extract raw waveforms from H5 file is `extract_event_traces.py`. The
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

Receiver function calculation is a very simple straightforward process using the `generate_rf.py` program.
It reads H5 file with waveform segments and generates LQT receiver functions without moveout.
The parameters to specify are:

 - H5 input file name
 - H5 output file name
 - filtering parameters
 - stations to exclude if any
 - sampling rate to interpolate

## Filtering and removing bad RF

There are number of ways to select good results and remove bad ones. 
Currently there are three methods implemented:
  - grouping by similarity - grouping method based on calculation of Euclidean distances and clustering by similarity
    (aca machine learning approach)
  - coherence - finding the coherent signals (in frequency domain) relative to median. Consequently, moveout should be
    applied to use this technique.
  - knive - analysing the change of RMS relative to median. Noisy stations will give higher input. Moveout should be
    applied to use this technique.

In the last two approaches, RFs are compared to the median and, consequently, the moveout should be applied to built one.
Currently it works with LQT convention that can be changed within the program. Existing libraries allow one to tailor
the process and more options can be easily added.

The actual script for generating QA'd RFs is `rf_quality_filter.py`.

Following parameters should be specified within program:

- H5 Input file name
- H5 Output file name

Although there are number of methods coded, the current version outputs only the data grouped by similarity.
Please feel free to change and modify this code. Output RF are cut -5 to 60 seconds relative to arrival time.

## Analysing RF using vespagrams

Vespagrams are images of receiver functions ordered by distance or slowness. The character (shapes) of conversion distribution
could tell if they are direct conversions such as Pms or its multiples such as PpPms or PpSms.
The good reading about application of vespagrams is *Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885*.

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


## Inversion

TBD


## Models folder

TBD

## Legacy code

The `legacy` folder is for legacy code which is retained temporarily for reference purposes, and whose methods
may in future be rolled into the existing RF workflow.
