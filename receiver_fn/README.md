# Receiver Function (RF) tools


## Data preparation

First step in data preparation is to extract segments of waveforms stored in the H5 waveform file containing all deployment data.
It is based on set of registered teleseismic earthquakes extracted from IRIS website. 
The main program to extract waveforms from H5 file is **prepare_rf_data.py** 
Following parameters must be set up in the program body:

 - File name of the H5 file that contains all waveform data
 - Start and end dates of the deployment
 - Station inventory XML file
 - Output file name

Firs, program will create catalogue of the events from 15 to 90 degrees of distance irelative to centre of Australia and quit. It should be run again directly or using `qsub` command.
Latter is preferable because data extraction is very time consuming process. If there is catalogue with specified dates it will continue to extraction without quitting.


## Receiver function calculation

It is very simple strightforward process using **generate_rf.py** program. It reads H5 file with waveform segments and generate LQT receiver functions without moveout.
The parameters to specify are:

 - H5 input file name
 - H5 output file name
 - filtering parameters
 - stations to exclude if any
 - sampling rate to interpolate

There is a work-around to preserve H5 header information due to the obspy library bug. It was reported and will be fixed soon.

The parallelized version of the same program is **generate_rf_parallel.py** 


## Filtering and removing bad RF

There are number of ways to select good results and remove bad ones. 
Currently there are three methods implemented:

  - grouping by similarity - grouping method based on calculation of euclidean distances and clustering by similarity ( aca machine learning approach)
  - coherence - finding the coherent signals (in frequency domain) relative to median. Consequently, moveout should be applied to use this technique
  - knive - analysing the change of RMS relative to median. Noisy stations will give higher input. Moveout should be applied to use this technique

In the last two approaches RFs are compared to the median and , consequently, the moveout should be applied to built one.
Currently it works with LQT convention that can be changed within the program. Existing libraries allow to tailor process and more options can be easily added.

The example of how methods work can be seen running **rf_smart_tune_up.py** code. It uses `rf_test.dat` file to illustrate the process.

The actual program that works with RF is **rf_smart_bin.py**

Following parameters should be specified within program:

- H5 Input file name
- H5 Output file name

Although there are number of methods coded the current version outputs only the data grouped by similarity. Please feel free to change and modify this code.
Output RF are cut -5 to 60 seconds relative to arrival time.

## Analysing RF using vespagrams

Vespagrams are images of receiver functions ordered by distance or sloweness. The character (shapes) of conversion distribution
could tell if they are direct conversions such as Pms or its multiples such as PpPms or PpSms.
The good reading about application of vespagrams is * Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885 *

The tool to plot vespagrams is **rf_slow_stack.py** and it takes set of RF stored in H5 file (python rf package). 
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

The **plot_map.py** draws a station map and piercing points of the receiver functions. It is highly adaptable plotting routine with use of high-resolution topography and automatic adjustment of the station name positions. It requires installation of `https://github.com/Phlya/adjustText` for better results.
