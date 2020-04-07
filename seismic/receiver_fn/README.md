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

The work in this module utilizes the RF framework developed in the Python [`rf` library
developed by Tom Eulenfeld](https://doi.org/10.21105/joss.01808).


## Event Trace Extraction

The first step in data preparation is to extract segments of waveforms stored in the H5 waveform
file containing all deployment data.
It is based on a set of registered teleseismic earthquakes extracted from the ISC web service.
The main program to extract raw waveforms from H5 file is `extract_event_traces.py`. The
following parameters must be supplied to the script:

 - File name of the H5 file that contains all waveform data
 - Event catalog to create or re-use
 - Station inventory XML file
 - Output file name

For details of usage, see `python extract_event_traces.py --help`.


## Receiver Function Calculation

Receiver functions are generated using the `generate_rf.py` script. It reads the event traces
produced by `extract_event_traces.py` receiver functions without moveout. The receiver function
rotation method is controlled with the `--rotation-type` parameter, and the deconvolution
technique used to produce the RF is controlled with the `--deconv-domain` parameter.

The mandatory parameters to specify are the input and output file names.  There are various
other options for tuning output range and filtering, and computational performance.
See `generate_rf.py --help` for details.

The default parameters are somewhat tuned to BH*
instrument response waveforms, and should be explored during a RF studty.
A high performance computing environment is recommended for large scale bulk RF generation.

During RF generation, some removal of problematic traces occurs for problems that would be
fatal for RF generation, as follows:
1. Remove traces with invalid inclination value in the metadata (e.g. NaN).
2. Remove traces in which not all ZNE channels are present for a given event.
3. Remove traces in which the number of samples in the Z, N and E channels are not the same.
4. Remove traces in which NaN is found in the data.
5. Remove traces with zero variance in the data.
6. Remove traces for which an exception is caught during RF generation.
7. Remove traces with empty channel data or channels missing after deconvolution.


## RF Quality Filtering

Usually the RFs produced by component rotation and deconvolution are not all useful for stacking,
and further quality metrics and quality filtering must be considered. The are a family of approaches
available and this is an area of active research. Consequently, rather than applying a fixed set of
RF quality filters, the approach taken here is to directly apply only the most simple and trusted
of filters to remove some RFS (such as when a signal is all zero), and then compute a range of
quality metrics and store them in the RF metadata.  This way, the actual filtering of the RFs by
quality criteria is deferred until later visualization and analysis, allowing greater choice and
flexibility in the analysis stage. For example, the quality metadata stored in each RF is used by the reporting
script `bulk_rf_report.py` (see below) to offer two fairly well established means of filtering RF
quality, and the user has a choice of which ones to apply when producing the report. This deferred
approach to RF quality filtering is also intended to readily support new machine learning approaches
to RF quality classification.

The script for applying RF quality metrics is `rf_quality_filter.py`.
The following filters are applied directly and immediately in this script to remove spurious RFs:
1. Remove traces containing NaNs
2. Remove traces containing 1 or less sample in the trace
3. Remove traces whose number of samples is not as expected. The expected number is the statistical _mode_
   of the number of samples across all traces in the stream.

Then the following additional metrics are computed on each RF and stored in the RF trace metadata. See source code
documentation for additional details:
1. Signal-to-noise ratio
2. Spectral entropy
3. Max amplitude and RMS amplitude of the RF component, and the same from the Z-component stored in the R and T components
4. Base 10 logarithms of the max and RMS amplitude metrics
5. Spectral and histogram statistics
6. Ratios of statistics before and after onset
7. Statistical moments of the RF frequency spectrum
8. Coherence of signals in the frequency domain
9. Grouping by similarity (DBSCAN method)

As one can guess by the diversity of quality metrics computed, this is an area of active research.

To use script `rf_quality_filter.py`, the following parameters should be provided:
* H5 Input file name
* H5 Output file name

For additional options, see `rf_quality_filter.py --help`.


## Visualization of RFs

Automatic report generation on a full set of RF results is performed using the `bulk_rf_report.py` script.
This script is provided with an input H5 file containing quality-filtered RFs (e.g. the output from
`rf_quality_filter.py` script) as well as an output file name, and then generates a PDF report containing
a standard set of visualisations of the RFs:

* Pinwheel diagram showing RFs source directions around the full azimuthal range
* Stacked RFs for R-component
* Stacked RFs for T-component
* _H-k_ stacking plot based on the stacked R-components

Then using the `bulk_rf_report.py` script, there are RF quality filter options that can be applied.
Use of at least one of these options is highly recommended to remove spurious RFs.
See `bulk_rf_report.py --help` for quality filtering options.


## Analysing RFs using Vespagrams

Vespagrams are images of receiver functions ordered by distance or slowness. The character (shape) of the conversion
distribution can indicate if they are direct conversions, such as Pms, or its multiples, such as PpPms or PpSms.
Recommended reading about the application of vespagrams is *Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885*.

The tool to plot vespagrams is `rf_plot_vespagram.py` which takes as input a set of RFs stored in H5 format.
See `rf_plot_vespagram.py --help`. Currently the plotting is set from -5 to 20 seconds time window and
5 to 9 s/degrees in slowness. The filtering is bandpass from 0.03 to 0.5 Hz. The parameters of the data processing
are saved within the PDF file _document properties_.


## Plotting Station Maps and Ray Piercing Points

The legacy script `legacy/plot_map.py` draws a station map and piercing points of the receiver functions. This script
is slated to be refactored. This script generates plots with high-resolution topography and automatic adjustment of
the station label positions using the `https://github.com/Phlya/adjustText` library.

Note that general network and station plotting can be found in the `seismic.inventory` module.


## Extraction of RFs for earth model inversion

RFs can be exported from H5 file to special format for RF inversion using the `rf_inversion_export.py` script.
This script takes an input H5 file (which may contain RFs prior to quality filtering since `rf_inversion_export.py`
has some filtering support) and a specified output folder. In the output folder, the script outputs one `.dat`
file per channel per station containing the RF trace amplitude, with moveout, as a function of time at specific
sampling rate of 6.25 Hz.  These `.dat` files are then ready for ingestion into RF inversion Fortran code.


## Inversion

Experimental codes for RF inversion to 1D earth model.


## Models folder

This folder contains miscellaneous earth model files used for RF analysis, such as custom 1D velocity models.


## Legacy code

The `legacy` folder is for legacy code which is retained temporarily for reference purposes, and whose methods
may in future be rolled into the existing RF workflow.
