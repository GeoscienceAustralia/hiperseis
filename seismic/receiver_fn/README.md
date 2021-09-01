# Receiver Function (RF) tools

## Dependencies
The workflow requires MPI (mpi4py) and parallel HDF5 (h5py) capabilities on the NCI.
Installation instructions for NCI (Gadi ) are as follows:

### Load system modules:
  1. `module purge` highly recommeneded to purge all modules to begin with
  2. `module load pbs` 
  3. `module load python3-as-python`
  4. `module load openmpi/2.1.6-mt`
  5. `module load hdf5/1.10.5p`

### Upgrade pip

 1. `pip3.6 install pip==21.1.2 --user`

### numpy

 1. `pip3.6 install numpy==1.18.5 --user`

### Install mpi4py that uses the correct OpenMPI libs
  1. `MPICC=/apps/openmpi/2.1.6-mt/bin/mpicc pip3 install mpi4py --user` Note that we use `pip3`, the system-provided pip for python 3.6

### Build Parallel H5PY

1. `pip3.6 install cython==0.29.22 --user`
2. `git clone --single-branch --branch 2.10.0.gadi_tweaks https://github.com/rh-downunder/h5py.git` Pull a branch (based on version 2.10.0) from h5py repository from github fork of h5py, adapted for Gadi, for purpose of custom build
3. `cd h5py`
4. `CC=mpicc python setup.py configure --mpi --hdf5=/apps/hdf5/1.10.5p/` Configure with mpi enabled
5. `python setup.py build` Build h5py
6. `python setup.py install --user` Install in user space

### Setup standard packages
  1. `pip3.6 install obspy==1.1.0 --user`
  2. `pip3.6 install click==7.1.2 --user `
  3. `pip3.6 install netCDF4==1.4.0 --user`
  4. `pip3.6 install pyasdf==0.5.1 --user`
  5. `pip3.6 install ordered_set ujson psutil --user`

### Initialize PhasePapy
Additionally, the workflow requires the PhasePapy package (a submodule of this repository) to be
initialized, if not already.
```
cd hiperseis
git submodule update --init --recursive
```

## Workflow

The RF workflow is based on the following main steps:

 - Data preparation: extract raw waveforms in time windows associated with selected seismic
   events from ASDF database for the deployment
 - Calculate Receiver Functions: rotation (ZRT/QLT) and deconvolution
 - Quality filtering on RFs: discard bad results by filtering
 - Perform analysis such as 1D Earth structure inversion, H-k stacking, etc.

There are a number of manuals and literature about RF. You can refer to examples provided by
[Charles Ammon](http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/rftn01.html), manuals
distributed with [python RF libraries](https://rf.readthedocs.io/en/latest/) or [chapter 4.1 of
Tom Richter's dissertation](http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf).

The work in this module utilizes the RF framework developed in the Python [`rf` library
developed by Tom Eulenfeld](https://doi.org/10.21105/joss.01808).


## Event Trace Extraction

The first step in data preparation is to extract relevant waveform data for a set of 
registered teleseismic earthquakes obtained from the ISC web service.
The main program to extract raw waveforms from an FDSN web-service or the 
FederatedASDFDataset infrastructure is `extract_event_traces.py`. The
following parameters must be supplied to the script:

 - Waveform source (URL to an FDSN web-service or path to a text file to initialize a \
   FederatedASDFDataset)
 - Event catalog to create or re-use
 - Output file name

All waveforms are processed through an automatic (AIC-based) picker and the following metadata 
are added to all traces where a pick is found:

```
'arrival_time' # actual arrival time of a P-wave
'slope_ratio' # quality metric for arrival -- typically, values>5 indicate \
                robust arrivals
```
The P-arrival quality metric (`slope_ratio`) can be particularly useful for filtering out 
low-quality waveforms in downstream processing steps.

For usage details, see output of `python extract_event_traces.py --help` below:

```
Usage: extract_event_traces.py [OPTIONS]

Options:
  --inventory-file FILE           Optional path to input inventory file
                                  corresponding to waveform source provided
                                  through, --waveform-database. Note that this
                                  parameter is required only when the waveform
                                  source is not a definition file for a
                                  FederatedASDFDataSet, in which case, the
                                  relevant inventory is extracted internally.

  --network-list TEXT             A space-separated list of networks (within
                                  quotes) to process.  [default: *]

  --station-list TEXT             A space-separated list of stations (within
                                  quotes) to process.  [default: *]

  --waveform-database TEXT        Location of waveform source database from
                                  which to extract traces. May be a recognized
                                  service provider from
                                  obspy.clients.fdsn.header.URL_MAPPINGS (e.g.
                                  "ISC"), an actual URL (e.g.
                                  "http://auspass.edu.au") or a file path. If
                                  detected as a URL, the obspy client
                                  get_waveform function will be used to
                                  retrieve waveforms from web service.
                                  Otherwise, if detected as a valid file path,
                                  then it must be the path to a definition
                                  file for a FederatedASDFDataSet, e.g. "/g/da
                                  ta/ha3/Passive/SHARED_DATA/Index/asdf_files.
                                  txt".  [required]

  --event-catalog-file FILE       Path to event catalog file, e.g.
                                  "catalog_7X_for_rf.xml". If file already
                                  exists, it will be loaded, otherwise it will
                                  be created by querying the ISC web service.
                                  Note that for traceability, start and end
                                  times will be appended to file name.
                                  [required]

  --event-trace-datafile FILE     Path to output file, e.g.
                                  "7X_event_waveforms.h5". Note that for
                                  traceability, start and end datetimes will
                                  be appended to file name.  [required]

  --start-time TEXT               Start datetime in ISO 8601 format, e.g.
                                  "2009-06-16T03:42:00". If empty, will be
                                  inferred from the inventory file.  [default:
                                  ]

  --end-time TEXT                 End datetime in ISO 8601 format, e.g.
                                  "2011-04-01T23:18:49". If empty, will be
                                  inferred from the inventory file.  [default:
                                  ]

  --taup-model TEXT               Theoretical tau-p Earth model to use for
                                  Trace stats computation. Other
                                  possibilities, such as ak135, are documented
                                  here: https://docs.obspy.org/packages/obspy.
                                  taup.html  [default: iasp91]

  --distance-range <FLOAT FLOAT>...
                                  Range of teleseismic distances (in degrees)
                                  to sample relative to the mean lat,lon
                                  location  [default: 30.0, 90.0]

  --magnitude-range <FLOAT FLOAT>...
                                  Range of seismic event magnitudes to sample
                                  from the event catalog.  [default: 5.5, 7.0]

  --catalog-only                  If set, only generate catalog file and exit.
                                  Used for preparing input file on HPC systems
                                  with no internet access.  [default: False]

  --help                          Show this message and exit.
```

### Typical Usage:

Download events for OA stations:
```
python extract_event_traces.py --waveform-database asdf_files.txt \
--start-time 2011-01-01T00:00:00 --end-time 2011-01-05T00:00:00 \
--event-catalog-file OA_events.xml --event-trace-datafile OA_waveforms.h5 \
--network-list "OA" --catalog-only
```

Generate OA waveforms:
```
python extract_event_traces.py --waveform-database asdf_files.txt \
--event-catalog-file OA_events.xml --event-trace-datafile OA_waveforms.h5 \
--network-list "OA"
```

## Receiver Function Calculation / Correction

The `generate_rf.py` script reads the event traces  produced by `extract_event_traces.py` 
and generates receiver functions. The script parallelizes RF computation over stations 
through MPI.  

For usage details, see output of `python generate_rf.py --help` below:

```
Usage: generate_rf.py [OPTIONS] INPUT_FILE OUTPUT_FILE
  
  INPUT_FILE : Input waveforms in H5 format
               (output of extract_event_traces.py)

  OUTPUT_FILE : Output H5 file name

Options:
  --network-list TEXT  A space-separated list of networks (within quotes) to
                       process.  [default: *]

  --station-list TEXT  A space-separated list of stations (within quotes) to
                       process.  [default: *]

  --config-file FILE   Run configuration file in JSON format
  --only-corrections   Compute and apply corrections for stations listed under
                       'correction' in the input json config file -- all other
                       stations are ignored. Note that preexisting data (for
                       relevant channels, if present) are deleted before
                       saving the corrections  [default: False]

  --help               Show this message and exit.
```

Run configuration files consist of 3 sub-dictionaries. One named `filtering` for
input stream filtering settings, one named `processing` for RF processing
settings, and one named `correction` for rotating / swapping / negating channel
data for one or more named stations with potential orientation discrepancies.
Each of these sub-dicts is described below:

```
        "filtering":  # Filtering settings
        {
          "resample_rate": float # Resampling rate in Hz
          "taper_limit": float   # Fraction of signal to taper at end, between 0 and 0.5
          "filter_band": (float, float) # Filter pass band (Hz). Not required for freq-domain deconvolution.
          "channel_pattern": # Ordered list of preferred channels, e.g. 'HH*,BH*',
                             # where channel selection is ambiguous.
          "baz_range": (float, float) or [(float, float), ...] # Discrete ranges of source back azimuth to use (degrees).
              # Each value must be between 0 and 360. May be a pair or a list of pairs for multiple ranges.
        }

        "processing":  # RF processing settings
        {
          "custom_preproc":
          {
            "import": 'import custom symbols',  # statement to import required symbols
            "func": 'preproc functor'  # expression to get handle to custom preprocessing functor
            "args": {}  # additional kwargs to pass to func
          }
          "trim_start_time": float # Trace trim start time in sec, relative to onset
          "trim_end_time": float # Trace trim end time in sec, relative to onset
          "rotation_type": str # Choice of ['zrt', 'lqt']. Rotational coordinate system
                               # for aligning ZNE trace components with incident wave direction
          "deconv_domain": str # Choice of ['time', 'freq', 'iter']. Whether to perform deconvolution
                               # in time or freq domain, or iterative technique
          "gauss_width": float # Gaussian freq domain filter width. Only required for freq-domain deconvolution
          "water_level": float # Water-level for freq domain spectrum. Only required for freq-domain deconvolution
          "spiking": float # Spiking factor (noise suppression), only required for time-domain deconvolution
          "normalize": bool # Whether to normalize RF amplitude
        }

        "correction": # corrections to be applied to data for named stations prior to RF computation
        {
          "plot_dir": str # path to folder where plots related to cahnnel rotations are to be saved
          "swap_ne": list # list of NET.STA.LOC for which N and E channels are to be swapped, e.g ["OA.BL27."],
          "rotate": list # list of NET.STA.LOC that are to be rotated to maximize P-arrival energy on \
                           the primary RF component, e.g ["OA.BL27."]
          "negate": list # list of NET.STA.LOC.CHA that are to be negated, e.g ["OA.BL27..HHZ"]
        }       
```
Default values for parameters in "filtering" and "processing", above, are drawn from the list below:

```
     DEFAULT_RESAMPLE_RATE_HZ = 20.0
     DEFAULT_FILTER_BAND_HZ = (0.02, 1.00)
     DEFAULT_TAPER_LIMIT = 0.05
     DEFAULT_TRIM_START_TIME_SEC = -50.0
     DEFAULT_TRIM_END_TIME_SEC = 150.0
     DEFAULT_ROTATION_TYPE = 'zrt'   # from ['zrt', 'lqt']
     DEFAULT_DECONV_DOMAIN = 'time'  # from ['time', 'freq', 'iter']
     DEFAULT_GAUSS_WIDTH = 1.0
     DEFAULT_WATER_LEVEL = 0.01
     DEFAULT_SPIKING = 0.5
     RAW_RESAMPLE_RATE_HZ = 20.0
     BANDPASS_FILTER_ORDER = 2        
```

The default parameters are somewhat tuned to BH* instrument response waveforms, and should be 
explored during a RF study. 

During RF generation, some problematic traces are removed that would be
fatal for RF generation, as follows:
1. Remove traces with invalid inclination value in the metadata (e.g. NaN).
2. Remove traces in which not all ZNE channels are present for a given event.
3. Remove traces in which the number of samples in the Z, N and E channels are not the same.
4. Remove traces in which NaN is found in the data.
5. Remove traces with zero variance in the data.
6. Remove traces for which an exception is caught during RF generation.
7. Remove traces with empty channel data or channels missing after deconvolution.

### Typical Usage:

Compute RFs in parallel:

```
mpirun -np NUM_PROCS python generate_rf.py OA_waveforms.h5  OA_rf.h5 \
--config-file config_rfs.json
```

config_rfs.json:
```
{
  "filtering": {
    "resample_rate": 10.0,
    "taper_limit": 0.05,
    "filter_band": [0.02, 1.0]
  },
  "processing": {
      "rotation_type": "ZRT",
      "deconv_domain": "iter",
      "normalize": true
  }
}
```

#### Correcting RFs

A small percentage of stations in some networks (e.g. OA) suffer from channel orientation 
problems. Those stations can be corrected through a combination of channel rotation, negation 
or swapping (N/E). Typically, one would generate the RFs and visualize them as described in 
[Visualization of RFs](#visualization-of-rfs) -- note that the
[RF Quality Filtering](#rf-quality-filtering) step can be skipped for the corrections and 
completed afterwards. A pdf report generated as described in
[Visualization of RFs](#visualization-of-rfs) helps identify stations with orientation problems.
Once a list of problematic stations is compiled, users need to run `generate_rf.py` with the same
input/output parameters, but with a `correction` block added to the json configuration file, as 
follows:

config_rfs.json:
```
{
  "filtering": {
    "resample_rate": 10.0,
    "taper_limit": 0.05,
    "filter_band": [0.02, 1.0]
  },
  "processing": {
      "rotation_type": "ZRT",
      "deconv_domain": "time",
      "normalize": true
  }
  "correction":{
      "plot_dir":"path/to/plots",
      "swap_ne":["OA.BL27.", "OA.BM31."]
      "rotate":["OA.BL27.", "OA.BM31."]
      "negate":["OA.CA21.", "OA.CB31."]
  }
}
```
Note that the correction step should be run with the `--only-corrections` flag, as follows:

```
mpirun -np NUM_PROCS python generate_rf.py OA_waveforms.h5  OA_rf.h5 \
--config-file config_rfs.json --only-corrections
```
The above will ensure that only those stations in the `correction` block are processed. Plots 
showing effects of channel rotation (if specified) are saved in the `path/to/plots` folder, 
along with a json file listing corrections to back-azimuths. Note that when back-azimuths 
are updated, the original back-azimuths are stored in the trace-header as `orig_back_azimuth`.

Correcting orientation problems is an iterative process that may require different combinations 
of rotation, channel swapping and negation for problematic stations -- consequently, the 
`correction` block may need to be updated, RFs generated and visualized, iteratively, to correct 
all problematic stations. To speed up the correction process, one can use time-domain
deconvolution, which is much faster than iterative deconvolution -- once optimal corrections 
are found for all problematic stations, one can then revert to iterative deconvolution 
to regenerate RFs for problematic stations.

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

Finally, when saving the traces with quality metrics in the metadata, the Z-component is dropped to
reduce file size, so only the R- and T-components are saved.  A diversity of
quality metrics are computed because this is an area of active research.

To use script `rf_quality_filter.py`, the following parameters should be provided:
* H5 Input file name
* H5 Output file name

For additional options, see `rf_quality_filter.py --help`.


## Visualization of RFs

Automatic report generation on a full set of RF results is performed using the MPI-parallelized 
`bulk_rf_report.py` script.
This script needs to be provided with an input H5 file containing RFs generated in 
[Receiver Function Calculation](#receiver-function-calculation) or that in 
[RF Quality Filtering](#rf-quality-filtering) and an output file name -- it then generates 
a PDF report containing a standard set of visualisations of the RFs:

* Pinwheel diagram showing RFs source directions around the full azimuthal range
* Stacked RFs for R-component
* Stacked RFs for T-component
* H-k stacking plot based on the stacked R-components
* Estimates of sediment thickness if de-reverberation was applied on a given station

Quality filters can be used to filter out spurious RFs, as shown in usage details from 
`bulk_rf_report.py --help` below:

```
Usage: bulk_rf_report.py [OPTIONS] INPUT_FILE OUTPUT_FILE

  INPUT_FILE : Input RFs in H5 format
               (output of generate_rf.py or rf_quality_filter.py)

  OUTPUT_FILE : Output pdf file name

Options:
  --network-list TEXT             A space-separated list of networks (within
                                  quotes) to process.  [default: *]

  --station-list TEXT             A space-separated list of stations (within
                                  quotes) to process.  [default: *]

  --event-mask-folder DIRECTORY   Folder containing event masks to use to
                                  filter traces. Such masks are generated
                                  using rf_handpick_tool

  --apply-amplitude-filter        Apply RF amplitude filtering to the RFs. The
                                  default filtering logic includes: Signal SNR
                                  >= 2.0 RMS amplitude of signal < 0.2 Maximum
                                  amplitude of signal < 1.0  [default: False]

  --apply-similarity-filter       Apply RF similarity filtering to the RFs.
                                  [default: False]

  --min-slope-ratio FLOAT         Apply filtering to the RFs based on the
                                  "slope_ratio" metric that indicates
                                  robustnessof P-arrival. Typically, a minimum
                                  slope-ratio of 5 is able to pick out strong
                                  arrivals. The default value of -1 does not
                                  apply this filter  [default: -1]

  --hk-weights <FLOAT FLOAT FLOAT>...
                                  Weightings per arrival multiple for H-k
                                  stacking  [default: 0.5, 0.4, 0.1]

  --hk-solution-labels [global|local|none]
                                  Method of labeling automatically selected
                                  solutions on H-k stack plots. global: find
                                  and label global maximum, local: find and
                                  label up to 3 local maxima after clustering,
                                  none: do not label any solutions on H-k
                                  stack plot.  [default: global]

  --depth-colour-range <FLOAT FLOAT>...
                                  The range of depth values from which to
                                  choose the maximum hk_stack value for
                                  plotting purposes. Note that this parameter
                                  has no effect on the computation of the
                                  hk_stack.  [default: 20, 70]

  --hk-hpf-freq FLOAT             If present, cutoff frequency for high pass
                                  filter to use prior to generating H-k
                                  stacking plot.

  --help                          Show this message and exit.
```

Note that, if `bulk_rf_report.py` is run on RFs produced in 
[Receiver Function Calculation](#receiver-function-calculation), `--min-slope-ratio` is 
the only quality-filter available. Otherwise, if run on the output generated in 
[RF Quality Filtering](#rf-quality-filtering), all quality-filters described above 
are available for use. Use of at least one of these options is highly recommended to 
remove spurious RFs.

### Typical Usage
Generate pdf report in parallel:

```
mpirun -np NUM_PROCS python bulk_rf_report.py OA_rf.h5 OA_rf_report.pdf \
--hk-hpf-freq 0.1 --min-slope-ratio 5 --station-list 'CJ27 BV22'
```

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


## Models folder

This folder contains miscellaneous earth model files used for RF analysis, such as custom 1D velocity models.


## Legacy code

The `legacy` folder is for legacy code which is retained temporarily for reference purposes, and whose methods
may in future be rolled into the existing RF workflow.
