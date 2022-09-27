# Body-wave Tools

## [Dependencies](../README.md#Dependencies)

## Workflow

The Body-wave workflow comprises the following steps:

1. Automatically pick waveforms for P- and S- arrivals based on a curated [catalogue](../catalogue)
2. Relocate local earthquakes and redefine phases of automatic arrivals utilizing Source-Specific-Station-Term corrections 
3. Cluster the output of 2 in order to balance uneven data-coverage

## 1. Automatic Arrival Picking

The MPI-parallel picker is launched via `pick.py`, which fetches waveforms from a FederatedASDF_Dataset source and carries out arrival detection. 
Run `pick.py` with `-h` for more details:
```
Usage: pick.py [OPTIONS] ASDF_SOURCE CSV_CATALOG_FILE OUTPUT_PATH

  ASDF_SOURCE: Text file containing a list of paths to ASDF files
  CSV_CATALOG_FILE: Path to catalog file in csv format

  OUTPUT_PATH: Output folder

Options:
  --min-magnitude INTEGER  Minimum magnitude of event; all events are analyzed
                           by default  [default: -1]

  --max-amplitude FLOAT    Maximum amplitude in waveforms analyzed; default is
                           1e8. This parameter is useful for filtering out
                           waveform snippets with spurious spikes  [default:
                           100000000.0]

  --network-list TEXT      A space-separated list of networks to process.
                           [default: *]

  --station-list TEXT      A space-separated list of stations to process.
                           [default: *]

  --restart                Restart job  [default: False]
  --save-quality-plots     Save plots of quality estimates  [default: (True)]
  -h, --help               Show this message and exit.
```
The output of `pick.py` comprises two text files: `p_combined.txt` and `s_combined.txt` which are then processed in step 2.

## 2. SSST Event Relocation and Phase Redefinition

The MPI-parallel SSST script launched via `ssst_relocate.py`. It takes the original catalog, used as input in step 1, and the outputs of step 2 and a config file as 
input. 
Run `ssst_relocate.py` with `-h` for more details:
```
Usage: ssst_relocate.py [OPTIONS] CATALOG_CSV CONFIG OUTPUT_FILE_NAME

  CATALOG_CSV: catalog in csv format CONFIG: config file in json format
  OUTPUT_FILE_NAME: name of output file

Options:
  --automatic-picks-p FILE   Automatic P picks in txt format (output of
                             pick.py).

  --automatic-picks-s FILE   Automatic S picks in txt format (output of
                             pick.py).

  --phases TEXT              A space-separated list of phases (case-sensitive
                             and within quotes) to process. Note that phases
                             not in the list are dropped.  [default: P Pg Pb
                             Pn S Sg Sb Sn]

  --p-residual-cutoff FLOAT  P-arrivals with abs(residual) exceeding this
                             cutoff are dropped.  [default: 5.0]

  --s-residual-cutoff FLOAT  S-arrivals with abs(residual) exceeding this
                             cutoff are dropped.  [default: 10.0]

  --min-slope-ratio INTEGER  Automatic arrivals with quality_measure_slope
                             less than this value are excluded from SST and
                             SSST computations and the event relocation
                             procedure -- the phases of such arrivals are
                             redefined nonetheless.  [default: 5]

  --ssst-niter INTEGER       Number of ssst iterations  [default: 5]
  --ball-radius-km INTEGER   Radius of sphere in km used for calculating ssst
                             corrections  [default: 55]

  -h, --help                 Show this message and exit.
```

The input json config file should contain three entries for each event-source: 1. whether events should be relocated or held fixed, 
2. whether preexisting arrivals in the original catalog should have their phases redefined or held fixed, 3. whether automatically picked 
arrivals should have their phases redefined or held fixed. A typical config file should look as follows:
```
{"ISC":{
    "events": "fixed",
    "preexisting_arrivals": "fixed",
    "automatic_arrivals": "redefine"
},
"USGS":{
    "events": "fixed",
    "preexisting_arrivals": "fixed",
    "automatic_arrivals": "redefine"
},
"EHB":{
    "events": "fixed",
    "preexisting_arrivals": "fixed",
    "automatic_arrivals": "redefine"
},
"GA":{
    "events": "relocate",
    "preexisting_arrivals": "redefine",
    "automatic_arrivals": "redefine"
},
"GG":{
    "events": "relocate",
    "preexisting_arrivals": "redefine",
    "automatic_arrivals": "redefine"
}}
```
The output of `ssst_relocate.py` comprises: 1. a PDF file containing detailed plots for the evolution of the distribution of travel-time residuals 
through all SSST-iterations, 2. An HDF5 file with complete snapshots of SSST-results at each iteration. The PDF report forms the basis for choosing 
appropriate parameters for rerunning the workflow for fine-tuning results.

## 3. Cluster Arrivals to Balance Uneven Data-Coverage

The clustering procedure is launched via `cluster_arrivals.py`. It takes the output of step 2, the HDF5 file, and a grid config file delineating the 
region to be clustered as input. The purpose of this script is to: i) define a coarse, global 3D mesh, embedded within which is a finer regional mesh 
ii) take all rays emanating from a source-cell and arriving at a station-cell and coalesce them into a single ray, identified by the median travel-time.

Run `cluster_arrivals.py` with `-h` for more details:
```
Usage: cluster_arrivals.py [OPTIONS] INPUT_H5 PARAM_FN OUTPUT_FILE_NAME

  INPUT_H5: hdf5 input (output of ssst_relocate.py) PARAM_FN: Grid
  parameterization file OUTPUT_FILE_NAME: name of output file

Options:
  --phases TEXT                   A space-separated list of phases (case-
                                  sensitive and within quotes) to cluster.
                                  Note that phases not in the list are
                                  dropped.  [default: P Pg Pb S Sg Sb]

  --min-slope-ratio INTEGER       Automatic arrivals with
                                  quality_measure_slope less than this value
                                  are discarded prior to the clustering step
                                  [default: 5]

  --p-residual-cutoff FLOAT       P-arrivals with abs(residual) exceeding this
                                  cutoff are dropped.  [default: 5.0]

  --s-residual-cutoff FLOAT       S-arrivals with abs(residual) exceeding this
                                  cutoff are dropped.  [default: 10.0]

  --outer-depth-refinement-factor INTEGER
                                  Depth-refinement factor for outer grid. Must
                                  be a power of 2.  [default: 1]

  --inner-depth-refinement-factor INTEGER
                                  Depth-refinement factor for inner grid. Must
                                  be a power of 2.  [default: 1]

  --only-paired-s                 S arrivals not accompanied by corresponding
                                  P arrivals are dropped before clustering.
                                  [default: False]

  -h, --help                      Show this message and exit.
```

The grid config file has the following format:
```
global_ncx global_ncy global_ncd dx dy
d1
d2
.
.
d_{global_ncd+1}
lon_min lon_max lat_min lat_max regional_ncx regional_ncy regional_ncd
d1
d2
.
.
d_{regional_ncd+1}
```
The first six lines above describe a coarse global mesh as follows:

`global_ncx`: number of longitudinal cells

`global_ncy`: number of latitudinal cells

`global_ncd`: number of cells in depth

`dx`: longitudinal resolution (`dx ≡ 360 / global_ncx`)

`dy`: latitudinal resolution (`dy ≡ 180 / global_ncy`)

`d1, d2, .. d_{global_ncd+1}`: a total of `global_ncd+1` depth nodes starting from `d1=0`, in km

The remaining six lines describe a regional mesh as follows:
`lon_min lon_max`: minimum/maximum longitude of regional mesh

`lat_min lat_max`: minimum/maximum latitude of regional mesh

`regional_ncx`: number of longitudinal cells

`regional_ncy`: number of latitudinal cells

`regional_ncd`: number of cells in depth

`d1, d2, .. d_{regional_ncd+1}`: a total of `regional_ncd+1` depth nodes starting from `d1=0`, in km

A typical grid config file for the Australasian region is as follows:
```
72 36 16 5. 5.
0.
110.
280.
410.
660.
840.
1020.
1250.
1400.
1600.
1850.
2050.
2250.
2450.
2600.
2750.
2889.
100. 158. -54. -9. 58 45 29
0.
10.
30.
50.
75.
100.
125.
150.
175.
200.
225.
250.
275.
300.
330.
370.
410.
460.
510.
560.
610.
660.
710.
810.
910.
1010.
1110.
1250.
1400.
1600.
```
Output files from the clustering step have the following columns:

`source_block`: cell id of earthquake hypocenter. Note that cells in the regional mesh are numbered [1 - NR], where NR is the total
number of cells in the regional mesh. Cells in the global mesh are numbered [NR+1 - NR+NG], where NG is the total number of 
cells in the global mesh

`station_block`: cell id of recording station

`residual`: travel-time residual (s)

`eorigin_ts`: event origin timestamp

`elon`: event longitude

`elat`: event latitude

`edepth_km`: event depth (km)

`slon`: station longitude

`slat`: station latitude

`selev_km`: station elevation (km)

`observed_tt`: observed travel-time (s)

`ecdist`: epicentral distance (degrees)

`phase`: seismic phase

`phase_type`: 1 for P and 2 for S
