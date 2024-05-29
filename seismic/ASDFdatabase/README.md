# Data-curation workflows

The ASDFDatabase module features a number of data-curation workflows described below.

## Waveform-Analytics

`waveform_analytics.py` generates detailed reports on data coverage and quality. It 
can be run in two different modes: `mseed` and `asdf` -- the former requires a folder 
containing miniseed data, while the latter requires a text file containing paths 
to ASDF files. 

Usage details for `mseed` mode:

```commandline
python waveform_analytics.py mseed -h
Usage: waveform_analytics.py mseed [OPTIONS] MSEED_FOLDER MSEED_PATTERN
                                   INSTRUMENT_RESPONSE SAMPLING_RATE
                                   OUTPUT_FOLDER

  MSEED_FOLDER: Path to folder containing mseed files

  MSEED_PATTERN: File pattern to be used to capture files pertaining to
  specific channels.                Note that pattern must be specified
  within quotes.

  INSTRUMENT_RESPONSE: Path to inventory containing instrument response in
  StationXML or .resp format

  SAMPLING_RATE: Sampling rate used to record the mssed files OUTPUT_FOLDER:
  Path to output folder

Options:
  --start-date TEXT  Start date in UTC format for processing data
  --end-date TEXT    End date in UTC format for processing data
  --nproc INTEGER    Number of parallel processes use. Default is to use all
                     available cores.  [default: -1]

  -h, --help         Show this message and exit.

```

Usage details for `asdf` mode:

```commandline
python waveform_analytics.py asdf -h
Usage: waveform_analytics.py asdf [OPTIONS] ASDF_SOURCE NETWORK STATION
                                  LOCATION CHANNEL INSTRUMENT_RESPONSE
                                  SAMPLING_RATE OUTPUT_FOLDER

  ASDF_SOURCE: Path to text file containing paths to ASDF files

  NETWORK: network code 
  STATION: station code 
  CHANNEL: channel code
  INSTRUMENT_RESPONSE: Path to inventory containing instrument response in
  StationXML or .resp format

  SAMPLING_RATE: Sampling rate used to record the mssed files OUTPUT_FOLDER:
  Path to output folder

Options:
  --start-date TEXT  Start date in UTC format for processing data
  --end-date TEXT    End date in UTC format for processing data
  -h, --help         Show this message and exit.

```

