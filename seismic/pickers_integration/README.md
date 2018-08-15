## How to use the pick extraction code

1. Create an input folder
2. Copy over the input event catalog to the input folder
3. Create an output folder
4. Invoke the command below:

```
    $ python3 extract_picks.py --help
    Usage: extract_picks.py [OPTIONS] SC3_HOST YEAR
    
      This script accomplishes the task of picking the p-phase and s-phase
      arrival times for given set of events and waveform data sets that have
      been ingested as miniseed data into a Seiscomp3 (SC3 for short) server
      SDS archive.
    
    Arguments:
    
      SC3_HOST    : the SC3 FDSN server where the permanent station waveform
      data is hosted
    
      YEAR        : the year to which the input event catalog corresponds
    
    Options:
      --input_folder TEXT   The input folder where ISC event catalog xml input
                            files are kept.
      --output_folder TEXT  The output folder where generated xml files are
                            generated.
      --parallel            Needs to be run in conjunction with MPI.(DO NOT USE.
                            NOT TESTED.)
      --help                Show this message and exit.
```

## The algorithm for pick extraction is as below

```
1. Query the FDSN client to download the permanent station inventory to be picked
2. Iterate over all the input catalog events:
     Iterate over all stations of the downloaded inventory
       For this station event combination, get the theoretical travel time
       Query the FDSN server for the waveform window around (origin-time + theoretical travel time)
       Apply picking procedure on the returned waveform time-series data
       If pick found:
         Calculate SNR for the extracted pick
       return Pick
```
The output events generated from this will be in the output folder.
