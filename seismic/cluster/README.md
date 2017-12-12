## How to generate inputs for 3D Travel Time inversion

Install `passive-seismic` software.

This will create an executable `cluster`. Use `help` on how to use the 
`cluster` to generate input files for 3D travel time inversion.

    $ cluster --help
    Usage: cluster [OPTIONS] COMMAND [ARGS]...

    Options:
      -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
                                      Level of logging
      --help                          Show this message and exit.
    
    Commands:
      gather
      match   Match source and station blocks and output...
      sort    Sort based on the source and station block...

    
## Three steps to travel time input files

There are three steps to the input file generation in the following sequence:

    cluster gather
    cluster sort
    cluster match
    cluster zone
    
Use help on each command like this on how to use each step:

    $ cluster gather --help
    
    Usage: cluster gather [OPTIONS] EVENTS_DIR
    
      Gather all source-station block pairs for all events in a directory.
    
    Options:
      -o, --output_file TEXT          output arrivals file basename
      -x, --nx INTEGER                number of segments from 0 to 360 degrees for
                                      longitude
      -y, --ny INTEGER                number of segments from 0 to 180 degrees for
                                      latitude
      -z, --dz FLOAT                  unit segment length of depth in meters
      -w, --wave_type [P S|Pn Sn|Pg Sg|p s]
                                      Wave type pair to generate inversion inputs
      --help                          Show this message and exit.


## Example workflow
    
    # gather
    cluster gather /g/data/ha3/sudipta/event_xmls -w "P S"
    
    # sort
    cluster sort outfile_P.csv 5. -s sorted_P.csv
    cluster sort outfile_S.csv 10. -s sorted_S.csv

    # match
    cluster match sorted_P.csv sorted_S.csv -p matched_P.csv -s matched_S.csv

    # zones
    cluster zone '0 -50.0 100 190' matched_P.csv -r region_P.csv -g global_P.csv
    cluster zone '0 -50.0 100 190' matched_S.csv -r region_S.csv -g global_S.csv
    
## PBS scripts

When gather a large number of arrivals, batch jobs can be utilised making use
 of the steps above. Example PBS scripts can be found in [here for python3.4](../../hpc/cluster.sh) and 
[here for python3.6](../../hpc/cluster36.sh). Instructions on how to setup hpc 
environement on `raijn` like supercomputer can be found [here for python3.4](../../hpc/README.rst) and 
[here for python3.6](../../hpc/READMEPY36.sh).
