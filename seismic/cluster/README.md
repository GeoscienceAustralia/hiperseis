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


### `Cluster gather`
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

### `cluster sort`

    $ cluster sort --help
    
    Usage: cluster sort [OPTIONS] OUTPUT_FILE RESIDUAL_CUTOFF
    
      Sort and filter the arrivals.
    
      Sort based on the source and station block number. There are two stages of
      filtering: 1. Filter based on the time residual 2. Filter based on median
      of observed travel time.
    
      If there are multiple source and station block combinations, we keep the
      row corresponding to the median observed travel time (observed_tt).
    
      :param output_file: output file from the gather stage :param sorted_file:
      str, optional     optional sorted output file path. Default: sorted.csv.
      :param residual_cutoff: float     residual seconds above which arrivals
      are rejected. :return: None
    
    Options:
      -s, --sorted_file FILENAME  output sorted and filter file.
      --help                      Show this message and exit.

### `cluster match`

    $ cluster match --help
    Usage: cluster match [OPTIONS] P_FILE S_FILE
    
      Match source and station blocks and output files with matched source and
      station blocks.
    
      :param p_file: str     path to sorted P arrivals :param s_file: str
      path to sorted S arrivals :param matched_p_file: str, optional     output
      p arrivals file with matched p and s arrivals source-block
      combinations :param matched_s_file: str, optional     output s arrivals
      file with matched p and s arrivals source-block     combinations
    
      :return:None
    
    Options:
      -p, --matched_p_file FILENAME  output matched p file.
      -s, --matched_s_file FILENAME  output matched s file.
      --help                         Show this message and exit.
 

### `cluster zone`
    $ cluster zone --help
    Usage: cluster zone [OPTIONS] str 'upperlat, bottomlat, leftlon, rightlon'
                        cluster_matched_file
    
      `zone'ing the arrivals into three regions.
    
    Options:
      -r, --region_file FILENAME      region file name.
      -g, --global_file FILENAME      global file name.
      -s, --grid_size FLOAT           grid size in degrees in the region. If grid
                                      size is provided, cross region file will be
                                      created.
      -c, --cross_region_file FILENAME
                                      cross region file name.
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

When gathering a large number of arrivals, batch jobs can be utilised making use
 of the steps above. Example PBS scripts can be found [here for python3.4](../../hpc/cluster.sh) and 
[here for python3.6](../../hpc/cluster36.sh). Instructions on how to setup hpc 
environment on `raijin` like supercomputer can be found [here for python3.4](../../hpc/README.rst) and 
[here for python3.6](../../hpc/READMEPY36.sh).
