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
    
Use help on each command like this on how to use each step:

    $ cluster gather --help
    
    Usage: cluster gather [OPTIONS] EVENTS_DIR

    Options:
      -o, --output_file TEXT          output arrivals file basename
      -x, --nx INTEGER                number of segments from 0 to 360
                                      degrees for longitude
      -y, --ny INTEGER                number of segments from 0 to 180
                                      degrees for latitude
      -z, --dz FLOAT                  unit segment length of depth in
                                      meters
      -w, --wave_type [P S|Pn Sn|Pg Sg|p s]
                                      Wave type pair to generate
                                      inversion inputs
      --help                          Show this message and exit.


## Example workflow

    cluster gather tests/mocks/events/engdahl_sample
