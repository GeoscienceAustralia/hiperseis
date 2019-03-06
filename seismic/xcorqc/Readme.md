# Cross-Correlation Workflow

The cross-correlation functionality now works with ASDF files -- for voluminous temporary stations data, a SeisDS object (based on a JSON database) can be utilized to speed up data access. An example, 7G_example.py, shows how data from temporary stations and reference stations can be cross-correlated over a 6 month period, stacked over 10 day intervals. Cross-correlation results are written out as NetCDF-4 files for each station-pair.

Interrogating cross-correlation results requires interactive visualization capabilities. [Panoply], a freely available cross-platform tool, which is also available on NCI VDIs, can be used to visualize the results interactively. A quick intro to [Panoply] is available [here].

## Fetching Permanent Stations Data

File `client_data.py` can be used to fetch permanent stations data, for a given time-range, at given location (lat, lon). Permanent station locations for AU can be found on the [FDSN website].

## Setting up NCI Raijin environment for running cross-correlation

Things to address
- who these instructions are for
- gotchas and dependencies in setting up
- the requirements of the x-correlation code
- the limitations on raijin python libraries
  - Python 2 only
- uses MPI for massively parallel execution
- requires custom build of h5py library
- requires use of 'modules' installed on raijin
  - modules interact with Python third party libraries and can modify user PATH and PYTHONPATH variables
- Python environment setup requires third party libraries to be installed
  - custom build of h5py
  - use system `pip`, DO NOT upgrade to pip version 19
- order of operations is critical
- in general virtual environments are preferred, but not well supported in py2.7 and this
  solution uses the python user space (`pip --user`), but not virtual environments for the
  sake of simplicity
- this is may not be the only solution, but it is currently the only known working solution
- list of libraries required: h5py, obspy, pyasdf, click, netCDF4 and their dependencies

[Panoply]:https://www.giss.nasa.gov/tools/panoply/
[here]:http://www.meteor.iastate.edu/classes/mt452/EdGCM/Documentation/EdGCM_Panoply.pdf
[FDSN website]:http://www.fdsn.org/networks/detail/AU/

