# Cross-Correlation Workflow

The cross-correlation functionality now works with ASDF files -- for voluminous temporary stations data, a SeisDS object (based on a JSON database) can be utilized to speed up data access. An example, 7G_example.py, shows how data from temporary stations and reference stations can be cross-correlated over a 6 month period, stacked over 10 day intervals. Cross-correlation results are written out as NetCDF-4 files for each station-pair.

Interrogating cross-correlation results requires interactive visualization capabilities. [Panoply], a freely available cross-platform tool, which is also available on NCI VDIs, can be used to visualize the results interactively. A quick intro to [Panoply] is available [here].

## Fetching Permanent Stations Data

Client_data.py can be used to fetch permanent stations data, for a given time-range, at given location (lat, lon). Permanent station locations for AU can be found on the [FDSN website].

[Panoply]:https://www.giss.nasa.gov/tools/panoply/
[here]:http://www.meteor.iastate.edu/classes/mt452/EdGCM/Documentation/EdGCM_Panoply.pdf
[FDSN website]:http://www.fdsn.org/networks/detail/AU/

