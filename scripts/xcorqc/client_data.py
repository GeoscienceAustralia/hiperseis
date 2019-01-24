from seismic.ASDFdatabase.ClientUtils import Client2ASDF
import os

# =========================== User Input Required =========================== #

# Output filename
fn = '/g/data/ha3/rakib/_ANU/7G(2013-2015)/refData/stka.6m.h5'

# Station Location (lat, lon)
# Station names and their coordinates can be found in:
# http://www.fdsn.org/networks/detail/AU/
loc = (-31.8769, 141.5952) # e.g. inka = (-27.741, 140.746), stka = (-31.8769, 141.5952)

# Buffer around location in degrees. Note that a larger buffer may net additional
# nearby stations
buffer = 0.0001

# Time-range (start, end)
timeRange = ("2015-01-01T00:00:00", "2015-06-01T00:00:00")

# =========================== End of User Input =============================#



# =========================== Begin Processing ============================= #

c = Client2ASDF()

c.queryByBBoxInterval(fn, [loc[1] - buffer, loc[1] + buffer,
                           loc[0] - buffer, loc[0] + buffer],
                      timeRange, 1, verbose=True)

# test that the file exists and was written
assert os.path.isfile(fn), "ASDF file not written"
