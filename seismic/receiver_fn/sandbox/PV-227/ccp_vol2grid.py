"""
Blurb

Using cartesian geometry as an approximation to the curved surface of Earth.

Requires:
- pyepsg
"""

import os

import numpy as np
from scipy.spatial import cKDTree
import cartopy as cp

# from seismic.units_utils import KM_PER_DEG

# Load point data whose x-y coordinates are (longitude, latitude)
infile = 'ccp_line_data_sample.csv'
data_vol = np.loadtxt(infile, delimiter=',')
# z here is any scalar spatial data metric (e.g. depth to Moho).
# Keep z in column vector format.
z = data_vol[:, 2][:, np.newaxis]
# xy coordinates in lon/lat
xy = data_vol[:, 0:2]

# Convert lon/lat to cartesian so we can compute distances between stations on equal basis

# Create KDTree using projected coordinate system
EPSG_CODE = 3577  # Would prefer 7845 if it worked, since it is based on GDA2020.
map_proj = cp.crs.epsg(EPSG_CODE)  # Requires internet connection
wgs_84 = cp.crs.Geodetic()
xy_map = map_proj.transform_points(wgs_84, xy[:, 0], xy[:, 1])
xy_map = xy_map[:, :2]  # Ignore z. xy here are in metres.
tree = cKDTree(xy_map)

bb_min = tree.mins
bb_max = tree.maxes
span = bb_max - bb_min
spacing_m = 10.0*1000  # km to metres
# Produce linear spacing in cartesian projection
x_n = int(np.ceil(span[0]/spacing_m))
y_n = int(np.ceil(span[1]/spacing_m))
x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], x_n),
                             np.linspace(bb_min[1], bb_max[1], y_n))
grid = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))

# Proximity weighting based on http://dx.doi.org/10.1098/rspa.2019.0352 with only
# one dataset
sigma = 25.0*1000  # km to metres
sig_sq = np.square(sigma)
max_dist = 3.5*sigma
dist_matrix = tree.sparse_distance_matrix(cKDTree(grid), max_distance=max_dist).tocsr()
dist_matrix = dist_matrix.transpose()
dist_sq_matrix = -dist_matrix.power(2)/sig_sq  # Note the sign negation

# TRICKY point: We compute the exponential terms as exp(x) - 1 rather than exp(x),
# so that we can remain in sparse matrix format until the matrix multiplication.
# After that, we will compensate for the missing +1 term in the sum.
dist_weighting_m1 = dist_sq_matrix.expm1()
mask = (dist_matrix != 0)
denom = dist_weighting_m1.sum(axis=1) + mask.sum(axis=1)
numer = (dist_weighting_m1 + mask)*z
z_interp = numer/denom
# Compute pointwise variance
numer = (dist_weighting_m1 + mask)*np.square(z)
S = np.sqrt(numer/denom - np.square(z_interp))

# Convert back from local projection to global lon/lat coordinates so that downstream
# consumers can use whatever coordinate systems and projections they want.
grid_lonlat = wgs_84.transform_points(map_proj, grid[:,0], grid[:,1])[:, :2]

# Collect data and write to file
data_vol_gridded = np.hstack((grid_lonlat, z_interp, S))
outfile = os.path.splitext(infile)[0] + '_gridded_epsg3577.csv'
np.savetxt(outfile, data_vol_gridded, fmt=['%.6f', '%.6f', '%.1f', '%.1f'], delimiter=',',
           header='Lon,Lat,Depth,Stddev')

# See https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html
# for how to use projection and transform settings with matplotlib.
# Since this script outputs in lon/lat, it is in CRS WGS84 which is default for cartopy.crs.Geodetic()
