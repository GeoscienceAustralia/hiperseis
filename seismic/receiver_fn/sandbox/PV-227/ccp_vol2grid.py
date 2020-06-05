"""
Blurb

Using cartesian geometry as an approximation to the curved surface of Earth.
"""

import os

import numpy as np
from scipy.spatial import cKDTree

from seismic.units_utils import KM_PER_DEG


infile = 'ccp_line_data_sample.csv'
data_vol = np.loadtxt(infile, delimiter=',')
xy = data_vol[:, 0:2]
tree = cKDTree(xy)
# Keep z in column vector format
z = data_vol[:, 2][:, np.newaxis]

bb_min = tree.mins
bb_max = tree.maxes
span = bb_max - bb_min
spacing_km = 10.0
spacing_deg = spacing_km/KM_PER_DEG
x_n = int(np.ceil(span[0]/spacing_deg))
y_n = int(np.ceil(span[1]/spacing_deg))
x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], x_n),
                             np.linspace(bb_min[1], bb_max[1], y_n))
grid = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))

# Proximity weighting based on http://dx.doi.org/10.1098/rspa.2019.0352 with only
# one dataset
sigma = 25.0/KM_PER_DEG
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

data_vol_gridded = np.hstack((grid, z_interp))
outfile = os.path.splitext(infile)[0] + '_gridded.csv'
np.savetxt(outfile, data_vol_gridded, fmt=['%.6f', '%.6f', '%.1f'], delimiter=',',
           header='Lon,Lat,Depth')
