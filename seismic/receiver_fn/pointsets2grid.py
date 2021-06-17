#!/usr/bin/env python
"""
Blend data from multiple methods together onto a common grid.

Multiple methods with per-method settings are passed in using a JSON 
configuration file.

See `config_moho_workflow_example.json` and the wiki page at
https://github.com/GeoscienceAustralia/hiperseis/wiki/Blending-and-Plotting-Point-Datasets
for an explanation of the config.

Bounds are optional. If provided, the interpolation grid will be limited
to this extent. If not provided, the interpolation grid is bounded by
the maximum extent of the aggregate datasets.

The gridded data will be written to output directory as 'moho_grid.csv'.
If output directory is not provided, the current working directory will 
be used.

The gridded gradient will be written to output directory as 
'moho_gradient.csv'.

Reference:
B. L. N. Kennett 2019, "Areal parameter estimates from multiple datasets",
Proc. R. Soc. A. 475:20190352, http://dx.doi.org/10.1098/rspa.2019.0352

Requires:

- pyepsg
- cartopy
"""
import os

import numpy as np
from scipy.spatial.distance import cdist

from seismic.receiver_fn.moho_config import DIST_METRIC
from collections import defaultdict

DEFAULT_CUTOFF = 3.6

def _grid(bb_min, bb_max, spacing):
    if bb_min[0] >= bb_max[0]:
        raise ValueError(f"Bounds xmin {bb_min[0]} must be less than bounds xmax {bb_max[1]}")
    if bb_min[1] >= bb_max[1]:
        raise ValueError(f"Bounds ymin {bb_min[1]} must be less than bounds ymax {bb_max[1]}")
    span = bb_max - bb_min
    n_x = int(np.ceil(span[0]/spacing)) + 1
    n_y = int(np.ceil(span[1]/spacing)) + 1
    x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], n_x),
                                 np.linspace(bb_min[1], bb_max[1], n_y))
    return n_x, x_grid, n_y, y_grid


def _bounds(xy_mins, xy_maxs, bounds=None):
    if bounds is not None:
        l, b, r, t = bounds
        bb_min = np.array((l, b))
        bb_max = np.array((r, t))
    else:
        bb_min = np.min(xy_mins, axis=0)
        bb_max = np.max(xy_maxs, axis=0)
    return bb_min, bb_max


def make_grid(params):
    """
    Run multi point dataset weighted averaging over Gaussian interpolation functions
    to produce aggregate dataset.
    Source data coordinates are assumed to be in lon/lat (WGS84).
    First 2 rows of output file contain grid dimensions, followed by CSV data.
    """
    print("Generating Moho grid from point data")
    xy_mins = []
    xy_maxs = []
    for data in params.method_datasets:
        if(len(data.sta)==0): continue
        print(f"Building grid for '{data.name}'")
        pt_data = np.array((data.lon, data.lat, data.val)).T
        xy_map = pt_data[:, :2]
        xy_mins.append(xy_map.min(axis=0))
        xy_maxs.append(xy_map.max(axis=0))

    bb_min, bb_max = _bounds(xy_mins, xy_maxs, params.bounds)
    n_x, x_grid, n_y, y_grid = _grid(bb_min, bb_max, params.grid_interval)
    grid_map = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))
    denom_agg = np.zeros(grid_map.shape[0], dtype=float)[:, np.newaxis]
    z_agg = np.zeros_like(denom_agg)
    s_agg = np.zeros_like(denom_agg)
    print("Distance matrix calculation may take several minutes for large datasets")
    for data in params.method_datasets:
        if(len(data.sta)==0): continue
        print(f"Processing '{data.name}'")
        sigma = data.scale_length 
        sig_sq = np.square(sigma)
        cutoff = DEFAULT_CUTOFF
        max_dist = cutoff*sigma
        pt_data = np.array((data.lon, data.lat, data.val)).T
        xy_map = pt_data[:, :2]
        print(f"{xy_map.shape[0]} samples")
        print(f"Calculating distance matrix...")
        kwargs = {'max_dist': max_dist}
        dist_matrix = cdist(grid_map, xy_map, metric=DIST_METRIC, **kwargs)
        dist_weighting = None
        if(params.interp_function == 'gaussian'):
            dist_weighting = np.exp(-np.power(dist_matrix/sigma, 2.))
        else:
            dist_weighting = np.exp(-dist_matrix/sigma)
        # end if

        # We return nan instead of 0 from distance functions if beyond max distance, because
        # theoretically a sample exactly on the grid point will have a dist of 0, and the weighting
        # will be exp(0) == 1. After calculating exp, set NaNs to 0 so samples beyond max distance
        # have a weight of 0, rather than NaN which will propagate.
        dist_weighting[np.isnan(dist_weighting)] = 0
        z = (pt_data[:, 2][:, np.newaxis])
        zw = data.total_weight[:, np.newaxis]

        denom = np.matmul(dist_weighting, zw)
        denom_agg += denom

        z_numer = np.matmul(dist_weighting, (z * zw))
        s_numer = np.matmul(dist_weighting, (z*z * zw))
        z_agg += z_numer
        s_agg += s_numer
    # end for

    # Generate NaN where data doesn't exist
    prior_settings = np.seterr(divide='ignore', invalid='ignore')
    # Convert weights < 0.02 to NaN (from B.K's code)
    denom_agg = np.where(denom_agg < params.weight_cutoff, np.nan, denom_agg)
    # Get depth Z and std S for each grid cell
    Z = z_agg/denom_agg
    S = np.sqrt(s_agg/denom_agg - np.square(Z))
    np.seterr(**prior_settings)
    
    # Calculate gradient
    Z_2d = Z.reshape((n_y, n_x))
    v, u = np.gradient(Z_2d)
    gradient = np.array((u.flatten(), v.flatten())).T

    # Collect data and write to file
    if not os.path.exists(params.output_dir):
        os.mkdir(params.output_dir)

    depth_gridded = np.hstack((grid_map, Z, S))
    with open(params.grid_data, mode='w') as f:
        np.savetxt(f, [n_x, n_y], fmt='%d')
        np.savetxt(f, depth_gridded, fmt=['%.6f', '%.6f', '%.2f', '%.2f'], delimiter=',',
                   header='Lon,Lat,Depth,Stddev')

    gradient_gridded = np.hstack((grid_map, gradient))
    with open(params.grad_data, mode='w') as f:
        np.savetxt(f, [n_x, n_y], fmt='%d')
        np.savetxt(f, gradient_gridded, fmt=['%.6f', '%.6f', '%.6f', '%.6f'], delimiter=',',
                   header='Lon,Lat,U,V')

    print(f"Complete, results saved to '{params.grid_data}' and '{params.grad_data}'")
