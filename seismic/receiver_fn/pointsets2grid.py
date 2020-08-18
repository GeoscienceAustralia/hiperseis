#!/usr/bin/env python
"""
Blend data from multiple methods together onto a common grid.

Multiple methods with per-method settings are passed in using a JSON 
configuration file, illustrated by the following example:

.. code-block:: JSON

    {
        "methods":
        [
            {
                "name": "CCP1",
                "csv_file": "sandbox/PV-227/ccp_line_data_sample1.csv",
                "weighting": 1.0,
                "scale_length_degrees": 0.2,
                "enable_sample_weighting": false
            },
            {
                "name": "CCP2",
                "csv_file": "sandbox/PV-227/ccp_line_data_sample2.csv",
                "weighting": 0.5,
                "scale_length_degrees": 0.2,
                "scale_length_cutoff": 3.6,
                "enable_sample_weighting": false
            }
        ],
        "output_spacing_degrees": 0.25,
        "bounds": [130.0, -20.0, 140.0, -17.0],
        "output_dir": "./ccp_moho_output",
    }

The CSV data is of format:
    LON LAT DEPTH [SAMPLE_WEIGHT]

Sample weighting can be switched on or off per method.
If enabled, requires the data file to have a fourth column containg
sample weights.

Bounds are optional. If provided, the interpolation grid will be limited
to this extent. If not provided, the interpolation grid is bounded by
the maximum extent of the aggregate datasets.

The gridded data will be placed in the output directory as 
'moho_grid.csv'. If output directory is not provided, the current 
working directory will be used.

The gridded gradient will be in output directory as 'moho_gradient.csv'.

Reference:
B. L. N. Kennett 2019, "Areal parameter estimates from multiple datasets",
Proc. R. Soc. A. 475:20190352, http://dx.doi.org/10.1098/rspa.2019.0352

Requires:

- pyepsg
- cartopy
"""
import os
import json
import math

import click
import numpy as np
from scipy.spatial.distance import cdist
import cartopy as cp

try:
    import kennett_dist
    K_DIST = True
except ImportError:
    print("Kennett distance function not found. Using haversine.")
    K_DIST = False

DEFAULT_CUTOFF = 3.6  # dimensionless factor applied to scale length

def _kennett_dist(s, r, **kwargs):
    """
    Spherical distance. Requires 'kennett_dist' f2py module.
    Copyright B.L.N Kennett 1978, R.S.E.S A.N.U
    """
    clons, clats = s
    clonr, clatr = r
    max_dist = kwargs.get('max_dist', np.inf)
    # Distance in degrees
    delta, _, _ = kennett_dist.ydiz(clats, clons, clatr, clonr)
    if delta > max_dist:
        return np.nan
    else:
        return delta


def _haversine(s, r, **kwargs):
    """
    Haversine distance.
    """
    lon1, lat1 = s
    lon2, lat2 = r
    max_dist = kwargs.get('max_dist', np.inf)
    R = 6372.8

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))

    delta_degrees = c * 180.0/math.pi

    if delta_degrees > max_dist:
        return np.nan
    else:
        return delta_degrees


DIST_METRIC = _kennett_dist if K_DIST else _haversine


def make_grid(config_file):
    """
    Run multi point dataset weighted averaging over Gaussian interpolation functions
    to produce aggregate dataset.
    Source data coordinates are assumed to be in lon/lat (WGS84).
    First 2 rows of output file contain grid dimensions, followed by CSV data.

    Example usage::

        python pointsets2grid.py --config-file config_pts2grid_example.json \
            --output-file test_pts2grid.csv

    :param config_file: Input filename of the job configuration in JSON format
    :param output_file: Name of output file
    :return: None
    """
    print("Generating Moho grid from point data")
    # Load config
    with open(config_file, mode='r') as f:
        job_config = json.load(f)

    try:
        methods = job_config['methods']
    except KeyError as e:
        raise Exception("No methods provided") from e

    bounds = job_config.get('bounds')
    if bounds is None:
        bb_min = np.array([np.inf, np.inf])
        bb_max = np.array([-np.inf, -np.inf])
    else:
        l, b, r, t = bounds
        bb_min = np.array((l, b))
        bb_max = np.array((r, t))

    for method_params in methods:
        method = method_params['name']
        pt_dataset = np.loadtxt(method_params['csv_file'], delimiter=',')
        method_params['pt_data'] = pt_dataset
        # If bounds not provided, find extent of datasets
        if bounds is None:
            xy_map = pt_dataset[:, :2]
            xy_min = xy_map.min(axis=0)
            xy_max = xy_map.max(axis=0)
            bb_min = np.minimum(bb_min, xy_min)
            bb_max = np.maximum(bb_max, xy_max)
        # Add sample weights to filedict if provided
        if method_params['enable_sample_weighting']:
            try:
                method_params['sample_weights'] = (pt_dataset[:, 3][:, np.newaxis]).T
            except IndexError as e:
                raise Exception(
                    f"Sample weights enabled for '{method}' but could not be found in "
                      "csv file. Ensure fourth column exists and contains sample weights.") from e
    
    span = bb_max - bb_min
    spacing = job_config["output_spacing_degrees"]
    n_x = int(np.ceil(span[0]/spacing)) + 1
    n_y = int(np.ceil(span[1]/spacing)) + 1
    x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], n_x),
                                 np.linspace(bb_min[1], bb_max[1], n_y))
    grid_map = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))
    denom_agg = np.zeros(grid_map.shape[0], dtype=float)[:, np.newaxis]
    z_agg = np.zeros_like(denom_agg)
    s_agg = np.zeros_like(denom_agg)
    print("Distance matrix calculation may take several minutes for large datasets")
    for method_params in methods:
        method = method_params['name']
        print(f"Processing '{method}'")
        sigma = method_params["scale_length_degrees"]
        sig_sq = np.square(sigma)
        cutoff = method_params.get("scale_length_cutoff", DEFAULT_CUTOFF)
        max_dist = cutoff*sigma
        pt_data = method_params["pt_data"]
        xy_map = pt_data[:, :2]
        print(f"{xy_map.shape[0]} samples")
        print(f"Calculating distance matrix...")
        kwargs = {'max_dist': max_dist}
        dist_matrix = cdist(grid_map, xy_map, metric=DIST_METRIC, **kwargs)
        dist_sq_matrix = (dist_matrix/sigma)**2
        dist_weighting = np.exp(-dist_sq_matrix)
        # We return nan instead of 0 from distance functions if beyond max distance, because 
        # theoretically a sample exactly on the grid point will have a dist of 0, and the weighting
        # will be exp(0) == 1. After calculating exp, set NaNs to 0 so samples beyond max distance
        # have a weight of 0, rather than NaN which will propagate.
        dist_weighting[np.isnan(dist_weighting)] = 0
        z = (pt_data[:, 2][:, np.newaxis]).T
        zw = method_params.get('sample_weights') if method_params.get('sample_weights') is not None \
                else np.ones_like(z)
        w = method_params["weighting"]
        zw = zw * w
        denom = ((dist_weighting) * zw).sum(axis=1)[:, np.newaxis]
        z_numer = ((dist_weighting) * (z * zw)).sum(axis=1)[:, np.newaxis]
        s_numer = ((dist_weighting)*(np.square(z) * zw)).sum(axis=1)[:, np.newaxis]
        denom_agg += denom
        z_agg += z_numer
        s_agg += s_numer
    
    # Generate NaN where data doesn't exist
    prior_settings = np.seterr(divide='ignore', invalid='ignore')
    # Convert weights < 0.02 to NaN (from B.K's code)
    denom_agg = np.where(denom_agg < 0.02, np.nan, denom_agg)
    # Get depth Z and std S for each grid cell
    Z = z_agg/denom_agg
    S = s_agg/denom_agg
    np.seterr(**prior_settings)
    
    # Calculate gradient
    Z_2d = Z.reshape((n_y, n_x))
    v, u = np.gradient(Z_2d)
    gradient = np.array((u.flatten(), v.flatten())).T

    # Collect data and write to file
    output_dir = job_config.get('output_dir', os.getcwd())
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    depth_gridded = np.hstack((grid_map, Z, S))
    grid_output_file = os.path.join(output_dir, 'moho_grid.csv')
    with open(grid_output_file, mode='w') as f:
        np.savetxt(f, [n_x, n_y], fmt='%d')
        np.savetxt(f, depth_gridded, fmt=['%.6f', '%.6f', '%.2f', '%.2f'], delimiter=',',
                   header='Lon,Lat,Depth,Stddev')

    gradient_gridded = np.hstack((grid_map, gradient))
    gradient_output_file = os.path.join(output_dir, 'moho_gradient.csv')
    with open(gradient_output_file, mode='w') as f:
        np.savetxt(f, [n_x, n_y], fmt='%d')
        np.savetxt(f, gradient_gridded, fmt=['%.6f', '%.6f', '%.6f', '%.6f'], delimiter=',',
                   header='Lon,Lat,U,V')

    print(f"Complete, results saved to '{grid_output_file}' and '{gradient_output_file}'")


@click.command()
@click.option('--config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    make_grid(config_file)


if __name__ == '__main__':
    main()
