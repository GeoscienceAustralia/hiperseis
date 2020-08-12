#!/usr/bin/env python
"""
Blend multiple datasets on different lon/lat point sets together onto a common grid.

Uses weighted Gaussian interpolation method of Kennett, but simplified to
ignore the individual point weighting.

Multiple datasets with per-dataset settings are passed in using a JSON configuration
file with layout illustrated by the following example::

    {
        "source_files":
        {
            "file1":
            {
                "weighting": 2.5,
                "scale_length_km": 10.0,
                "scale_length_cutoff": 3.5
            },
            "file2":
            {
                "weighting": 0.5,
                "scale_length_km": 20.0
            },
            "file3":
            {
                "weighting": 1.0,
                "scale_length_km": 10.0,
                "scale_length_cutoff": 3
            }
        },
        "proj_crs_code": 3577,
        "output_spacing_km": 10.0
    }

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
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import cartopy as cp

try:
    import kennett_dist
    DIST_METRIC=_kennett_dist
except ImportError:
    DIST_METRIC=_haversine

DEFAULT_CRS_CODE = 3395  # WGS84 World Mercator projection
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
        return 0.
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
        return 0.
    else:
        return delta_degrees


@click.command()
@click.option('--config-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True)
def main(config_file, output_file):
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
    # Explicitly ensure output file extension is .csv
    outfile_base, ext = os.path.splitext(output_file)
    if ext.lower() != '.csv':
        output_file = outfile_base + '.csv'

    # Load config
    with open(config_file, mode='r') as f:
        job_config = json.load(f)
    src_files = job_config['source_files']
    bounds = job_config.get('bounds')
    if not src_files:
        print("No source files provided, exiting")
        return

    if bounds is None:
        bb_min = np.array([np.inf, np.inf])
        bb_max = np.array([-np.inf, -np.inf])
    else:
        l, b, r, t = bounds
        bb_min = np.array((l, b))
        bb_max = np.array((r, t))

    for fname, filedict in src_files.items():
        pt_dataset = np.loadtxt(fname, delimiter=',')
        filedict['pt_data'] = pt_dataset
        # If bounds not provided, find extent of datasets
        if bounds is None:
            xy_map = pt_dataset[:, :2]
            xy_min = xy_map.min(axis=0)
            xy_max = xy_map.max(axis=0)
            bb_min = np.minimum(bb_min, xy_min)
            bb_max = np.maximum(bb_max, xy_max)
        # Add sample weights to filedict if provided
        if filedict.get('enable_sample_weighting', False):
            filedict['sample_weights'] = pt_dataset[:, 3]
    
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
    for _fname, filedict in src_files.items():
        sigma = filedict["scale_length_degrees"]
        sig_sq = np.square(sigma)
        cutoff = filedict.get("scale_length_cutoff", DEFAULT_CUTOFF)
        max_dist = cutoff*sigma
        pt_data = filedict["pt_data"]
        xy_map = pt_data[:, :2]
        kwargs = {'max_dist': max_dist}
        dist_matrix = cdist(grid_map, xy_map, metric=DIST_METRIC, **kwargs)
        dist_sq_matrix = (dist_matrix/sigma)**2
        dist_weighting_m1 = np.exp(-dist_sq_matrix)
        dist_weighting_m1 = np.where(dist_weighting_m1 == 1, 0, dist_weighting_m1)
        z = pt_data[:, 2][:, np.newaxis]
        zw = filedict.get('sample_weights')
        w = filedict["weighting"]
        if zw is not None:
            zw = zw[:, np.newaxis]
            denom = ((dist_weighting_m1) * zw.T * w).sum(axis=1)[:, np.newaxis]
            z_numer = ((dist_weighting_m1) * (z * zw).T * w).sum(axis=1)[:, np.newaxis]
        else:
            denom = (dist_weighting_m1 * w).sum(axis=1)[:, np.newaxis]
            z_numer = ((dist_weighting_m1) * z.T * w).sum(axis=1)[:, np.newaxis]
        if zw is not None:
            s_numer = ((dist_weighting_m1)*(np.square(z) * zw).T * w).sum(axis=1)[:, np.newaxis]
        else:
            s_numer = ((dist_weighting_m1)*np.square(z).T * w).sum(axis=1)[:, np.newaxis]
        denom_agg += denom
        z_agg += z_numer
        s_agg += s_numer
    
    # Generate NaN where data doesn't exist
    prior_settings = np.seterr(divide='ignore', invalid='ignore')
    # Convert weights < 0.02 to NaN (from B.K's code)
    denom_agg = np.where(denom_agg < 0.02, np.nan, denom_agg)
    Z = z_agg/denom_agg
    S = s_agg/denom_agg
    np.seterr(**prior_settings)
    
    # Collect data and write to file
    data_agg_gridded = np.hstack((grid_map, Z, S))
    with open(output_file, mode='w') as f:
        np.savetxt(f, [n_x, n_y], fmt='%d')
        np.savetxt(f, data_agg_gridded, fmt=['%.6f', '%.6f', '%.2f', '%.2f'], delimiter=',',
                   header='Lon,Lat,Depth,Stddev')
        # To load this data back in:
        # >>> with open('test_pts2grid.csv', 'r') as f:
        # ...   nx = int(f.readline())
        # ...   ny = int(f.readline())
        # ...   ds = np.loadtxt(f, delimiter=',')
        # >>> # For each column in ds, reshape((ny, nx))
    # end with

    # See https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html
    # for how to use projection and transform settings with matplotlib.
    # Since this script outputs in lon/lat, it is in CRS WGS84 which is default for cartopy.crs.Geodetic()

# end func


if __name__ == '__main__':
    main()
# end if
