#!/usr/bin/env python
"""
Andrew's original version of pointsets2grid.py that projects the lat/lon
points onto a flat projection so the distance matrix can be calculated
using Euclidean distance instead of spherical distance.

This makes it extremely fast, but it doesn't produce the same output
as B. Kennett's FORTRAN implementation.

This version also multiples each sample by the dataset weight. 
Previously the mean depth estimates from each dataset for each grid 
point were then weighted against the dataset weight, so the weighted
average of the dataset means was being taken.

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
        "bounds": [131.0, -36.0, 146.0, -27.0],
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

import click
import numpy as np
from scipy.spatial import cKDTree
import cartopy as cp


DEFAULT_CRS_CODE = 3395  # WGS84 World Mercator projection
DEFAULT_CUTOFF = 3.6  # dimensionless factor applied to scale length


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
    # end if

    # Load config
    with open(config_file, mode='r') as f:
        job_config = json.load(f)
    # end with
    epsg_code = job_config.get("proj_crs_code", DEFAULT_CRS_CODE)
    src_files = job_config["source_files"]
    if not src_files:
        return

    wgs_84 = cp.crs.Geodetic()
    map_proj = cp.crs.epsg(epsg_code)  # Requires internet connection
    assert map_proj.proj4_params['units'].lower() == 'm', 'Require projected length units of metres'

    bounds = job_config.get('bounds')
    if bounds is None:
        bb_min = np.array([np.inf, np.inf])
        bb_max = np.array([-np.inf, -np.inf])
    else:
        l, b, r, t = bounds
        bb_min = np.array(map_proj.transform_point(l, b, wgs_84))
        bb_max = np.array(map_proj.transform_point(r, t, wgs_84))

    for fname, filedict in src_files.items():
        pt_dataset = np.loadtxt(fname, delimiter=',')
        # Convert xy coordinates in lon/lat to map eastings, northings
        xyz_map = map_proj.transform_points(wgs_84, pt_dataset[:, 0], pt_dataset[:, 1])
        xyz_map[:, 2] = pt_dataset[:, 2]
        filedict['proj_data'] = xyz_map
        if bounds is None:
            xy_map = pt_dataset[:, :2]
            xy_min = xy_map.min(axis=0)
            xy_max = xy_map.max(axis=0)
            bb_min = np.minimum(bb_min, xy_min)
            bb_max = np.maximum(bb_max, xy_max)
        # Add sample weights to filedict if provided
        if filedict.get('enable_sample_weighting', False):
            filedict['sample_weights'] = pt_dataset[:, 3]

    # Generate regular grid in projected coordinate system covering the area of interest.
    # All the input datasets will be interpolated onto this grid.
    span = bb_max - bb_min
    spacing_m = job_config["output_spacing_km"]*1000  # km to metres
    n_x = int(np.ceil(span[0]/spacing_m))
    n_y = int(np.ceil(span[1]/spacing_m))
    x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], n_x),
                                 np.linspace(bb_min[1], bb_max[1], n_y))
    map_grid = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))
    map_tree = cKDTree(map_grid)

    denom_agg = np.zeros(map_grid.shape[0], dtype=float)[:, np.newaxis]
    z_agg = np.zeros_like(denom_agg)
    s_agg = np.zeros_like(denom_agg)
    # Loop over source datasets and interpolate each to map grid
    for _fname, filedict in src_files.items():
        sigma = filedict["scale_length_km"]*1000  # km to metres
        sig_sq = np.square(sigma)
        cutoff = filedict.get("scale_length_cutoff", DEFAULT_CUTOFF)
        max_dist = cutoff*sigma
        proj_data = filedict["proj_data"]
        xy_map = proj_data[:, :2]
        src_tree = cKDTree(xy_map)
        dist_matrix = src_tree.sparse_distance_matrix(map_tree, max_distance=max_dist).tocsr()
        dist_matrix = dist_matrix.transpose()
        dist_sq_matrix = -dist_matrix.power(2)/sig_sq  # Note the sign negation
        # Keep z in column vector format.
        z = proj_data[:, 2][:, np.newaxis]
        zw = filedict.get('sample_weights')
        w = filedict['weighting']
        # Compute interpolation of z onto target grid
        # TRICKY point: We compute the exponential terms as exp(x) - 1 rather than exp(x),
        # so that we can remain in sparse matrix format until the matrix multiplication.
        # After that, we will compensate for the missing +1 term in the sum.
        dist_weighting_m1 = dist_sq_matrix.expm1()
        mask = (dist_matrix != 0)
        if zw is not None:
            zw = zw[:, np.newaxis]
            denom = (dist_weighting_m1 + mask) * (zw * w)
            z_numer = (dist_weighting_m1 + mask) * (z * zw * w)
        else:
            denom = (dist_weighting_m1.sum(axis=1) + mask.sum(axis=1)) * w
            z_numer = (dist_weighting_m1 + mask) * (z * w)
        # Compute pointwise variance
        if zw is not None:
            s_numer = (dist_weighting_m1 + mask) * (np.square(z) * zw * w)
        else:
            s_numer = (dist_weighting_m1 + mask) * (np.square(z) * w)
        denom_agg += denom
        z_agg += z_numer
        s_agg += s_numer

    # Generate NaN where data doesn't exist
    prior_settings = np.seterr(divide='ignore', invalid='ignore')
    # Convert weights < 0.02 to NaN (from B.K's code)
    # Disabled when using flat projection - weight scalse are much smaller
    # denom_agg = np.where(denom_agg < 0.02, np.nan, denom_agg)
    Z = z_agg/denom_agg
    S = s_agg/denom_agg
    np.seterr(**prior_settings)
 
    # Convert back from local projection to global lon/lat coordinates so that downstream
    # consumers can use whatever coordinate systems and projections they want.
    grid_lonlat = wgs_84.transform_points(map_proj, map_grid[:, 0], map_grid[:, 1])[:, :2]

    # Collect data and write to file
    data_agg_gridded = np.hstack((grid_lonlat, Z, S))
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
