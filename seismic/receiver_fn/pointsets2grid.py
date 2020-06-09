"""
Blend multiple datasets on different lon/lat point sets together onto a common grid.

Uses weighted Gaussian interpolation method of Kennett, but simplified to
ignore the individual point weighting.

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
DEFAULT_CUTOFF = 3.5  # dimensionless factor applied to scale length


"""
Example config file:
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
"""


@click.command()
@click.option('--config-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True)
def main(config_file, output_file):
    """
    Source data coordinates are assumed to be in lon/lat (WGS84).

    :param config_file:
    :param output_file:
    :return:
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

    # Figure out bounds of aggregate dataset in projected CRS, load data,
    # and project loaded source data into projected CRS.
    bb_min = np.array([np.inf, np.inf])
    bb_max = np.array([-np.inf, -np.inf])
    wgs_84 = cp.crs.Geodetic()
    map_proj = cp.crs.epsg(epsg_code)  # Requires internet connection
    for fname, filedict in src_files.items():
        pt_dataset = np.loadtxt(fname, delimiter=',')
        # Convert xy coordinates in lon/lat to map eastings, northings
        xyz_map = map_proj.transform_points(wgs_84, pt_dataset[:, 0], pt_dataset[:, 1])
        xyz_map[:, 3] = pt_dataset[:, 3]
        filedict['proj_data'] = xyz_map
        xy_map = xyz_map[:, :2]
        xy_min = xy_map.min(axis=0)
        xy_max = xy_map.max(axis=0)
        bb_min = np.minimum(bb_min, xy_min)
        bb_max = np.minimum(bb_max, xy_max)
    # end for

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
        # Compute interpolation of z onto target grid
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
        # Store result
        filedict["z_interp"] = z_interp
        filedict["z_uncertainty"] = S
    # end for

    # Produce weighted mean of all source datasets
    pass
    exit(0)

    # # Load point data whose x-y coordinates are (longitude, latitude)
    # infile = 'ccp_line_data_sample.csv'
    # data_vol = np.loadtxt(infile, delimiter=',')
    # # z here is any scalar spatial data metric (e.g. depth to Moho).
    # # Keep z in column vector format.
    # z = data_vol[:, 2][:, np.newaxis]
    # # xy coordinates in lon/lat
    # xy = data_vol[:, 0:2]

    # # Convert lon/lat to cartesian so we can compute distances between stations on equal basis

    # # Create KDTree using projected coordinate system
    # map_proj = cp.crs.epsg(EPSG_CODE)  # Requires internet connection
    # xy_map = map_proj.transform_points(wgs_84, xy[:, 0], xy[:, 1])
    # xy_map = xy_map[:, :2]  # Ignore z. xy here are in metres.
    # tree = cKDTree(xy_map)
    #
    # bb_min = tree.mins
    # bb_max = tree.maxes
    # span = bb_max - bb_min
    # spacing_m = 10.0*1000  # km to metres
    # # Produce linear spacing in cartesian projection
    # x_n = int(np.ceil(span[0]/spacing_m))
    # y_n = int(np.ceil(span[1]/spacing_m))
    # x_grid, y_grid = np.meshgrid(np.linspace(bb_min[0], bb_max[0], x_n),
    #                              np.linspace(bb_min[1], bb_max[1], y_n))
    # grid = np.hstack((x_grid.flatten()[:, np.newaxis], y_grid.flatten()[:, np.newaxis]))

    # # Proximity weighting based on http://dx.doi.org/10.1098/rspa.2019.0352 with only
    # # one dataset
    # sigma = 25.0*1000  # km to metres
    # sig_sq = np.square(sigma)
    # max_dist = 3.5*sigma
    # dist_matrix = tree.sparse_distance_matrix(cKDTree(grid), max_distance=max_dist).tocsr()
    # dist_matrix = dist_matrix.transpose()
    # dist_sq_matrix = -dist_matrix.power(2)/sig_sq  # Note the sign negation
    #
    # # TRICKY point: We compute the exponential terms as exp(x) - 1 rather than exp(x),
    # # so that we can remain in sparse matrix format until the matrix multiplication.
    # # After that, we will compensate for the missing +1 term in the sum.
    # dist_weighting_m1 = dist_sq_matrix.expm1()
    # mask = (dist_matrix != 0)
    # denom = dist_weighting_m1.sum(axis=1) + mask.sum(axis=1)
    # numer = (dist_weighting_m1 + mask)*z
    # z_interp = numer/denom
    # # Compute pointwise variance
    # numer = (dist_weighting_m1 + mask)*np.square(z)
    # S = np.sqrt(numer/denom - np.square(z_interp))

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

# end func


if __name__ == '__main__':
    main()
# end if
