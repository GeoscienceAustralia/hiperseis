#!/usr/bin/env python
"""Ticket PV-116. Collect all inversion results and put them together
into one array amenable to 3D mapping.
"""

import sys
import os
import re
import glob

import click

import numpy as np
from scipy.interpolate import interp1d

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

# File to look for in case output folder containing solution
SOLUTION_FILE = 'Posterior.out'


@click.command()
@click.argument('input-folder', type=click.Path(exists=True, file_okay=False), required=True)
@click.argument('output-file', type=click.Path(exists=False, file_okay=True), required=True)
@click.option('--folder-mask', type=str, required=True,
              help='Regular expression mask used for identifying solution folders to scrape.'
                   'e.g "OA_B[S-V]*_OUT" (use quotes to prevent shell from expanding wildcard)')
@click.option('--station-database', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Location of station database used to generate FederatedASDFDataSet. '
                   'Provides station location metadata. e.g. "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt".')
@click.option('--max-depth-km', type=float, required=True, help='Max depth in km')
@click.option('--depth-levels', type=int, required=True, help='Number of levels in depth discretization')
def aggregate(input_folder, output_file, folder_mask, station_database,
              max_depth_km, depth_levels):
    """
    Scrape together all the trans-D inversion solutions and collect into volumetric dataset.

    :param input_folder: Folder containing solutions to scrape together
    :type input_folder: str or Path
    :param output_file: Output file (must not exist already)
    :type output_file: str or Path (pdf extension expected)
    """

    # Open station database from which to get station lat,lon coordinates
    station_location_db = FederatedASDFDataSet(station_database).unique_coordinates

    # Process folders in alphanumerical order
    folders = sorted(glob.glob(os.path.join(input_folder, folder_mask)))

    # regex pattern for matching case strings containing network, station and channel codes
    case_pattern = '^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)'
    matcher = re.compile(case_pattern)

    # Container for storing Vs as a function of depth for each station.
    station_profiles = []

    # Loop over folders one at a time
    for f in folders:
        # If it is a folder and has a solution file in it.
        if os.path.isdir(f) and os.path.isfile(os.path.join(f, SOLUTION_FILE)):
            _, case_folder = os.path.split(f)
            case_meta = matcher.match(case_folder)
            # Extract network, station and channel metadata from folder name
            net = case_meta.group(1)
            sta = case_meta.group(2)
            cha = case_meta.group(3)
            station_id = '.'.join([net, sta, cha])

            soln_file = os.path.join(f, SOLUTION_FILE)
            # station_coords are in lon,lat order
            station_coords = station_location_db['.'.join([net, sta])]

            print(station_id, station_coords)

            # Open solution file and collect relevant fields
            with open(soln_file, 'r') as posterior:
                post_dat = posterior.readlines()
            # end with
            _0, depth_discretization, depth_max = post_dat[0].strip('\n').split(None)
            depth_discretization = int(depth_discretization)
            depth_max = float(depth_max)
            z_range = depth_max*(np.arange(depth_discretization) + 0.5)/depth_discretization

            Vs_min, Vs_max, vel_discretization, _width = post_dat[1].strip('\n').split(None)
            vel_discretization = int(vel_discretization)
            Vs_min, Vs_max = float(Vs_min), float(Vs_max)
            vel_range = Vs_min + (Vs_max - Vs_min)*(np.arange(vel_discretization) + 0.5)/vel_discretization
            # Each row of posterior_distribution corresponds to a discrete depth. At each depth,
            # we have a velocity probability distribution based on MCMC sampling.
            posterior_distribution = np.reshape(np.array([float(x.strip('\n')) for x in post_dat[2:]]),
                                                (depth_discretization, vel_discretization))

            # Normalize the distribution at each depth.
            post = posterior_distribution/np.expand_dims(np.sum(posterior_distribution, axis=-1), -1)
            assert np.allclose(np.sum(post, axis=-1), 1)

            # Compute mean at each depth, reducing the 2D posterior to 1D
            # velocity as a function of depth.
            vel_mean = np.dot(post, vel_range)

            # Create 4-column 2D matrix storing results for this station.
            xy_range = np.array([[station_coords[0], station_coords[1]]]*depth_levels)
            interpolator = interp1d(z_range, vel_mean, kind='cubic')
            z_interp = max_depth_km*(np.arange(depth_levels) + 0.5)/depth_levels
            vel_interp = interpolator(z_interp)
            data_all = np.column_stack([xy_range, z_interp, vel_interp])
            station_profiles.append(data_all)

        # end if
    # end for

    data_all = np.vstack(station_profiles)
    np.save(output_file, data_all, allow_pickle=False)
    print('Saved {} size array to {}'.format(data_all.shape, output_file))
# end func


@click.command()
@click.argument('data-filename', type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True)
def test_plot(data_filename):
    """
    Load scraped dataset and generate test plot based on it.

    :param data_filename: Name of dataset file to load. Run 'aggregate' commmand
        to generate dataset.
    """
    import cartopy as cp
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator

    # Load data file and plot slices on geographic plot
    print("Loading data")
    data_all = np.load(data_filename)

    xy = data_all[:, 0:2]
    bb_min = np.nanmin(xy, axis=0)
    bb_max = np.nanmax(xy, axis=0)
    print("Bounding box: {} to {}".format(bb_min, bb_max))
    xy_min = bb_min - 0.5
    xy_max = bb_max + 0.5
    xy_max[1] += 1.0  # Increase view to the north

    # Generate slices
    xyz = data_all[:, 0:3]
    # depth_range = np.linspace(0.0, 60.0, 61)
    print("Gridding data")
    xy_res = 51
    xy_range = np.linspace(xy_min, xy_max, xy_res)
    xyz_grid = np.meshgrid(xy_range.T[0], xy_range.T[1], 30.0)
    slice_30 = griddata(xyz, data_all[:, 3], np.array(xyz_grid).reshape((3, xy_res*xy_res)).T, rescale=True)
    slice_30 = slice_30.reshape((xy_res, xy_res))

    # Plot Vs contours at 30 km.
    xlocator = FixedLocator(np.arange(np.floor(xy_min[0] - 0.5), np.ceil(xy_max[0] + 0.5)))
    ylocator = FixedLocator(np.arange(np.floor(xy_min[1] - 0.5), np.ceil(xy_max[1] + 0.5)))

    map_proj = cp.crs.PlateCarree()
    resolution = '110m'

    print("Plotting")
    fig = plt.figure(figsize=(16, 9))
    ax = plt.subplot(1, 1, 1, projection=map_proj)
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])
    gridliner = ax.gridlines(draw_labels=True, linestyle=':', xlocs=xlocator, ylocs=ylocator)

    ax.add_feature(cp.feature.COASTLINE.with_scale(resolution))
    ax.add_feature(cp.feature.OCEAN.with_scale(resolution))
    ax.add_feature(cp.feature.LAND.with_scale(resolution))
    ax.text(-0.08, 0.5, 'Latitude (deg)', rotation='vertical', transform=ax.transAxes, rotation_mode='anchor',
            va='bottom', ha='center')
    ax.text(0.5, -0.14, 'Longitude (deg)', rotation='horizontal', transform=ax.transAxes, rotation_mode='anchor',
            va='bottom', ha='center')

    # Add contours
    xg, yg  = np.squeeze(xyz_grid[0]), np.squeeze(xyz_grid[1])
    ctr = plt.contourf(xg, yg, slice_30, cmap='copper_r')
    cb = plt.colorbar(ctr, pad=0.1, shrink=0.9, fraction=0.05)
    cb.set_label('$V_s$ (km/s)')

    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')

    plt.title('OA $V_s$ at 20 km depth', y=1.05)

    out_name = os.path.splitext(data_filename)[0] + '.png'
    plt.savefig(out_name, dpi=300)
    plt.close()

# end func


@click.group()
def main():
    pass
# end if


if __name__ == '__main__':
    main.add_command(aggregate)
    main.add_command(test_plot, name='plot')
    main()
# end if
