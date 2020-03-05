#!/usr/bin/env python
"""Ticket PV-116. Collect all inversion results and put them together
into one array amenable to 3D mapping.
"""

import os
import re
import glob

import click

import numpy as np

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
def main(input_folder, output_file, folder_mask, station_database):
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

            Vs_min, Vs_max, vel_discretization, width = post_dat[1].strip('\n').split(None)
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
            xy_range = np.array([[station_coords[0], station_coords[1]]]*depth_discretization)
            data_all = np.column_stack([xy_range, z_range, vel_mean])
            station_profiles.append(data_all)

        # end if
    # end for

    data_all = np.vstack(station_profiles)
    np.save(output_file, data_all, allow_pickle=False)
    print('Saved {} size array to {}'.format(data_all.shape, output_file))

# end func

if __name__ == '__main__':
    main()
# end if
