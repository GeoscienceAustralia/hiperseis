#!/usr/bin/env python
"""
Load inventory file from hdf5 format and export station location plots to graphics files.

Usage::
    python generate_network_plots.py inventory_20190206.h5

where inventory_20190206.h5 here is an example file name.  The output folder for the
graphics files is inferred from the inventory file name, which would be "inventory_20190206"
in this example.
"""

import os
import click

import pandas as pd

from seismic.inventory.plotting import save_network_local_plots, save_station_local_plots

try:
    import tqdm
    show_progress = True
except ImportError:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")
# end try


@click.command()
@click.option('--plot-type', type=click.Choice(['network', 'station']), default='network',
              help='Whether to plot on map per network or one map per station')
@click.argument('inventory-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(inventory_h5_file, plot_type):
    """Generate map plots of station locations for an inventory of stations.

    :param inventory_h5_file: Name of input HDF5 file containing a station inventory
    :type inventory_h5_file: str or Path
    :param plot_type: Type of plots to make - per network or per station
    :type plot_type: str
    """
    # Set up file names
    print("Loading file {0}".format(inventory_h5_file))
    basename, _ = os.path.splitext(inventory_h5_file)
    _, basename = os.path.split(basename)
    folder_name = "plots." + basename
    print("Saving plots to {0}".format(folder_name))

    # Read inventory from hdf5 file
    db = pd.read_hdf(inventory_h5_file, mode='r', key='inventory')

    if show_progress:
        pbar = tqdm.tqdm(total=len(db), ascii=True)
        progressor = pbar.update
    else:
        progressor = None
    # end if

    # Save network plots to files
    plot_functor = {'network': save_network_local_plots, 'station': save_station_local_plots}
    plot_functor[plot_type](db, folder_name, progressor=progressor)
    if show_progress:
        pbar.close()
    # end if
# end main

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
# end if
