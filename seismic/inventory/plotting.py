#!/usr/bin/env python
"""Helper functions for plotting FDSN Network objects.
"""

import os
import sys
from collections import defaultdict

from obspy.core.inventory import Inventory
import matplotlib.pyplot as plt

from seismic.inventory.pdconvert import dataframe_to_network

# pylint: disable=invalid-name

if sys.version_info[0] < 3:
    import pathlib2 as pathlib  # pylint: disable=import-error
else:
    import pathlib  # pylint: disable=import-error
# end if

no_instruments = defaultdict(lambda: None)


def save_network_local_plots(df, plot_folder, progressor=None, include_stations_list=True):
    """
    Save visual map plot per network, saved to file netcode.png.

    :param df: Dataframe of station records to save.
    :type df: pandas.DataFrame conforming to seismic.inventory.table_format.TABLE_SCHEMA
    :param plot_folder: Name of output folder
    :type plot_folder: str
    :param progressor: Progress bar functor to receive progress updates, defaults to None
    :type progressor: Callable object receiving incremental update on progress, optional
    :param include_stations_list: If True, also export stationtxt file alongside each png file, defaults to True
    :type include_stations_list: bool, optional
    """
    dest_path = os.path.join(plot_folder, "networks")
    pathlib.Path(dest_path).mkdir(parents=True, exist_ok=True)
    failed = []
    for netcode, data in df.groupby('NetworkCode'):
        net = dataframe_to_network(netcode, data, no_instruments)
        plot_fname = os.path.join(dest_path, netcode + ".png")
        try:
            fig = net.plot(projection="local", resolution="l", outfile=plot_fname, continent_fill_color="#e0e0e0",
                           water_fill_color="#d0d0ff", color="#c08080")
            plt.close(fig)
        except Exception:
            failed.append(netcode)
        # end try

        if include_stations_list:
            inv = Inventory(networks=[net], source='EHB')
            inv_fname = os.path.join(dest_path, netcode + ".txt")
            inv.write(inv_fname, format="stationtxt")

        if progressor:
            progressor(len(data))
    # end for

    if failed:
        print("FAILED plotting on the following networks:")
        print("\n".join(failed))
    else:
        print("SUCCESS!")
    # end if
# end func


def save_station_local_plots(df, plot_folder, progressor=None, include_stations_list=True):
    """
    Save visual map plot per station, saved to file netcode.stationcode.png.

    :param df: Dataframe of station records to save.
    :type df: pandas.DataFrame conforming to seismic.inventory.table_format.TABLE_SCHEMA
    :param plot_folder: Name of output folder
    :type plot_folder: str
    :param progressor: Progress bar functor to receive progress updates, defaults to None
    :type progressor: Callable object receiving incremental update on progress, optional
    :param include_stations_list: If True, also export stationtxt file alongside each png file, defaults to True
    :type include_stations_list: bool, optional
    """
    dest_path = os.path.join(plot_folder, "stations")
    pathlib.Path(dest_path).mkdir(parents=True, exist_ok=True)
    failed = []
    for (netcode, statcode), data in df.groupby(['NetworkCode', 'StationCode']):
        net = dataframe_to_network(netcode, data, no_instruments)
        station_name = ".".join([netcode, statcode])
        plot_fname = os.path.join(dest_path, station_name + ".png")
        try:
            fig = net.plot(projection="local", resolution="l", outfile=plot_fname, continent_fill_color="#e0e0e0",
                           water_fill_color="#d0d0ff", color="#c08080")
            plt.close(fig)
        except Exception:
            failed.append(station_name)
        # end try

        if include_stations_list:
            inv = Inventory(networks=[net], source='EHB')
            inv_fname = os.path.join(dest_path, station_name + ".txt")
            inv.write(inv_fname, format="stationtxt")

        if progressor:
            progressor(len(data))
    # end for

    if failed:
        print("FAILED plotting on the following stations:")
        print("\n".join(failed))
    else:
        print("SUCCESS!")
    # end if
# end func
