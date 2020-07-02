#!/usr/bin/env python
"""Split a station inventory file into a separate file per network.
"""

from __future__ import print_function

import os
import sys
from collections import defaultdict

import click
from obspy.core.inventory import Inventory

from seismic.inventory.inventory_util import load_station_xml
from seismic.inventory.fdsnxml_convert import sc3_conversion_available, toSc3ml

PY2 = sys.version_info[0] < 3  # pylint: disable=invalid-name

if PY2:
    import pathlib2 as pathlib  # pylint: disable=import-error
else:
    import pathlib  # pylint: disable=import-error

try:
    import tqdm
    show_progress = True
except ImportError:  # pragma: no cover
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")


def split_inventory_by_network(obspy_inv, output_folder, validate=False):
    """Export a station XML file per network for each network in given obspy Inventory.

    :param obspy_inv: Obspy Inventory containing the networks to export to file.
    :type obspy_inv: obspy.core.inventory.inventory.Inventory
    :param output_folder: Folder in which to output the per-network XML files. Will be created if doesn't yet exist.
    :type output_folder: str or pathlib.Path
    :param validate: Whether to validate the station data on write, defaults to False
    :type validate: bool, optional
    """
    pathlib.Path(output_folder).mkdir(exist_ok=True)

    if show_progress:
        pbar = tqdm.tqdm(total=len(obspy_inv.networks), ascii=True)
        std_print = pbar.write
    else:
        std_print = print

    # Since duplicate network codes can occur, we ensure that output file names are unique by keeping an instance
    # count for the occurrences of each network, and appending this to the file name.
    network_count = defaultdict(int)
    for network in obspy_inv:
        if show_progress:
            pbar.update()
            pbar.set_description("Network {}".format(network.code))
        net_inv = Inventory(networks=[network], source=obspy_inv.source)
        fname = "network_{}_{}.xml".format(network.code, network_count[network.code])
        network_count[network.code] += 1
        try:
            net_inv.write(os.path.join(output_folder, fname), format="stationxml", validate=validate)
        except Exception as e:
            std_print(str(e))
            std_print("FAILED writing file {0} for network {1}, continuing".format(fname, network.code))
            continue
    if show_progress:
        pbar.close()


def inventory_split(inv_file, output_folder, sc3ml=False):
    """Split an inventory file (station XML) into separate inventory file per network, stored in folder output_folder.

    :param inv_file: Station inventory to be split. This file is not changed by the script.
    :type inv_file: str or path
    :param output_folder: Folder where per-network station xml files will be generated.
    :type output_folder: str or path to folder
    :param sc3ml: If True, try to convert output files to sc3ml format if possible using seiscomp3. Defaults to False.
    :type sc3ml: bool, optional
    """
    print("Loading file {}".format(inv_file))
    stn_inventory = load_station_xml(inv_file)
    print("Splitting inventory into folder {}".format(output_folder))
    split_inventory_by_network(stn_inventory, output_folder)
    # Can be large in memory, so release as soon as no longer needed.
    del stn_inventory

    # Perform sc3ml conversion if requested and if supported.
    if sc3ml and sc3_conversion_available():
        sc3ml_output_folder = output_folder + "_sc3ml"
        try:
            toSc3ml(output_folder, sc3ml_output_folder)
        except OSError as e:
            print("WARNING: Unable to convert to sc3ml!")
            print(str(e))
            if os.path.isdir(sc3ml_output_folder):
                import uuid
                print("         Renaming {} to avoid accidental use!".format(sc3ml_output_folder))
                stashed_name = sc3ml_output_folder + ".BAD." + str(uuid.uuid4())[-8:]
                os.rename(sc3ml_output_folder, stashed_name)


@click.command()
@click.option('--inv-file', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True,
              help='Station inventory to be split. This file is not changed by the script.')
@click.option('--output-folder', type=click.Path(exists=False, dir_okay=True), required=True,
              help='Folder where per-network station xml files will be generated.')
@click.option('--sc3ml', is_flag=True, default=False, show_default=True,
              help='Try to convert output files to sc3ml format if possible using seiscomp3.')
def main(inv_file, output_folder, sc3ml=False):  # pragma: no cover
    """Main function.
    """
    inventory_split(inv_file, output_folder, sc3ml)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
