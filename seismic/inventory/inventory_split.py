#!/usr/bin/env python
"""Split a station inventory file into a separate file per network.
"""

import os
import sys

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
except ImportError:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")


def split_inventory_by_network(obspy_inv, output_folder, sc3ml_convert=False):

    pathlib.Path(output_folder).mkdir(exist_ok=True)

    if show_progress:
        pbar = tqdm.tqdm(total=len(obspy_inv.networks), ascii=True)
        std_print = pbar.write
    else:
        std_print = print

    for network in obspy_inv:
        pbar.update()
        net_inv = Inventory(networks=[network], source=obspy_inv.source)
        fname = "network_{}.xml".format(network.code)
        try:
            net_inv.write(os.path.join(output_folder, fname), format="stationxml", validate=True)
        except Exception as e:
            std_print(e)
            std_print("FAILED writing file {0} for network {1}, continuing".format(fname, network.code))
            continue
    if show_progress:
        pbar.close()

    # Perform sc3ml conversion if requested and if supported.
    if sc3ml_convert and sc3_conversion_available():
        sc3ml_output_folder = output_folder + "_sc3ml"
        try:
            toSc3ml(output_folder, sc3ml_output_folder)
        except OSError:
            print("WARNING: Unable to convert to sc3ml!")
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
@click.option('--sc3ml', is_flag=True, help='Try to convert output files to sc3ml format if possible using seiscomp3.')
def main(inv_file, output_folder, sc3ml=False):
    print("Loading file {}".format(inv_file))
    stn_inventory = load_station_xml(inv_file)
    print("Splitting inventory into folder {}".format(output_folder))
    split_inventory_by_network(stn_inventory, output_folder, sc3ml_convert=sc3ml)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
