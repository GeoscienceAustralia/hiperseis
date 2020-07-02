#! /usr/bin/env python
"""
Add GPS clock time correction csv_data into inventory file to get a modified station xml file

CreationDate:
    24/02/2020

Developer:
    fei.zhang@ga.gov.au

"""

import os
import sys

from obspy import read_inventory
from obspy.core.util import AttribDict


def get_csv_correction_data(path_csvfile):
    """
    Read in the csv data from an input file, get the network_code, station_code, csv_data. Format::

        $ head 7D.DE43_clock_correction.csv
        net,sta,date,clock_correction
        7D,DE43,2012-11-27,1.0398489013215846
        7D,DE43,2012-11-28,0.9408504322549281
        7D,DE43,2012-11-29,0.8418519631882714
        7D,DE43,2012-11-30,0.7428534941216148
        7D,DE43,2012-12-01,0.6438550250549583

    :param path_csvfile: input csv file in /g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/
    :return: (network_code, station_code, csv_data)
    """

    with open(path_csvfile, "r") as csvfid:
        all_csv = csvfid.read()

    line2 = all_csv.split('\n')[1]
    # print(line2)

    my_items = line2.split(",")

    network_code = my_items[0].strip()  # network_code = "7D"
    station_code = my_items[1].strip()  # station_code = "DE43"

    return (network_code, station_code, all_csv)


def add_gpscorrection_into_stationxml(csv_file, input_xml, out_xml=None):
    """
    Read in the correction CSV data from a file, get the station metadata node from input_xml file,
    then add the CSV data into the station xml node to write into out_xml

    :param csv_file: input csv file with correction data
    :param input_xml: input original stationXML file which contains the metadata for the network and station of csv_file
    :param out_xml:  Directory of the output xml file
    :return: full path of the output xml file
    """

    ns = "https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0"

    (net, sta, csv_data) = get_csv_correction_data(csv_file)

    # path2_myxml = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"
    my_inv = read_inventory(input_xml, format='STATIONXML')

    # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.select.html#obspy.core.inventory.inventory.Inventory.select

    selected_inv = my_inv.select(network=net, station=sta)

    # print(selected_inv)

    my_tag = AttribDict()
    my_tag.namespace = ns
    my_tag.value = csv_data

    selected_inv.networks[0].stations[0].extra = AttribDict()
    selected_inv.networks[0].stations[0].extra.gpsclockcorrection = my_tag

    stationxml_with_csv = '%s.%s_station_inv_modified.xml' % (net, sta)

    if out_xml is not None and os.path.isdir(out_xml):
        stationxml_with_csv = os.path.join(out_xml, stationxml_with_csv)

    selected_inv.write(stationxml_with_csv, format='STATIONXML',
                       nsmap={'GeoscienceAustralia': 'https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0'})

    # my_inv.write('modified_inventory.xml', format='STATIONXML')

    return stationxml_with_csv


def extract_csvdata(path2xml):
    """
    Read the station xml file and extract the csv data to be parsed by pandas

    :param path2xml: path_to_stationxml
    :return: csv_str
    """
    import io
    import pandas as pd

    new_inv = read_inventory(path2xml, format='STATIONXML')

    csv_str = new_inv.networks[0].stations[0].extra.gpsclockcorrection.value
    # print(csv_str)
    # print(type(csv_str))

    df_clock_correction = pd.read_csv(io.StringIO(csv_str))

    print(df_clock_correction.head())

    return csv_str


# ----------------------------------------------------------------------------------------------------------------
# Quick test code
# Example How to run:
# python add_time_corrections.py  /g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/7D.CZ40_clock_correction.csv
#                                ../../tests/testdata/7D_2012_2013.xml      ~/tmpdir/
# ----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    USAGE = "python %s csv_file inventory_file [out_dir]" % sys.argv[0]

    if len(sys.argv) < 3:
        print(USAGE)
        sys.exit(1)
    else:
        time_correction_csvfile = sys.argv[1]
        my_inventory = sys.argv[2]

    if len(sys.argv) >= 4:
        out_dir = sys.argv[3]
    else:
        out_dir = None

    output_xml = add_gpscorrection_into_stationxml(time_correction_csvfile, my_inventory, out_xml=out_dir)

    # Optional test to extract the CSV data and make a pandas dataframe object for future use.
    # csvstr = extract_csvdata(output_xml)
