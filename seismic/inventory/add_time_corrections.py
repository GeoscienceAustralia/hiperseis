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
from obspy import Inventory
from obspy.core.inventory import Network
from obspy.core.util import AttribDict


def get_csv_correction_data(path_csvfile):
    """
    Read in the csv data from an input file, get the network_code, station_code, csv_data

    Args:
        path_csvfile: input csv file in /g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/
        $ head 7D.DE43_clock_correction.csv
        net,sta,date,clock_correction
        7D,DE43,2012-11-27,1.0398489013215846
        7D,DE43,2012-11-28,0.9408504322549281
        7D,DE43,2012-11-29,0.8418519631882714
        7D,DE43,2012-11-30,0.7428534941216148
        7D,DE43,2012-12-01,0.6438550250549583

    Returns: (network_code, station_code, csv_data)

    """
    with open(path_csvfile, "r") as csvfid:
        all_lines = csvfid.readlines()

    print("The csv file length", len(all_lines))
    print("the first line", all_lines[0])
    print("the last  line", all_lines[-1])

    line2 = all_lines[1]
    print("line2=", line2)

    my_items = line2.split(",")

    network_code = my_items[0]  # network_code = "7D"
    station_code = my_items[1]  # station_code = "DE43"

    return (network_code, station_code, all_lines)


def add_gpscorrection_into_stationxml(csv_file, input_xml, out_xml=None):
    """
    read in the correction CSV data from a file, get the station metadata node from input_xml file, 
    then add the CSV data into the station xml node to write into out_xml
    Args:
        csv_file: input csv file with correction data
        input_xml: input original stationXML file which contains the metadata for the network and station of csv_file
        out_xml:  Directory of the output xml file

    Returns:
        full path of the output xml file

    """

    ns = "https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0"

    (net, sta, csv_data) = get_csv_correction_data(csv_file)

    # path2_myxml = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"
    my_inv = read_inventory(input_xml)

    # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.select.html#obspy.core.inventory.inventory.Inventory.select

    selected_inv = my_inv.select(station=sta)

    print(selected_inv)

    my_tag = AttribDict()
    my_tag.namespace = ns
    my_tag.value = csv_data

    selected_inv.networks[0].stations[0].extra = AttribDict()
    selected_inv.networks[0].stations[0].extra.gpsclockcorrection = my_tag

    stationxml_with_csv = 'modified_inventory_%s.xml' % sta

    if out_xml is not None and os.path.isdir(out_xml):
        stationxml_with_csv = os.path.join(out_xml, stationxml_with_csv)

    selected_inv.write(stationxml_with_csv, format='STATIONXML',
                       nsmap={'GeoscienceAustralia': 'https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0'})

    # my_inv.write('modified_inventory.xml', format='STATIONXML')

    return stationxml_with_csv


def extract_csvdata(path2xml):
    """
    Read and extract the csv data from an inventory file

    Args:
        path2xml: path_to_stationxml

    Returns:
        csv_str

    """

    new_inv = read_inventory(path2xml)

    csv_str = new_inv.networks[0].stations[0].extra.gpsclockcorrection.value

    return csv_str

# ----------------------------------------------------------------------------------------------------------------
# Quick test code
# Example How to run:
# python add_time_corrections.py  /g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/7D.CZ40_clock_correction.csv
#                                ../../tests/testdata/7D_2012_2013.xml      ~/tmpdir/
# ----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    USAGE = "python %s csv_file inventory_file [out_dir}" %sys.argv[0]

    # time_correction_csvfile = "/home/feizhang/Githubz/hiperseis/tests/testdata/corrections/7D.DE43_clock_correction.csv"
    # time_correction_csvfile = "/home/feizhang/Githubz/hiperseis/tests/testdata/corrections/7D.CZ40_clock_correction.csv"
    time_correction_csvfile = "/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/7D.CZ40_clock_correction.csv"

    my_inventory = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"
    my_inventory = "../../tests/testdata/7D_2012_2013.xml"

    if len(sys.argv)<3:
        print ("USAGE")
    else:
        time_correction_csvfile = sys.argv[1]
        my_inventory = sys.argv[2]

    if len(sys.argv) >= 4:
        out_dir = sys.argv[3]
    else:
        out_dir = None

    output_xml = add_gpscorrection_into_stationxml(time_correction_csvfile, my_inventory, out_xml=out_dir)

    csvstr = extract_csvdata(output_xml)

    print(csvstr)
    print(type(csvstr))
