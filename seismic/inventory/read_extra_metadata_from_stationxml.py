#! /usr/bin/env python
"""
Description:
    Example python script to extract extra GA's metadata

CreationDate:   2020-11-24
Developer:      fei.zhang@ga.gov.au

"""

import os, sys
import obspy
from obspy.core import UTCDateTime

def check_invenory(in_stationxml, output_dir="/tmp/"):
    """
    Check the obspy IO-handling of inventory xml file
    :return:
    """
    big_inv = obspy.read_inventory(in_stationxml)  # read in fresh

    filename_new_xml = "check1_" + os.path.basename( in_stationxml)

    new_big_xml = os.path.join( output_dir, filename_new_xml)

    print("New XML file = ", new_big_xml)

    # Write out the inventory into a new file to see if they are identical?
    # big_inv.write(new_big_xml, format="stationxml", validate=True) #without nsmap

    GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"
    big_inv.write(new_big_xml, format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                          validate=True)


    # ObsPy Problem:  Missing some startDate for some stations
    # Why the obspy api missed some start_date for stations?

    # grep "<Station" OA_stations2.xml
    # grep "<Station" OA_stations_2017-2018_new.xml
    #     <Station code="BS24" startDate="2017-09-27T05:16:47.000000Z" endDate="2018-07-04T04:38:50.000000Z">
    #     <Station code="BS25" endDate="2018-06-27T11:44:37.000000Z">
    #     <Station code="BS26" startDate="2017-09-27T00:49:10.000000Z" endDate="2018-01-12T04:02:13.000000Z">
    #     <Station code="BS27" startDate="2017-09-26T00:42:49.000000Z" endDate="2018-07-03T03:04:17.000000Z">

    for net in big_inv.networks:
        print("The Network Code", net.code)
        number_of_stations = len(net.stations)
        print("The total number of stations =", number_of_stations)

    # identify None start or end Date
    for i in range(number_of_stations):
        if net.stations[i].start_date is None or net.stations[i].end_date is None:
            print(net.stations[i].code, net.stations[i].start_date, net.stations[i].end_date)

    return

def read_response(resp_file):

    #resp_obj = obspy.read_inventory('/run/media/felixw/workdisc_1/6D6-Trillium-250sps.resp', format='RESP')
    inv_obj = obspy.read_inventory(resp_file, format='RESP')

    print (type(inv_obj))

    # write out to check: inv_obj.write("test_resp.xml", format="stationxml")

    #<Channel code="BHZ" startDate="2015-01-01T00:00:00.000000Z" locationCode="">
    datetime = UTCDateTime("2015-01-01T00:00:00")
    resp_obj = inv_obj.get_response("XX.XX..BHZ", datetime)

    print(resp_obj)

    return  resp_obj

# =============================================
# Section for quick run of this script's functions
# ---------------------------------------------
if __name__ == "__main__":
    # call main function
    #check_invenory(sys.argv[1])

    resp=read_response("/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/Equipments/RESP.COMPACT120S.MINIMUS.txt")

    print("The returned object type: ", type(resp))

