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
import json
import pandas as pd

def get_extra_metadata(input_station_xml, net, sta, sta_index=0):
    """
    use obspy to read in an input stationxml file that may include GA extra metadata in JSON format.
    Then get the details of the extra metadata
    Ref: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/notebooks/read_parse_stationxml_with_json.ipynb

    Args:
        input_station_xml: a station XML file
        net:  network code
        sta; station code
        sta_index: 0,1,2,.. in case of multiple same-coded stations, configured differently over non-overlapping dates.

        Returns: Json-formatted extra metadata
    """
    _inv = obspy.read_inventory(input_station_xml, format='STATIONXML')

    stn_meta = _inv.select(network=net, station=sta)  # station code name eg, "CF28"

    selected_stations = stn_meta.networks[0].stations

    if sta_index >= len(selected_stations):
        print("There are %s stations with the code %s" % (len(selected_stations), sta))
        print("The sta_index is out-of-the-range", sta_index)
        return None

    a_station = selected_stations[sta_index]  # more than one stations due to start stop time differences

    show_station_prop(a_station)

    extra_meta = a_station.extra.GAMetadata.value

    print(type(extra_meta))

    return extra_meta


def show_station_prop(astation):
    print("********** Show a station's properties *********** ")

    print("The station name code", astation.code)
    print("The station creation and termination date", astation.creation_date, astation.termination_date)

    print("The station start and end date", astation.start_date, astation.end_date)

    print(astation.equipments)

    return

def get_gps_correction(extra_meta):

    mdata = json.loads(extra_meta)

    print(type(mdata))  # a dictionry
    print(mdata.keys())

    # print(type(mdata['GPS_CORRECTION']))  # a list of dict
    # for corr in mdata['GPS_CORRECTION']:
    #     print(corr["date"], corr["clock_correction"], "seconds")

    df_clock_correction = pd.DataFrame(mdata['GPS_CORRECTION'])

    return df_clock_correction

def get_orientation_correction(extra_meta):
    mdata = json.loads(extra_meta)

    print(mdata["ORIENT_CORRECTION"])

    return mdata["ORIENT_CORRECTION"]['azimuth_correction']


def check_obspy_functions_on_inventoryxml(in_stationxml, output_dir="/tmp/"):
    """
    Check the obspy IO-handling of inventory xml file:
    read-in existing stationXML file, then write out a new xml.
    Check if "identical".
    It was found due to obspy version difference, an issue is observed: new stationXML missing some startDate for some stations!!
    :return:
    """
    big_inv = obspy.read_inventory(in_stationxml)  # read in fresh

    filename_new_xml = "check1_" + os.path.basename(in_stationxml)

    new_big_xml = os.path.join(output_dir, filename_new_xml)

    print("New XML file = ", new_big_xml)

    # Write out the inventory into a new file to see if they are identical?
    # big_inv.write(new_big_xml, format="stationxml", validate=True) #without nsmap

    GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"
    big_inv.write(new_big_xml, format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                  validate=True)

    # ObsPy Problem:  Missing some startDate for some stations
    # Why the obspy missed some start_date for stations?

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


def read_response_from_txtfile(resp_file):
    """
    Use obspy to read in Equipments/RESP.COMPACT120S.MINIMUS.txt
    and write it out in stationxml format.
    Args:
        resp_file: a standard response file like Equipments/RESP.COMPACT120S.MINIMUS.txt

    Returns: Obspy RESPONSE object

    """

    # resp_obj = obspy.read_inventory('/run/media/felixw/workdisc_1/6D6-Trillium-250sps.resp', format='RESP')
    inv_obj = obspy.read_inventory(resp_file, format='RESP')

    print(type(inv_obj))

    # write out to check:

    inv_obj.write("test_resp.xml", format="stationxml")

    # <Channel code="BHZ" startDate="2015-01-01T00:00:00.000000Z" locationCode="">
    datetime = UTCDateTime("2015-01-01T00:00:00")
    resp_obj = inv_obj.get_response("XX.XX..BHZ", datetime)
    # this resp_obj can be used in stationXML file writing.

    print(resp_obj)

    return resp_obj


# =============================================
# Section for quick run of this script's functions
# (hiperseispy37) fxz547@vdi-n25 /g/data/ha3/fxz547/Githubz/hiperseis (develop)
# $ python seismic/inventory/read_extra_metadata_from_stationxml.py /g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/NewInventoryXML/OA_stations_2017-2018_with_GAextra_responses.xml
# ---------------------------------------------
if __name__ == "__main__":
    # call main function
    # check_obspy_functions_on_inventoryxml(sys.argv[1])

    station_code = "CF28"
    extra_meta = get_extra_metadata(sys.argv[1], net="OA", sta= station_code, sta_index=0)
    orc = get_orientation_correction(extra_meta)

    pds = get_gps_correction(extra_meta)

    # query the pandas dataframe (CSV) to obtain correction for a specific date
    adate = "2018-06-06"
    query4date = pds.loc[pds["date"] == adate].clock_correction
    print(query4date.values)

    extra_meta = get_extra_metadata(sys.argv[1],  net="OA", sta= station_code,  sta_index=1)
    orc = get_orientation_correction(extra_meta)

    pds = get_gps_correction(extra_meta)
    # query the pandas dataframe (CSV) to obtain correction for a specific date
    adate = "2018-08-08"
    query4date = pds.loc[pds["date"] == adate].clock_correction
    print(query4date.values)


    extra_meta = get_extra_metadata(sys.argv[1],  net="OA", sta= station_code,  sta_index=2)
    orc = get_orientation_correction(extra_meta)
    print("Orientation Correcton: ", orc)

    pds = get_gps_correction(extra_meta)
    # query the pandas dataframe (CSV) to obtain correction for a specific date
    query4date = pds.loc[pds["date"] == adate].clock_correction
    print(query4date.values)


    extra_meta = get_extra_metadata(sys.argv[1], net="OA", sta=station_code, sta_index=3)  # could be out of range

    # Response txt file
    # resp=read_response_from_txtfile("/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/Equipments/RESP.COMPACT120S.MINIMUS.txt")
    # print("The returned object type: ", type(resp))
