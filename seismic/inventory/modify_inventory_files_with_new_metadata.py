#!/usr/bin/env python
# coding: utf-8
"""
Use Obspy to modify inventory XML files to include new metadata
(See also the notebook sandbox/modify_inventory_files_with_new_metadata.ipynb)

Fei Zhang
created on 2020-09-23

* modify station xml to inlcude new GA-metadata
* split a big xml file into multiple single-station xml files

"""

import os
import sys

import obspy
from obspy.core.inventory import Inventory
from obspy.core.util import AttribDict

sys.path.append("/home/fzhang/Githubz/hiperseis")
sys.path.append("/g/data/ha3/fxz547/Githubz/hiperseis/")

from seismic.inventory.extract_equipments_from_csv import EquipmentExtractor

from seismic.inventory.station_metadata_extra import StationMetadataExtra, get_csv_correction_data, get_orientation_corr

METADB_DIR = "/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata"

ORIG_INVENTORY_FILE = os.path.join(METADB_DIR, "SrcInventoryXML", "OA_stations_2017-2018.xml")
# "/Datasets/InventoryXml/OA_stations_2017-2018.xml"

# output
OUTPUT_DIR = os.path.join(METADB_DIR, "Output_Dir")  # "/Datasets/InventoryXml/OA_stations_2017-2018"

if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)


def split_inventory(big_inv):
    """
    Split the big inventory station XML file into multiple one-station.xml file, according the station CODE
    Args:
        big_inv: obspy inventory obj

    Returns:

    """

    for a_net in big_inv.networks:
        for a_sta in a_net.stations:
            print(a_net.code, a_sta.code)  # This print 328 (OA, sta) pairs - they are NOT unique pairs!!!

            a_inv = big_inv.select(network=a_net.code, station=a_sta.code)  # .copy()

            # modify station metadata
            #         my_sensor=obspy.core.inventory.util.Equipment(type="Sensor", description="Nanometrics Trillium Compact 120s",serial_number="004940")
            #         my_digitizer = obspy.core.inventory.util.Equipment(type="Digitizer", description="Guralp Minimus",serial_number="MIN-A456")
            #         a_sta.equipments = [my_sensor, my_digitizer]

            sta_file_name = "%s_%s_station.xml" % (a_net.code, a_sta.code)

            outxml = os.path.join(OUTPUT_DIR, sta_file_name)
            a_inv.write(outxml, format="stationxml", validate=True)  # nsmap={'GeoscienceAustralia': GA_NameSpace})
            # 119 OA_*_station.xml file written. Some files written 3 times!!
            # fzhang@zubuntu1804 /Datasets/InventoryXml/OA_stations_2017-2018
            #  $ ls *station.xml| wc
            #     119     119    2380


def write_new_version_inventory(big_inv):
    # Construct a new inventory object of networks.
    # This will use new obspy version and new attributes for Inventory

    filename_newv_xml = os.path.basename(ORIG_INVENTORY_FILE).replace(".xml", "_new_version.xml")

    new_version_xml = os.path.join(OUTPUT_DIR, filename_newv_xml)

    print("New XML file = ", new_version_xml)

    inv2 = Inventory(
        # We'll add networks later.
        networks=[],
        # The source should be the id whoever create the file.
        source="Geoscience Australia EFTF AusArray")

    for a_net in big_inv.networks:
        # Re-write each network of the big inventory into the new object inv2
        # inv2.networks = []
        inv2.networks.append(a_net)
        inv2.write(new_version_xml, format="stationxml", validate=True)  # every Station got equipment

        print("New Version XML file = ", new_version_xml, len(a_net.stations))

    # The new network inventory has new <Module>ObsPy 1.2.1 and new <Source>


def modify_invenory(big_inv):
    """
    Modify the existing station XML files inclding new metadata:
    - add sensor digitizer
    - add extra metadata: GPS correction
    - add extra metadata: Orientation correction
    Args:
        big_inv:

    Returns:

    """

    # Construct a new inventory object of networks.
    # This will use new obspy version and new attributes:
    inv2 = Inventory(
        # We'll add networks later.
        networks=[],
        # The source should be the id whoever create the file.
        source="Geoscience Australia EFTF AusArray PST")

    GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"

    # Original station inventory XML file to be modified
    # /g/data/ha3/Passive/_AusArray/OA/ASDF_cleaned/OA_stations_2017-2018.xml
    # in_station_xml_file = "/Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml"

    # extra metadata info file(s) to be read and formatted into JSON
    in_csv_file = os.path.join(METADB_DIR, "OA_All_GPSClock_Correction.csv")
    # "/Datasets/GPS_ClockCorr/OA.CE28_clock_correction.csv" # "./OA.CF28_clock_correction.csv"

    in_json_file = os.path.join(METADB_DIR, "OrientationCorr", "OA_ori_error_estimates.json")
    # "/Datasets/Orientation_Correction_json/OA_ori_error_estimates.json"
    # output dir for modified station inventory xml files
    out_dir = OUTPUT_DIR
    # "/home/fzhang/tmpdir"

    net, sta, csv_data = get_csv_correction_data(in_csv_file)
    net_sta, oricorr_json_data = get_orientation_corr(in_json_file)

    my_equip_obj = EquipmentExtractor(csvfile=os.path.join(METADB_DIR, "Equipments/FieldSiteVisitLive.csv"))

    for a_net in big_inv.networks:

        print("The number of station-nodes in the network =", len(a_net.stations))

        for a_sta in a_net.stations:
            # print(a_net.code, a_sta.code)  # this contains 328 pairs, but they are NOT unique, station code may repeat.

            a_inv = big_inv.select(network=a_net.code, station=a_sta.code)  # .copy appears to have no effect here

            # print (a_sta.code, " stations has %s channels"%len(a_sta))

            _sensors = my_equip_obj.get_sensors(a_net.code, a_sta.code)
            if len(_sensors) > 0:
                sensor_desc = _sensors[0].get("Description")
                sensor_sernumb = _sensors[0].get("SerNumber")
            else:
                print("%s %s  No sensors !" % (a_net.code, a_sta.code))
                sensor_desc = "NA Sensor for (%s,%s)" % (a_net.code, a_sta.code)
                sensor_sernumb = "NA"

            _digitizers = my_equip_obj.get_digitizer(a_net.code, a_sta.code)
            if len(_digitizers) > 0:
                dig_desc = _digitizers[0].get("Description")
                dig_sernumb = _digitizers[0].get("SerNumber")
            else:
                print("%s %s  No digitizers !" % (a_net.code, a_sta.code))
                dig_desc = "NA Digitizer for (%s,%s)" % (a_net.code, a_sta.code)
                dig_sernumb = "NA"

            # modify station metadata
            my_sensor = obspy.core.inventory.util.Equipment(type="Sensor", description=sensor_desc,
                                                            serial_number=sensor_sernumb)

            # my_digitizer = obspy.core.inventory.util.Equipment(type="Digitizer", description="Guralp Minimus",serial_number="MIN-A456")
            my_digitizer = obspy.core.inventory.util.Equipment(type="Digitizer", description=dig_desc,
                                                               serial_number=dig_sernumb)

            a_sta.equipments = [my_sensor, my_digitizer]

            # get station start_ end_date and split csv_data
            start_dt = a_sta.start_date
            end_dt = a_sta.end_date

            ajson = StationMetadataExtra(a_net.code, a_sta.code, start_datetime=start_dt, end_datetime=end_dt)

            # generate/format extra metadata from inputs
            mpdf = ajson.add_gps_correction_from_csv(csv_data)

            # updated the ajson object with more metadata, such as orientation corr
            ajson.add_orientation_correction(oricorr_json_data)

            ajson.write_metadata2json(
                os.path.join(out_dir, "%s.%s_%s_extra_metadata.json" % (a_net.code, a_sta.code, str(start_dt))))

            # Now, ready to write the ajson obj into new xml file
            mformat = "JSON"

            my_tag = AttribDict()
            my_tag.namespace = GA_NameSpace

            my_tag.value = ajson.make_json_string()  # store all the extra metadata into a json string.

            a_sta.extra = AttribDict()
            a_sta.extra.GAMetadata = my_tag

            # prepare to write out a modified xml file
            stationxml_with_extra = '%s.%s_station_metadata_%s.xml' % (a_net.code, a_sta.code, mformat)

            if out_dir is not None and os.path.isdir(out_dir):
                stationxml_with_extra = os.path.join(out_dir, stationxml_with_extra)

            a_inv.write(stationxml_with_extra, format='STATIONXML',
                        nsmap={'GeoscienceAustralia': GA_NameSpace})
        # Problem:
        #         sta_file_name2 = "%s_%s_station2.xml"%(a_net.code, a_sta.code)
        #         # OA_CE28 was written 3-times!!!!!! due to multiple (OA,CE28)-station-nodes
        #         There will be 119 xml files written in this loop of 328 items. However, the final results missed 119 equipments!!
        #         outxml2 = os.path.join(OUTPUT_DIR, sta_file_name2)

        #         inv2.networks = a_inv.networks

        #         inv2.write(outxml2,format="stationxml", validate=True) # nsmap={'GeoscienceAustralia': GA_NameSpace})

        # After the modification of ALL the station objects,
        # write the big inventory in new object inv2
        inv2.networks = []
        inv2.networks.append(a_net)
        inv2.write(a_net.code + "_stations2.xml", format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                   validate=True)  # every Station got equipment

        # and the original write out again to check if it has been modified
        big_inv.write(a_net.code + "_stations_post.xml", format="stationxml",
                      nsmap={'GeoscienceAustralia': GA_NameSpace},
                      validate=True)  # also has the Sensors etc


if __name__ == "__main__":

    print(obspy.__version__)

    big_inv = obspy.read_inventory(ORIG_INVENTORY_FILE)

    os.path.basename(ORIG_INVENTORY_FILE)

    filename_new_xml = os.path.basename(ORIG_INVENTORY_FILE).replace(".xml", "_new.xml")

    new_big_xml = os.path.join(OUTPUT_DIR, filename_new_xml)

    print("New XML file = ", new_big_xml)

    # Write the inventory into a new file to see if they are identical?
    big_inv.write(new_big_xml, format="stationxml", validate=True)
    # Problem:  Missing a lot of startDate for some stations
    # Why the obspy read_invnetory() missed some start_date for stations?

    # grep "<Station" OA_stations2.xml
    # grep "<Station" OA_stations_2017-2018_new.xml
    #     <Station code="BS24" startDate="2017-09-27T05:16:47.000000Z" endDate="2018-07-04T04:38:50.000000Z">
    #     <Station code="BS25" endDate="2018-06-27T11:44:37.000000Z">
    #     <Station code="BS26" startDate="2017-09-27T00:49:10.000000Z" endDate="2018-01-12T04:02:13.000000Z">
    #     <Station code="BS27" startDate="2017-09-26T00:42:49.000000Z" endDate="2018-07-03T03:04:17.000000Z">

    for net in big_inv.networks:
        print("The Network Code", net.code)
        number_of_stations = len(net.stations)
        print("The total number of stations/nodes =", number_of_stations)

        for i in range(number_of_stations):
            print(net.stations[i].code, net.stations[i].start_date, net.stations[i].end_date)

    modify_invenory(big_inv)
