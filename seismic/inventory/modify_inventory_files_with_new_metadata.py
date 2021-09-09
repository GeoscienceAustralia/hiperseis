#!/usr/bin/env python
# coding: utf-8
"""
Use Obspy to modify inventory XML files to include new metadata
(See also the notebook sandbox/modify_inventory_files_with_new_metadata.ipynb)

Fei Zhang
Created on 2020-09-23

"""

import os
import sys
import obspy
from obspy.core.inventory import Inventory
from obspy.core.util import AttribDict
from obspy.core import UTCDateTime

from seismic.inventory.extract_equipments_from_csv import EquipmentExtractor
from seismic.inventory.station_metadata_extra import StationMetadataExtra, get_csv_correction_data, get_orientation_corr

GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"  # global name space for XML

class InvXML_Modifier:
    def __init__(self, input_inv_xml, output_dir_path):
        self.input_xml = input_inv_xml  # ORIG_INVENTORY File
        self.output_dir = output_dir_path

        self.inv_obj = obspy.read_inventory( self.input_xml )  # read in inventory object

    def split_inventory(self, outdir = None):
        """
        Split the big inventory station XML file into multiple one-station.xml file, according the station_codes
        Args:
            big_inv: obspy inventory obj
        Returns:
        """

        big_inv = self.inv_obj

        if outdir is None:
            my_output_dir = self.output_dir
        else:
            my_output_dir = outdir

        for a_net in big_inv.networks:
            for a_sta in a_net.stations:
                print(a_net.code, a_sta.code)  # This print 328 (OA, sta) pairs - they are NOT unique pairs!!!

                a_inv = big_inv.select(network=a_net.code, station=a_sta.code)  # .copy()

                # modify station metadata
                # my_sensor=obspy.core.inventory.util.Equipment(type="Sensor", description="Nanometrics Trillium Compact 120s",serial_number="004940")
                # my_digitizer = obspy.core.inventory.util.Equipment(type="Digitizer", description="Guralp Minimus",serial_number="MIN-A456")
                # a_sta.equipments = [my_sensor, my_digitizer]

                sta_file_name = "%s_%s_station.xml" % (a_net.code, a_sta.code)

                outxml = os.path.join(my_output_dir, sta_file_name)
                a_inv.write(outxml, format="stationxml", validate=True)  # nsmap={'GeoscienceAustralia': GA_NameSpace})
                # 119 OA_*_station.xml file written. Some files written 3 times!!
                # fzhang@zubuntu1804 /Datasets/InventoryXml/OA_stations_2017-2018
                #  $ ls *station.xml| wc
                #     119     119    2380

        return my_output_dir

    def write_new_version_inventory(self, new_xml_file = None):
        """
        Construct a new inventory object of networks, making use of new obspy version and new attributes for Inventory
        re-write the input xml file in to new xml file

                # The new network inventory has new <Module>ObsPy 1.2.1 and new <Source>
        :return: path of new xml file
        """

        if new_xml_file is None:
            filename_newv_xml = "new_version_"+os.path.basename(self.input_xml)
            new_version_xml = os.path.join(self.output_dir, filename_newv_xml)
        else:
            new_version_xml = new_xml_file

        print("Target New XML file = ", new_version_xml)

        inv2 = Inventory(
            # We'll add networks later.
            networks=[],
            # The source should be the id whoever create the file.
            source="Geoscience Australia EFTF AusArray")

        for a_net in self.inv_obj.networks:
            # Re-write each network of the big inventory into the new object inv2
            # inv2.networks = []
            inv2.networks.append(a_net)
            print("The network %s has %s stations."%(a_net.code, len(a_net.stations)))

        inv2.write(new_version_xml, format="stationxml", validate=True)  # every Station got equipment

        return new_version_xml

    def modify_invenory(self, gps_clock_corr_csv=None, orient_corr_json = None, equipment_csv= None ):
        """
        Modify the existing station XML files to include new metadata:
        - add equipment sensor digitizer
        - add extra metadata: GPS correction
        - add extra metadata: Orientation correction
        Args:

        Returns: the final station_xml file modified with new metadata: inv2_xml_file

        """

        # Construct a new inventory object of networks.
        # This will use new obspy version and new attributes:
        inv2 = Inventory(
            # We'll add networks later.
            networks=[],
            # The source should be the id whoever create the file.
            source="Geoscience Australia EFTF AusArray PST")


        # output dir for modified station inventory xml files
        out_dir = self.output_dir  # "/home/fzhang/tmpdir"

        net, sta, csv_data = get_csv_correction_data(gps_clock_corr_csv)
        net_sta, oricorr_json_data = get_orientation_corr(orient_corr_json)
        my_equip_obj = EquipmentExtractor(csvfile=equipment_csv)

        big_inv = self.inv_obj

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
                    # sensor_desc = "NA Sensor for (%s,%s)" % (a_net.code, a_sta.code)
                    sensor_desc = "Nanometrics Trillium Compact 120s"
                    sensor_sernumb = "N/A"

                _digitizers = my_equip_obj.get_digitizer(a_net.code, a_sta.code)
                if len(_digitizers) > 0:
                    dig_desc = _digitizers[0].get("Description")
                    dig_sernumb = _digitizers[0].get("SerNumber")
                else:
                    print("%s %s  No digitizers !" % (a_net.code, a_sta.code))
                    #dig_desc = "NA Digitizer for (%s,%s)" % (a_net.code, a_sta.code)
                    dig_desc = "Guralp Minimus"
                    dig_sernumb = "N/A"

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
            # sta_file_name2 = "%s_%s_station2.xml"%(a_net.code, a_sta.code)
            # # OA_CE28 was written 3-times!!!!!! due to multiple (OA,CE28)-station-nodes
            # There will be 119 xml files written in this loop of 328 items. However, the final results missed 119 equipments!!
            # outxml2 = os.path.join(OUTPUT_DIR, sta_file_name2)
            #
            # inv2.networks = a_inv.networks
            #
            # inv2.write(outxml2,format="stationxml", validate=True) # nsmap={'GeoscienceAustralia': GA_NameSpace})

            # After the modification of ALL the station objects,
            # write the big inventory in new object inv2
            inv2.networks = []
            inv2.networks.append(a_net)
            inv2_xml_file = os.path.join(out_dir, a_net.code + "_stations2.xml")
            inv2.write(inv2_xml_file, format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                       validate=True)  # every Station got equipment

            # Add responses:
            resp_obj = read_response()
            self.add_response_into_stationxml(inv2, resp_obj)

            # and the original write out again to check what has been modified?
            post_orig =  os.path.join(out_dir, a_net.code + "_stations_post_orig.xml")
            big_inv.write(post_orig, format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                          validate=True)  # also has the Sensors etc

            return inv2_xml_file

    def add_response_into_stationxml(self, invent_obj, resp_obj, out_xml_file="station_inventory_with_responses.xml"):
        """
        write the resp_obj into each channel of the invent_obj
        Args:
            invent_obj: station_xml inventory object
            resp_obj: response Object

        Returns: update_invent_obj
        """

        for anet in invent_obj.networks:
            for astation in anet.stations:
                for acha in astation.channels:
                    acha.response=resp_obj

        invent_obj.write(out_xml_file, format="stationxml", nsmap={'GeoscienceAustralia': GA_NameSpace},
                       validate=True)

        return out_xml_file


    def check_invenory(self):
        """
        Check the obspy IO-handling of inventory xml file
        :return:
        """
        big_inv = obspy.read_inventory(self.input_xml)  # read in fresh

        filename_new_xml = "check_new_" + os.path.basename(self.input_xml)

        new_big_xml = os.path.join(self.output_dir, filename_new_xml)

        print("New XML file = ", new_big_xml)

        # Write out the inventory into a new file to see if they are identical?
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
            print("The total number of stations-nodes =", number_of_stations)

        # identify None start or end Date
        for i in range(number_of_stations):
            if net.stations[i].start_date is None or net.stations[i].end_date is None:
                print(net.stations[i].code, net.stations[i].start_date, net.stations[i].end_date)

        return "True-False"


def read_response(resp_file ="/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/Equipments/RESP.COMPACT120S.MINIMUS.txt" ):
    """
    read_response from a standard txt file that obspy can understand
    Args:
        resp_file: path2 the response.txt file

    Returns: Response Obj

    """

    inv_obj = obspy.read_inventory(resp_file, format='RESP')

    # write out to check:
    # print (type(inv_obj))
    # inv_obj.write("test_resp.xml", format="stationxml")

    # get the Response object
    # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.get_response.html
    datetime = UTCDateTime("2015-01-01T00:00:00")
    resp_obj = inv_obj.get_response("XX.XX..BHZ", datetime)

    print("The returned object type: ", type(resp_obj), resp_obj)

    return resp_obj

#---------------------------------------------------------------------
# Modify some input file paths if necessary 
# source gadi_env.sh
# python seismic/inventory/modify_inventory_files_with_new_metadata.py
# TODO: Refactor the input file paths as commandline input
#---------------------------------------------------------------------

if __name__ == "__main__":

    print("The Obspy Version is ", obspy.__version__)

    METADB_DIR = "/Datasets/Station_Extra_Metadata" #/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata"
    METADB_DIR = "/g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata"

    # ORIG_INVENTORY_FILE = os.path.join(METADB_DIR, "SrcInventoryXML", "7D_2012_2013.xml")
    # ORIG_INVENTORY_FILE = os.path.join(METADB_DIR, "SrcInventoryXML", "7X_2009_2011_ASDF.xml")
    # "/Datasets/InventoryXml/OA_stations_2017-2018.xml"
    ORIG_INVENTORY_FILE = os.path.join(METADB_DIR, "SrcInventoryXML", "OA_stations_2017-2018.xml")

    # extra metadata info file(s) to be read and formatted into JSON
    # "/Datasets/GPS_ClockCorr/OA.CE28_clock_correction.csv" # "./OA.CF28_clock_correction.csv"
    in_csv_file = os.path.join(METADB_DIR, "Ody_ClockCorr2021Sept.csv") # All_GPSClock_Correction.csv")

    # "/Datasets/Orientation_Correction_json/OA_ori_error_estimates.json"
    in_json_file = os.path.join(METADB_DIR, "OrientationCorr", "OA_ori_error_estimates.json")

    in_equip_csv = os.path.join(METADB_DIR, "Equipments/FieldSiteVisitLive.csv")

    # output
    OUTPUT_DIR = os.path.join(METADB_DIR, "Output_Dir")  # "/Datasets/InventoryXml/OA_stations_2017-2018"

    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    my_obj = InvXML_Modifier(ORIG_INVENTORY_FILE, OUTPUT_DIR)

    my_obj.check_invenory()
    my_obj.split_inventory()
    my_obj.write_new_version_inventory()
    
    # provide the right input data
    my_obj.modify_invenory(gps_clock_corr_csv=in_csv_file, orient_corr_json=in_json_file,equipment_csv=in_equip_csv)
