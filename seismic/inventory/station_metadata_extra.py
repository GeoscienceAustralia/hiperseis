#! /usr/bin/env python
"""
Json Class Model for GA Extra Metadata

References:
    https://gajira.atlassian.net/browse/PV-324
    https://gajira.atlassian.net/browse/PV-312

CreationDate:   07/08/2020
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     07/08/2020   FZ
    LastUpdate:

"""
import os
import json
import pandas as pd
import obspy
from obspy.core.util import AttribDict


class StationMetadataExtra:  # CapWords naming convention.
    """
    Class Model for GA Extra Metadata
    """

    def __init__(self, network_code, station_code, start_datetime=None, end_datetime=None):
        """
         A class to encapsulate extra metadata items which not defined in FDSN stationXML https://www.fdsn.org/xml/station/
        :param network_code:
        :param station_code:
        :param start_datetime:
        :param end_datetime:
        """
        self.net = network_code
        self.sta = station_code
        self.start = start_datetime
        self.end = end_datetime

        self.mdata = {"network": self.net,
                      "station": self.sta
                      }  # initially dict

    def make_json_string(self):
        """
        format this metadata object (dict) into a json string
        :return: a json string
        """
        json_str = self.mdata
        return json.dumps(json_str, indent=2)

    def add_gps_corrections(self, csvfile):
        """
        read from the csvfile, make a gps correction list for this net sta
        :param csvfile:
        :return: json_sub_node, mwhich was added into self.mdata
        """
        pdf = pd.DataFrame(pd.read_csv(csvfile, sep=",", header=0, index_col=False))

        # to json object
        gps_corr = pdf.to_json(orient="records", date_format="epoch", double_precision=10,
                               force_ascii=True, date_unit="ms", default_handler=None, indent=4)

        # add the json object as string to the self.mdata
        self.mdata.update({"GPS_CORRECTION": json.loads(gps_corr)})

        return gps_corr

    def add_orientation_correction(self, infile):
        """
        add json element for orientation_correction from input data file infile
        :param infile: an input file of certain format, containing the required orientation_correction
        :return: json_sub_node, which was added into self.mdata
        """
        return True

    def add_extra_metata_to_stationxml(self, in_station_xml_file, out_dir=None):
        """
        Add the current json-string (rep extra metadata) into a given station XML file
        :param in_station_xml_file:  path to station_xml
        :return: new_station_xml
        """
        # read in

        ns = "https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.1"

        # path2_myxml = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"
        my_inv = obspy.read_inventory(in_station_xml_file, format='STATIONXML')

        # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.select.html#obspy.core.inventory.inventory.Inventory.select

        selected_inv = my_inv.select(network=self.net, station=self.sta)

        # print(selected_inv)

        my_tag = AttribDict()
        my_tag.namespace = ns
        my_tag.value = self.make_json_string()

        selected_inv.networks[0].stations[0].extra = AttribDict()
        selected_inv.networks[0].stations[0].extra.GA_META = my_tag

        new_stationxml = '%s.%s_station_inv_modified.xml' % (self.net, self.sta)

        if out_dir is not None and os.path.isdir(out_dir):
            new_stationxml_with_json = os.path.join(out_dir, new_stationxml)
        else:
            new_stationxml_with_json = new_stationxml

        selected_inv.write(new_stationxml_with_json, format='STATIONXML', nsmap={'GeoscienceAustralia': ns})

        # my_inv.write('modified_inventory.xml', format='STATIONXML')

        return new_stationxml_with_json

    def validate_json_schema(self):
        """
        validate the current json by https://github.com/Julian/jsonschema
        :return: boolean
        """

        return True


# ===============================================================================
if __name__ == "__main__":
    obj = StationMetadataExtra("OA", "CF28")
    print(obj.make_json_string(), type(obj.make_json_string()))

    gpsc = obj.add_gps_corrections("./OA.CF28_clock_correction.csv")

    # console print(json.dumps(myInst.mdata,indent=2))

    with open("./extra_mdata.dev.json", "w") as f:
        json.dump(obj.mdata, f, indent=2)

    obj.add_extra_metata_to_stationxml("/Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml")
