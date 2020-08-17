#! /usr/bin/env python
"""
Class Model for GA Extra Metadata

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
import sys
import io
import json
import pandas as pd
import obspy
from obspy.core.util import AttribDict
from obspy.core import UTCDateTime

from seismic.inventory import add_time_corrections


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
                      }  # initially dict for net.sta

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

        # select the network and sta
        pdf2 = pdf.loc[(pdf['net'] == self.net) & (pdf['sta'] == self.sta)].copy()

        # drop two columns inplace pdf2 itself will be changed, otherwise will return a new df
        pdf2.drop(['net', 'sta'], axis=1, inplace=True)
        print("The shapes = ", pdf.shape, pdf2.shape)

        # to json object
        gps_corr = pdf2.to_json(orient="records", date_format="epoch", double_precision=10,
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


    def get_metadata_from_csv(self, csv_data):  # ,net, sta, start_dt, end_dt):
        """
        Select the csv rows according  net, sta, start_dt, end_dt

        Args:
            csv_data: a string of CSV
            start_dt: obspy.core.utcdatetime.UTCDateTime
            end_dt:  obspy.core.utcdatetime.UTCDateTime

        Returns: a subset of csv in pandas df, selected according to (net, sta, start_dt, end_dt)
        """
        pdf = pd.read_csv(io.StringIO(csv_data))
        # pdf.insert(4,"utcdate", UTCDateTime(0))
        # print(pdf.head())
        pdf["utcdate"] = pdf.apply(lambda row: UTCDateTime(row.date), axis=1)

        _crit = (pdf['net'] == self.net) & (pdf['sta'] == self.sta) & (
                pdf['utcdate'] > self.start) & (pdf['utcdate'] < self.end)
        pdf2 = pdf.loc[_crit].copy()  # use copy() to fix "SettingWithCopyWarning"

        # drop columns inplace pdf2 itself will be changed, otherwise will return a new df
        pdf2.drop(['net', 'sta', 'utcdate'], axis=1, inplace=True)
        print("The shapes = ", pdf.shape, pdf2.shape)

        # to json object
        gps_corr = pdf2.to_json(orient="records", date_format="epoch", double_precision=10,
                                force_ascii=True, date_unit="ms", default_handler=None, indent=2)

        # add the json object as string to the self.mdata
        self.mdata.update({"GPS_CORRECTION": json.loads(gps_corr)})

        # print(pdf2.head())
        # return pdf2.to_csv(index=False)
        return pdf2

    def metadata2json(self, metadta_json_file):
        """
        write out the current metatdaya object into a file
        :param metadta_json_file:
        :return:
        """
        with open( metadta_json_file, "w") as f:
            json.dump(self.mdata, f, indent=2)

        return metadta_json_file

    def validate_json_schema(self):
        """
        validate the current json by https://github.com/Julian/jsonschema
        :return: boolean
        """

        return True

# ======================================================================================================================
# Example How to Run:
# python station_metadata_extra.py OA.CF28_clock_correction.csv /Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml /tmp
# =====================================================================================================================
if __name__ == "__main__":

    USAGE = "python %s gps_clock_corr_csv station_inventory_xml [out_dir]" % sys.argv[0]

    GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"

    # extra metadata info file(s) to be read and formatted into JSON
    in_csv_file = "./OA.CF28_clock_correction.csv"

    # Original station inventory XML file to be modified
    # /g/data/ha3/Passive/_AusArray/OA/ASDF_cleaned/OA_stations_2017-2018.xml
    in_station_xml_file = "/Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml"

    # output dir for modified station inventory xml files
    out_dir = "/tmp/"

    if len(sys.argv) < 3:
        print(USAGE)
        sys.exit(1)
    else:
        in_csv_file = sys.argv[1]
        in_station_xml_file = sys.argv[2]

    if len(sys.argv) >= 4:
        out_dir = sys.argv[3]
    else:
        out_dir = None

    # get the metadata and it's associated network.station
    (net, sta, csv_data) = add_time_corrections.get_csv_correction_data(in_csv_file)

    _inv = obspy.read_inventory(in_station_xml_file, format='STATIONXML')
    selected_inv = _inv.select(network=net, station=sta)

    # selected_inv may be 0,1, 2, multiple stations, each have a start_date end_date
    station_list = selected_inv.networks[0].stations
    if station_list is None or len(station_list) == 0:  # no further process for this dummy case
        print("(network station) %s.%s NOT present in the XML file %s" % (net, sta, in_station_xml_file))
        sys.exits(1)

    # redefine the selected_inv
    for a_station in station_list:  # loop over all Stations
        # get station star end date and split csv_data
        start_dt = a_station.start_date
        end_dt = a_station.end_date

        ajson = StationMetadataExtra(net, sta, start_datetime=start_dt, end_datetime=end_dt)

        # generate/format extra metadata from inputs
        mpdf = ajson.get_metadata_from_csv(csv_data)
        # updated the ajson object with more metadata, such as orientation corr

        ajson.metadata2json(os.path.join(out_dir,"%s.%s_%s_extra_metadata.json"%(net,sta,str(start_dt))))

        # Now, ready to write the ajson obj into new xml file
        mformat = "JSON"

        my_tag = AttribDict()
        my_tag.namespace = GA_NameSpace

        my_tag.value = ajson.make_json_string()  # store all the extra metadata into a json string.

        a_station.extra = AttribDict()
        a_station.extra.GAMetadata = my_tag

    # prepare to write out a modified xml file

    stationxml_with_extra = '%s.%s_station_metadata_%s.xml' % (net, sta, mformat)

    if out_dir is not None and os.path.isdir(out_dir):
        stationxml_with_extra = os.path.join(out_dir, stationxml_with_extra)

    selected_inv.write(stationxml_with_extra, format='STATIONXML',
                       nsmap={'GeoscienceAustralia': GA_NameSpace})
