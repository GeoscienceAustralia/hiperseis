#! /usr/bin/env python
"""
Class Model for GA Extra Metadata to be included into station inventory XML files

CreationDate:   07/08/2020
Developer:      fei.zhang@ga.gov.au

"""
import os
import sys
import io
import json
import pandas as pd
import obspy
from obspy.core.util import AttribDict
from obspy.core import UTCDateTime


class StationMetadataExtra:  # CapWords naming convention.
    """
    Class Model for GA Extra Metadata
    """

    def __init__(self, network_code, station_code, start_datetime=None, end_datetime=None):
        """
         A class to encapsulate extra metadata items which are not defined in the standard of FDSN
         https://www.fdsn.org/xml/station/

        :param network_code:
        :param station_code:
        :param start_datetime:
        :param end_datetime:
        """
        self.net = network_code
        self.sta = station_code
        self.start = start_datetime
        self.end = end_datetime

        # A dictionary, skeleton JSON model, for a given network.station
        self.mdata = {
            "network": self.net,
            "station": self.sta,

            # GPS clock corrrections
            "GPS_CORRECTION": [],

            # Orientation correction "azimuth_correction" applicable for datetime interval (start_dt, end_dt)
            "ORIENT_CORRECTION": {
                # "start_dt": "2017-11-07T09:07:34.930000Z",
                # "end_dt": "2018-08-23T03:52:29.528000Z",
                # "azimuth_correction": -5.0
            },
        }

    def make_json_string(self):
        """
        format this metadata object (dict) into a json string
        :return: a json string
        """
        json_str = self.mdata
        return json.dumps(json_str, indent=2)


    def add_gps_correction_from_csv(self, csv_data):  # ,net, sta, start_dt, end_dt):
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

        # to json object, ignore default variables
        # gps_corr = pdf2.to_json(orient="records", date_format="epoch", double_precision=10,
        #                         force_ascii=True, date_unit="ms", default_handler=None, indent=2)
        gps_corr = pdf2.to_json(orient="records", indent=2)

        # add the json object as string to the self.mdata
        self.mdata.update({"GPS_CORRECTION": json.loads(gps_corr)})

        # print(pdf2.head())
        # return pdf2.to_csv(index=False)
        return pdf2


    def add_orientation_correction(self, input_json_file):
        """
        add json element for orientation_correction from input_json_file
        :param input_json_file: an input file of certain format, containing the required orientation_correction
        :return: json_sub_node, which was added into self.mdata
        """

        # update the "ORIENT_CORRECTION":
        #self.mdata.update({"ORIENT_CORRECTION":: json.loads(orient_corr)})
        return True


    def write_metadata2json(self, metadta_json_file):
        """
        write out the current metatdaya object into a file
        :param metadta_json_file:
        :return:
        """
        with open(metadta_json_file, "w") as f:
            json.dump(self.mdata, f, indent=2)

        return metadta_json_file

    def validate_json_schema(self):
        """
        validate the current json by https://github.com/Julian/jsonschema
        :return: boolean
        """

        return True


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

# ======================================================================================================================
# Example How to Run:
# python station_metadata_extra.py OA.CF28_clock_correction.csv /Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml ~/tmpdir
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
    (net, sta, csv_data) = get_csv_correction_data(in_csv_file)

    _inv = obspy.read_inventory(in_station_xml_file, format='STATIONXML')
    selected_inv = _inv.select(network=net, station=sta)

    # selected_inv may be 0,1, 2, multiple stations, each have a start_date end_date
    station_list = selected_inv.networks[0].stations
    if station_list is None or len(station_list) == 0:  # no further process for this dummy case
        print("(network station) %s.%s NOT present in the XML file %s" % (net, sta, in_station_xml_file))
        sys.exits(2)

    # redefine the selected_inv
    for a_station in station_list:  # loop over all Stations
        # get station star end date and split csv_data
        start_dt = a_station.start_date
        end_dt = a_station.end_date

        ajson = StationMetadataExtra(net, sta, start_datetime=start_dt, end_datetime=end_dt)

        # generate/format extra metadata from inputs
        mpdf = ajson.add_gps_correction_from_csv(csv_data)

        # updated the ajson object with more metadata, such as orientation corr
        ajson.add_orientation_correction("input_json_file")

        ajson.write_metadata2json(os.path.join(out_dir, "%s.%s_%s_extra_metadata.json" % (net, sta, str(start_dt))))

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
