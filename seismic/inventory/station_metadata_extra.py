#! /usr/bin/env python
"""
Class Model for GA Extra Metadata, which will be included into station inventory XML file as EXTRA tagged elements

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

    def add_gps_correction_from_csv(self, csv_data=None):  # ,net, sta, start_dt, end_dt):
        """
        Select the csv rows according  net, sta, start_dt, end_dt

        Args:
            csv_data: a string of CSV
            start_dt: obspy.core.utcdatetime.UTCDateTime
            end_dt:  obspy.core.utcdatetime.UTCDateTime

        Returns: a subset of csv in pandas df, selected according to (net, sta, start_dt, end_dt)
        """
        if csv_data is None:
            return None  # No correction metadata to be added

        pdf = pd.read_csv(io.StringIO(csv_data))
        # pdf.insert(4,"utcdate", UTCDateTime(0))
        # print(pdf.head())
        pdf["utcdate"] = pdf.apply(lambda row: UTCDateTime(row.date), axis=1)

        _crit = (pdf['net'] == self.net) & (pdf['sta'] == self.sta) & (
                pdf['utcdate'] > self.start) & (pdf['utcdate'] < self.end)
        pdf2 = pdf.loc[_crit].copy()  # use copy() to fix "SettingWithCopyWarning"

        # drop columns inplace pdf2 itself will be changed, otherwise will return a new df
        pdf2.drop(['net', 'sta', 'utcdate'], axis=1, inplace=True)
        print(self.net, self.sta, " shapes = ", pdf.shape, pdf2.shape)

        # to json object, ignore default variables
        # gps_corr = pdf2.to_json(orient="records", date_format="epoch", double_precision=10,
        #                         force_ascii=True, date_unit="ms", default_handler=None, indent=2)
        gps_corr = pdf2.to_json(orient="records", indent=2)

        # add the json object as string to the self.mdata
        self.mdata.update({"GPS_CORRECTION": json.loads(gps_corr)})

        # print(pdf2.head())
        # return pdf2.to_csv(index=False)
        return pdf2

    def add_orientation_correction(self, jason_data_list):
        """
        add orientation_correction from a list of json corrections (dictionary)
        The self.mdata is updated if orientatoin correction found

        :param jason_data_list, a list of dictinaries like
        [
        {
            "network":"OA",
            "station":"BS24",
            "start_dt": "2017-11-07T09:07:34.930000Z",
            "end_dt":   "2018-08-23T03:52:29.528000Z",
            "azimuth_correction": -5.0
        },
        ...
        ]
        :return: Adict with key net.sta deleted. 
        """

        for orcorr in jason_data_list: # search the right net.sta
            if (orcorr.get("network") == self.net and orcorr.get("station") == self.sta):  # get wouldn't have KeyError, return None if no key found
                #print(orcorr, type(orcorr))

                orcorr2=orcorr.copy() # make a copy to avoid alter the original dict obj. or copy()
                orcorr2.pop("network")
                orcorr2.pop("station") 

                # update the "ORIENT_CORRECTION":
                self.mdata.update({"ORIENT_CORRECTION": orcorr2})
                return orcorr2
            else:
                pass # not matching network.station, continue the loop search for the net.sta

        return None

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


def get_orientation_corr(input_json_file):
    """
    read inpu json file to get the orientation correction metadata, rehash if necessary
    :param input_json_file: a json file by AndrewMedlin original program
    :return:
    """
    with open(input_json_file) as f:
        medata = json.load(f)

    # print(type(medata))

    # print(medata.keys())  # dict_keys(['7X.MA01', '7X.MA11', '7X.MA12', '...)]

    net_sta_list = []
    orient_corr_list = []
    for netsta in medata.keys():
        # netsta = '7X.MA01'
        network, station = netsta.split('.')

        net_sta_list.append((network, station))

        # print(medata[netsta])

        try:
            az_corr = medata[netsta]['azimuth_correction']
        except KeyError:
            print("KeyError Warning: No 'azimuth_correction' data for ", netsta)

        az_corr = medata[netsta].get('azimuth_correction')
        start = medata[netsta]['date_range'][0]
        end = medata[netsta]['date_range'][1]
        # ['2009-09-10T02:57:07.160000Z', '2010-06-05T05:35:00.623400Z'] <class 'list'>
        # print(data[netsta]['azimuth_correction'], type(data[netsta]['azimuth_correction'])) # 38.0 <class 'float'>

        # print(start, end, az_corr)
        # Make a dictionary:
        dicorr = {"network": network, "station": station, "start_dt": start, "end_dt": end,
                  "azimuth_correction": az_corr}
        orient_corr_list.append(dicorr)

    return (net_sta_list, orient_corr_list)


# ======================================================================================================================
# Example How to Run:
# (hiperseis) fzhang@zubuntu1804 ~/Githubz/hiperseis
# python seismic/inventory/station_metadata_extra.py /Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml
#   /Datasets/corrections/OA.CF28_clock_correction.csv /Datasets/Orientation_Correction_json/OA_ori_error_estimates.json ~/tmpdir/
# 
# python seismic/inventory/station_metadata_extra.py tests/testdata/7X_2009_2011_ASDF.xml /g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/GPSClockCorr/7X.MA43_clock_correction.csv  /g/data/ha3/Passive/SHARED_DATA/Inventory/Station_Extra_Metadata/OrientationCorr/7X_ori_error_estimates.json  ~/tmpdir/
# =====================================================================================================================
if __name__ == "__main__":

    USAGE = "python %s station_inventory_xml gps_clock_corr_csv ori_error_estimates_json_file [out_dir]" % sys.argv[0]

    GA_NameSpace = "https://github.com/GeoscienceAustralia/hiperseis"

    # Original station inventory XML file to be modified
    # /g/data/ha3/Passive/_AusArray/OA/ASDF_cleaned/OA_stations_2017-2018.xml
    in_station_xml_file = "/Datasets/StationXML_with_time_corrections2/OA.CF28_station_inv_modified.xml"

    # extra metadata info file(s) to be read and formatted into JSON
    in_csv_file = "./OA.CF28_clock_correction.csv"

    # output dir for modified station inventory xml files
    out_dir = "~/tmpdir"

    if len(sys.argv) < 4:
        print(USAGE)
        sys.exit(1)
    else:
        in_station_xml_file = sys.argv[1]
        in_csv_file = sys.argv[2]
        in_json_file = sys.argv[3]

    if len(sys.argv) >= 5:
        out_dir = sys.argv[4]
    else:
        out_dir = None

    network_station_pairs=[]

    # get the metadata and it's associated network.station
    if os.path.exists(in_csv_file):
        (net, sta, csv_data) = get_csv_correction_data(in_csv_file)
        network_station_pairs.append((net, sta))
    else:
        csv_data = None 

    if os.path.exists(in_json_file):
        (net_sta, oricorr_json_data) = get_orientation_corr(in_json_file)
        network_station_pairs = network_station_pairs + net_sta
    else:
        oricorr_json_data = []  # None is not iterable

    print(network_station_pairs, len(network_station_pairs))

    # read in the initial station XML
    inv_obj = obspy.read_inventory(in_station_xml_file, format='STATIONXML')

    for (net, sta) in network_station_pairs:
        selected_inv = inv_obj.select(network=net, station=sta)

        # selected_inv may be 0,1, 2, multiple stations, each have a start_date end_date
        if selected_inv is None or len(selected_inv.networks) == 0:
            pass  # nothing to do, no network 
        else:
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
                ajson.add_orientation_correction(oricorr_json_data)

                ajson.write_metadata2json(
                    os.path.join(out_dir, "%s.%s_%s_extra_metadata.json" % (net, sta, str(start_dt))))

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
