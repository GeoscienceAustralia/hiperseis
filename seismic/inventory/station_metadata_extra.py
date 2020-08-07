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


import json
import pandas as pd

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

        self.mdata = { "network": self.net,
                       "station": self.sta
                     } # initially dict

    def validate_json_schema(self):
        """
        validate the current json by https://github.com/Julian/jsonschema
        :return: boolean
        """

        return True

    def make_json_string(self):
        """
        format this metadata object (dict) into a json string
        :return:
        """
        json_str = self.mdata
        return json.dumps(json_str, indent=4)

    def add_gps_corrections(self, csvfile):
        """
        read from the csvfile, make a gps correction list for this net sta
        :param csvfile:
        :return: json_sub_node, mwhich was added into self.mdata
        """
        pdf = pd.DataFrame(pd.read_csv(csvfile, sep=",", header=0, index_col=False))
        mygps_corr = pdf.to_json(orient="records", date_format="epoch", double_precision=10,
                    force_ascii=True, date_unit="ms", default_handler=None, indent=4)

        self.mdata.update({"GPS_CORRECTION": json.loads(mygps_corr)})

        return mygps_corr

    def add_orientation_correction(self, infile):
        """
        add json element for orientation_correction from input data file infile
        :param infile: an input file of certain format, containing the required orientation_correction
        :return: json_sub_node, which was added into self.mdata
        """

# ===============================================================================
if __name__ == "__main__":
    myInst= StationMetadataExtra("OA","CF28")
    print(myInst.make_json_string(), type(myInst.make_json_string()))

    gpsc = myInst.add_gps_corrections("./OA.CF28_clock_correction.csv")

    # console print(json.dumps(myInst.mdata,indent=2))

    with open("./extra_mdata.dev.json", "w") as f:
        json.dump(myInst.mdata, f, indent=2)



