#! /usr/bin/env python
"""
Read input csv file or other formats, generate json formatted output file

Sample input metadata:./OA.CF28_clock_correction.csv

Output Json (Schema) string look like following.
Sort by date and group the metadata by date ranges related to the given station.

meta_dict = {
    "network_code":"OA",
    "station_code":"CF28",

  "gps_clock_corrections": [
    {
      "date": "2018-01-04",
      "seconds": -1.3375814425470127
    },
    {
      "date": "2018-01-05",
      "seconds": -1.110449564656099
    },
    {
      "date": "2018-01-06",
      "seconds": -0.9032476255118933
    }
  ],
      "orient_correction": {
        "start_dt": "2017-11-07T09:07:34.930000Z",
        "end_dt":   "2018-08-23T03:52:29.528000Z",
        "azimuth_correction": -5.0
    }
}
References:
    https://gajira.atlassian.net/browse/PV-324
    https://gajira.atlassian.net/browse/PV-312

CreationDate:   06/08/20
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     6/08/20   FZ
    LastUpdate:

"""
# https://www.google.com/search?client=firefox-b-e&q=python+covert+csv+to+json
import os

import pandas as pd

def main(input_csv, output_json ):
    """
    define my main function
    :return:
    """
    print("starting main()")

    pdf = pd.DataFrame(pd.read_csv(input_csv, sep=",", header=0, index_col=False))
    pdf.to_json(output_json, orient="records", date_format="epoch", double_precision=10,
                     force_ascii=True, date_unit="ms", default_handler=None, indent=4)
    # Section to define functions or class

    return


# =============================================
# Section for quick test of this script
# ---------------------------------------------
if __name__ == "__main__":
    # call main function
    main("./OA.CF28_clock_correction.csv", "output.json")