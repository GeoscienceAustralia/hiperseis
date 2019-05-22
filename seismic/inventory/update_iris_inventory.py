#!/usr/bin/env python
"""Automatically update IRIS-ALL.xml file from IRIS web portal.

   Outout file is saved as FDSN station xml.
   Script also generates human readable form as IRIS-ALL.txt.

   Example usages:
   ---------------

   `python update_iris_inventory.py`

   `python update_iris_inventory.py -o outfile.xml`

   `python update_iris_inventory.py --netmask=U* --statmask=K*`

   `python update_iris_inventory.py --netmask=UW,LO --output outfile.xml`
"""

import os
import sys
import argparse
import time
import re

import requests
from seismic.inventory.iris_query import form_channel_request_url, set_text_encoding

default_output_file = "IRIS-ALL.xml"

# Dictionary of known illegal XML elements that obspy will reject if ingested. All such substrings
# are replaced with their replacement prior to ingestion into obspy.
KNOWN_ILLEGAL_ELEMENTS_PATTERNS = {  # pylint: disable=invalid-name
    r"<Azimuth>360\.[1-9]\d*</Azimuth>": "",
    r"<Azimuth>36[1-9](\.\d*)?</Azimuth>": "",
    r"<Azimuth>3[7-9]\d(\.\d*)?</Azimuth>": "",
    r"<Azimuth>-90</Azimuth>": "",
    r"<Latitude>-90.878944</Latitude>": r"<Latitude>-90</Latitude>",
}


def cleanup(tmp_filename):
    """
    Helper function to clean up temporary file on disk.

    :param tmp_filename: File name to clean up
    :type tmp_filename: str
    """
    try:
        os.remove(tmp_filename)
    except OSError:
        print("WARNING: Failed to remove temporary file " + tmp_filename)


def update_iris_station_xml(req, output_file, options=None):
    """
    Pull the latest IRIS complete station inventory (down to station level, not including
    instrument responses) from IRIS web service and save to file in FDSN station xml format.

    :param req: Request object to use for URI query
    :type req: Object conforming to interface of 'requests' library
    :param output_file: Destination file to generate
    :type output_file: str
    :param options: Filtering options for network, station and channel codes, defaults to None
    :param options: Python dict of key-values pairs matching command line options, optional
    """
    iris_url = form_channel_request_url() if options is None else form_channel_request_url(**options)
    # Download latest IRIS station database as FDSN station xml.
    try:
        print("Requesting data from server...")
        iris = req.get(iris_url)
        set_text_encoding(iris)
    except req.exceptions.RequestException:
        print("FAILED to retrieve URL content at " + iris_url)
        return

    # Repair errors with IRIS data
    print("Correcting known data errors...")
    iris_fixed = repair_iris_metadata(iris)
    # Close the query to free resources
    iris.close()

    with open(output_file, 'w') as f:
        f.write(iris_fixed)

    # Create human-readable text form of the IRIS station inventory (Pandas stringified table)
    output_txt = os.path.splitext(output_file)[0] + ".txt"
    regenerate_human_readable(iris_fixed, output_txt)


def repair_iris_metadata(iris):
    """Perform text subtitutions to fix known errors in station xml returned from IRIS.

    :param iris: Response to IRIS query request containing response text
    :type iris: requests.models.Response
    :return: The text from the response with known faulty data substituted with fixed data.
    :rtype: str (Python 3) or unicode (Python 2)
    """

    def repair_match(match):
        for pattern, replacement in KNOWN_ILLEGAL_ELEMENTS_PATTERNS.items():
            if re.match(pattern, match.group()):
                return replacement
        return ""

    # This code repairs illegal data matching KNOWN_ILLEGAL_ELEMENTS_PATTERNS from iris.text
    matcher = re.compile("|".join(list(KNOWN_ILLEGAL_ELEMENTS_PATTERNS.keys())))
    iris_text_fixed = matcher.sub(repair_match, iris.text)

    return iris_text_fixed


def regenerate_human_readable(iris_data, outfile):
    """
    Generate human readable, tabular version of the IRIS database.

    :param iris_data: String containing result string returned from IRIS query (without data errors).
    :type iris_data: str
    :param outfile: Output text file name
    :type outfile: str
    """

    print("Generating human readable version...")
    from seismic.inventory.pdconvert import inventory2Dataframe
    import pandas as pd
    from obspy import read_inventory

    if sys.version_info[0] < 3:
        from cStringIO import StringIO as sio  # pylint: disable=import-error
    else:
        from io import BytesIO as sio

    iris_str = iris_data.encode('utf-8')
    print("  Ingesting query response into obspy...")
    obspy_input = sio(iris_str)
    try:
        station_inv = read_inventory(obspy_input)
    except:
        dumpfile = 'fdsn_stn_inv_dump.xml'
        print("FAILED ingesting server response into obspy, dumping server response string to " + dumpfile)
        with open(dumpfile, 'w') as f:
            f.write(iris_str.decode('utf-8'))
        raise

    print("  Converting to dataframe...")
    inv_df = inventory2Dataframe(station_inv)

    with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 1000):
        print("  Converting to tabular text file " + outfile)
        inv_str = str(inv_df)
        with open(outfile, "w") as f:
            f.write(inv_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--netmask", help="Filter mask to apply to network codes, "
                        "e.g. U* to get all network codes starting with \"U\"")
    parser.add_argument("-s", "--statmask", help="Filter mask to apply to station codes. "
                        "Filter strings should not include quotation marks.")
    parser.add_argument("-c", "--chanmask", help="Filter mask to apply to channel codes.")
    parser.add_argument("-o", "--output", help="Name of output file.", default=default_output_file)
    args = vars(parser.parse_args())
    filter_args = {k: v for k, v in args.items() if v is not None and k in ["netmask", "statmask", "chanmask"]}
    output_filename = args['output']
    print("Destination file: " + output_filename)
    time.sleep(1)

    if filter_args:
        update_iris_station_xml(requests, output_filename, filter_args)
    else:
        update_iris_station_xml(requests, output_filename)
