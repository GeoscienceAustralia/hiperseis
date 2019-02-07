#!/usr/bin/env python
"""Automatically update IRIS-ALL.xml file from IRIS web portal.

   IRIS-ALL.xml file is saved as Seiscomp3 station xml.
   Script also generates human readable form as IRIS-ALL.txt.
   Requires seiscomp3 to be available in the system path.

   Example usages:
   ---------------
   python update_iris_inventory.py
   python update_iris_inventory.py -o outfile.xml
   python update_iris_inventory.py --netmask=U* --statmask=K*
   python update_iris_inventory.py --netmask=UW,LO --output outfile.xml
"""

import os
import sys

import requests as req
import tempfile as tmp
import subprocess
import argparse
import time
import re
from .iris_query import formChannelRequestUrl, setTextEncoding
from .fdsnxml_convert import toSc3ml


DEFAULT_OUTPUT_FILE = "IRIS-ALL.xml"

# Tuple of known illegal XML elements that obspy will reject if ingested. All such substrings
# are removed prior to ingestion into obspy.
KNOWN_ILLEGAL_ELEMENTS_PATTERNS = (
    r"<Azimuth>360\.[1-9]\d*</Azimuth>",
    r"<Azimuth>36[1-9](\.\d*)?</Azimuth>",
    r"<Azimuth>3[7-9]\d(\.\d*)?</Azimuth>",
)


def cleanup(tmp_filename):
    """
    Helper function to clean up temporary file on disk.

    :param tmp_filename: File name to clean up
    :type tmp_filename: str
    """
    try:
        os.remove(tmp_filename)
    except:
        print("WARNING: Failed to remove temporary file " + tmp_filename)


def updateIrisStationXml(output_file, options=None):
    """
    Pull the latest IRIS complete station inventory (down to station level, not including
    instrument responses) from IRIS web service and save to file in sc3ml format.

    :param output_file: Destination sc3ml file to generate
    :type output_file: str
    :param options: Filtering options for network, station and channel codes, defaults to None
    :param options: Python dict of key-values pairs matching command line options, optional
    """
    iris_url = formChannelRequestUrl() if options is None else formChannelRequestUrl(**options)
    # Download latest IRIS station database as FDSN station xml.
    try:
        print("Requesting data from server...")
        iris = req.get(iris_url)
        setTextEncoding(iris)
    except:
        print("FAILED to retrieve URL content at " + iris_url)
        return

    try:
        # Since there are not available Python bindings for the conversion, we have to dump FDSN stxml to file
        # first, then convert to seiscomp3 stxml inventory using system call.
        ifile = tmp.NamedTemporaryFile("w", delete=False)
        ifile.write(iris.text)
        ifile.close()
        # Convert using helper module
        print("Converting to SC3 format...")
        toSc3ml(ifile.name, output_file)
        print("Successfully updated file " + output_file)
    except:
        cleanup(ifile.name)
        raise

    cleanup(ifile.name)

    # Create human-readable text form of the IRIS station inventory (Pandas stringified table)
    output_txt = os.path.splitext(output_file)[0] + ".txt"
    regenerateHumanReadable(iris, output_txt)


def regenerateHumanReadable(iris, outfile):
    """
    Generate human readable, tabular version of the IRIS database.

    :param iris: Query response object containing result string returned from IRIS query.
    :type iris: requests.Response
    :param outfile: Output text file name
    :type outfile: str
    """

    print("Generating human readable version...")
    from pdconvert import inventory2Dataframe
    import pandas as pd
    from obspy import read_inventory

    if sys.version_info[0] < 3:
        import cStringIO as sio  # pylint: disable=import-error
    else:
        import io as sio

    print("  Ingesting query response into obspy...")
    # This code removes all instances of the substrings matching KNOWN_ILLEGAL_ELEMENTS_PATTERNS from iris.text,
    # then encodes it as UTF-8.
    matcher = re.compile("|".join(KNOWN_ILLEGAL_ELEMENTS_PATTERNS))
    iris_str = matcher.sub("", iris.text).encode('utf-8')
    obspy_input = sio.BytesIO(iris_str)
    try:
        station_inv = read_inventory(obspy_input)
    except:
        DUMPFILE = 'fdsn_stn_inv_dump.xml'
        print("FAILED ingesting server response into obspy, dumping server response string to " + DUMPFILE)
        with open(DUMPFILE, 'w') as f:
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
    parser.add_argument("-n", "--netmask", help="Filter mask to apply to network codes, e.g. U* to get all network codes starting with \"U\"")
    parser.add_argument("-s", "--statmask", help="Filter mask to apply to station codes. Filter strings should not include quotation marks.")
    parser.add_argument("-c", "--chanmask", help="Filter mask to apply to channel codes.")
    parser.add_argument("-o", "--output", help="Name of output file.", default=DEFAULT_OUTPUT_FILE)
    args = vars(parser.parse_args())
    filter_args = {k: v for k, v in args.items() if v is not None and k != "output"}
    output_filename = args['output']
    print("Destination file: " + output_filename)
    time.sleep(1)

    if filter_args:
        updateIrisStationXml(output_filename, filter_args)
    else:
        updateIrisStationXml(output_filename)
