#!/usr/bin/env python
"""Automatically update IRIS-ALL.xml file from IRIS web portal.

   IRIS-ALL.xml file is saved as Seiscomp3 station xml.
   Script also generates human readable form as IRIS-ALL.txt.
   Requires seiscomp3 to be available in the system path.

   Example usages:
   ---------------

   `python update_iris_inventory.py`

   `python update_iris_inventory.py -o outfile.xml`

   `python update_iris_inventory.py --netmask=U* --statmask=K*`

   `python update_iris_inventory.py --netmask=UW,LO --output outfile.xml`
"""

import os
import sys

import requests as req
import tempfile as tmp
import subprocess
import argparse
import time
import re
from seismic.inventory.iris_query import formChannelRequestUrl, setTextEncoding
from seismic.inventory.fdsnxml_convert import toSc3ml


DEFAULT_OUTPUT_FILE = "IRIS-ALL.xml"
OUTPUT_FORMATS = ["SC3", "FDSN"]

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
    except:
        print("WARNING: Failed to remove temporary file " + tmp_filename)


def update_iris_station_xml(output_file, output_format, sc3_compatible,  options=None):
    """
    Pull the latest IRIS complete station inventory (down to station level, not including
    instrument responses) from IRIS web service and save to file in chosen format.

    :param output_file: Destination file to generate
    :type output_file: str
    :param output_format: Destination file output format
    :type output_format: str
    :param options: Filtering options for network, station and channel codes, defaults to None
    :param options: Python dict of key-values pairs matching command line options, optional
    """
    assert output_format in OUTPUT_FORMATS
    iris_url = formChannelRequestUrl() if options is None else formChannelRequestUrl(**options)
    # Download latest IRIS station database as FDSN station xml.
    try:
        print("Requesting data from server...")
        iris = req.get(iris_url)
        setTextEncoding(iris)
    except:
        print("FAILED to retrieve URL content at " + iris_url)
        return

    # Repair errors with IRIS data
    sc3_format = (output_format == 'SC3')
    print("Correcting known data errors...")
    iris_fixed = repair_iris_metadata(iris, sc3_compatible=(sc3_format or sc3_compatible))
    iris.close()

    if sc3_format:
        try:
            # Since there are not available Python bindings for the conversion, we have to dump FDSN stxml to file
            # first, then convert to seiscomp3 stxml inventory using system call.
            ifile = tmp.NamedTemporaryFile("w", delete=False)
            ifile.write(iris_fixed)
            ifile.close()
            # Convert using helper module
            print("Converting to SC3 format...")
            toSc3ml(ifile.name, output_file)
            print("Successfully updated file " + output_file)
        except:
            cleanup(ifile.name)
            raise
        cleanup(ifile.name)
    elif output_format == 'FDSN':
        with open(output_file, 'w') as f:
            f.write(iris_fixed)
    else:
        assert False, "Unknown output format {}".format(output_format)

    # Create human-readable text form of the IRIS station inventory (Pandas stringified table)
    output_txt = os.path.splitext(output_file)[0] + ".txt"
    regenerate_human_readable(iris_fixed, output_txt)


def repair_iris_metadata(iris, sc3_compatible=False):
    """Perform text subtitutions to fix known errors in station xml returned from IRIS.
       If sc3_compatible:
       * Add a non-empty nominal instrument response where missing.
       * Replace empty station and channel dates with dates based on parent.

    :param iris: [description]
    :type iris: [type]
    :param sc3_compatible: [description]
    :type sc3_compatible: [type]
    :return: [description]
    :rtype: [type]
    """

    def repair_match(match):
        for pattern, replacement in KNOWN_ILLEGAL_ELEMENTS_PATTERNS.items():
            if re.match(pattern, match.group()):
                return replacement
        return ""

    # This code repairs illegal data matching KNOWN_ILLEGAL_ELEMENTS_PATTERNS from iris.text
    matcher = re.compile("|".join(list(KNOWN_ILLEGAL_ELEMENTS_PATTERNS.keys())))
    iris_text_fixed = matcher.sub(repair_match, iris.text)

    if sc3_compatible:
        # TODO: functionalize the response and date cleanup code from fdsnxml_convert.py and apply here.
        pass

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
    parser.add_argument("-n", "--netmask", help="Filter mask to apply to network codes, "
                        "e.g. U* to get all network codes starting with \"U\"")
    parser.add_argument("-s", "--statmask", help="Filter mask to apply to station codes. "
                        "Filter strings should not include quotation marks.")
    parser.add_argument("-c", "--chanmask", help="Filter mask to apply to channel codes.")
    parser.add_argument("-o", "--output", help="Name of output file.", default=DEFAULT_OUTPUT_FILE)
    parser.add_argument("-f", "--format", default="SC3", help="Output format, from {}".format(OUTPUT_FORMATS))
    parser.add_argument("--sc3-compatible", action="store_true",
                        help="Make sure output is patched for SC3 compatibility (regardless of output format")
    args = vars(parser.parse_args())
    filter_args = {k: v for k, v in args.items() if v is not None and k in ["netmask", "statmask", "chanmask"]}
    output_filename = args['output']
    output_format = args['format']
    sc3_compat = args['sc3_compatible']
    if output_format not in OUTPUT_FORMATS:
        print("Allowed output formats: {}".format(OUTPUT_FORMATS))
        exit(0)
    print("Destination file: " + output_filename + " (" + output_format + " format)")
    time.sleep(1)

    if filter_args:
        update_iris_station_xml(output_filename, output_format, sc3_compat, filter_args)
    else:
        update_iris_station_xml(output_filename, output_format, sc3_compat)
