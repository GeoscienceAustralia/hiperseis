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


DEFAULT_OUTPUT_FILE = "IRIS-ALL.xml"

# Tuple of known illegal XML elements that obspy will reject if ingested. All such substrings
# are removed prior to ingestion into obspy.
KNOWN_ILLEGAL_ELEMENTS_PATTERNS = (
   r"<Azimuth>360\.[1-9]\d*</Azimuth>",
   r"<Azimuth>36[1-9](\.\d*)?</Azimuth>",
   r"<Azimuth>3[7-9]\d(\.\d*)?</Azimuth>",
                          )


def formRequestUrl(netmask="*", statmask="*", chanmask="*"):
   """Form request URL to download station inventory in stationxml format, down to channel level,
      with the given filters applied to network codes, station codes and channel codes.

      Hardwired to exclude restricted channels and exclude comments to reduce file size.
   """
   return "http://service.iris.edu/fdsnws/station/1/query?net=" + netmask + \
          "&sta=" + statmask + \
          "&cha=" + chanmask + \
          "&level=channel&format=xml&includerestricted=false&includecomments=false&nodata=404"


def cleanup(tmp_filename):
   try:
      os.remove(tmp_filename)
   except:
      print("WARNING: Failed to remove temporary file " + tmp_filename)


def setTextEncoding(resp):
   """For the given response object, set its encoding from the contents of the text returned from server.
   """
   encoding_pattern = r"^<\?xml .* encoding=\"(.+)\""
   matcher = re.compile(encoding_pattern)
   first_line = resp.text.split('\n', 1)[0]
   match = matcher.search(first_line)
   assert match
   encoding = match[1]
   resp.encoding = encoding
   assert resp.encoding is not None


def updateIrisStationXml(output_file, options=None):
   iris_url = formRequestUrl() if options is None else formRequestUrl(**options)
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
      # Convert using system call
      cmd = "fdsnxml2inv --quiet --formatted " + ifile.name + " " + output_file
      print("Converting to SC3 format...")
      if sys.version_info[0] < 3:
         subprocess.check_call(cmd.split())
      else:
         subprocess.check_call(cmd.split(), timeout=3600)
      print("Successfully updated file " + output_file)
   except:
      cleanup(ifile.name)
      raise

   cleanup(ifile.name)

   # Create human-readable text form of the IRIS station inventory (Pandas stringified table)
   output_txt = os.path.splitext(output_file)[0] + ".txt"
   regenerateHumanReadable(iris, output_txt)


def regenerateHumanReadable(iris, outfile):
   print("Generating human readable version...")
   from pdconvert import inventory2Dataframe
   import pandas as pd
   from obspy import read_inventory
   import functools
   if sys.version_info[0] < 3:
      import cStringIO as sio
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
      DUMPFILE = 'fdsn_stn_inv_dump.txt'
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
