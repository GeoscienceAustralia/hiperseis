#!/usr/bin/env python
"""Automatically update IRIS-ALL.xml file from IRIS web portal.

   IRIS-ALL.xml file is saved as Seiscomp3 station xml.
   
   Script also generates human readable form as IRIS-ALL.txt.

   Requires seiscomp3 to be available in the system path.
"""

import os
import sys

import requests as req
import tempfile as tmp
import subprocess
import argparse

OUTPUT_FILE = "IRIS-ALL.xml"
OUTPUT_TXT = os.path.splitext(OUTPUT_FILE)[0] + ".txt"


def formRequestUrl(netmask="*", statmask="*", locmask="*", chanmask="*"):
   """Form request URL to download station inventory in stationxml format, down to channel level,
      with the given filters applied to network codes, station codes, locations and channel codes.
   """
   return "http://service.iris.edu/fdsnws/station/1/query?net=" + netmask + \
          "&sta=" + statmask + \
          "&loc=" + locmask + \
          "&cha=" + chanmask + \
          "&level=channel&format=xml&includerestricted=false&includecomments=false&nodata=404"


def cleanup(tmp_filename):
   try:
      os.remove(tmp_filename)
   except:
      print("WARNING: Failed to remove temporary file " + tmp_filename)

def updateIrisStationXml(options=None):
   iris_url = formRequestUrl() if options is None else formRequestUrl(**options)
   # Download latest IRIS station database as FDSN station xml.
   try:
      print("Requesting data from server...")
      iris = req.get(iris_url)
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
      cmd = "fdsnxml2inv --quiet --formatted " + ifile.name + " " + OUTPUT_FILE
      print("Converting to SC3 format...")
      subprocess.check_call(cmd.split(), timeout=3600)
      print("Successfully updated file " + OUTPUT_FILE)
   except:
      cleanup(ifile.name)
      raise

   cleanup(ifile.name)

   # Create human-readable text form of the IRIS station inventory (Pandas stringified table)
   print("Generating human readable version...")
   regenerateHumanReadable(OUTPUT_FILE, OUTPUT_TXT)


def regenerateHumanReadable(infile, outfile):
   from pdconvert import inventory2Dataframe
   import pandas as pd
   inv_df = inventory2Dataframe(infile)
   with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 1000):
      inv_str = str(inv_df)
      with open(outfile, "w") as f:
         f.write(inv_str)


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("-n", "--netmask", help="Filter mask to apply to network codes, e.g. U* to get all network codes starting with \"U\"")
   parser.add_argument("-s", "--statmask", help="Filter mask to apply to station codes. Filter strings should not include quotation marks.")
   parser.add_argument("-l", "--locmask", help="Filter mask to apply to station locations.")
   parser.add_argument("-c", "--chanmask", help="Filter mask to apply to channel codes.")
   args = vars(parser.parse_args())
   args = {k: v for k, v in args.items() if v is not None}

   if args:
      updateIrisStationXml(args)
   else:
      updateIrisStationXml()
