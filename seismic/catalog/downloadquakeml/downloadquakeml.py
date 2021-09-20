#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Updated on Tue Aug  3 10:34:00 2021
@author: sheecegardezi
"""
import argparse
import sys
import logging
import os
from pathlib import Path
from eatws.pipeline import download_eatws_event_data


def validate_directory(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def main():
    parser = argparse.ArgumentParser(description='Download Seismologic Data')

    parser.add_argument("-a", "--agency", nargs='+', default="GA", required=True)
    parser.add_argument("-s", "--start_time", required=True)
    parser.add_argument("-e", "--end_time", required=True)
    parser.add_argument("-n", "--min_magnitude", default=0)
    parser.add_argument("-x", "--max_magnitude", default=9.94)
    parser.add_argument("-o", "--output", type=validate_directory, help="output dir for QuakeML files", required=True)
    parser.add_argument("-l", "--log", type=str.upper, default="INFO", dest="logLevel",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Set the logging level",
                        required=False)

    # Debugging
    sys.argv = [
        'C:/Users/sheec/Desktop/Project/genquakeml/genquakeml/__main__.py',
        '-o', '/home/ssr/Desktop/Work/hiperseis/seismic/catalog/downloadquakeml/data',
        '-s', '2021-05-24T09:27:00',
        '-e', '2021-06-24T09:27:00',
        '-a', "GA"
    ]

    args = parser.parse_args()

    logging.info("Running GenQuakeML")
    if args.logLevel:
        logging.basicConfig(level=getattr(logging, args.logLevel))

    for agency in args.agency:
        if agency == "GA":
            download_eatws_event_data(start_time=args.start_time, end_time=args.end_time,
                                      min_magnitude=args.min_magnitude, max_magnitude=args.max_magnitude,
                                      output_directory=args.output)
        else:
            raise Exception("Data acquisition not implemented for: " + agency)


if __name__ == '__main__':
    main()
