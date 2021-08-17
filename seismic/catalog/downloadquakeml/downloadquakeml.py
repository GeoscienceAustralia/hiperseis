#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Updated on Tue Aug  3 10:34:00 2021
@author: sheecegardezi
"""
import argparse
import sys
import logging
import os
from eatws.pipeline import download_eatws_event_data


def validate_directory(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def main():

    parser = argparse.ArgumentParser(description='Download Seismologic Data')

    parser.add_argument("-a", "--agency", nargs='+', default="GA", required=True)
    parser.add_argument("-s", "--start", required=True)
    parser.add_argument("-e", "--end", required=True)
    parser.add_argument("-o", "--output", type=validate_directory, help="output dir for QuakeML files", required=True)
    parser.add_argument("-l", "--log", type=str.upper, default="INFO", dest="logLevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Set the logging level", required=False)

    # Debugging
    sys.argv = [
        'C:/Users/sheec/Desktop/Project/genquakeml/genquakeml/__main__.py',
        '-o', 'C:/Users/sheec/Desktop/Project/downloadquakeml/data',
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
            download_eatws_event_data(start_time=args.start, end_time=args.end, output_directory=args.output)
        else:
            raise Exception("Data acquisition not implemented for: "+agency)


if __name__ == '__main__':
    main()
