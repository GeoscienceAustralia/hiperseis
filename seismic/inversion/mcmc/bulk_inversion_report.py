#!/usr/bin/env python
# coding: utf-8
"""Produce PDF report of network stations showing RF inversion results
"""

import os
import re
import glob
import logging
import click

import tqdm.auto as tqdm

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import seismic.receiver_fn.inversion.plot_inversion as plot_inversion

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

# File to look for in case output folder by which we assume the case has completed running
RUN_COMPLETE_INDICATOR = 'Posterior.out'
# Expected extension appended to case id to identify case output folder
CASE_FOLDER_EXTENSION = '_OUT'


@click.command()
@click.argument('input-folder', type=click.Path(exists=True, file_okay=False), required=True)
@click.argument('output-file', type=click.Path(exists=False, file_okay=True), required=True)
@click.option('--file-mask', type=str, required=True,
              help='Regular expression filename mask used for identifying input .dat files'
                   'e.g "OA_B[S-V]*.dat" (use quotes to prevent shell from expanding wildcard)')
def main(input_folder, output_file, file_mask):
    '''
    :param output_file: Output file (must not exist already)
    :type output_file: str or Path (pdf extension expected)
    '''
    # Input .dat file are expected to be named according to following pattern
    # so that network and station ID can be extracted:
    # {net}_{sta}_{cha}_rf.dat

    files = sorted(glob.glob(file_mask))

    case_pattern = '^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)'
    matcher = re.compile(case_pattern)
    with PdfPages(output_file) as pdf:
        pbar = tqdm.tqdm(total=len(files))
        for f in files:
            pbar.update()
            case_name, ext = os.path.splitext(f)
            case_folder = case_name + CASE_FOLDER_EXTENSION
            if os.path.isdir(case_folder) and os.path.isfile(os.path.join(case_folder, RUN_COMPLETE_INDICATOR)):
                case_meta = matcher.match(case_name)
                net = case_meta.group(1)
                sta = case_meta.group(2)
                cha = case_meta.group(3)
                station = '.'.join([net, sta, cha])

                try:
                    plot_inversion.plot_bodin_inversion(pdf, case_folder, f, station=station)
                except OSError as e:
                    log.error('Error {} processing case {}'.format(str(e), case_name))
                # end try
            else:
                pbar.write('Skipping {}, output not found'.format(case_name))
            # end if
        # end for
        pbar.close()
    # end with

# end main


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()  # pylint: disable=no-value-for-parameter
# end if
