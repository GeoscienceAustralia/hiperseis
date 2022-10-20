#!/bin/env python
"""
Description:
    Find data-gaps in ASDF archives, by network.

References:

CreationDate:   10/01/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     10/01/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
from obspy import UTCDateTime
import click
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('output-filename', required=True,
                type=click.Path(dir_okay=False))
@click.option('--network', type=str, default=None, show_default=True,
                help="Network name")
@click.option('--station', type=str, default=None, show_default=True,
              help="Station name")
@click.option('--location', type=str, default=None, show_default=True,
              help="Location name")
@click.option('--channel', type=str, default=None, show_default=True,
              help="Channel name")
@click.option('--start-date', default=None,
              help="Start date-time in UTC format. If specified, 'end-date' must also be specified; " 
                   "otherwise this parameter is ignored.",
              type=str, show_default=True)
@click.option('--end-date', default=None,
              help="End date-time in UTC format. If specified, 'start-date' must also be specified; " 
                   "otherwise this parameter is ignored.",
              type=str, show_default=True)
@click.option('--min-gap-length', default=86400, type=float, show_default=True,
              help="Minimum length of gaps in seconds to report")
def process(asdf_source, output_filename, network, station, location, channel, start_date, end_date, min_gap_length):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    OUTPUT_FILENAME: Output file-name\n

    Example usage:

    """

    try:
        start_date = UTCDateTime(start_date).timestamp if start_date else None
        end_date   = UTCDateTime(end_date).timestamp if end_date else None
    except Exception as e:
        print(str(e))
        assert 0, 'Invalid input'
    # end try

    assert min_gap_length > 0, '--min-gap-length must be > 0'

    ds = FederatedASDFDataSet(asdf_source)
    gaps = ds.find_gaps(network, station, location, channel, start_date, end_date, min_gap_length)

    # Dump results to file
    fh = open(output_filename, 'w+')
    for i in np.arange(len(gaps)):
        row = gaps[i]
        fh.write('{} {} {} {} {} {}\n'.format(row['net'], row['sta'],
                                              row['loc'] if len(row['loc']) else '--',
                                              row['cha'], UTCDateTime(row['st']),
                                              UTCDateTime(row['et'])))
    # end for
    fh.close()
# end func

if (__name__ == '__main__'):
    process()
# end if
