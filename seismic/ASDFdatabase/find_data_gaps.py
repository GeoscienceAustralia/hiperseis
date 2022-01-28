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

import os, sys

from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from tqdm import tqdm
import click

def dump_gaps(asdf_source, network, start_date, end_date, min_gap_length, output_filename):
    ds = FederatedASDFDataSet(asdf_source)
    conn = ds.fds.conn

    query = 'select st, et, net, sta, loc, cha from wdb where net="{}"'.format(network)
    if(start_date and end_date):
        query += ' and st>={} and et<={}'.format(start_date, end_date)
    # end if
    query += ' order by st, et'

    rows = conn.execute(query).fetchall()
    rows = np.array(rows, dtype=[('st', 'float'), ('et', 'float'),
                                ('net', 'object'),
                                ('sta', 'object'),
                                ('loc', 'object'),
                                ('cha', 'object')])


    tree = lambda: defaultdict(tree)
    nested_dict = tree()
    for i in tqdm(np.arange(rows.shape[0])):
        net = rows['net'][i]
        sta = rows['sta'][i]
        loc = rows['loc'][i]
        cha = rows['cha'][i]
        st = rows['st'][i]
        et = rows['et'][i]

        if (type(nested_dict[net][sta][loc][cha]) == defaultdict):
            nested_dict[net][sta][loc][cha] = []
        # end if

        nested_dict[net][sta][loc][cha].append([st, et])
    # end for

    fh = open(output_filename, 'w+')
    for net in nested_dict.keys():
        for sta in nested_dict[net].keys():
            for loc in nested_dict[net][sta].keys():
                for cha in nested_dict[net][sta][loc].keys():
                    arr = nested_dict[net][sta][loc][cha]
                    if (len(arr)):
                        arr = np.array(arr)

                        st = arr[:, 0]
                        et = arr[:, 1]
                        assert np.allclose(np.array(sorted(st)), st), 'Start-times array not sorted!'
                        gaps = np.argwhere((st[1:] - et[:-1]) >= min_gap_length)

                        if (len(gaps)):
                            for i, idx in enumerate(gaps):
                                idx = idx[0]

                                fh.write('{} {} {} {} {} {}\n'.format(net, sta,
                                                                      loc if len(loc) else '--',
                                                                      cha, UTCDateTime(et[idx]),
                                                                      UTCDateTime(st[idx + 1])))
                            # end for
                        # end if
                    # end if
                # end for
                # break
            # end for
        # end for
    # end for
    fh.close()
# end func


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('network-name', required=True,
                type=str)
@click.argument('output-filename', required=True,
                type=click.Path(dir_okay=False))
@click.option('--start-date', default=None,
              help="Start date-time in UTC format. If specified, 'end-date' must also be specified; " 
                   "otherwise this parameter is ignored.",
              type=str)
@click.option('--end-date', default=None,
              help="End date-time in UTC format. If specified, 'start-date' must also be specified; " 
                   "otherwise this parameter is ignored.",
              type=str)
@click.option('--min-length', default=86400, type=float,
              help="Minimum length of gaps in seconds to report")
def process(asdf_source, network_name, output_filename, start_date, end_date, min_length):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    NETWORK_NAME: Network name \n
    OUTPUT_FILENAME: Output file-name\n

    Example usage:

    """

    try:
        start_date = UTCDateTime(start_date).timestamp if start_date else None
        end_date   = UTCDateTime(end_date).timestamp if end_date else None
        length     = int(min_length)
    except Exception as e:
        print(str(e))
        assert 0, 'Invalid input'
    # end try

    assert min_length > 0, '--min-length must be > 0'

    dump_gaps(asdf_source, network_name, start_date, end_date, min_length, output_filename)
# end func

if (__name__ == '__main__'):
    process()
# end if