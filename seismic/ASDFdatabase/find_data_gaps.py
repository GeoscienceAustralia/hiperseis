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
from collections import defaultdict
from tqdm import tqdm


def find_gaps(asdf_source, network=None, station=None, location=None,
              channel=None, start_date_ts=None, end_date_ts=None,
              min_gap_length=86400):
    ds = FederatedASDFDataSet(asdf_source)
    conn = ds.fds.conn

    clause_added = 0
    query = 'select net, sta, loc, cha, st, et from wdb '
    if(network or station or location or channel or (start_date_ts and end_date_ts)): query += " where "

    if(network):
        query += ' net="{}" '.format(network)
        clause_added += 1
    # end if

    if(station):
        if(clause_added): query += ' and sta="{}" '.format(station)
        else: query += ' sta="{}" '.format(station)
        clause_added += 1
    # end if

    if(location):
        if(clause_added): query += ' and loc="{}" '.format(location)
        else: query += ' loc="{}" '.format(location)
        clause_added += 1
    # end if

    if(channel):
        if(clause_added): query += ' and cha="{}" '.format(channel)
        else: query += ' cha="{}" '.format(channel)
        clause_added += 1
    # end if

    if(start_date_ts and end_date_ts):
        if(clause_added): query += ' and st>={} and et<={}'.format(start_date_ts, end_date_ts)
        else: query += ' st>={} and et<={}'.format(start_date_ts, end_date_ts)
    # end if
    query += ' order by st, et'

    rows = conn.execute(query).fetchall()

    array_dtype = [('net', 'U10'), ('sta', 'U10'),
                   ('loc', 'U10'), ('cha', 'U10'),
                   ('st', 'float'), ('et', 'float')]
    rows = np.array(rows, dtype=array_dtype)

    # Process rows
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

    result = []
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

                                result.append((net, sta, loc, cha, et[idx], st[idx + 1]))
                            # end for
                        # end if
                    # end if
                # end for
                # break
            # end for
        # end for
    # end for
    result = np.array(result, dtype=array_dtype)

    return result
# end func

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

    gaps = find_gaps(asdf_source, network, station, location, channel, start_date, end_date, min_gap_length)

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