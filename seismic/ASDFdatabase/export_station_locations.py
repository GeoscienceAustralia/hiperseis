#!/bin/env python
"""
Description:
    Export station locations in kml format.

References:

CreationDate:   16/06/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     16/06/23   RH
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
@click.option('--output-format', type=click.Choice(['csv', 'kml']),
              default='csv', show_default=True,
              help='Output format')
@click.option('--network', type=str, default=None, show_default=True,
              help="Network name")
@click.option('--station', type=str, default=None, show_default=True,
              help="Station name")
@click.option('--location', type=str, default=None, show_default=True,
              help="Location name")
@click.option('--channel', type=str, default=None, show_default=True,
              help="Channel name")
@click.option('--start-date', default='1900-01-01',
              help="Start date-time in UTC format",
              type=str, show_default=True)
@click.option('--end-date', default='2100-01-01',
              help="End date-time in UTC format.",
              type=str, show_default=True)
def process(asdf_source, output_filename, output_format, network, station, location, channel,
            start_date, end_date):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    OUTPUT_FILENAME: Output file-name\n

    Example usage:

    """

    def write_kml(rows, ofn):
        if('.kml' not in ofn): ofn += '.kml'
        fh = open(ofn, "w")
        fh.write("""<?xml version="1.0" encoding="UTF-8"?> 
                    <kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2"> \
                 \n""")
        fh.write("<Document>\n")
        for row in rows:
            fh.write("\t<Placemark>\n")
            fh.write("\t\t<name>" + '.'.join([row[0], row[1]]) + "</name>\n")
            fh.write("\t\t<Point>")
            fh.write("\t\t\t<coordinates>" + ','.join([str(row[2]), str(row[3])]) + "</coordinates>")
            fh.write("\t\t</Point>\n")
            fh.write("\t</Placemark>\n")
        # end for
        fh.write("</Document>\n")
        fh.write("</kml>\n")
        fh.close()
    # end func

    def write_csv(rows, ofn):
        if ('.csv' not in ofn): ofn += '.csv'
        np.savetxt(ofn, rows, header='net, sta, lon, lat',
                   fmt='%s, %s, %f, %f')
    # end func

    start_date_ts = None
    end_date_ts = None
    try:
        start_date_ts = UTCDateTime(start_date).timestamp if start_date else None
        end_date_ts   = UTCDateTime(end_date).timestamp if end_date else None
    except Exception as e:
        print(str(e))
        assert 0, 'Invalid input'
    # end try

    ds = FederatedASDFDataSet(asdf_source)

    query = 'select ns.net, ns.sta, ns.lon, ns.lat from netsta as ns, wdb as wdb ' \
            'where ns.net=wdb.net and ns.sta=wdb.sta '

    if (network):
        query += ' and ns.net="{}" '.format(network)
    # end if

    if (station):
        query += ' and ns.sta="{}" '.format(station)
    # end if

    if (location):
        query += ' and wdb.loc="{}" '.format(location)
    # end if

    if (channel):
        query += ' and wdb.cha="{}" '.format(channel)
    # end if

    if (start_date_ts and end_date_ts):
        query += ' and wdb.st>={} and wdb.et<={}'.format(start_date_ts, end_date_ts)
    # end if

    query += ' group by ns.net, ns.sta'

    # fetch rows
    rows = ds.fds.conn.execute(query).fetchall()
    array_dtype = [('net', 'U10'), ('sta', 'U10'),
                   ('lon', 'float'), ('lat', 'float')]
    rows = np.array(rows, dtype=array_dtype)

    if(output_format == 'kml'): write_kml(rows, output_filename)
    else: write_csv(rows, output_filename)
# end func

if (__name__ == '__main__'):
    process()
# end if
