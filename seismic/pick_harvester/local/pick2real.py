#!/bin/env python
"""
Description:
    Exports the outputs of pick_eqt.py in a form that REAL can ingest

    REAL: Rapid Earthquake Association and Location,
          Seismol. Res. Lett., 90.6, 2276-2284, 2019, https://doi.org/10.1785/0220190052
References:

CreationDate:   20/07/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     20/07/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os, shutil

from collections import defaultdict
import numpy as np
from obspy import UTCDateTime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import click
from tqdm import tqdm
import pandas as pd
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
def export_picks(fds: FederatedASDFDataSet, startDate: UTCDateTime, endDate: UTCDateTime,
                 picks: defaultdict(list), outputFolder:str, applyClockCorrection: bool=False):

    stationCache = defaultdict(list)
    step = 3600 * 24 # day
    # outputCache is keyed by day -> net.sta.cha -> phase
    outputCache = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for phase, arrivals in picks.items():
        current = startDate.timestamp
        logger.info('Processing {} phases..'.format(phase))
        pbar = tqdm(total = int((endDate.timestamp - current) / step),
                    desc='Processing daily data:')
        while(current + step <= endDate.timestamp):
            pbar.update()
            imask = (arrivals['timestamp'] >= current) & \
                    (arrivals['timestamp'] <= (current + step))
            dayArrivals = arrivals[imask]
            dayKey = UTCDateTime(current).strftime('%Y%m%d')
            for i in np.arange(len(dayArrivals)):
                lon, lat, elev = None, None, None
                net, sta, loc, cha = dayArrivals.iloc[i, 0:4]
                nsc = '{}.{}.{}'.format(net, sta, cha)
                if(nsc in stationCache.keys()):
                    lon, lat, elev = stationCache[nsc]
                else:
                    result = fds.get_stations(current,
                                              current + step,
                                              net, sta, loc, cha)
                    if(len(result)):
                        _, _, _, _, lon, lat, elev = result[0][:]
                        stationCache[nsc] = [lon, lat, elev/1e3]
                    # end if
                # end if

                if(lon): # station exists
                    clockCorr = 0
                    if(applyClockCorrection):
                        result = fds.fds._get_correction(net, sta, loc,
                                                         UTCDateTime(current),
                                                         UTCDateTime(current+step))
                        if(result): clockCorr = result[1]
                    # end if

                    row = dayArrivals.iloc[i]
                    # timestamp is referenced to the current day
                    # probability is from 0-100
                    # amplitude is set to 0 as we do not seek to compute earthquake magnitude
                    outputCache[dayKey][nsc][phase].append([row['timestamp'] - current - clockCorr,
                                                            row['probability'] * 100,
                                                            0])
                # end if
            # end for

            #print(len(arrivals), len(dayArrivals))
            current += step
        # wend
    # end for

    # output data to expected folder-structure
    for dk, dv in outputCache.items():
        for nsck, nscv in outputCache[dk].items():
            net, sta, cha = nsck.split('.')
            for pk, pv in outputCache[dk][nsck].items():
                path = os.path.join(outputFolder, dk)
                if(not os.path.exists(path)): os.mkdir(path)
                fn = os.path.join(path, '{}.{}.{}.txt'.format(net, sta, pk))
                with open(fn, 'w+') as fh:
                    for row in pv:
                        fh.write('{} {} {}\n'.format(row[0], row[1], row[2]))
                    # end for
                # end with
            # end for
        # end for
    # end for

    # output stations.db
    fn = os.path.join(outputFolder, 'station.dat')
    with open(fn, 'w+') as fh:
        for nsck, val in stationCache.items():
            net, sta, cha = nsck.split('.')
            lon, lat, elev = val
            fh.write('{} {} {} {} {} {}\n'.format(lon, lat, net, sta, cha, elev))
        # end for
    # end with
    #print(outputCache)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--p-arrivals', type=click.Path(exists=True),
              default=None, show_default=True,
              help="P-arrivals text file, as output by pick_eqt.py")
@click.option('--s-arrivals', type=click.Path(exists=True),
              default=None, show_default=True,
              help="S-arrivals text file, as output by pick_eqt.py")
@click.option('--start-date', type=str, default='1900-01-01', show_default=True,
              help="Date to generate outputs from")
@click.option('--end-date', type=str, default='2100-01-01', show_default=True,
              help="Date to generate outputs until")
@click.option('--apply-clock-correction', is_flag=True, default=False,
              help='Apply GPS time corrections')
def process(asdf_source, output_folder, p_arrivals, s_arrivals,
            start_date, end_date, apply_clock_correction):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    OUTPUT_FOLDER: Output folder \n
    """

    if(apply_clock_correction):
        val = os.environ.get('GPS_CLOCK_CORRECTION')
        if(val != '1'):
            assert 0, 'Environment variable "GPS_CLOCK_CORRECTION" not set to 1. Aborting..'
        # end if
    # end if

    fds = None
    if(start_date and len(start_date) != 10): assert 0, \
        'invalid start-date; must be specified as YYYY-MM-DD. Aborting..'
    if(end_date and len(end_date) != 10): assert 0, \
        'invalid end-date; must be specified as YYYY-MM-DD. Aborting..'

    try:
        start_date = UTCDateTime(start_date)
        end_date = UTCDateTime(end_date)
    except Exception as e:
        print(str(e))
        assert 0, 'Invalid start- or end-date. Aborting..'
    # end try

    try:
        fds = FederatedASDFDataSet(asdf_source)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to instantiate FederatedASDFDataSet with input file: {}'.format(asdf_source)
    # end try

    # load P and S arrivals
    picks = defaultdict(list)

    dataStartTimes = []
    dataEndTimes = []
    if(p_arrivals):
        logger.info('Loading P arrivals from {}..'.format(p_arrivals))
        df = pd.read_csv(p_arrivals)
        df = df.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
        df.columns = df.columns.str.replace('#', '')
        df.columns = df.columns.str.replace(' ', '')
        picks['P'] = df
        dataStartTimes.append(UTCDateTime(np.min(df['timestamp'])))
        dataEndTimes.append(UTCDateTime(np.max(df['timestamp'])))
    # end if

    if(s_arrivals):
        logger.info('Loading S arrivals from {}..'.format(s_arrivals))
        df = pd.read_csv(s_arrivals)
        df = df.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)
        df.columns = df.columns.str.replace('#', '')
        df.columns = df.columns.str.replace(' ', '')
        picks['S'] = df
        dataStartTimes.append(UTCDateTime(np.min(df['timestamp'])))
        dataEndTimes.append(UTCDateTime(np.max(df['timestamp'])))
    # end if

    # check for consistency
    for k, df in picks.items():
        assert np.all(df.columns == ['net', 'sta', 'loc', 'cha', 'timestamp', 'probability']), \
               'Unexpected columns found in {} arrivals file. '.format(k)
    # end for

    if(start_date < np.min(dataStartTimes)):
        start_date = np.min(dataStartTimes)
        start_date = UTCDateTime(start_date.year, start_date.month, start_date.day)
    # end if
    if(end_date > np.max(dataEndTimes)):
        end_date = np.max(dataEndTimes)
        end_date = UTCDateTime(end_date.year, end_date.month, end_date.day)
        end_date += 24*3600
    # end if

    export_picks(fds, start_date, end_date, picks, output_folder, apply_clock_correction)
# end func

if (__name__ == '__main__'):
    process()
# end if
