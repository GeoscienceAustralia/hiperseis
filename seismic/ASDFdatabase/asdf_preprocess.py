#!/usr/env python
"""
Description:
    Reads waveforms from an ASDF file, optionally applies instrument response correction,
    resamples and outputs them to another ASDF file. This preprocessing is crucial for
    large-scale studies involving > 10000 Green's Functions, e.g. in ambient noise
    tomography. This approach significantly reduces IO bottlenecks and computational
    costs associated with having to apply instrument response corrections on data from
    a given station in alternative workflows.

References:

CreationDate:   18/07/19

Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     18/07/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import click
import os
from mpi4py import MPI
import pyasdf
from tqdm import tqdm
from obspy import UTCDateTime, read_inventory, Inventory, Stream
from collections import defaultdict
from obspy.core.util.misc import get_window_times
import gc
from obspy.core.util.misc import limit_numpy_fft_cache

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def getStationInventory(master_inventory, inventory_cache, netsta):
    netstaInv = None
    if (master_inventory):
        if (inventory_cache is None): inventory_cache = defaultdict(list)
        net, sta = netsta.split('.')

        if (isinstance(inventory_cache[netsta], Inventory)):
            netstaInv = inventory_cache[netsta]
        else:
            inv = master_inventory.select(network=net, station=sta)
            if(len(inv.networks)):
                inventory_cache[netsta] = inv
                netstaInv = inv
            # end if
        # end if
    # end if

    return netstaInv, inventory_cache
# end func

def create_station_asdf(input_asdf, output_folder, resample_rate,
                        instrument_response_inventory, instrument_response_output, water_level):
    # mpi attributes
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    # input asdf file
    ids = pyasdf.ASDFDataSet(input_asdf, mode='r')

    # get stations
    stations = ids.get_all_coordinates().keys()

    # local work-load
    stations = split_list(stations, nproc)[rank]

    # read inventory
    stationInvCache = None
    # read inventory
    inv = None
    try:
        inv = read_inventory(instrument_response_inventory)
    except Exception as e:
        print (e)
        raise RuntimeError('Failed to read inventory: %s' % (instrument_response_inventory))
    # end try

    for s in stations:
        # output asdf file
        ofn = os.path.join(output_folder,
                           os.path.splitext(os.path.basename(input_asdf))[0] + '.%s.h5'%(s))

        if (os.path.exists(ofn)): os.remove(ofn)
        ods = pyasdf.ASDFDataSet(ofn, mode='w', mpi=False, compression='gzip-3')

        sta = ids.waveforms[s]
        for tag in tqdm(sta.list(), desc='Rank %d, Station %s:'%(rank, s)):
            # get response object
            sinv, stationInvCache = getStationInventory(inv, stationInvCache, s)

            st = sta[tag]
            dayst = Stream()
            for tr in st:
                start_time = tr.stats.starttime
                offset = (UTCDateTime(year=start_time.year, month=start_time.month,
                                      day=start_time.day) - start_time)
                for wtr in tr.slide(3600*24, 3600*24, offset=offset, include_partial_windows=True):
                    wtr = wtr.copy()
                    dayst += wtr
                # end for
            # end for
            gc.collect() # force cleanup fft-related internal caches

            # remove response
            if(sinv):
                for tr in dayst:
                    limit_numpy_fft_cache(max_size_in_mb_per_cache=10)
                    try:
                        tr.remove_response(sinv, output=instrument_response_output.upper(),
                                           water_level=water_level)
                    except Exception as e:
                        print (e)
                    # end try
                    gc.collect()
                # end for
            # end if

            # detrend and taper
            taper_length = 20.0  # seconds
            for tr in dayst:
                if tr.stats.npts < 4 * taper_length * tr.stats.sampling_rate:
                    dayst.remove(tr)
                else:
                    tr.detrend(type="demean")
                    tr.detrend(type="linear")
                    tr.taper(max_percentage=None, max_length=1.0)
                # end if
            # end for
            gc.collect()

            # apply low-pass filter and create day traces
            for tr in dayst:
                tr.filter('lowpass', freq=resample_rate * 0.5, corners=6, zerophase=True)
                tr.interpolate(resample_rate, method='weighted_average_slopes')
            # end for
            gc.collect()

            # add traces
            for tr in dayst:
                try:
                    ods.add_waveforms(tr, tag='raw_recording')
                except Exception as e:
                    print (e)
                    print (tr)
                # end try
            # end for
            #break
        # end for
        gc.collect()

        ods.add_stationxml(ids.waveforms[s].StationXML)

        print ('Closing asdf file..')
        del ods

        #break
    # end for
    del ids
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-asdf', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--resample-rate', default=10, show_default=True,
              help="Resample rate in Hz; default is 10 Hz")
@click.option('--instrument-response-inventory', default=None,
              type=click.Path('r'),
              help="FDSNxml inventory containing instrument response information. Note that when this parameter is provided "
                   ", instrument response corrections are automatically applied for matching stations with response "
                   "information.")
@click.option('--instrument-response-output',
              type=click.Choice(['vel', 'disp']),
              default='vel', help="Output of instrument response correction; must be either 'vel' (default) for velocity"
                                  " or 'disp' for displacement. Note, this parameter has no effect if instrument response"
                                  " correction is not performed.")
@click.option('--water-level', default=50., help="Water-level in dB to limit amplification during instrument response correction"
                                                 "to a certain cut-off value. Note, this parameter has no effect if instrument"
                                                 "response correction is not performed.")
def process(input_asdf, output_folder, resample_rate, instrument_response_inventory,
            instrument_response_output, water_level):

    create_station_asdf(input_asdf, output_folder, resample_rate, instrument_response_inventory,
                        instrument_response_output, water_level)
# end func

if (__name__ == '__main__'):
    process()