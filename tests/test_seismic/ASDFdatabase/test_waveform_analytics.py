#!/bin/env python
"""
Description:
    Tests various aspects of the functionality in waveform_analytics.py

References:

CreationDate:   08/04/24
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     08/04/24   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import tempfile
import numpy as np
from obspy.core import Trace, Stream
from obspy.signal.spectral_estimation import PPSD
from obspy import UTCDateTime
from seismic.ASDFdatabase.analytics.station_analytics import StationAnalytics, PeriodDetailsCache
from seismic.inventory.response import ResponseFactory
from shutil import rmtree
from scipy.interpolate import interp1d
import sys

path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
path = os.path.join(path, 'seismic/inventory/')
resp_file = os.path.join(path, 'RESP.COMPACT120S.MINIMUS.txt')
tempdir = tempfile.mkdtemp()

def test_fast_psd():
    network = 'AA'
    station = 'BB'
    location = ''
    channel = 'BHZ'
    sampling_rate = 10
    scaling_factor = 1

    def get_time_range_func(net, sta, loc, cha):
        return UTCDateTime('2007-01-01'), UTCDateTime('2007-01-02')
    # end func

    def get_waveforms_func(net, sta, loc, cha, st, et):
        np.random.seed(10)
        t = Trace(data=np.random.random(86400 * sampling_rate) * scaling_factor,
              header={'network':net, 'station':sta, 'location':loc,
                      'channel':cha, 'sampling_rate':sampling_rate,
                      'starttime':UTCDateTime('2007-01-01')})
        return Stream([t])
    # end func

    # create a flat response
    rf = ResponseFactory()
    rf.createFromPAZ('flat_response', 'LAPLACE (RADIANS/SECOND)', normFactor=1,
                     normFreq=1,
                     stageGain=1,
                     stageGainFreq=1,
                     poles=[0 + 1j],
                     zeros=[0 + 1j])
    resp = rf.getResponse('flat_response')

    pdc = PeriodDetailsCache()
    pd = pdc.get_period_details(sampling_rate)
    # generate spectrum for fast_ppsd
    sa = StationAnalytics(get_time_range_func,
                          get_waveforms_func,
                          resp,
                          tempdir,
                          None,
                          start_time = None,
                          end_time = None,
                          nproc=1)

    sa.analyse_data(network, station, location, channel)
    # create an interpolation object for the spectrum from fast_ppsd
    fast_ppsd_io = None
    for start_time, end_time in zip(sa.st_list, sa.et_list):
        output_fn_stem = '{}.{}.{}.{}.{}.{}.{}'.format(network,
                                                    station,
                                                    location,
                                                    channel,
                                                    sampling_rate,
                                                    start_time,
                                                    end_time)

        output_fn_stem = output_fn_stem.replace(':', '__')
        output_fn_npz = os.path.join(sa.temp_folder, output_fn_stem + '.npz')

        results = np.load(output_fn_npz)
        spec = results['sparse_spec']

        fast_ppsd_io = interp1d(pd.sparse_periods, spec)
    # end for

    # generate spectrum using obspy ppsdf, with a flat response
    stream = get_waveforms_func(network, station, location, channel, None, None)
    ppsd = PPSD(stream[0].stats, {'sensitivity': 1.0,
                                  'gain': 1.0,
                                  'poles': [0 + 1j],
                                  'zeros': [0 + 1j]},
                db_bins=(-150, 50, 1.))
    ppsd.add(stream)

    periods, spec = ppsd.get_mean()
    ppsd_io = interp1d(periods, spec)

    # check conformity of spectra
    common_periods = np.linspace(np.max([pd.sparse_periods[0], periods[0]]),
                             np.min([pd.sparse_periods[-1], periods[-1]]), 100)

    corr = np.corrcoef(fast_ppsd_io(common_periods), ppsd_io(common_periods))[0,1]
    assert corr > 0.9
# end func
