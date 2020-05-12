#!/usr/bin/env python
# coding: utf-8
"""
Prototype script for applying the de-reverberation filter and sedimentary time
delay correction described in:

Yu, Y., J. Song, K. H. Liu, and S. S. Gao (2015), Determining crustal structure beneath seismic stations overlying a low-velocity sedimentary layer using receiver functions, J. Geophys. Res. Solid Earth, 120, 3208–3218, doi:10.1002/2014JB011610. http://dx.doi.org/10.1002/2014JB011610

Applies to R-channel only, so output file will have only R-channel.
"""

import logging

import click
import numpy as np
import scipy.signal as signal

import rf

import seismic.receiver_fn.rf_util as rf_util


def _reverb_removal_filter(offset, amplitude):
    reverb_removal_kernel = np.zeros(offset + 1)
    reverb_removal_kernel[0] = 1
    reverb_removal_kernel[-1] = amplitude
    return reverb_removal_kernel
# end func


def _estimate_Pbs_timedelta(tr):
    """Estimate time delta to PbS phase arrival relative to theoretical onset
    by finding first peak in time window about onset.
    """
    # TODO: Change to find first peak rather than highest peak
    lead_time = tr.stats.onset - tr.stats.starttime
    times_rel = tr.times() - lead_time
    mask = np.array((times_rel >= -0.5) & (times_rel <= 5.0))
    loc = np.nonzero(mask & (tr.data == np.max(tr.data[mask])))
    return tr.times()[loc[0][0]] - lead_time
# end func


def _dereverb_yu_autocorr(rf_stream):

    # Generate reverb removal filter

    x = np.array([tr.data for tr in rf_stream])

    autocorr = np.apply_along_axis(lambda _x: np.correlate(_x, _x, 'full')/np.dot(
        _x, _x), axis=1, arr=x)

    ac_causal = autocorr[:, autocorr.shape[1]//2:]

    ac_extrema = np.apply_along_axis(lambda _x: signal.argrelextrema(_x, np.less)[0][0],
                                 axis=1, arr=ac_causal)

    de_reverb = [_reverb_removal_filter(ac_extrema[i], -ac_causal[i, ac_extrema[i]])
                 for i in range(len(x))]

    # Apply reverb removal filter
    x_norev = [np.convolve(_x, de_reverb[i], 'full') for i, _x in enumerate(x)]
    x_norev = [_x[:(-ac_extrema[i] + 1)] for i, _x in enumerate(x_norev)]

    traces = []
    for i, tr in enumerate(rf_stream):
        new_tr = tr.copy()
        new_tr.data = x_norev[i]
        traces.append(new_tr)
    # end for
    rf_stream_noreverb = rf.RFStream(traces)

    # Compute time offsets allowed for in the H-k stacking code. These time offsets
    # are from Yu's paper.
    for i, tr in enumerate(rf_stream_noreverb):
        pbs_delta = _estimate_Pbs_timedelta(tr)
        tt_delta = ac_extrema[i]/tr.stats.sampling_rate
        tr.stats.sediment = {'t1_offset': pbs_delta,
                             't2_offset': tt_delta - pbs_delta,
                             't3_offset': tt_delta}
    # end for

    return rf_stream_noreverb
# end func


@click.command()
@click.option('--network', type=str)
@click.option('--station', type=str)
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False))
@click.argument('output-file', type=click.Path(exists=False, dir_okay=False))
def main(input_file, output_file, network, station):
    """
    Apply dereverberation to the R-component of receiver function.
    This is not a fully fledged sediment compensation method, and
    convention stacking applied after this may not reflect true depths
    to signal features.

    :param input_file: RF file to load in obspyh5 format with rf indexing
    :param output_file: RF file to write in obspyh5 format with rf indexing
    :param network: Network code of data to load from input file
    :param station: Station code of data to load from input file
    """
    # Data loading
    data_all = rf_util.read_h5_rf(input_file, network, station)

    # Select R-channel data
    db = rf_util.rf_to_dict(data_all)
    station_db = db[station]

    rf_type = 'ZRT'
    channel = rf_util.choose_rf_source_channel(rf_type, station_db)
    logging.warning("Selected channel: {}".format(channel))
    channel_data = station_db[channel]

    rf_stream = rf.RFStream(channel_data).sort(['back_azimuth'])
    rf_stream_noreverb = _dereverb_yu_autocorr(rf_stream)
    rf_stream_noreverb.write(output_file, format='h5')
# end func


if __name__ == '__main__':
    main()
# end if
