#!/usr/bin/env python
"""Utility stream processing functions.
"""

import copy

import numpy as np


def zne_order(tr):
    """Channel ordering sort key function

    :param tr: Trace whose ordinal is to be determined.
    :type tr: obspy.Trace or RFTrace
    :return: Numeric index indicating ZNE sort order of traces in a stream
    """
    trace_ordering = {'Z': 0, 'N': 1, 'E': 2}
    component = tr.stats.channel[-1].upper() if tr.stats.channel else False
    if component and component in trace_ordering:
        return trace_ordering[component]
    else:
        return 3
# end func


def zrt_order(tr):
    """Channel ordering sort key function

    :param tr: Trace whose ordinal is to be determined.
    :type tr: obspy.Trace or RFTrace
    :return: Numeric index indicating ZRT sort order of traces in a stream
    """
    trace_ordering = {'Z': 0, 'R': 1, 'T': 2}
    component = tr.stats.channel[-1].upper() if tr.stats.channel else False
    if component and component in trace_ordering:
        return trace_ordering[component]
    else:
        return 3
# end func


def sinc_resampling(t, y, t_new):
    """Resample signal y for known times t onto new times t_new.
    Sampling rates do not need to match and time windows do not need
    to overlap. t_new should not have a lower sampling rate than t.

    :param t: numpy.array of times
    :param y: numpy.array of sample values
    :param t_new: numpy.array of new times to interpolate onto
    :return: numpy.array of new interpolated sample values
    """
    dt = np.mean(np.diff(t))
    Ts, T = np.meshgrid(t_new, t)
    y_new = np.matmul(np.sinc((Ts.T - T.T) / dt), y)
    return y_new
# end func


def back_azimuth_filter(baz, baz_range):
    """Check if back azimuth `baz` is within range. Inputs must be in the range [0, 360] degrees.

    :param baz_range: Pair of angles in degrees.
    :type baz_range: List or array of 2 floats, min and max back azimuth
    :return: True if baz is within baz_range, False otherwise.
    """
    assert not np.any(np.isinf(np.array(baz_range)))
    assert not np.any(np.isnan(np.array(baz_range)))
    assert 0 <= baz <= 360
    assert 0 <= baz_range[0] <= 360
    assert 0 <= baz_range[1] <= 360
    baz_range = list(copy.copy(baz_range))
    while baz_range[0] > baz_range[1]:
        baz_range[0] -= 360
    return ((baz_range[0] <= baz <= baz_range[1]) or
            (baz_range[0] <= baz - 360 <= baz_range[1]) or
            (baz_range[0] <= baz + 360 <= baz_range[1]))
# end func


def swap_ne_channels(_event_id, stream):
    """
    Swap N and E channels on a stream. Changes the input stream.

    :param _event_id: Ignored
    :param stream: Stream whose N and E channels are to be swapped
    :return: Stream with channel swapping applied
    """
    stream_n = stream.select(component='N')
    stream_e = stream.select(component='E')
    data_n = copy.deepcopy(stream_n[0].data)
    stream_n[0].data = stream_e[0].data
    stream_e[0].data = data_n
    return stream
# end func


def negate_channel(_event_id, stream, channel):
    trace_selected = stream.select(component=channel)[0]
    trace_selected.data = -trace_selected.data
    return stream
# end func


def correct_back_azimuth(_event_id, stream, baz_correction):
    for tr in stream:
        # Each station may have a custom correction. Expect that all such
        # possible corrections are represented in baz_correction argument.
        stats = tr.stats
        sta = stats.station
        if isinstance(baz_correction, (float, int)):
            correction = baz_correction
        else:
            assert isinstance(baz_correction, dict)
            correction = baz_correction[sta]
        # end if
        baz = stats.back_azimuth
        baz += correction
        while baz < 0:
            baz += 360
        while baz >= 360:
            baz -= 360
        stats.back_azimuth = baz
    # end for
    return stream
# end func
