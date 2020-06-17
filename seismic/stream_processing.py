#!/usr/bin/env python
"""Utility stream processing functions.
"""

import copy
import functools
import numbers
import json

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

    :param baz: Value to check
    :type baz: int or float
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
    """
    Negate the data in the given channel of the stream

    :param _event_id: Ignored
    :param stream: Stream containing channel to flip
    :param channel: Single character string indicating which component to flip
    :return: Stream with channel data negated
    """
    trace_selected = stream.select(component=channel)[0]
    trace_selected.data = -trace_selected.data
    return stream
# end func


@functools.singledispatch
def scalarize(_obj, _stats):
    raise NotImplementedError('Cannot convert generic object to scalar')
# end func


@scalarize.register(numbers.Number)
def _(val, _stats):
    return float(val)
# end func


@scalarize.register(dict)
def _(d, stats):
    net = stats.network
    sta = stats.station
    key = '.'.join([net, sta])
    record = d.get(key)
    if record is None:
        return 0.0
    else:
        return record.get('azimuth_correction', 0.0)
    # end if
# end func


@scalarize.register(str)
def _(filename, stats):
    @functools.lru_cache()
    def load_correction_database(_fname):
        with open(_fname, 'r') as _f:
            _db = json.load(_f)
        # end with
        return _db
    # end func
    # Load function uses lru_cache, effectively equivalent to runonce, so this function
    # can be performant when called at high frequency.
    db = load_correction_database(filename)
    return scalarize(db, stats)
# end func


def correct_back_azimuth(_event_id, stream, baz_correction):
    """
    Apply modification to the back azimuth value in the stream stats

    :param _event_id: Ignored
    :param stream: Stream to which correction is applied
    :param baz_correction: Any object with a registered `scalarize` function for
        generating an angle correction for a trace in degrees. E.g. could be a
        numeric value, a dictionary of correction values, or a file produced by
        script `analyze_station_orientations.py`
    :return: Stream with modified back azimuth
    """
    for tr in stream:
        # Each station may have a custom correction. Expect that all such
        # possible corrections are represented in baz_correction argument.
        stats = tr.stats
        correction = scalarize(baz_correction, stats)
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


def assert_homogenous_stream(stream, funcname):
    """
    Verify that the given stream does not contain mixture of stations or channels/components.

    :param stream: Stream containing one or more traces
    :type stream: obspy.core.stream.Stream or rf.RFStream
    """
    # Check station and channel uniqueness. It is not sensible to expect RF similarity for
    # different stations or channels.
    if not stream:
        return
    # end if
    expected_station = stream[0].stats.station
    expected_channel = stream[0].stats.channel
    assert np.all(np.array([(tr.stats.station == expected_station) for tr in stream])), \
        'Mixed station data incompatible with function {}'.format(funcname)
    assert np.all(np.array([(tr.stats.channel == expected_channel) for tr in stream])), \
        'Mixed channel data incompatible with function {}'.format(funcname)
# end func
