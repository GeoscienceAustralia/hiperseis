#!/usr/bin/env python
"""Utility stream processing functions.
"""

import copy
import functools
import logging
import numbers
import json
import os
from collections import defaultdict
from obspy.signal.filter import lowpass
from obspy.core import Stream, Trace

import numpy as np

# pylint: disable=invalid-name

def zerophase_resample(item, resample_hz):
    def _zerophase_resample(trc, resample_hz):
        if(np.ma.is_masked(trc.data)): raise TypeError

        if(resample_hz < trc.stats.sampling_rate):
            trc.data = lowpass(trc.data, resample_hz / 2., trc.stats.sampling_rate, corners=2, zerophase=True)
        # end if

        trc.resample(resample_hz, no_filter=True)
    # end func

    if(isinstance(item, Stream)):
        for trc in item: _zerophase_resample(trc, resample_hz)
    elif(isinstance(item, Trace)): _zerophase_resample(item, resample_hz)
    else: raise TypeError
# end func

def zne_order(tr):
    """Channel ordering sort key function for ZNE ordering

    :param tr: Trace whose ordinal is to be determined.
    :type tr: obspy.Trace or rf.RFTrace
    :return: Numeric index indicating ZNE sort order of traces in a stream
    :rtype: int
    """
    trace_ordering = {'Z': 0, 'N': 1, 'E': 2}
    component = tr.stats.channel[-1].upper() if tr.stats.channel else False
    if component and component in trace_ordering:
        return trace_ordering[component]
    else:
        return 3
# end func


def zrt_order(tr):
    """Channel ordering sort key function for ZRT ordering

    :param tr: Trace whose ordinal is to be determined.
    :type tr: obspy.Trace or rf.RFTrace
    :return: Numeric index indicating ZRT sort order of traces in a stream
    :rtype: int
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

    :param t: 1D array of times
    :type t: numpy.array
    :param y: 1D array of sample values
    :type y: numpy.array
    :param t_new: 1D array of new times to interpolate onto
    :type t_new: numpy.array
    :return: 1D array of new interpolated sample values
    :rtype: numpy.array
    """
    dt = np.mean(np.diff(t))
    Ts, T = np.meshgrid(t_new, t)
    y_new = np.matmul(np.sinc((Ts.T - T.T) / dt), y)
    return y_new
# end func


def back_azimuth_filter(back_azi, back_azi_range):
    """Check if back azimuth `back_azi` is within range. Inputs must be in the range [0, 360] degrees.

    :param back_azi: Value to check
    :type back_azi: int or float
    :param back_azi_range: Pair of angles in degrees.
    :type back_azi_range: List or array of 2 floats, min and max back azimuth
    :return: True if back_azi is within back_azi_range, False otherwise.
    :rtype: bool
    """
    assert not np.any(np.isinf(np.array(back_azi_range)))
    assert not np.any(np.isnan(np.array(back_azi_range)))
    assert 0 <= back_azi <= 360
    assert 0 <= back_azi_range[0] <= 360
    assert 0 <= back_azi_range[1] <= 360
    back_azi_range = list(copy.copy(back_azi_range))
    while back_azi_range[0] > back_azi_range[1]:
        back_azi_range[0] -= 360
    return ((back_azi_range[0] <= back_azi <= back_azi_range[1]) or
            (back_azi_range[0] <= back_azi - 360 <= back_azi_range[1]) or
            (back_azi_range[0] <= back_azi + 360 <= back_azi_range[1]))
# end func


def swap_ne_channels(_event_id, stream):
    """
    Swap N and E channels on a stream. Changes the input stream.

    :param _event_id: Ignored
    :param stream: Stream whose N and E channels are to be swapped
    :type stream: obspy.Stream
    :return: Stream with channel swapping applied
    :rtype: obspy.Stream
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
    :type stream: obspy.Stream
    :param channel: Single character string indicating which component to flip
    :type channel: str
    :return: Stream with channel data negated
    :rtype: obspy.Stream
    """
    trace_selected = stream.select(component=channel)[0]
    trace_selected.data = -trace_selected.data
    return stream
# end func


@functools.singledispatch
def scalarize(_obj, _stats):
    """Fallback scalarize function for non-specialized type
    """
    raise NotImplementedError('Cannot convert generic object to scalar')
# end func


@scalarize.register(numbers.Number)
def _(val, _stats):
    """Scalarize function for numeric types

    :rtype: float
    """
    return float(val)
# end func


@scalarize.register(dict)
def _(d, stats):
    """Scalarize function for dict type, *assumed to be a container of obspy Stats*,
    to azimuth correction scalar.

    :rtype: float
    """
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

def recompute_inclinations(_event_id, stream):
    """
    :param _event_id: ignored
    :param stream: obspy.Stream or rf.RFSTream

    Rotates a copy of the input stream (ZNE) to ZRT and recomputes incidence angles
    by minimizing the energy on the L component through a grid search. The recomputed
    incidence angles are assigned to those in the input stream.
    @param stream:
    @return:
    """
    if(len(stream) == 0): return stream

    assert stream.traces[0].stats.phase == 'S', \
        'Expecting S waveforms, but found {}. Aborting..'.format(stream.traces[0].stats.phase)

    zrt = stream.copy().rotate('NE->RT')

    zrt_trimmed = zrt.copy().trim2(-10, 10)
    zr = np.array([zrt_trimmed[0].data,
                   zrt_trimmed[1].data])
    angles = np.linspace(0, np.pi / 2, 91)
    lq_ratios = []
    for a in angles:
        rot2D = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        lq = np.dot(rot2D, zr)
        power = np.sum(np.power(lq, 2.), axis=1)
        lq_ratios.append(power[0] / power[1])
    # end for
    lq_ratios = np.array(lq_ratios)
    inc = np.degrees(angles[np.argmin(lq_ratios)])

    for t in stream: t.stats.inclination = inc
    return
# end func


def correct_back_azimuth(_event_id, stream, baz_correction):
    """
    Apply modification to the back azimuth value in the stream stats

    :param _event_id: Ignored
    :param stream: Stream to which correction is applied
    :type stream: obspy.Stream or rf.RFSTream
    :param baz_correction: Any object with a registered `scalarize` function for
        generating an angle correction for a trace in degrees. E.g. could be a
        numeric value, a dictionary of correction values, or a file produced by
        script `rf_station_orientations.py`
    :return: Stream with modified back azimuth
    :rtype: Same as type(stream)
    """
    for tr in stream:
        # Each station may have a custom correction. Expect that all such
        # possible corrections are represented in baz_correction argument.
        stats = tr.stats
        correction = scalarize(baz_correction, stats)
        back_azi = stats.back_azimuth
        back_azi += correction
        while back_azi < 0:
            back_azi += 360
        while back_azi >= 360:
            back_azi -= 360

        stats.update({'orig_back_azimuth':stats.back_azimuth})
        stats.back_azimuth = back_azi
    # end for
    return stream
# end func

def assert_homogenous_stream(stream, funcname):
    """
    Verify that the given stream does not contain mixture of stations or channels/components.

    :param stream: Stream containing one or more traces
    :type stream: obspy.Stream or rf.RFStream
    :return: None
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
