#!/usr/bin/env python
"""Utility stream processing functions.
"""

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
