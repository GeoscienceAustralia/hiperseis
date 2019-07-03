#!/usr/bin/env python
"""Utility functions to help with RF processing and analysis.
"""

import numpy as np
from scipy.signal import hilbert

KM_PER_DEG = 111.1949 # pylint: disable=invalid-name


def phase_weights(stream):
    """Phase weighting takes all the traces in a stream and computes a weighting for each sample in the
       stream between 0 and 1. The weighting represents how consistent is the phase angle of the signal
       at the same point in the time series across all streams.

       See https://doi.org/10.1111/j.1365-246X.1997.tb05664.x

    :param stream: Stream containing one or more traces from which phase coherence weightings will be generated.
    :type stream: obspy.core.trace.Trace
    :return: Array of normalized weighting factors with same length as traces in stream.
    :rtype: numpy.array
    """
    traces = np.array([tr.data for tr in stream])
    # Hilbert transform to separate complex amplitude from complex phase.
    analytic = hilbert(traces)
    # The complex part of the hilbert transform contains sine of the phase angle.
    # numpy.angle extracts the angle in radians from the imaginary component.
    angle = np.angle(analytic)
    # Using just the phase angle component (throwing away the amplitude) generate complex number
    # representing just the phase angle of the signal at each sample.
    i_phase = np.exp(1j * angle)
    # Compute the mean of all the complex phases. If they're random, due to summing in the complex
    # domain they will tend to cancel out. If they're not random, they will tend to sum coherently and
    # generated a large stacked complex amplitude.
    tphase = np.abs(np.mean(i_phase, axis=0))
    # Return normalized result against max amplitude, so that the most coherent part of the signal
    # has a scaling of 1.
    return tphase/np.max(tphase)


def find_rf_group_ids(stream):
    """For the given stream, which is expected to have an rf_group attribute in its traces' metadata, determine
       the unique set of group ids that the traces contain.

    :param stream: Stream containing traces with rf_group ids associated with them.
    :type stream: obspy.core.trace.Trace
    :return: Set of rf_group ids found in the traces
    :rtype: set(int)
    """
    # AttributeError may be raised here if rf_group attribute does not exist in the stats.
    group_ids = set((trace.stats.rf_group for trace in stream))
    return group_ids
