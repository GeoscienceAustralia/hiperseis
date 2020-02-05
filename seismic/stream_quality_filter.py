#!/usr/bin/env python
"""
Helper functions for curating and quality controlling stream objects.
"""

import numpy as np


def curate_stream3c(ev_id, stream3c, logger=None):
    """
    Apply quality curation criteria to a stream. Modifies the stream in-place if required. Traces in stream
    must be in ZNE order.

    The following checks are made. If any of these checks fails, the function returns False:
    * Inclination value is not NaN
    * The stream has 3 channels
    * Each trace in the stream has the same number of samples
    * None of the traces have any NaN values in the time series data
    * None of the traces have zero variance

    The following cleanups are attempted on the stream:
    * All 3 traces have the same time range

    :param ev_id: The event id
    :type ev_id: int or str
    :param stream3c: Stream with 3 components of trace data
    :type stream3c: obspy.Stream
    :param logger: Logger in which to log messages
    :type logger: logger.Logger object
    :return: True if checks pass, False otherwise
    """

    # Apply essential sanity checks before trying to compute RFs.
    for tr in stream3c:
        if np.isnan(tr.stats.inclination):
            if logger:
                logger.warning("Invalid inclination found in stream {} (skipping):\n{}".format(ev_id, stream3c))
            return False
    # end for

    if len(stream3c) != 3:
        if logger:
            logger.warning("Unexpected number of channels in stream {} (skipping):\n{}".format(ev_id, stream3c))
        return False
    # end if

    # Strongly assert expected ordering of traces. This must be respected so that
    # RF normalization works properly.
    assert stream3c.traces[0].stats.channel[-1] == 'Z', stream3c.traces[0].stats.channel
    assert stream3c.traces[1].stats.channel[-1] == 'N', stream3c.traces[1].stats.channel
    assert stream3c.traces[2].stats.channel[-1] == 'E', stream3c.traces[2].stats.channel

    # If traces have inconsistent time ranges, clip to time range during which they overlap. Not guaranteed
    # to make time ranges consistent due to possible inconsistent sampling times across channels.
    start_times = np.array([tr.stats.starttime for tr in stream3c])
    end_times = np.array([tr.stats.endtime for tr in stream3c])
    if not (np.all(start_times == start_times[0]) and np.all(end_times == end_times[0])):
        clip_start_time = np.max(start_times)
        clip_end_time = np.min(end_times)
        stream3c.trim(clip_start_time, clip_end_time)
    # end if

    if len(stream3c[0]) != len(stream3c[1]) or len(stream3c[0]) != len(stream3c[2]):
        if logger:
            logger.warning("Channels in stream {} have different lengths, cannot generate RF (skipping):\n{}"
                           .format(ev_id, stream3c))
        return False
    # end if

    for tr in stream3c:
        # Each tr here is one component.
        # Check for any NaNs in any component. Discard such traces, as we don't want any NaNs
        # propagating through workflow.
        if np.any(np.isnan(tr.data)):
            if logger:
                logger.warning("NaN detected in trace {} of stream {} (skipping):\n{}"
                               .format(tr.stats.channel, ev_id, stream3c))
            return False
        # end if
        # Check for all zeros or all same value in any component. This is infeasible for an operational
        # station, and has been observed as a failure mode in practice.
        if np.std(tr.data) == 0:
            if logger:
                logger.warning("Invariant data detected in trace {} of stream {} (skipping):\n{}"
                               .format(tr.stats.channel, ev_id, stream3c))
            return False
        # end if
    # end for

    return True
# end func
