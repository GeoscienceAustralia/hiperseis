from collections import defaultdict
import logging
import numpy as np
from scipy.signal import correlate
import rf
# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

def has_reverberations(cha_data, dt_max=0.2):
    """
    Checks if the data shows signs of the presence of reverberations

    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: Iterable(rf.RFTrace)
    :param dt_max: if the median temporal offset between an RF peak and the onset time > dt_max,
                   this function returns true
    :type dt_max: float

    return: Bool
    """

    dt_array = []
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime

        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt_array.append(relative_time[loc])
    # end for
    dt_array = np.array(dt_array)

    if(np.median(dt_array) > dt_max): return True

    return False
# end func

def apply_reverberation_filter(cha_data):
    result_stream = []
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime

        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt = relative_time[loc]

        data = correlate(tr.data, tr.data, mode='full')
        data /= np.max(data)
        data = data[len(data) // 2:]

        r0 = -(np.min(data))
        Dt = np.argmin(data) * 1. / tr.stats.sampling_rate

        tr_copy = tr.copy()
        resonance_filter = np.zeros(len(tr.data))
        resonance_filter[0] = 1
        resonance_filter[int(Dt * tr.stats.sampling_rate)] = r0
        tr_copy.data = np.convolve(tr_copy.data, resonance_filter, mode='full')
        tr_copy.data = tr_copy.data[:len(tr_copy.data) // 2 + 1]

        assert tr.data.shape == tr_copy.data.shape, 'Input/output length mismatch detected in ' \
                                                    'reverberation removal routine'

        tr_copy.stats.update({'t1_offset':dt,
                              't2_offset':Dt - dt,
                              't3_offset':Dt})

        result_stream.append(tr_copy)
    # end for

    return rf.RFStream(result_stream)
# end func
