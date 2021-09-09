from collections import defaultdict
import logging
import copy
import os
import numpy as np
from scipy import signal
from scipy.signal import hilbert, correlate

import obspy
import rf
import h5py

from seismic.stream_processing import assert_homogenous_stream
from seismic.receiver_fn.rf_network_dict import NetworkRFDict

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

class Corrections:
    def __init__(self, config_correction, hdf_keys):
        """
        Validate correction-config block
        """
        accepted_keys = ['rotate', 'negate', 'swap_ne', 'plot_dir']

        for k in config_correction.keys():
            if (k not in accepted_keys):
                assert 0, 'Unsupported key {} found in correction-block. ' \
                          'Supported keys are [{}]'.format(k, ','.join(accepted_keys))
            # end if
        # end for

        self._swap_ne_list = config_correction.setdefault('swap_ne', [])
        self._rotate_list = config_correction.setdefault('rotate', [])
        self._negate_list = config_correction.setdefault('negate', [])
        self._plot_dir = config_correction.setdefault('plot_dir', None)
        self._corrections = defaultdict(lambda: None)

        if (type(self._swap_ne_list) != list): assert 0, 'config:correction: Expected a list of net.sta.loc for key: swap_ne'
        if (type(self._rotate_list) != list): assert 0, 'config:correction: Expected a list of net.sta.loc for key: rotate'
        if (type(self._negate_list) != list): assert 0, 'config:correction: Expected a list of net.sta.loc.cha for key: negate'

        # NE channel swaps
        try:
            for item in self._swap_ne_list:
                net, sta, loc = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['swap_ne']): self._corrections['swap_ne'] = defaultdict(bool)
                        self._corrections['swap_ne'][hdf_item] = True
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:swap_ne block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc item was found in config:correction:rotate block'
        # end try

        # Rotations
        try:
            for item in self._rotate_list:
                net, sta, loc = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['rotate']): self._corrections['rotate'] = defaultdict(bool)
                        self._corrections['rotate'][hdf_item] = True
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:rotate block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc item was found in config:correction:rotate block'
        # end try

        # Negations
        try:
            for item in self._negate_list:
                net, sta, loc, cha = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['negate']): self._corrections['negate'] = defaultdict(list)
                        self._corrections['negate'][hdf_item].append = cha[-1]
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:negate block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc.cha item was found in config:correction:negate block'
        # end try
    # end func

    def needsChannelSwap(self, netstaloc):
        if(self._corrections['swap_ne']):
            return self._corrections['swap_ne'][netstaloc]
        else:
            return False
        # end if
    # end func

    def needsRotation(self, netstaloc):
        if(self._corrections['rotate']):
            return self._corrections['rotate'][netstaloc]
        else:
            return False
        # end if
    # end func

    def needsNegation(self, netstaloc):
        if(self._corrections['negate']):
            result = bool(len(self._corrections['negate'][netstaloc]))

            if(result): return self._corrections['negate'][netstaloc]
            else: return False
        else:
            return False
    # end func

    def needsCorrections(self, netstaloc):
        return self.needsChannelSwap(netstaloc) or \
               self.needsRotation(netstaloc) or \
               self.needsNegation(netstaloc)
    # end func

    @property
    def plot_dir(self):
        return self._plot_dir
    # end func
# end class

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
