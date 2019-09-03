#!/usr/bin/env python
"""Helper functions for producing synthetic pseudo-Receiver function traces
"""

import numpy as np
from scipy import signal

# pylint: disable=invalid-name

def generate_synth_rf(arrival_times, arrival_amplitudes, fs_hz=100.0, window_sec=(-10, 30), f_cutoff_hz=2.0):
    """Simple generator of synthetic R component receiver function with pulses at given arrival times.

    :param arrival_times: Iterable of arrival times as numerical values in seconds
    :type arrival_times: iterable of float
    :param arrival_amplitudes: Iterable of arrival amplitudes
    :type arrival_amplitudes: iterable of float
    :param fs_hz: Sampling rate (Hz) of output signal, defaults to 100.0
    :type fs_hz: float, optional
    :param window_sec: Time window over which to create signal (sec), defaults to (-10, 30)
    :type window_sec: tuple, optional
    :param f_cutoff_hz: Cutoff frequency (Hz) for low-pass filtering to generate realistic result, defaults to 2.0
    :type f_cutoff_hz: float, optional
    :return: Array of times and corresponding signal amplitudes
    :rtype: numpy.array, numpy.array
    """
    # Compute array of time values and indexes of arrivals
    duration = window_sec[1] - window_sec[0]
    N = int(fs_hz*duration)
    times = np.linspace(window_sec[0], window_sec[1], N)
    arrivals_index = np.round((np.array(arrival_times) - times[0])*fs_hz).astype(int)

    # Generate kernel of delta functions at specified arrival times
    kern = np.zeros_like(times)
    kern[arrivals_index] = np.array(arrival_amplitudes)  # pylint: disable=unsupported-assignment-operation

    # Filter to pass low frequencies
    waveform = signal.butter(4, f_cutoff_hz/fs_hz)
    signal_filt = signal.filtfilt(waveform[0], waveform[1], kern)

    # Normalize signal so max positive amplitude is 1.
    signal_filt = signal_filt/np.max(signal_filt)

    return times, signal_filt
# end func
