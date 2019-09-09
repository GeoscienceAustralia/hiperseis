#!/usr/bin/env python
"""Helper functions for producing synthetic pseudo-Receiver function traces
"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d

import obspy
import rf

import seismic.receiver_fn.rf_util as rf_util

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


def synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, ds, log=None, include_t3=False):
    """Synthesize RF R-component data set over range of inclinations and distances
    and get result as a rf.RFStream instance.

    :param H: Moho depth (km)
    :type H: float
    :param V_p: P body wave velocity in uppermost layer
    :type V_p: float
    :param V_s: S body wave velocity in uppermost layer
    :type V_s: float
    :param inclinations: Array of inclinations for which to create RFs
    :type inclinations: numpy.array(float)
    :param distances: Array of teleseismic distances corresponding to inclinations
    :type distances: numpy.array(float)
    :param ds: Final sampling rate (Hz) for the downsampled output signal
    :type ds: float
    :param log: Logger to send output to, defaults to None
    :type log: logger, optional
    :param include_t3: If True, include the third expected multiple PpSs+PsPs
    :type include_t3: bool, optional
    :return: Stream containing synthetic RFs
    :rtype: rf.RFStream
    """
    assert len(inclinations) == len(distances), "Must provide 1:1 inclination and distance pairs"

    k = V_p/V_s
    traces = []
    for i, inc_deg in enumerate(inclinations):
        theta_p = np.deg2rad(inc_deg)
        p = np.sin(theta_p)/V_p

        t1 = H*(np.sqrt((k*k/V_p/V_p) - p*p) - np.sqrt(1.0/V_p/V_p - p*p))
        t2 = H*(np.sqrt((k*k/V_p/V_p) - p*p) + np.sqrt(1.0/V_p/V_p - p*p))
        arrivals = [t1, t2]
        if include_t3:
            t3 = t1 + t2
            arrivals.append(t3)
        if log is not None:
            log.info("Inclination {:3g} arrival times: {}".format(inc_deg, arrivals))

        arrivals = [0] + arrivals
        amplitudes = [1, 0.5, 0.4]
        if include_t3:
            amplitudes.append(-0.3)
        window = (-5.0, 50.0)  # sec
        fs = 100.0  # Hz
        _, synth_signal = generate_synth_rf(arrivals, amplitudes, fs_hz=fs, window_sec=window)

        now = obspy.UTCDateTime.now()
        dt = float(window[1] - window[0])
        end = now + dt
        onset = now - window[0]
        header = {'network': 'SY', 'station': 'TST', 'location': 'GA', 'channel': 'HHR', 'sampling_rate': fs,
                  'starttime': now, 'endtime': end, 'onset': onset,
                  'station_latitude': -19.0, 'station_longitude': 137.0,  # arbitrary (approx location of OA deployment)
                  'slowness': p*rf_util.KM_PER_DEG, 'inclination': inc_deg,
                  'back_azimuth': 0, 'distance': float(distances[i])}
        tr = rf.rfstream.RFTrace(data=synth_signal, header=header)
        tr = tr.decimate(int(np.round(fs/ds)), no_filter=True)
        traces.append(tr)
    # end for

    stream = rf.RFStream(traces)

    return stream
# end func


def convert_inclination_to_distance(inclinations, model="iasp91", nominal_source_depth_km=10.0):
    """Helper function to convert range of inclinations to teleseismic distance in degrees.

    :param inclinations: Array of inclination angles in degrees
    :type inclinations: numpy.array(float)
    :param model: Name of model to use for ray tracing, defaults to "iasp91"
    :type model: str, optional
    :param nominal_source_depth_km: Assumed depth of source events, defaults to 10.0
    :type nominal_source_depth_km: float, optional
    :return: Array of teleseismic distances in degrees corresponding to input inclination angles.
    :rtype: numpy.array(float)
    """
    # Generate function mapping ray parameter to teleseismic distance.
    # The distances are not strictly required for H-k stacking, but rf behaves better when they are there.
    ts_distance = np.linspace(25, 95, 71)
    inc = np.zeros_like(ts_distance)
    model = obspy.taup.TauPyModel(model=model)
    source_depth_km = nominal_source_depth_km
    for i, d in enumerate(ts_distance):
        ray = model.get_ray_paths(source_depth_km, d, phase_list=['P'])
        inc[i] = ray[0].incident_angle  # pylint: disable=unsupported-assignment-operation
    # end for

    # Create interpolator to convert inclination to teleseismic distance
    interp_dist = interp1d(inc, ts_distance, bounds_error=False,
                           fill_value=(np.max(ts_distance), np.min(ts_distance)))

    # Loop over valid range of inclinations and generate synthetic RFs
    distances = interp_dist(inclinations)

    return distances
# end func
