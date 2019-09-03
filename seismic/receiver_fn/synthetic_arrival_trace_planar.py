#!/usr/bin/env python
"""Use a basic planar, 2-layer model of only the crust and the Moho to generate
synthetic arrival traces for known model characteristics. Intended to be used
for model validation.
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt

import obspy
import rf

import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking

# pylint: disable=invalid-name

def generate_synth_rf(arrival_times, arrival_amplitudes, fs_hz=100.0, window_sec=(-10, 30)):
    duration = window_sec[1] - window_sec[0]
    N = int(fs_hz*duration)
    times = np.linspace(window_sec[0], window_sec[1], N)
    arrivals_index = np.round((np.array(arrival_times) - times[0])*fs_hz).astype(int)
    kern = np.zeros_like(times)
    kern[arrivals_index] = np.array(arrival_amplitudes)  # pylint: disable=unsupported-assignment-operation
    # kern = np.ones(int(0.5*fs_hz))
    waveform = scipy.signal.butter(4, 2.0/fs_hz)
    signal_filt = scipy.signal.filtfilt(waveform[0], waveform[1], kern)
    signal_filt = signal_filt/np.max(signal_filt)
    return times, signal_filt


def main():
    H = 40  # km

    V_p = 6.4  # km/s, top (crust) layer
    V_s = 4.0  # km/s, top (crust) layer

    k = V_p/V_s
    print("H = {:3g}".format(H))
    print("k = {:3g}".format(k))

    inc_deg = 15.0
    theta_p = np.deg2rad(inc_deg)
    p = np.sin(theta_p)/V_p

    t1 = H*(np.sqrt((k*k/V_p/V_p) - p*p) - np.sqrt(1.0/V_p/V_p - p*p))
    t2 = H*(np.sqrt((k*k/V_p/V_p) - p*p) + np.sqrt(1.0/V_p/V_p - p*p))
    print("Arrival times: {}".format([t1, t2]))

    arrivals = [0, t1, t2]
    amplitudes = [1, 0.5, 0.4]
    window = (-10.0, 50.0)  # sec
    fs = 100.0
    _, signal = generate_synth_rf(arrivals, amplitudes, fs_hz=fs, window_sec=window)

    now = obspy.UTCDateTime.now()
    dt = float(window[1] - window[0])
    end = now + dt
    onset = now - window[0]
    header = {'network': 'OA', 'station': 'TST', 'location': '0M', 'channel': 'HHR', 'sampling_rate': fs,
              'starttime': now, 'endtime': end, 'onset': onset,
              'station_latitude': -19.0, 'station_longitude': 137.0,
              'slowness': p*rf_util.KM_PER_DEG, 'inclination': inc_deg,
              'back_azimuth': 0, 'distance': 45}
    tr = rf.rfstream.RFTrace(data=signal, header=header)
    tr = tr.decimate(10, no_filter=True)
    tr_list = [tr, tr, tr]
    # stream = rf.RFStream(tr_list)
    # _ = rf_plot_utils.plot_rf_stack(stream)
    # plt.show()

    station_db = {'HHR': tr_list}

    k, h, hk = rf_stacking.compute_hk_stack(station_db, 'HHR', include_t3=False)
    rf_plot_utils.plot_hk_stack(k, h, hk[0], save_file="c:\\temp\\Ps.png", show=False)
    rf_plot_utils.plot_hk_stack(k, h, hk[1], save_file="c:\\temp\\Ppps.png", show=False)
    stack = rf_stacking.compute_weighted_stack(hk, weighting=(0.5, 0.5))
    rf_plot_utils.plot_hk_stack(k, h, stack, save_file="c:\\temp\\synth_stack_0.5_0.5.png", show=False)
    # stack = rf_stacking.compute_weighted_stack(hk, weighting=(0.8, 0.2))
    # rf_plot_utils.plot_hk_stack(k, h, stack, save_file="c:\\temp\\synth_stack_0.8_0.2.png", show=False)

if __name__ == "__main__":
    main()
