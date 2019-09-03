#!/usr/bin/env python
"""Use a basic planar, 2-layer model of only the crust and the Moho to generate
synthetic arrival traces for known model characteristics. Intended to be used
for model validation.
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import obspy
import rf

import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking
from seismic.receiver_fn.rf_synthetic import generate_synth_rf

# pylint: disable=invalid-name

def main():
    H = 50  # km

    V_p = 6.4  # km/s, top (crust) layer
    V_s = 4.0  # km/s, top (crust) layer

    k = V_p/V_s
    print("H = {:3g}".format(H))
    print("k = {:3g}".format(k))

    # Generate function mapping ray parameter to teleseismic distance
    ts_distance = np.linspace(25, 95, 71)
    inc = np.zeros_like(ts_distance)
    model = obspy.taup.TauPyModel(model="iasp91")
    source_depth_km = 10.0
    for i, d in enumerate(ts_distance):
        ray = model.get_ray_paths(source_depth_km, d, phase_list=['P'])
        inc[i] = ray[0].incident_angle  # pylint: disable=unsupported-assignment-operation
    # end for
    plt.plot(inc, ts_distance, '-+')
    plt.xlabel('P-arrival inclination (deg)')
    plt.ylabel('Teleseismic distance (deg)')
    plt.title("Relationship between incident angle and teleseismic distance\n(source depth = {:2g} km)".format(
        source_depth_km))
    plt.grid("#80808080", linestyle=":")
    plt.savefig("c:\\temp\\inclination_dist.png", dpi=300)

    # Create interpolator to convert inclination to teleseismic distance
    interp_dist = interp1d(inc, ts_distance, bounds_error=False, fill_value=(np.max(ts_distance), np.min(ts_distance)))

    # Loop over valid range of inclinations and generate synthetic RFs
    inclinations = np.linspace(np.min(inc), np.max(inc), 10)
    # inclinations = np.array([1.0]*10)
    traces = []
    for inc_deg in inclinations:
        theta_p = np.deg2rad(inc_deg)
        p = np.sin(theta_p)/V_p

        t1 = H*(np.sqrt((k*k/V_p/V_p) - p*p) - np.sqrt(1.0/V_p/V_p - p*p))
        t2 = H*(np.sqrt((k*k/V_p/V_p) - p*p) + np.sqrt(1.0/V_p/V_p - p*p))
        print("Inclination {:3g} arrival times: {}".format(inc_deg, [t1, t2]))

        arrivals = [0, t1, t2]
        amplitudes = [1, 0.5, 0.4]
        window = (-5.0, 50.0)  # sec
        fs = 100.0
        _, signal = generate_synth_rf(arrivals, amplitudes, fs_hz=fs, window_sec=window)

        now = obspy.UTCDateTime.now()
        dt = float(window[1] - window[0])
        end = now + dt
        onset = now - window[0]
        header = {'network': 'SY', 'station': 'TST', 'location': 'GA', 'channel': 'HHR', 'sampling_rate': fs,
                  'starttime': now, 'endtime': end, 'onset': onset,
                  'station_latitude': -19.0, 'station_longitude': 137.0,
                  'slowness': p*rf_util.KM_PER_DEG, 'inclination': inc_deg,
                  'back_azimuth': 0, 'distance': float(interp_dist(inc_deg))}
        tr = rf.rfstream.RFTrace(data=signal, header=header)
        tr = tr.decimate(10, no_filter=True)
        traces.append(tr)
    # end for

    stream = rf.RFStream(traces)
    _ = rf_plot_utils.plot_rf_stack(stream, time_window=(-10, 30), save_file="c:\\temp\\synth_rf.png", dpi=300)

    station_db = {'HHR': traces}

    k, h, hk = rf_stacking.compute_hk_stack(station_db, 'HHR',
                                            h_range=np.linspace(40.0, 60.0, 101),
                                            k_range=np.linspace(1.5, 1.8, 151),
                                            include_t3=False)
    rf_plot_utils.plot_hk_stack(k, h, hk[0], title="Synthetic Ps component",
                                save_file="c:\\temp\\Ps.png", show=False)
    rf_plot_utils.plot_hk_stack(k, h, hk[1], title="Synthetic Ppps component",
                                save_file="c:\\temp\\Ppps.png", show=False)
    w = (0.5, 0.5)
    w_str = "w={:2g},{:2g}".format(*w)
    stack = rf_stacking.compute_weighted_stack(hk, weighting=w)
    rf_plot_utils.plot_hk_stack(k, h, stack, num=len(traces), title="Synthetic H-k stack ({})".format(w_str),
                                save_file="c:\\temp\\synth_stack_{}.png".format(w_str), show=False)
    # stack = rf_stacking.compute_weighted_stack(hk, weighting=(0.8, 0.2))
    # rf_plot_utils.plot_hk_stack(k, h, stack, save_file="c:\\temp\\synth_stack_0.8_0.2.png", show=False)

    max_loc = np.unravel_index(np.argmax(stack), stack.shape)
    h_max = h[max_loc[0], 0]
    k_max = k[0, max_loc[1]]
    print("Numerical solution (H, k) = ({:3g}, {:3g})".format(h_max, k_max))

if __name__ == "__main__":
    main()
