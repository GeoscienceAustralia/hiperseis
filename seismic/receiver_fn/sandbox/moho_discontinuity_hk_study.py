#!/usr/bin/env python
"""Study behavior of H-k stacking in presence of Moho discontinuity. Ticket PST-433.
"""

import os
import logging

import numpy as np
# import obspy
import matplotlib.pyplot as plt

import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking
from seismic.receiver_fn.rf_synthetic import synthesize_rf_dataset, convert_inclination_to_distance

# pylint: disable=invalid-name,logging-format-interpolation

logging.basicConfig()


def main():
    output_folder = 'moho_discontinuity'
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    # end if

    # Generate data files corresponding to two different Moho depths
    V_s = 3.8  # km/s
    V_p = 6.4  # km/s
    H_0 = 40  # km
    k = V_p/V_s
    log.info("H_0 = {:.3f}".format(H_0))
    log.info("k = {:.3f}".format(k))

    # Loop over valid range of inclinations and generate synthetic RFs
    num_traces = 15
    inclinations = np.linspace(14.0, 27.0, num_traces)
    distances = convert_inclination_to_distance(inclinations)
    final_sample_rate_hz = 10.0
    # stream_0 = synthesize_rf_dataset(H_0, V_p, V_s, inclinations, distances, final_sample_rate_hz, log, baz=0)
    stream_0 = synthesize_rf_dataset(H_0, V_p, V_s, inclinations, distances, final_sample_rate_hz, log,
                                    #  amplitudes=[1, 0.5, 0.2], baz=0)
                                     amplitudes=[1, 0.4, 0.3], baz=0)
    stream_0.write(os.path.join(output_folder, "synth_rf_data_H={}.h5".format(H_0)), format='h5')

    H_1 = 50
    # stream_1 = synthesize_rf_dataset(H_1, V_p, V_s, inclinations, distances, final_sample_rate_hz, log, baz=90)
    stream_1 = synthesize_rf_dataset(H_1, V_p, V_s, inclinations, distances, final_sample_rate_hz, log,
                                    #  amplitudes=[1, 0.4, 0.3], baz=90)
                                     amplitudes=[1, 0.5, 0.2], baz=90)
    stream_1.write(os.path.join(output_folder, "synth_rf_data_H={}.h5".format(H_1)), format='h5')

    # Plot each dataset separately, then aggregated together
    rf_fig_0 = rf_plot_utils.plot_rf_stack(stream_0, time_window=(-5, 30))
    plt.title("Exact H = {:.3f}, k = {:.3f}".format(H_0, k))
    rf_fig_0.savefig(os.path.join(output_folder, "synth_rf_H={}.png".format(H_0)), dpi=300)
    rf_fig_1 = rf_plot_utils.plot_rf_stack(stream_1, time_window=(-5, 30))
    plt.title("Exact H = {:.3f}, k = {:.3f}".format(H_1, k))
    rf_fig_1.savefig(os.path.join(output_folder, "synth_rf_H={}.png".format(H_1)), dpi=300)
    stream_all = stream_0 + stream_1
    rf_fig_all = rf_plot_utils.plot_rf_stack(stream_all, time_window=(-5, 30))
    plt.title("Exact H = {:.3f}+{:.3f}, k = {:.3f}".format(H_0, H_1, k))
    rf_fig_all.savefig(os.path.join(output_folder, "synth_rf_H={}+{}.png".format(H_0, H_1)), dpi=300)

    # Run H-k stacking on synthetic data
    station_db = {'HHR': stream_all}
    k_grid, h_grid, hk = rf_stacking.compute_hk_stack(station_db['HHR'], include_t3=False, root_order=2)
    w = (0.5, 0.5)
    w_str = "w={:1.2f},{:1.2f}".format(*w)
    stack = rf_stacking.compute_weighted_stack(hk, weighting=w)
    hk_fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, stack, num=len(stream_all),
                                         title="Synthetic H-k stack ({})".format(w_str))
    maxima = rf_stacking.find_local_hk_maxima(k_grid, h_grid, stack, min_rel_value=0.9)
    log.info("Numerical solutions:{}".format(maxima))
    plt.text(k, H_0 + 1.5, "Solution $H_0$ = {:.3f}, k = {:.3f}".format(H_0, k),
             color="#ffffff", fontsize=16, horizontalalignment='left')
    plt.text(k, H_1 + 1.5, "Solution $H_1$ = {:.3f}, k = {:.3f}".format(H_1, k),
             color="#ffffff", fontsize=16, horizontalalignment='left')
    for h_max, k_max, _, _, _ in maxima:
        plt.scatter(k_max, h_max, marker='+', c="#000000", s=20)
    # end for
    # Save result
    hk_fig.savefig(os.path.join(output_folder, "synth_stack_{}.png".format(w_str)), dpi=300)

# end main


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    main()
