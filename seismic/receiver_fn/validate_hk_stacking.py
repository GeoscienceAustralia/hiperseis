#!/usr/bin/env python
"""Use a basic planar, 2-layer model of only the crust and the Moho to generate
synthetic arrival traces for known model characteristics. Intended to be used
for model validation.
"""

import os
import logging

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import obspy

import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking
from seismic.receiver_fn.rf_synthetic import synthesize_rf_dataset

# pylint: disable=invalid-name

logging.basicConfig()


def main():

    save_dist_plot = False
    output_folder = 'hk_validation'
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    # end if

    n = 0
    V_s = 3.8  # km/s, top (crust) layer

    H_range = np.linspace(30.0, 60.0, 4)  # km
    V_p_range = np.linspace(5.7, 7.22, 5)  # km/s
    include_t3 = True

    for H in H_range:
        for V_p in V_p_range:
            k = V_p/V_s
            log.info("H = {:.3f}".format(H))
            log.info("k = {:.3f}".format(k))

            # Generate function mapping ray parameter to teleseismic distance.
            # The distances are not strictly required for H-k stacking, but rf behaves better when they are there.
            ts_distance = np.linspace(25, 95, 71)
            inc = np.zeros_like(ts_distance)
            model = obspy.taup.TauPyModel(model="iasp91")
            source_depth_km = 10.0
            for i, d in enumerate(ts_distance):
                ray = model.get_ray_paths(source_depth_km, d, phase_list=['P'])
                inc[i] = ray[0].incident_angle  # pylint: disable=unsupported-assignment-operation
            # end for
            if save_dist_plot:
                plt.plot(inc, ts_distance, '-+')
                plt.xlabel('P-arrival inclination (deg)')
                plt.ylabel('Teleseismic distance (deg)')
                plt.title("Relationship between incident angle and teleseismic distance\n"
                          "(source depth = {:2.0f} km)".format(source_depth_km))
                plt.grid("#80808080", linestyle=":")
                plt.savefig(os.path.join(output_folder, "inclination_dist.png"), dpi=300)
            # end if

            # Create interpolator to convert inclination to teleseismic distance
            interp_dist = interp1d(inc, ts_distance, bounds_error=False,
                                   fill_value=(np.max(ts_distance), np.min(ts_distance)))

            # Loop over valid range of inclinations and generate synthetic RFs
            num_traces = 20
            inclinations = np.linspace(np.min(inc), np.max(inc), num_traces)
            distances = interp_dist(inclinations)
            final_sample_rate_hz = 10.0
            stream = synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, final_sample_rate_hz, log, include_t3)
            stream.write(os.path.join(output_folder, "synth_rf_data.h5"), format='h5')
            # Plot synthetic RFs
            rf_fig = rf_plot_utils.plot_rf_stack(stream, time_window=(-5, 30))
            plt.title("Exact H = {:.3f}, k = {:.3f}".format(H, k))
            rf_fig.savefig(os.path.join(output_folder, "synth_rf_{:03d}.png".format(n)), dpi=300)

            # Run H-k stacking on synthetic data
            station_db = {'HHR': stream}

            k_grid, h_grid, hk = rf_stacking.compute_hk_stack(station_db, 'HHR', include_t3=include_t3)
            fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk[0], title="Synthetic Ps component")
            fig.savefig(os.path.join(output_folder, "Ps_{:03d}.png".format(n)), dpi=300)
            plt.close()
            fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk[1], title="Synthetic PpPs component")
            fig.savefig(os.path.join(output_folder, "PpPs_{:03d}.png".format(n)), dpi=300)
            plt.close()
            if include_t3:
                fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk[2], title="Synthetic PpSs+PsPs component")
                fig.savefig(os.path.join(output_folder, "PpSs+PsPs_{:03d}.png".format(n)), dpi=300)
                plt.close()
                w = (0.34, 0.33, 0.33)
                w_str = "w={:1.2f},{:1.2f},{:1.2f}".format(*w)
            else:
                w = (0.5, 0.5)
                w_str = "w={:1.2f},{:1.2f}".format(*w)
            # end if
            stack = rf_stacking.compute_weighted_stack(hk, weighting=w)
            hk_fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, stack, num=len(stream),
                                                 title="Synthetic H-k stack ({})".format(w_str))
            plt.text(1.45, 65.0, "Expected H = {:.3f}, k = {:.3f}".format(H, k), color="#ffffff",
                     fontsize=16, fontweight='bold')
            # Determine H-k location of numerical maximum in the solution space.
            h_max, k_max = rf_stacking.find_global_hk_maximum(k_grid, h_grid, stack)
            log.info("Numerical solution (H, k) = ({:.3f}, {:.3f})".format(h_max, k_max))
            plt.scatter(k_max, h_max, marker='+', c="#000000", s=20)
            if k_max >= 1.7:
                plt.text(k_max - 0.01, h_max + 1, "Solution H = {:.3f}, k = {:.3f}".format(h_max, k_max),
                         color="#ffffff", fontsize=16, horizontalalignment='right')
            else:
                plt.text(k_max + 0.01, h_max + 1, "Solution H = {:.3f}, k = {:.3f}".format(h_max, k_max),
                         color="#ffffff", fontsize=16)
            # end if

            # Save result
            hk_fig.savefig(os.path.join(output_folder, "synth_stack_{}_{:03d}.png".format(w_str, n)), dpi=300)

            n += 1
        # end for
    # end for

# end main


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    main()
