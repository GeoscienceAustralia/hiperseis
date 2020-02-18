#!/usr/bin/env python
# coding: utf-8

import logging

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

import numpy as np
from scipy import stats
import scipy.optimize as optimize

from tqdm.auto import tqdm
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.wavefield_continuation_tao import WfContinuationSuFluxComputer
from seismic.inversion.wavefield_decomp.model_properties import LayerProps
from seismic.stream_quality_filter import curate_stream3c
from seismic.receiver_fn.rf_util import compute_vertical_snr

# pylint: disable=invalid-name, logging-format-interpolation, too-many-arguments, too-many-statements, too-many-locals


def example_1():
    # Example 1: Computing energy for a single model proposition.
    crust_props = LayerProps(Vp_c, 3.7, rho_c, 35)
    single_layer_model = [crust_props]

    logging.info("Computing single point mean SU flux...")
    energy, energy_per_event, wf_mantle = flux_comp(mantle_props, single_layer_model)
    logging.info(energy)
# end func


def example_2(station):
    # Example 2: Computing energy across a parametric space of models.

    def plot_Esu_space(H, k, Esu, network, station, save=True, show=True):
        colmap = 'plasma'
        plt.figure(figsize=(16, 12))
        plt.contourf(k, H, Esu, levels=50, cmap=colmap)
        plt.colorbar()
        plt.contour(k, H, Esu, levels=10, colors='k', linewidths=1, antialiased=True)
        plt.xlabel('Crustal $\kappa$', fontsize=14)
        plt.ylabel('Crustal $H$ (km)', fontsize=14)
        plt.tick_params(right=True, labelright=True, which='both')
        plt.tick_params(top=True, labeltop=True, which='both')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.minorticks_on()
        plt.xlim(np.min(k), np.max(k))
        plt.ylim(np.min(H), np.max(H))
        plt.grid(linestyle=':', color="#80808080")
        plt.title('{}.{} $E_{SU}$ energy vs. crustal properties'.format(network, station), fontsize=20, y=1.05)
        if save:
            plt.savefig('example_{}.{}_crust_props.png'.format(network, station), dpi=300)
        if show:
            plt.show()
        plt.close()
    # end func

    logging.info("Computing 2D parametric space mean SU flux...")
    H_space = np.linspace(25, 45, 51)
    k_space = np.linspace(1.5, 2.1, 51)
    H, k, Esu = flux_comp.grid_search(mantle_props, [LayerProps(Vp_c, None, rho_c, None)], 0, H_space, k_space,
                                      flux_window=FLUX_WINDOW)

    plot_Esu_space(H, k, Esu, 'OA', station)
# end func


def example_3():
    # Example 3: Using a global energy minimization solver to find solution.

    logging.info("Computing optimal crust properties by SU flux minimization...")

    # Use single layer (crust only) model in this example.
    Vp = [Vp_c]
    rho = [rho_c]
    fixed_args = (flux_comp, mantle_props, Vp, rho, FLUX_WINDOW)
    H_initial = 40.0
    k_initial = np.mean((k_min, k_max))
    model_initial = np.array([H_initial, Vp_c/k_initial])
    H_min, H_max = (25.0, 45.0)
    bounds = optimize.Bounds([H_min, Vp_c/k_max], [H_max, Vp_c/k_min])

    # Find local minimum relative to initial guess.
    soln = optimize.minimize(objective_fn, model_initial, fixed_args, bounds=bounds)
    H_crust, Vs_crust = soln.x
    logging.info('Success = {}, Iterations = {}, Function evaluations = {}'.format(soln.success, soln.nit, soln.nfev))
    logging.info('Solution H_crust = {}, Vs_crust = {}, SU energy = {}'.format(H_crust, Vs_crust, soln.fun))
# end func


def example_4():
    # Demonstrate syntactic usage of scipy global optimizers:
    k_initial = np.mean((k_min, k_max))
    model_initial_poor = np.array([35, Vp_c/k_initial])

    Vp = [Vp_c]
    rho = [rho_c]
    fixed_args = (flux_comp, mantle_props, Vp, rho, FLUX_WINDOW)
    H_min, H_max = (25.0, 45.0)
    bounds = optimize.Bounds([H_min, Vp_c/k_max], [H_max, Vp_c/k_min])

    # - Basin hopping
    # logging.info('Trying basinhopping...')
    # soln_bh = optimize.basinhopping(objective_fn, model_initial_poor, T=0.3,
    #                                 minimizer_kwargs={'args': fixed_args, 'bounds': bounds})
    # logging.info('Result:\n{}'.format(soln_bh))

    # - Differential evolution
    logging.info('Trying differential_evolution...')
    soln_de = optimize.differential_evolution(objective_fn, bounds, fixed_args, workers=-1,
                                              popsize=25, tol=1.0e-3, mutation=(0.5, 1.2), recombination=0.5)
    logging.info('Result:\n{}'.format(soln_de))

    # - SHGO (VERY EXPENSIVE AND/OR not convergent)
    # logging.info('Trying shgo...')
    # soln_shgo = optimize.shgo(objective_fn, list(zip(bounds.lb, bounds.ub)), fixed_args,
    #                           options={'f_min': 0, 'f_tol': 1.0e-3, 'maxev': 10000})
    # logging.info('Result:\n{}'.format(soln_shgo))

    # - Dual annealing
    # logging.info('Trying dual_annealing...')
    # soln_da = optimize.dual_annealing(objective_fn, list(zip(bounds.lb, bounds.ub)), fixed_args, x0=model_initial_poor,
    #                                   initial_temp=2000.0, maxfun=10000)
    # logging.info('Result:\n{}'.format(soln_da))


    # - Differential evolution
    logging.info('Trying BRUTE FORCE...')
    H_brute, Vs_brute = optimize.brute(objective_fn, tuple(zip(bounds.lb, bounds.ub)), fixed_args, Ns=51, workers=-1)
    logging.info('Result:\n{}'.format((H_brute, Vs_brute)))

# end func


def example_5():
    # Example 5: Adding a sedimentary layer and directly using global minimizer

    # Assumed sediment property constants
    Vp_s = 2.1
    rho_s = 1.97

    Vp = [Vp_s, Vp_c]
    rho = [rho_s, rho_c]
    fixed_args = (flux_comp, mantle_props, Vp, rho, FLUX_WINDOW)
    H_initial = [1.0, 35.0]  # sediment, crust
    Vs_initial = [0.8, 3.4]  # sediment, crust
    model_initial_sed = np.array(zip(H_initial, Vs_initial))
    H_sed_min, H_sed_max = (0, 3.5)
    Vs_sed_min, Vs_sed_max = (0.3, 2.5)
    H_cru_min, H_cru_max = (20.0, 60.0)
    Vs_cru_min, Vs_cru_max = (Vp_c/k_max, Vp_c/k_min)
    bounds = optimize.Bounds([H_sed_min, Vs_sed_min, H_cru_min, Vs_cru_min],
                             [H_sed_max, Vs_sed_max, H_cru_max, Vs_cru_max])

    logging.info('Differential_evolution (sedimentary)...')
    soln_de = optimize.differential_evolution(objective_fn, bounds, fixed_args, workers=-1,
                                              popsize=25, tol=1.0e-3, mutation=(0.5, 1.2), recombination=0.5)
    logging.info('Result:\n{}'.format(soln_de))
# end func


# def job_caller(i, j, callable, mantle, earth_model, flux_window):
#     energy, _, _ = callable(mantle, earth_model, flux_window=flux_window)
#     return (i, j, energy)
# # end func


def objective_fn(model, callable, mantle, Vp, rho, flux_window):
    num_layers = len(model)//2
    earth_model = []
    for i in range(num_layers):
        earth_model.append(LayerProps(Vp[i], model[2*i + 1], rho[i], model[2*i]))
    # end for
    earth_model = np.array(earth_model)
    energy, _, _ = callable(mantle, earth_model, flux_window=flux_window)
    return energy
# end func


if __name__ == "__main__":

    def stream_snr_compute(stream):
        stream.taper(0.05)
        compute_vertical_snr(stream)
    # end func

    def amplitude_nominal(stream, max_amplitude):
        return ((np.max(np.abs(stream[0].data)) <= max_amplitude) and
                (np.max(np.abs(stream[1].data)) <= max_amplitude) and
                (np.max(np.abs(stream[2].data)) <= max_amplitude))
    # end func

    target_station = 'BT23'
    # target_station = 'CI23'
    logging.info("Loading input file...")
    src_file = (r"/g/data/ha3/am7399/shared/OA_RF_analysis/" +
                r"OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5")
    data_all = NetworkEventDataset(src_file, network='OA', station=target_station, location='0M')

    # Time window of original data to use for processing. All traces must have at least this extent
    # about the onset time.
    TIME_WINDOW = (-20, 50)
    # Narrower time window used for integration of energy flux
    FLUX_WINDOW = (-10, 20)
    # Cut window for selecting central wavelet
    CUT_WINDOW = (-5, 30)

    # -----------------------------------------------------------------------------
    # Apply windowing, filtering and QC to loaded dataset before passing to Tao's algorithm.
    logging.info("Cleaning input data...")


    # Trim streams to time window
    data_all.apply(lambda stream:
                   stream.trim(stream[0].stats.onset + TIME_WINDOW[0], stream[0].stats.onset + TIME_WINDOW[1]))

    # Apply curation to streams prior to rotation
    data_all.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

    # Rotate to ZRT coordinates
    data_all.apply(lambda stream: stream.rotate('NE->RT'))

    # Detrend the traces
    data_all.apply(lambda stream: stream.detrend('linear'))

    # Run high pass filter to remove high amplitude, low freq noise, if present.
    f_min = 0.05
    data_all.apply(lambda stream: stream.filter('highpass', freq=f_min, corners=2, zerophase=True))

    # Compute SNR of Z component to use as a quality metric
    data_all.apply(stream_snr_compute)

    # Filter by SNR
    data_all.curate(lambda _1, _2, stream: stream[0].stats.snr_prior >= 3.0)

    # It does not make sense to filter by similarity, since these are raw waveforms, not RFs,
    # and the waveform will be dominated by the source waveform which differs for each event.

    # Filter streams with incorrect number of traces
    discard = []
    for sta, ev_db in data_all.by_station():
        num_pts = np.array([tr.stats.npts for st in ev_db.values() for tr in st])
        expected_pts = stats.mode(num_pts)[0][0]
        for evid, stream in ev_db.items():
            if ((stream[0].stats.npts != expected_pts) or
                (stream[1].stats.npts != expected_pts) or
                (stream[2].stats.npts != expected_pts)):
                discard.append((sta, evid))
            # end if
        # end for
    # end for
    data_all.prune(discard)

    # Filter streams with spuriously high amplitude
    MAX_AMP = 10000
    data_all.curate(lambda _1, _2, stream: amplitude_nominal(stream, MAX_AMP))

    # -----------------------------------------------------------------------------
    # Pass cleaned up data set for test station to flux computer class.
    data_OA = data_all.station(target_station)
    fs_processing = 10.0  # Hz
    logging.info("Ingesting source data streams...")
    flux_comp = WfContinuationSuFluxComputer(data_OA.values(), fs_processing, TIME_WINDOW, CUT_WINDOW)

    # Define bulk properties of mantle (lowermost half-space)
    mantle_props = LayerProps(vp=8.0, vs=4.5, rho=3.3, thickness=np.Infinity)

    # Assumed crust property constants
    Vp_c = 6.4
    rho_c = 2.7
    k_min, k_max = (1.5, 2.1)

    # -----------------------------------------------------------------------------
    # Example 1: Computing energy for a single model proposition.
    example_1()

    # -----------------------------------------------------------------------------
    # Example 2: Computing energy across a parametric space of models.
    example_2(target_station)

    # -----------------------------------------------------------------------------
    # Example 3: Using a global energy minimization solver to find solution.
    example_3()

    # -----------------------------------------------------------------------------
    # Example 4: Demonstrate syntactic usage of scipy global optimizers.
    # example_4()

    # -----------------------------------------------------------------------------
    # Example 5: Adding a sedimentary layer and directly using global minimizer
    example_5()

# end if
