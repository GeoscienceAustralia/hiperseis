#!/usr/bin/env python
# coding: utf-8
"""
Batch execution interfaces for wavefield continuation methods and solvers.
"""

import json
import logging

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

import click
import numpy as np
import h5py
from scipy import stats
import scipy.optimize as optimize

from seismic.inversion.wavefield_decomp.model_properties import LayerProps
from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.wavefield_continuation_tao import WfContinuationSuFluxComputer
from seismic.stream_quality_filter import curate_stream3c
from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster


DEFAULT_TEMP = 1.0
DEFAULT_AR = 0.4


def run_mcmc(waveform_data, config, output_file=None):

    # Create flux computer
    flux_computer_opts = config["su_energy_opts"]
    logging.info('Flux computer options:\n{}'.format(json.dumps(flux_computer_opts, indent=4)))
    fs_processing = flux_computer_opts["sampling_rate"]
    time_window = flux_computer_opts["time_window"]
    cut_window = flux_computer_opts["cut_window"]
    logging.info("Ingesting source data streams...")
    flux_comp = WfContinuationSuFluxComputer(waveform_data, fs_processing, time_window, cut_window)

    # Create model
    mantle_config = config["mantle_properties"]
    mantle_props = LayerProps(mantle_config['Vp'], mantle_config['Vs'], mantle_config['rho'], np.Infinity)
    logging.info('Mantle properties: {}'.format(mantle_props))

    layers = config["layers"]
    for i, layer_config in enumerate(layers):
        logging.info('Layer {}: {}'.format(i, layer_config))
    # end for

    # Solve model
    solver_opts = config["solver"]
    logging.info('Solver options:\n{}'.format(json.dumps(solver_opts, indent=4)))
    flux_window = flux_computer_opts["flux_window"]
    Vp = [layer["Vp"] for layer in layers]
    rho= [layer["rho"] for layer in layers]
    fixed_args = (flux_comp, mantle_props, Vp, rho, flux_window)
    bounds_min = []
    bounds_max = []
    for layer in layers:
        if "k_range" in layer:
            Vp_layer = layer["Vp"]
            bounds_min.extend([layer["H_range"][0], Vp_layer/layer["k_range"][1]])
            bounds_max.extend([layer["H_range"][1], Vp_layer/layer["k_range"][0]])
        else:
            assert "Vs_range" in layer
            bounds_min.extend([layer["H_range"][0], layer["Vs_range"][0]])
            bounds_max.extend([layer["H_range"][1], layer["Vs_range"][1]])
        # end if
    # end for
    bounds = optimize.Bounds(np.array(bounds_min), np.array(bounds_max))
    logging.info('Running MCMC solver...')
    temp = solver_opts.get("temp", DEFAULT_TEMP)
    burnin = solver_opts["burnin"]
    max_iter = solver_opts["max_iter"]
    target_ar = solver_opts.get("target_ar", DEFAULT_AR)
    collect_samples = solver_opts.get("collect_samples", None)
    soln = optimize_minimize_mhmcmc_cluster(
        mcmc_solver_wrapper, bounds, fixed_args, T=temp, burnin=burnin, maxiter=max_iter, target_ar=target_ar,
        collect_samples=collect_samples, logger=logging, verbose=True)
    logging.info('Result:\n{}'.format(soln))

    # Save solution
    # if soln.success and output_file is not None:
    #     pass  # TBD: write output file in HDF5 and include configu dictionary.
    #     with h5py.File(output_file, 'w') as f:
    #         f['config'] = json.dumps(config)
    # end with

    # # end if

    return True

# end func


def mcmc_solver_wrapper(model, obj_fn, mantle, Vp, rho, flux_window):
    num_layers = len(model)//2
    earth_model = []
    for i in range(num_layers):
        earth_model.append(LayerProps(Vp[i], model[2*i + 1], rho[i], model[2*i]))
    # end for
    earth_model = np.array(earth_model)
    energy, _, _ = obj_fn(mantle, earth_model, flux_window=flux_window)
    return energy
# end func


def curate_seismograms(data_all, curation_opts):

    def stream_snr_compute(_stream):
        _stream.taper(0.05)
        compute_vertical_snr(_stream)
    # end func

    def amplitude_nominal(_stream, max_amplitude):
        return ((np.max(np.abs(_stream[0].data)) <= max_amplitude) and
                (np.max(np.abs(_stream[1].data)) <= max_amplitude) and
                (np.max(np.abs(_stream[2].data)) <= max_amplitude))
    # end func

    logging.info("Curating input data...")
    logging.info('Curation options:\n{}'.format(json.dumps(curation_opts, indent=4)))

    if "baz_range" in curation_opts:
        # Filter by back-azimuth
        baz_range = curation_opts["baz_range"]
        data_all.curate(lambda _1, _2, stream: baz_range[0] <= stream[0].stats.back_azimuth <= baz_range[1])
    # end if

    # Apply curation to streams prior to rotation
    data_all.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

    # Rotate to ZRT coordinates
    data_all.apply(lambda stream: stream.rotate('NE->RT'))

    # Detrend the traces
    data_all.apply(lambda stream: stream.detrend('linear'))

    if "freq_min" in curation_opts:
        # Run high pass filter to remove high amplitude, low freq noise, if present.
        f_min = curation_opts["freq_min"]
        data_all.apply(lambda stream: stream.filter('highpass', freq=f_min, corners=2, zerophase=True))
    # end if

    if "min_snr" in curation_opts:
        min_snr = curation_opts["min_snr"]
        # Compute SNR of Z component to use as a quality metric
        data_all.apply(stream_snr_compute)

        # Filter by SNR
        data_all.curate(lambda _1, _2, stream: stream[0].stats.snr_prior >= min_snr)
    # end if

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

    if "max_raw_amplitude" in curation_opts:
        # Filter streams with spuriously high amplitude
        max_amp = curation_opts["max_raw_amplitude"]
        data_all.curate(lambda _1, _2, stream: amplitude_nominal(stream, max_amp))
    # end if

# end func


@click.command()
@click.argument('config_file', type=click.File('r'), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--network', type=str, required=True)
@click.option('--station', type=str, required=True)
@click.option('--location', type=str, default='')
@click.option('--output-file', type=click.Path(dir_okay=False))
def main(config_file, waveform_file, network, station, location='', output_file=None):
    """Runner for analysis of single station. For multiple stations, run separate processes.

    The output file is in HDF5 format. The configuration details are added to the output file for traceability.

    :param config_file:
    :param waveform_file:
    :param output_file:
    :return: None
    """
    config = json.load(config_file)
    # logging.info("Config:\n{}".format(json.dumps(config, indent=4)))
    logging.info("Waveform source: {}".format(waveform_file))
    logging.info("Network.Station: {}.{}".format(network, station))
    logging.info("Destination file: {}".format(output_file))

    stype = config['solver']['type']
    if stype.lower() == 'mcmc':
        runner = run_mcmc
    else:
        logging.error("Unknown solver type: {}".format(stype))
        return
    # end if

    # Load input data
    logging.info('Ingesting waveform file {}'.format(waveform_file))
    waveform_data = NetworkEventDataset(waveform_file, network=network, station=station, location=location)
    config.update({"waveform_file": waveform_file})

    # Trim entire dataset to max time window required.
    time_window = config["su_energy_opts"]["time_window"]
    # Trim streams to time window
    waveform_data.apply(lambda stream:
                        stream.trim(stream[0].stats.onset + time_window[0], stream[0].stats.onset + time_window[1]))

    # Curate input data if curation options given
    if "curation_opts" in config:
        curation_opts = config["curation_opts"]
        if curation_opts:
            curate_seismograms(waveform_data, curation_opts)
        # end if
    # end if

    status = runner(waveform_data.station(station).values(), config, output_file)
    if not status:
        logging.error('Runner failed')
    # end if

# end func


if __name__ == '__main__':
    main()
# end if
