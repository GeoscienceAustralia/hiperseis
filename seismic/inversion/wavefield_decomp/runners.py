#!/usr/bin/env python
# coding: utf-8
"""
Batch execution interfaces for wavefield continuation methods and solvers.
"""

import os
import json
import logging
from datetime import datetime

import click
import numpy as np
from scipy import stats
import scipy.optimize as optimize
import h5py

from seismic.inversion.wavefield_decomp.model_properties import LayerProps
from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.wavefield_continuation_tao import WfContinuationSuFluxComputer
from seismic.stream_quality_filter import curate_stream3c
from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster


# Custom logging format to add timestamps to each output line.
LOG_FORMAT = {'fmt': '%(asctime)s %(levelname)-8s %(message)s',
              'datefmt': '%Y-%m-%d %H:%M:%S'}

# Default MCMC solver "temperature" parameter.
DEFAULT_TEMP = 1.0
# Default target acceptance rate for MCMC solver.
DEFAULT_AR = 0.4


def run_mcmc(waveform_data, config, logger):
    """
    Top level runner function for MCMC solver on SU flux minimization for given settings.

    :param waveform_data: Iterable container of obspy.Stream objects.
    :param config: Dict of job settings. See example files for fields and format of settings.
    :param logger: Log message receiver. Optional, pass None for no logging.
    :return: Solution object based on scipy.optimize.OptimizeResult
    """

    # Create flux computer
    flux_computer_opts = config["su_energy_opts"]
    logger.info('Flux computer options:\n{}'.format(json.dumps(flux_computer_opts, indent=4)))
    fs_processing = flux_computer_opts["sampling_rate"]
    time_window = flux_computer_opts["time_window"]
    cut_window = flux_computer_opts["cut_window"]
    if logger:
        logger.info("Ingesting source data streams...")
    flux_comp = WfContinuationSuFluxComputer(waveform_data, fs_processing, time_window, cut_window)

    # Create model
    mantle_config = config["mantle_properties"]
    mantle_props = LayerProps(mantle_config['Vp'], mantle_config['Vs'], mantle_config['rho'], np.Infinity)
    logger.info('Mantle properties: {}'.format(mantle_props))

    layers = config["layers"]
    for i, layer_config in enumerate(layers):
        logger.info('Layer {}: {}'.format(i, layer_config))
    # end for

    # Solve model
    solver_opts = config["solver"]
    logger.info('Solver options:\n{}'.format(json.dumps(solver_opts, indent=4)))
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
    logger.info('Running MCMC solver...')
    temp = solver_opts.get("temp", DEFAULT_TEMP)
    burnin = solver_opts["burnin"]
    max_iter = solver_opts["max_iter"]
    target_ar = solver_opts.get("target_ar", DEFAULT_AR)
    collect_samples = solver_opts.get("collect_samples", None)
    soln = optimize_minimize_mhmcmc_cluster(
        mcmc_solver_wrapper, bounds, fixed_args, T=temp, burnin=burnin, maxiter=max_iter, target_ar=target_ar,
        collect_samples=collect_samples, logger=logger, verbose=True)

    # TODO: Compute energy flux and wavefield at top of mantle per seismogram for each solution, and return with soln

    return soln

# end func


def mcmc_solver_wrapper(model, obj_fn, mantle, Vp, rho, flux_window):
    """
    Wrapper callable for passing to MCMC solver in scipy style, which unpacks inputs into
    vector variables for solver.

    :param model: Per-layer model values (the vector being solved for) as flat array of
        (H, Vs) value pairs ordered by layer.
    :param obj_fn: Callable to WfContinuationSuFluxComputer to compute SU flux.
    :param mantle: Mantle properties stored in class LayerProps instance.
    :param Vp: Array of Vp values ordered by layer.
    :param rho: Arroy of rho values ordered by layer.
    :param flux_window: Pair of floats indicating the time window over which to perform SU flux integration
    :return: Integrated SU flux energy at top of mantle
    """
    num_layers = len(model)//2
    earth_model = []
    for i in range(num_layers):
        earth_model.append(LayerProps(Vp[i], model[2*i + 1], rho[i], model[2*i]))
    # end for
    earth_model = np.array(earth_model)
    energy, _, _ = obj_fn(mantle, earth_model, flux_window=flux_window)
    return energy
# end func


def curate_seismograms(data_all, curation_opts, logger):
    """
    Curation function to remove bad data from streams.

    :param data_all: NetworkEventDataset containing seismograms to curate.
    :param curation_opts: Dict containing curation options.
    :param logger: Logger for emitting log messages
    :return: None, curation operates directly on data_all
    """

    def stream_snr_compute(_stream):
        _stream.taper(0.05)
        compute_vertical_snr(_stream)
    # end func

    def amplitude_nominal(_stream, max_amplitude):
        return ((np.max(np.abs(_stream[0].data)) <= max_amplitude) and
                (np.max(np.abs(_stream[1].data)) <= max_amplitude) and
                (np.max(np.abs(_stream[2].data)) <= max_amplitude))
    # end func

    logger.info("Curating input data...")
    logger.info('Curation options:\n{}'.format(json.dumps(curation_opts, indent=4)))

    # Apply curation to streams prior to rotation
    data_all.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

    if "baz_range" in curation_opts:
        # Filter by back-azimuth
        baz_range = curation_opts["baz_range"]
        data_all.curate(lambda _1, _2, stream: baz_range[0] <= stream[0].stats.back_azimuth <= baz_range[1])
    # end if

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


def save_mcmc_solution(soln_configs, input_file, output_file, job_timestamp, logger=None):
    """
    Save solution to HDF5 file.

    :param soln_configs: List of (solution, configuration) pairs to save (one per station).
    :param input_file: Name of input file. Saved to job node for traceability
    :param output_file: Name of output file to create
    :param job_timestamp: Job timestamp that will be used to generate the top level job group
    :param logger: [OPTIONAL] Log message destination
    :return: None
    """
    assert isinstance(soln_configs, list)
    assert isinstance(job_timestamp, str)
    # TODO: migrate this to member of a new class for encapsulating an MCMC solution

    # Version string for layout of data in a node containing MCMC solution
    FORMAT_VERSION = '0.2'

    # Convert timestamp to valid Python identifier
    job_timestamp = 'T' + job_timestamp.replace('-', '_').replace(' ', '__').replace(':', '').replace('.', '_')

    with h5py.File(output_file, 'a') as h5f:
        job_root = h5f.create_group(job_timestamp)
        job_root.attrs['input_file'] = input_file
        for soln, config in soln_configs:
            station_id = config.get("station_id")
            if not station_id:
                if logger:
                    logger.warning('Unidentified solution, no station id!')
                continue
            # end if
            if soln is None or not soln.success:
                if logger:
                    logger.warning('Solver failed for station {}'.format(station_id))
                continue
            # end if
            station_id = station_id.replace('.', '_')
            station_node = job_root.create_group(station_id)
            station_node.attrs['config'] = json.dumps(config)
            station_node.attrs['format_version'] = FORMAT_VERSION
            try:
                station_node['x'] = soln.x
                # Each cluster may be a different size, so we can't dump directly into h5py (which
                # requires hyper-rectangular array shape).
                cluster_node = station_node.create_group('clusters')
                for idx, cluster in enumerate(soln.clusters):
                    cluster_node[str(idx)] = cluster
                # end for
                station_node['bins'] = soln.bins
                station_node['distribution'] = soln.distribution
                station_node['acceptance_rate'] = soln.acceptance_rate
                station_node['success'] = soln.success
                station_node['status'] = soln.status
                station_node['message'] = soln.message
                station_node['fun'] = soln.fun
                station_node['jac'] = h5py.Empty('f') if soln.jac is None else soln.jac
                station_node['nfev'] = soln.nfev
                station_node['njev'] = soln.njev
                station_node['nit'] = soln.nit
                station_node['maxcv'] = h5py.Empty('f') if soln.maxcv is None else soln.maxcv
                station_node['samples'] = h5py.Empty('f') if soln.samples is None else soln.samples
                station_node['bounds'] = np.array([soln.bounds.lb, soln.bounds.ub])
                station_node['version'] = soln.version
            except TypeError as exc:
                if logger:
                    logger.error('Error saving station {} solution'.format(station_id))
                    logger.error(repr(exc))
                # end  if
            # end try
        # end for
    # end with

# end func


def run_station(config_file, waveform_file, network, station, location, logger):
    """Runner for analysis of single station. For multiple stations, set up config file to run batch
    job using mpi_job CLI.

    The output file is in HDF5 format. The configuration details are added to the output file for traceability.

    :param config_file: Config filename specifying job settings
    :param waveform_file: Event waveform source file for seismograms, generated using extract_event_traces.py script
    :param network: Network code of station to analyse
    :param station: Station code to analyse
    :param location: Location code of station to analyse. Can be '' (empty string) if not set.
    :return: Pair containing (solution, configuration) containers. Configuration will have additional traceability
        information.
    """
    with open(config_file, 'r') as cf:
        config = json.load(cf)
    # end with
    # logger.info("Config:\n{}".format(json.dumps(config, indent=4)))
    station_id = "{}.{}.{}".format(network, station, location)
    logger.info("Network.Station.Location: {}".format(station_id))
    config.update({"station_id": station_id})

    stype = config['solver']['type']
    if stype.lower() == 'mcmc':
        runner = run_mcmc
    else:
        logger.error("Unknown solver type: {}".format(stype))
        return (None, config)
    # end if

    # Load input data
    logger.info('Ingesting waveform file {}'.format(waveform_file))
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
            curate_seismograms(waveform_data, curation_opts, logger)
        # end if
    # end if

    soln = runner(waveform_data.station(station).values(), config, logger)

    return soln, config

# end func


@click.command()
@click.argument('config_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Event waveform source file for seismograms, generated using extract_event_traces.py script')
@click.option('--network', type=str, required=True, help='Network code for the station to load')
@click.option('--station', type=str, required=True, help='Station code whose data is to be loaded')
@click.option('--location', type=str, default='', show_default=True, help='Location code for the station')
@click.option('--output-file', type=click.Path(dir_okay=False),
              help='Name of the output file in which solutions should be saved')
def station_job(config_file, waveform_file, network, station, location='', output_file=None):
    """
    CLI dispatch function for single station. See help strings for option documentation.

    Example usage:
        python runners.py example_config.json \
        --waveform-file /g/data/ha3/am7399/shared/OA_RF_analysis/OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5 \
        --network OA --station CD23 --location 0M --output-file out_test.h5

    :param config_file: JSON file containing job configuration parameters.
    :type config_file: str or pathlib.Path
    :return: Integer status code
    """
    job_timestamp = str(datetime.now())

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    log_fmt = logging.Formatter(**LOG_FORMAT)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_fmt)
    logger.addHandler(console_handler)

    logger.info("Waveform source: {}".format(waveform_file))
    logger.info("Destination file: {}".format(output_file))

    soln_config = run_station(config_file, waveform_file, network, station, location, logger=logger)

    # Save solution
    save_mcmc_solution([soln_config], waveform_file, output_file, job_timestamp, logger)

    return 0
# end func


@click.command()
@click.argument('config_file', type=click.File('r'), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Event waveform source file for seismograms, generated using extract_event_traces.py script')
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output file in which solutions should be saved')
def mpi_job(config_file, waveform_file, output_file):
    """
    CLI dispatch function for MPI run over batch of stations. See help strings for option documentation.

    Example MPI usage:
        mpiexec -n 8 python runners.py example_batch.json \
        --waveform-file /g/data/ha3/am7399/shared/OA_RF_analysis/OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5 \
        --output-file test_batch_output.h5

    :param config_file: JSON file containing batch configuration parameters.
    :type config_file: str or pathlib.Path
    :return: Integer status code
    """

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()  # TODO: Use this to more intelligently split the jobs across available processors
    rank = comm.Get_rank()

    job_timestamp = comm.bcast(str(datetime.now()), root=0)

    logger = logging.getLogger(__name__ + str(rank))
    logger.setLevel(logging.INFO)
    log_fmt = logging.Formatter(**LOG_FORMAT)
    job_id = os.getenv('PBS_JOBID')
    if job_id is None:
        job_id = comm.bcast(os.getpid(), root=0)
    # end if
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_fmt)
    file_handler = logging.FileHandler('_'.join([job_id, str(rank)]) + '.log')
    file_handler.setFormatter(log_fmt)
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    logger.info("Waveform source: {}".format(waveform_file))

    batch_config = json.load(config_file)
    jobs = batch_config.get("jobs", None)
    if jobs is None:
        logger.error("Job list missing, file {} not a batch configuration.".format(config_file.name))
        return 1
    # end if

    if rank == 0:
        node_args = [(job_config_file, waveform_file) + tuple(job_id.split('.'))
                     for job_id, job_config_file in jobs.items()]
    else:
        node_args = None
    # end if

    node_args = comm.scatter(node_args, root=0)
    soln_config = run_station(*node_args, logger=logger)
    soln_config = comm.gather(soln_config, root=0)

    if rank == 0:
        save_mcmc_solution(soln_config, waveform_file, output_file, job_timestamp, logger)
    # end if

    return 0

# end func


def main():
    import sys

    if len(sys.argv) <= 6:
        status = mpi_job()
    else:
        status = station_job()
    # end if

    sys.exit(status)
# end func


if __name__ == '__main__':
    main()
# end if
