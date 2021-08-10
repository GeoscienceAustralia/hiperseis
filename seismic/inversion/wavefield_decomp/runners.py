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
import scipy.optimize as optimize
import h5py

from seismic.model_properties import LayerProps
from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.wavefield_continuation_tao import WfContinuationSuFluxComputer
from seismic.stream_quality_filter import curate_seismograms
from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster, DEFAULT_CLUSTER_EPS


# pylint: disable=invalid-name, logging-format-interpolation


# Custom logging format to add timestamps to each output line.
LOG_FORMAT = {'fmt': '%(asctime)s %(levelname)-8s %(message)s',
              'datefmt': '%Y-%m-%d %H:%M:%S'}

# Default MCMC solver "temperature" parameter.
DEFAULT_TEMP = 1.0
# Default target acceptance rate for MCMC solver.
DEFAULT_AR = 0.4
DEFAULT_ROTATION = 'NE->RT'


def run_mcmc(waveform_data, config, logger):
    """
    Top level runner function for MCMC solver on SU flux minimization for given settings.

    :param waveform_data: Collection of waveform data to process
    :type waveform_data: Iterable obspy.Stream objects
    :param config: Dict of job settings. See example files for fields and format of settings.
    :type config: dict
    :param logger: Log message receiver. Optional, pass None for no logging.
    :type logger: logging.Logger or NoneType
    :return: Solution object based on scipy.optimize.OptimizeResult
    :rtype: scipy.optimize.OptimizeResult with additional custom attributes
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
    rho = [layer["rho"] for layer in layers]
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
    cluster_eps = solver_opts.get("cluster_eps", DEFAULT_CLUSTER_EPS)
    N = solver_opts.get("max_solutions", 3)
    soln = optimize_minimize_mhmcmc_cluster(
        mcmc_solver_wrapper, bounds, fixed_args, T=temp, N=N, burnin=burnin, maxiter=max_iter, target_ar=target_ar,
        cluster_eps=cluster_eps, collect_samples=collect_samples, logger=logger)

    # Record number of independent events processed
    soln.num_input_seismograms = len(waveform_data)

    if soln.success:
        # Compute energy flux per seismogram for each solution, and return with soln
        Esu_per_x = []
        for i, _x in enumerate(soln.x):
            num_layers = len(layers)
            earth_model = []
            for j in range(num_layers):
                earth_model.append(LayerProps(Vp[j], _x[2*j + 1], rho[j], _x[2*j]))
            # end for
            earth_model = np.array(earth_model)
            energy, energy_per_event, _ = flux_comp(mantle_props, earth_model, flux_window=flux_window)
            # Check computed energy matchs what solver produced
            assert np.isclose(energy, soln.fun[i])
            Esu_per_x.append(energy_per_event)
        # end for
        # Upward S-wave energy at top of mantle per seismogram per solution point.
        soln.esu = Esu_per_x

        # Compute per-event seismograms at bottom of layers that flagged it, and return with soln.
        # For now we just re-use the original flux computer, but if we want to use a broader
        # cut_window here, we will need to create a new one.
        subsurface = {}
        for i, layer in enumerate(layers):
            if bool(layer.get("save_seismogram", False)):
                layer_name = layer["name"]
                base_seismograms = []
                for _x in soln.x:
                    earth_model = []
                    for j in range(i + 1):
                        earth_model.append(LayerProps(Vp[j], _x[2*j + 1], rho[j], _x[2*j]))
                    # end for
                    earth_model = np.array(earth_model)
                    layer_base_vel = flux_comp.propagate_to_base(earth_model)
                    base_seismograms.append(layer_base_vel)
                # end for
                subsurface[layer_name] = base_seismograms
            # end if
        # end for
        soln.subsurface = subsurface
    else:
        soln.esu = soln.subsurface = None
    # end if

    return soln
# end func


def mcmc_solver_wrapper(model, obj_fn, mantle, Vp, rho, flux_window):
    """
    Wrapper callable for passing to MCMC solver in scipy style, which unpacks inputs into
    vector variables for solver.

    :param model: Per-layer model values (the vector being solved for) as flat array of
        (H, Vs) value pairs ordered by layer.
    :type model: numpy.array
    :param obj_fn: Callable to WfContinuationSuFluxComputer to compute SU flux.
    :param mantle: Mantle properties stored in class LayerProps instance.
    :type mantle: seismic.model_properties.LayerProps
    :param Vp: Array of Vp values ordered by layer.
    :type Vp: numpy.array
    :param rho: Array of rho values ordered by layer.
    :type rho: numpy.array
    :param flux_window: Pair of floats indicating the time window over which to perform SU flux integration
    :type flux_window: (float, float)
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

def save_mcmc_solution(soln_configs, input_file, output_file, job_timestamp, job_tracking, logger=None):
    """
    Save solution to HDF5 file. In general soln_configs will be an ensemble of per-station solutions.

    :param soln_configs: List of (solution, configuration) pairs to save (one per station).
    :type soln_configs: list((solution, dict))
    :param input_file: Name of input file used for the job. Saved to job node for traceability
    :type input_file: str
    :param output_file: Name of output file to create
    :type output_file: str or pathlib.Path
    :param job_timestamp: Job timestamp that will be used to generate the top level job group
    :type job_timestamp: str
    :param job_tracking: Dict containing job identification information for traceability
    :type job_tracking: dict or dict-like
    :param logger: [OPTIONAL] Log message destination
    :type logger: logging.Logger
    :return: None
    """
    assert isinstance(soln_configs, list)
    assert isinstance(job_timestamp, str)
    # TODO: migrate this to member of a new class for encapsulating an MCMC solution

    def write_data_empty(dataset):
        """
        Write dataset that might be empty. If empty, return h5py.Empty('f').
        See also function `read_data_empty`.
        """
        return h5py.Empty('f') if dataset is None else dataset
    # end func

    def write_list_dataset(dest_node, list_of_data):
        """Write list to a datase node (dict), storing item sequence index so that it can
        be later restored back to an ordered list. See also function `read_list_dataset`.
        """
        for idx, ds in enumerate(list_of_data):
            dest_node[str(idx)] = ds
        # end for
    # end func

    # Version string for layout of data in a node containing MCMC solution
    FORMAT_VERSION = '0.6'

    # Convert timestamp to valid Python identifier
    job_timestamp = 'T' + job_timestamp.replace('-', '_').replace(' ', '__').replace(':', '').replace('.', '_')

    with h5py.File(output_file, 'a') as h5f:
        job_root = h5f.create_group(job_timestamp)
        job_root.attrs['input_file'] = input_file
        job_root.attrs['job_tracking'] = json.dumps(job_tracking)
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
                    if soln is not None and soln.get('message'):
                        logger.error('Station {} reported error: {}'.format(station_id, soln.message))
                    # end if
                continue
            # end if
            station_id = station_id.replace('.', '_')
            station_node = job_root.create_group(station_id)
            station_node.attrs['config'] = json.dumps(config)
            station_node.attrs['format_version'] = FORMAT_VERSION
            try:
                station_node['x'] = soln.x
                station_node['num_input_seismograms'] = soln.num_input_seismograms
                assert len(soln.x) == len(soln.clusters)
                # Each cluster may be a different size, so we can't dump directly into h5py (which
                # requires hyper-rectangular array shape).
                cluster_node = station_node.create_group('clusters')
                write_list_dataset(cluster_node, soln.clusters)
                cluster_energy_node = station_node.create_group('cluster_energy')
                write_list_dataset(cluster_energy_node, soln.cluster_funvals)
                per_event_energy_node = station_node.create_group('per_event_energy')
                write_list_dataset(per_event_energy_node, soln.esu)
                # Subsurface seismograms, e.g. at bottom of sedimentary layer
                subsurface_node = station_node.create_group('subsurface')
                for layer_name, layer_seismograms in soln.subsurface.items():
                    layer_node = subsurface_node.create_group(layer_name)
                    write_list_dataset(layer_node, layer_seismograms)
                # end for
                station_node['bins'] = soln.bins
                station_node['distribution'] = soln.distribution
                station_node['acceptance_rate'] = soln.acceptance_rate
                station_node['success'] = soln.success
                station_node['status'] = soln.status
                station_node['message'] = soln.message
                station_node['fun'] = soln.fun
                station_node['jac'] = write_data_empty(soln.jac)
                station_node['nfev'] = soln.nfev
                station_node['njev'] = soln.njev
                station_node['nit'] = soln.nit
                station_node['maxcv'] = write_data_empty(soln.maxcv)
                station_node['samples'] = write_data_empty(soln.samples)
                station_node['sample_energies'] = write_data_empty(soln.sample_funvals)
                station_node['bounds'] = np.array([soln.bounds.lb, soln.bounds.ub])
                station_node['version'] = soln.version
                station_node['rnd_seed'] = soln.rnd_seed
            except TypeError as exc:
                if logger:
                    logger.error('Error saving station {} solution'.format(station_id))
                    logger.error(repr(exc))
                # end  if
            # end try
        # end for
    # end with

# end func


def load_mcmc_solution(h5_file, job_timestamp=None, logger=None):
    """Load Monte Carlo Markov Chain solution from HDF5 file.

    :param h5_file: File from which to load solution
    :type h5_file: str or pathlib.Path
    :param job_timestamp: Timestamp of job whose solution is to be loaded
    :type job_timestamp: str or NoneType
    :param logger: Output logging instance
    :type logger: logging.Logger
    :return: (solution, job configuration), job timestamp
    :rtype: (solution, dict), str
    """
    assert isinstance(job_timestamp, (str, type(None)))
    # TODO: migrate this to member of a new class for encapsulating an MCMC solution

    def read_data_empty(dataset):
        """
        Read dataset that might be empty. If empty, return None.
        See also function `write_data_empty`.

        :param dataset: The h5py.Dataset node to read.
        :type dataset: h5py.Dataset
        :return: Dataset value or None
        :rtype: numpy.array or NoneType
        """
        if not dataset.shape:
            value = None
        else:
            value = dataset.value
        # end if
        return value
    # end func

    def read_list_dataset(source_node):
        """Read list from a datase node containing ordered collection of items.
        See also function `write_list_dataset`.
        """
        list_data = []
        for idx, ds in source_node.items():
            list_data.append((int(idx), ds.value))
        # end for
        # Sort clusters by idx, then throw away the idx values.
        list_data.sort(key=lambda i: i[0])
        return [d[1] for d in list_data]
    # end func

    soln_configs = []
    with h5py.File(h5_file, 'r') as h5f:
        while job_timestamp is None:
            timestamps = list(h5f.keys())
            if len(timestamps) > 1:
                for i, ts in enumerate(timestamps):
                    job_node = h5f[ts]
                    job_tracking = json.loads(job_node.attrs['job_tracking']) \
                        if 'job_tracking' in job_node.attrs else ''
                    if job_tracking:
                        job_tracking = '(' + ', '.join([': '.join([k, str(v)]) for k, v in job_tracking.items()]) + ')'
                    # end if
                    print('[{}]'.format(i), ts, job_tracking)
                # end for
                index = input('Choose dataset number to load: ')
                if index.isdigit() and (0 <= int(index) < len(timestamps)):
                    index = int(index)
                # end if
            else:
                index = 0
            # end if
            job_timestamp = timestamps[index] if isinstance(index, int) else None
        # end while

        job_root = h5f[job_timestamp]
        # source_data_file = job_root.attrs['input_file']
        for station_id, station_node in job_root.items():
            if logger:
                logger.info('Loading {}'.format(station_id.replace('_', '.')))
            # end if
            job_config = json.loads(station_node.attrs['config'])
            format_version = station_node.attrs['format_version']
            job_config.update({'format_version': format_version})
            if logger:
                logger.info('H5 storage format version: {}'.format(format_version))
            # end if

            try:
                soln = optimize.OptimizeResult()
                soln.x = station_node['x'].value
                soln.num_input_seismograms = station_node['num_input_seismograms'].value

                cluster_node = station_node['clusters']
                soln.clusters = read_list_dataset(cluster_node)

                assert len(soln.x) == len(soln.clusters)

                cluster_energy_node = station_node['cluster_energy']
                soln.cluster_funvals = read_list_dataset(cluster_energy_node)

                per_event_energy_node = station_node['per_event_energy']
                soln.esu = read_list_dataset(per_event_energy_node)

                # Subsurface seismograms
                subsurface_node = station_node['subsurface']
                subsurface = {}
                for layer_name, layer_node in subsurface_node.items():
                    subsurface[layer_name] = read_list_dataset(layer_node)
                # end for
                soln.subsurface = subsurface

                soln.bins = station_node['bins'].value
                soln.distribution = station_node['distribution'].value
                soln.acceptance_rate = station_node['acceptance_rate'].value
                soln.success = bool(station_node['success'].value)
                soln.status = int(station_node['status'].value)
                soln.message = station_node['message'].value
                soln.fun = station_node['fun'].value
                soln.jac = read_data_empty(station_node['jac'])
                soln.nfev = int(station_node['nfev'].value)
                soln.njev = int(station_node['njev'].value)
                soln.nit = int(station_node['nit'].value)
                soln.maxcv = read_data_empty(station_node['maxcv'])
                soln.samples = read_data_empty(station_node['samples'])
                soln.sample_funvals = read_data_empty(station_node['sample_energies'])
                bounds = station_node['bounds'].value
                soln.bounds = optimize.Bounds(bounds[0], bounds[1])
                soln.version = station_node['version'].value
                if 'rnd_seed' in station_node:
                    soln.rnd_seed = int(station_node['rnd_seed'].value)
                else:
                    soln.rnd_seed = None
                # end if

                soln_configs.append((soln, job_config))
            except TypeError as exc:
                if logger:
                    logger.error('Error loading station {} solution'.format(station_id))
                    logger.error(repr(exc))
            # end try
        # end for
    # end with

    return soln_configs, job_timestamp
# end func


def run_station(config_file, waveform_file, network, station, location, logger):
    """Runner for analysis of single station. For multiple stations, set up config file to run batch
    job using mpi_job CLI.

    The output file is in HDF5 format. The configuration details are added to the output file for traceability.

    :param config_file: Config filename specifying job settings
    :type config_file: dict
    :param waveform_file: Event waveform source file for seismograms, generated using `extract_event_traces.py` script
    :type waveform_file: str or pathlib.Path
    :param network: Network code of station to analyse
    :type network: str
    :param station: Station code to analyse
    :type station: str
    :param location: Location code of station to analyse. Can be '' (empty string) if not set.
    :type location: str
    :param logger: Output logging instance
    :type logger: logging.Logger
    :return: Pair containing (solution, configuration) containers. Configuration will have additional traceability
        information.
    :rtype: (solution, dict)
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

    try:
        # Ordering of seismograms important here, since storage of sequential values in solution
        # depend on it. Here the input seismograms are ordered by event ID.
        soln = runner(waveform_data.station(station).values(), config, logger)
    except Exception as e:
        logger.error('Runner failed on station {}'.format(station_id))
        logger.exception(e)
        soln = optimize.OptimizeResult()
        soln.success = False
        soln.message = str(e)
    # end try

    # Add ordered event IDs so source waveforms can be re-extraced later
    # from source file if necessary.
    try:
        ordered_event_ids = [st[0].stats.event_id for st in waveform_data.station(station).values()]
    except Exception as e:
        logger.error('Event ID collection failed on station {}'.format(station_id))
        logger.exception(e)
        ordered_event_ids = []
    # end try
    config.update({"event_ids": ordered_event_ids})

    return soln, config

# end func


@click.group()
def main():
    """Main dispatch function to click sub-commands.
    """
# end func


@main.command(name='single-job')
@click.argument('config_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Event waveform source file for seismograms, generated using extract_event_traces.py script')
@click.option('--network', type=str, required=True, help='Network code for the station to load')
@click.option('--station', type=str, required=True, help='Station code whose data is to be loaded')
@click.option('--location', type=str, default='', show_default=True, help='Location code for the station')
@click.option('--output-file', type=click.Path(dir_okay=False),
              help='Name of the output file in which solutions should be saved. May be an existing file, '
                   'as multiple solutions may be saved in a single file as they are indexed by job timestamp.')
def station_job(config_file, waveform_file, network, station, location='', output_file=None):
    """
    CLI dispatch function for single station. See help strings for option documentation.

    Example usage::

        python runners.py single-job example_config.json \
        --waveform-file /g/data/ha3/am7399/shared/OA_RF_analysis/OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5 \
        --network OA --station CD23 --location 0M --output-file out_test.h5

    :param config_file: JSON file containing job configuration parameters.
    :type config_file: str or pathlib.Path
    :param waveform_file: HDF5 file containing event waveforms, generated using `extract_event_traces.py` script
    :type waveform_file: str or pathlib.Path
    :param network: Network code of job to run
    :type network: str
    :param station: Station code of job to run
    :type station: str
    :param location: Location code of job to run
    :type location: str
    :param output_file: Name of the output file in which solutions should be saved
    :type output_file: str or pathlib.Path
    :return: Integer status code
    :rtype: int
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
    job_tracking = {'job_name': '.'.join([network, station, location])}
    job_id = os.getenv('PBS_JOBID')
    if job_id is not None:
        job_tracking.update({"job_id": job_id})
    # end if
    save_mcmc_solution([soln_config], waveform_file, output_file, job_timestamp, job_tracking, logger=logger)

    return 0
# end func


@main.command(name='batch-job')
@click.argument('config_file', type=click.File('r'), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Event waveform source file for seismograms, generated using extract_event_traces.py script')
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output file in which solutions should be saved. May be an existing file, '
                   'as multiple solutions may be saved in a single file as they are indexed by job timestamp.')
def mpi_job(config_file, waveform_file, output_file):
    """
    CLI dispatch function for MPI run over batch of stations. See help strings for option documentation.

    Example MPI usage::

        mpiexec -n 8 python runners.py batch-job example_batch.json \
        --waveform-file /g/data/ha3/am7399/shared/OA_RF_analysis/OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5 \
        --output-file test_batch_output.h5

    :param config_file: JSON file containing batch configuration parameters.
    :type config_file: str or pathlib.Path
    :param waveform_file: HDF5 file containing event waveforms, generated using `extract_event_traces.py` script
    :type waveform_file: str or pathlib.Path
    :param output_file: Name of the output file in which solutions should be saved
    :type output_file: str or pathlib.Path
    :return: Integer status code
    :rtype: int
    """

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    # nproc = comm.Get_size()  # TODO: Use this to more intelligently split the jobs across available processors
    rank = comm.Get_rank()

    job_timestamp = comm.bcast(str(datetime.now()), root=0)

    logger = logging.getLogger(__name__ + str(rank))
    logger.setLevel(logging.INFO)
    log_fmt = logging.Formatter(**LOG_FORMAT)
    job_id = os.getenv('PBS_JOBID')
    if job_id is None:
        job_id = str(comm.bcast(os.getpid(), root=0))
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
        job_tracking = {}
        # Label with job name from job queueing if available and no name given
        if "job_name" in batch_config:
            job_tracking.update({"job_name": batch_config["job_name"]})
        else:
            job_name = os.getenv('PBS_JOBNAME')
            if job_name is not None:
                job_tracking.update({"job_name": job_name})
            # end if
        # end if
        job_id = os.getenv('PBS_JOBID')
        if job_id is not None:
            job_tracking.update({"job_id": job_id})
        # end if

        save_mcmc_solution(soln_config, waveform_file, output_file, job_timestamp, job_tracking, logger=logger)
    # end if

    return 0

# end func


if __name__ == '__main__':
    main()
# end if
