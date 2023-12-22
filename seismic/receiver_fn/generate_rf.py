#!/usr/bin/env python
"""Generate Receiver Functions (RF) from collection of 3-channel seismic traces.
"""
import copy

from mpi4py import MPI
import logging

import numpy as np
import click
import os, json
import tqdm.auto as tqdm
from seismic.receiver_fn import rf_util
from seismic.receiver_fn.rf_config import RFConfig, Corrections
from seismic.receiver_fn.generate_rf_helper import transform_stream_to_rf
from seismic.rf_station_orientations import analyze_station_orientations
from seismic.network_event_dataset import NetworkEventDataset
from seismic.stream_processing import zne_order, negate_channel, \
    swap_ne_channels, correct_back_azimuth, recompute_inclinations
from seismic.stream_io import remove_group, get_obspyh5_index
from rf import RFStream
from collections import defaultdict
from seismic.misc import split_list

# pylint: disable=invalid-name, logging-format-interpolationa

logging.basicConfig()

def event_waveforms_to_rf(input_file: str, output_file: str, config: RFConfig,
                          network_list='*', station_list='*', only_corrections=False):
    """
    Main entry point for generating RFs from event traces.

    :param input_file: Event waveform source file for seismograms, generated using extract_event_traces.py script
    :type input_file: str or pathlib.Path
    :param output_file: Name of hdf5 file to produce containing RFs
    :type output_file: str or pathlib.Path
    :param config: Dictionary containing job configuration parameters
    :type config: dict
    :return: None
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_hdfkeys = None
    corrections = None

    config_filtering = config.config_filtering
    config_processing = config.config_processing
    config_correction = config.config_correction

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # set h5 root to appropriate waveform path
    h5_root = None
    rf_type = config.config_processing['rf_type']
    if(rf_type == 'prf'):
        h5_root = 'waveforms' # default is P for backward compatibility
    else:
        h5_root = 'waveforms/S'
    # end if

    if(rank == 0):
        logger.info("Processing source file {}".format(input_file))

        # retrieve all available hdf_keys
        proc_hdfkeys = None
        try:
            proc_hdfkeys = get_obspyh5_index(input_file, seeds_only=True, root=h5_root)
        except Exception as e:
            print('Failed to read {} with error {}. '
                  'Ensure file has data at key "{}". Aborting..'.format(input_file,
                                                                        e, h5_root))
        # ene try

        # trim stations to be processed based on the user-provided network- and station-list
        proc_hdfkeys = rf_util.trim_hdf_keys(proc_hdfkeys, network_list, station_list)

        if(only_corrections): # trim the hdf_keys if processing only corrections
            corrections = Corrections(config, proc_hdfkeys)
            proc_hdfkeys = [item for item in proc_hdfkeys if corrections.needsCorrections(item)]
        # end if
    # end if

    # broadcast workload to all procs
    proc_hdfkeys = comm.bcast(proc_hdfkeys, root=0)

    # load corrections-config
    corrections = Corrections(config, proc_hdfkeys)

    # split stations over all ranks
    proc_hdfkeys = split_list(proc_hdfkeys, nproc)

    pbar = tqdm.tqdm(total=len(proc_hdfkeys[rank]))
    proc_rf_stream = RFStream()
    baz_corrections = defaultdict(dict)
    for proc_hdfkey in proc_hdfkeys[rank]:
        nsl = proc_hdfkey  # network-station-location
        pbar.set_description("Rank {}: {}".format(rank, nsl))

        net, sta, loc = nsl.split('.')
        # note that ned contains a single station
        ned = NetworkEventDataset(input_file, network=net, station=sta, location=loc, root=h5_root)

        # corrections
        if (corrections.needsCorrections(proc_hdfkey)):
            # channel negation
            ch_list = corrections.needsNegation(proc_hdfkey)
            if(ch_list):
                for ch in ch_list:
                    logger.info('Rank {}: {}: Negating component {}'.format(rank, nsl, ch))
                    ned.apply(lambda st: negate_channel(None, st, ch))
                # end for
            # end if

            # channel swaps
            if(corrections.needsChannelSwap(proc_hdfkey)):
                logger.info('Rank {}: {}: Applying a NE channel swap '.format(rank, nsl))
                ned.apply(lambda st: swap_ne_channels(None, st))
            # end if

            # channel rotations through baz correction. Note that this is only
            # applicable for P-RFs
            if (corrections.needsRotation(proc_hdfkey) and rf_type == 'prf'):
                # TODO : expose the following parameters
                bazcorr_curation_opts = {"min_slope_ratio": 5,
                                         "min_snr": 2.0,
                                         "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}}
                bazcorr_config_filtering = {"resample_rate": 4.0,
                                            "taper_limit": 0.05,
                                            "filter_band": [0.01, 0.5]}

                result = analyze_station_orientations(copy.deepcopy(ned),
                                                      curation_opts=bazcorr_curation_opts,
                                                      config_filtering=bazcorr_config_filtering,
                                                      save_plots_path=corrections.plot_dir)

                baz_corrections.update({nsl: result[list(result.keys())[0]]})

                try:
                    logger.info('Rank {}: {}: Applying a baz correction '
                                'of {}'.format(rank,
                                               nsl,
                                               result[nsl]['azimuth_correction']))
                    ned.apply(lambda st: correct_back_azimuth(None, st, result[nsl]['azimuth_correction']))
                except:
                    logger.warning('Channel rotation failed for {}. Moving along..'.format(nsl))
                # end try
            # end if

            # recompute inclinations. Note that this is only applicable for S-RFs
            if (corrections.needsInclinationRecomputed(proc_hdfkey) and rf_type == 'srf'):
                logger.info('Rank {}: {}: Recomputing inclinations'.format(rank, nsl))
                ned.apply(lambda st: recompute_inclinations(None, st))
            # end if
        # end if

        status_list = []
        for sta, db_evid in ned.by_station(): # note that ned contains a single station
            for evid, stream in db_evid.items():
                rf_3ch = None
                try:
                    stream = RFStream(stream.traces)
                    stream.traces = sorted(stream.traces, key=zne_order)
                    # Strongly assert expected ordering of traces. This must be respected so that
                    # RF normalization works properly.
                    assert stream.traces[0].stats.channel[-1] == 'Z'
                    assert stream.traces[1].stats.channel[-1] == 'N'
                    assert stream.traces[2].stats.channel[-1] == 'E'

                    rf_3ch = transform_stream_to_rf(evid, stream, config)
                except Exception as e:
                    print(str(e) + ' in station {}, event {}'.format(nsl, evid))
                    status_list.append(False)
                    continue
                # end try

                if rf_3ch is None:
                    status_list.append(False)
                else:
                    status_list.append(True)
                    proc_rf_stream += rf_3ch
                # end if
            # end for
        # end for
        pbar.update()

        num_tasks = len(status_list)
        num_success = np.sum(status_list)
        num_rejected = num_tasks - num_success
        logger.info("{}: {}/{} streams returned valid RF, {} streams rejected".format(nsl, num_success, num_tasks, num_rejected))
    # end for
    pbar.close()

    # gather and output baz corrections if any for P RFs
    if(rf_type == 'prf'):
        baz_corrections = comm.gather(baz_corrections, root=0)
        if(rank == 0):
            if(len(baz_corrections) and corrections.plot_dir):
                output_dict = {}
                for d in baz_corrections:
                    if(len(d)): output_dict.update(d)
                # end for
                json.dump(output_dict, open(os.path.join(corrections.plot_dir,
                                                         'azimuth_corrections.json'), 'w'))
            # end if
        # end if
    # end if

    # serialize writing of rf_streams to disk
    for irank in np.arange(nproc):
        if(irank == rank):
            if(len(proc_rf_stream)):
                for hdf_key in proc_hdfkeys[rank]:
                    # remove existing traces if there are any
                    try:
                        remove_group(output_file, hdf_key, logger=logger)
                    except Exception as e:
                        logger.warning(str(e))
                    # end try

                    logger.info("Writing RF stream(s) for {} on rank {}...".format(hdf_key, rank))
                # end for

                proc_rf_stream.write(output_file, format='H5', mode='a')
            # end if
        # end if
        comm.Barrier()
    # end for

    if(rank == 0):
        print("Finishing...")
        print("generate_rf SUCCESS!")
    # end if
# end func

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) or a text file '
                                                  'with station names in each row, w/wo location codes.', type=str,
              show_default=True)
@click.option('--config-file', type=click.Path(dir_okay=False), default=None,
              show_default=True, help="Run configuration file in JSON format")
@click.option('--only-corrections', is_flag=True, default=False, show_default=True,
              help="Compute and apply corrections for stations listed under 'correction' in the "
                   "input json config file -- all other stations are ignored. Note that "
                   "preexisting data (for relevant channels, if present) are deleted before "
                   "saving the corrections")
def _main(input_file, output_file, network_list, station_list, config_file, only_corrections):
    """
    INPUT_FILE : Input waveforms in H5 format\n
                 (output of extract_event_traces.py)\n
    OUTPUT_FILE : Output H5 file name
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    config = RFConfig(config_file)

    # Dispatch call to worker function. See worker function for documentation.
    event_waveforms_to_rf(input_file, output_file, config, network_list=network_list, station_list=station_list,
                          only_corrections=only_corrections)
# end main

if __name__ == "__main__":
    _main()  # pylint: disable=no-value-for-parameter
# end if
