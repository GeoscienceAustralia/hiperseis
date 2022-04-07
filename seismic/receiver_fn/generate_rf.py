#!/usr/bin/env python
"""Generate Receiver Functions (RF) from collection of 3-channel seismic traces.
"""
import copy

from mpi4py import MPI
import logging
import json

import numpy as np
import click
import os
import tqdm.auto as tqdm
from seismic.receiver_fn.rf_corrections import Corrections
from seismic.receiver_fn import rf_util
from seismic.receiver_fn.generate_rf_helper import transform_stream_to_rf
from seismic.rf_station_orientations import analyze_station_orientations
from seismic.network_event_dataset import NetworkEventDataset
from seismic.stream_processing import zne_order, negate_channel, swap_ne_channels, correct_back_azimuth
from seismic.stream_io import remove_group, get_obspyh5_index
from rf import RFStream
from collections import defaultdict

# pylint: disable=invalid-name, logging-format-interpolationa

logging.basicConfig()

def event_waveforms_to_rf(input_file, output_file, config, network_list='*', station_list='*', only_corrections=False):
    """
    Main entry point for generating RFs from event traces.

    Config file consists of 3 sub-dictionaries. One named "filtering" for
    input stream filtering settings, one named "processing" for RF processing
    settings, and one named "correction" for rotating/swapping/negating channel
    data for one or more named stations with potential orientation discrepancies.
    Each of these sub-dicts is described below::

        "filtering":  # Filtering settings
        {
          "resample_rate": float # Resampling rate in Hz
          "taper_limit": float   # Fraction of signal to taper at end, between 0 and 0.5
          "filter_band": (float, float) # Filter pass band (Hz). Not required for freq-domain deconvolution.
          "channel_pattern": # Ordered list of preferred channels, e.g. 'HH*,BH*',
                             # where channel selection is ambiguous.
          "baz_range": (float, float) or [(float, float), ...] # Discrete ranges of source back azimuth to use (degrees).
              # Each value must be between 0 and 360. May be a pair or a list of pairs for multiple ranges.
        }

        "processing":  # RF processing settings
        {
          "custom_preproc":
          {
            "import": 'import custom symbols',  # statement to import required symbols
            "func": 'preproc functor'  # expression to get handle to custom preprocessing functor
            "args": {}  # additional kwargs to pass to func
          }
          "trim_start_time": float # Trace trim start time in sec, relative to onset
          "trim_end_time": float # Trace trim end time in sec, relative to onset
          "rotation_type": str # Choice of ['zrt', 'lqt']. Rotational coordinate system
                               # for aligning ZNE trace components with incident wave direction
          "deconv_domain": str # Choice of ['time', 'freq', 'iter']. Whether to perform deconvolution
                               # in time or freq domain, or iterative technique
          "gauss_width": float # Gaussian freq domain filter width. Only required for freq-domain deconvolution
          "water_level": float # Water-level for freq domain spectrum. Only required for freq-domain deconvolution
          "spiking": float # Spiking factor (noise suppression), only required for time-domain deconvolution
          "normalize": bool # Whether to normalize RF amplitude
        }

        "correction": # corrections to be applied to data for named stations prior to RF computation
        {
          "plot_dir": str # path to folder where plots related to orientation corrections are to be saved
          "swap_ne": list # list of NET.STA.LOC for which N and E channels are to be swapped, e.g ["OA.BL27."],
          "rotate": list # list of NET.STA.LOC that are to be rotated to maximize P-arrival energy on \
                           the primary RF component, e.g ["OA.BL27."]
          "negate": list # list of NET.STA.LOC.CHA that are to be negated, e.g ["OA.BL27..HHZ"]
        }

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

    config_filtering = config.setdefault("filtering", {})
    channel_pattern = config_filtering.get("channel_pattern")
    config_processing = config.setdefault("processing", {})
    config_correction = config.setdefault("correction", {})

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if channel_pattern is not None:
        channel_pattern = channel_pattern.strip()
        logger.info("Using channel matching pattern {}".format(channel_pattern))
    # end if

    if(rank == 0):
        logger.info("Processing source file {}".format(input_file))

        # retrieve all available hdf_keys
        proc_hdfkeys = get_obspyh5_index(input_file, seeds_only=True)

        # trim stations to be processed based on the user-provided network- and station-list
        proc_hdfkeys = rf_util.trim_hdf_keys(proc_hdfkeys, network_list, station_list)

        if(only_corrections): # trim the hdf_keys if processing only corrections
            corrections = Corrections(config_correction, proc_hdfkeys)
            proc_hdfkeys = [item for item in proc_hdfkeys if corrections.needsCorrections(item)]
        # end if
    # end if

    # broadcast workload to all procs
    proc_hdfkeys = comm.bcast(proc_hdfkeys, root=0)

    # load corrections-config
    corrections = Corrections(config_correction, proc_hdfkeys)

    # split stations over all ranks
    proc_hdfkeys = rf_util.split_list(proc_hdfkeys, nproc)

    pbar = tqdm.tqdm(total=len(proc_hdfkeys[rank]))
    proc_rf_stream = RFStream()
    baz_corrections = defaultdict(dict)
    for proc_hdfkey in proc_hdfkeys[rank]:
        nsl = proc_hdfkey  # network-station-location
        pbar.set_description("Rank {}: {}".format(rank, nsl))

        net, sta, loc = nsl.split('.')
        # note that ned contains a single station
        ned = NetworkEventDataset(input_file, network=net, station=sta, location=loc)

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

            # channel rotations through baz correction
            if(corrections.needsRotation(proc_hdfkey)):
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

                    rf_3ch = transform_stream_to_rf(evid, stream, config_filtering,
                                                    config_processing)
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

    # gather and output baz corrections
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
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
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
    if config_file is not None:
        with open(config_file, 'r') as cf:
            config = json.load(cf)
        # end with

        config_correction = config.setdefault("correction", {})
        if(not len(config_correction) and only_corrections):
            assert 0, 'A correction block is required in the config file for --only-corrections'
        # end if
    else:
        config = {}  # all default settings
    # end if

    # Dispatch call to worker function. See worker function for documentation.
    event_waveforms_to_rf(input_file, output_file, config, network_list=network_list, station_list=station_list,
                          only_corrections=only_corrections)
# end main

if __name__ == "__main__":
    _main()  # pylint: disable=no-value-for-parameter
# end if
