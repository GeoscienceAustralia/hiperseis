#!/usr/bin/env python
"""Generate Receiver Functions (RF) from collection of 3-channel seismic traces.
"""

import logging
from multiprocessing import Process, Manager
import json

import numpy as np
import click

# from rf import IterMultipleComponents

try:
    from joblib import Parallel, delayed
    parallel_available = True
except ImportError:
    parallel_available = False
# end try

from seismic.receiver_fn.rf_process_io import async_write
from seismic.receiver_fn.rf_h5_file_event_iterator import IterRfH5FileEvents
from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.receiver_fn.rf_deconvolution import rf_iter_deconv
from seismic.stream_quality_filter import curate_stream3c
from seismic.stream_processing import back_azimuth_filter


logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

DEFAULT_RESAMPLE_RATE_HZ = 20.0
DEFAULT_FILTER_BAND_HZ = (0.02, 1.00)
DEFAULT_TAPER_LIMIT = 0.05
DEFAULT_TRIM_START_TIME_SEC = -50.0
DEFAULT_TRIM_END_TIME_SEC = 150.0
DEFAULT_ROTATION_TYPE = 'zrt'   # from ['zrt', 'lqt']
DEFAULT_DECONV_DOMAIN = 'time'  # from ['time', 'freq', 'iter']
DEFAULT_GAUSS_WIDTH = 1.0
DEFAULT_WATER_LEVEL = 0.01
DEFAULT_SPIKING = 0.5
RAW_RESAMPLE_RATE_HZ = 20.0
BANDPASS_FILTER_ORDER = 2


def transform_stream_to_rf(oqueue, ev_id, stream3c, config_filtering,
                           config_processing, **kwargs):
    """Generate P-phase receiver functions for a single 3-channel stream.
    See documentation for function event_waveforms_to_rf for details of
    config dictionary contents.

    :param oqueue: Output queue where filtered streams are queued
    :type oqueue: multiprocessing.Manager.Queue
    :param ev_id: The event id
    :type ev_id: int or str
    :param stream3c: Stream with 3 components of trace data
    :type stream3c: rf.RFStream
    :param config_filtering: Dictionary containing stream filtering settings
    :type config_filtering: dict
    :param config_processing: Dictionary containing RF processing settings
    :type config_processing: dict
    :param kwargs: Keyword arguments that will be passed to filtering and deconvolution functions.
    :type kwargs: dict
    :return: True if stream containing receiver function is pushed into output queue
        oqueue, False otherwise
    :rtype: bool
    """

    resample_rate_hz = config_filtering.get("resample_rate", DEFAULT_RESAMPLE_RATE_HZ)
    filter_band_hz = config_filtering.get("filter_band", DEFAULT_FILTER_BAND_HZ)
    assert resample_rate_hz >= 2.0*filter_band_hz[1], "Too low sample rate will alias signal"
    taper_limit = config_filtering.get("taper_limit", DEFAULT_TAPER_LIMIT)
    baz_range = config_filtering.get("baz_range")

    trim_start_time_sec = config_processing.get("trim_start_time", DEFAULT_TRIM_START_TIME_SEC)
    trim_end_time_sec = config_processing.get("trim_end_time", DEFAULT_TRIM_END_TIME_SEC)
    rotation_type = config_processing.get("rotation_type", DEFAULT_ROTATION_TYPE)
    deconv_domain = config_processing.get("deconv_domain", DEFAULT_DECONV_DOMAIN)

    logger = logging.getLogger(__name__)
    logger.info("Event #{}".format(ev_id))

    assert deconv_domain.lower() in ['time', 'freq', 'iter']

    if not curate_stream3c(ev_id, stream3c, logger):
        return False

    if baz_range is not None:
        if not isinstance(baz_range[0], list):
            baz_range = [baz_range]
        baz = stream3c[0].stats.back_azimuth
        if not any([back_azimuth_filter(baz, tuple(b)) for b in baz_range]):
            return False
    # end if

    # Compute SNR of prior z-component after some low pass filtering.
    # Apply conservative anti-aliasing filter before downsampling raw signal. Cutoff at half
    # the Nyquist freq to make sure almost no high freq energy leaking through the filter can
    # alias down into the frequency bands of interest.
    stream_z = stream3c.select(component='Z').copy().filter('lowpass', freq=RAW_RESAMPLE_RATE_HZ/4.0,
                                                            corners=2, zerophase=True)
    # Since cutoff freq is well below Nyquist, we use a lower Lanczos kernel size (default is a=20).
    stream_z = stream_z.interpolate(RAW_RESAMPLE_RATE_HZ, method='lanczos', a=10)
    # Trim original trace to time window
    stream_z.trim2(trim_start_time_sec, trim_end_time_sec, reftime='onset')
    stream_z.detrend('linear')
    stream_z.taper(taper_limit, **kwargs)
    compute_vertical_snr(stream_z)

    rotation_type = rotation_type.lower()
    assert rotation_type in ['zrt', 'lqt']
    if rotation_type == 'zrt':
        rf_rotation = 'NE->RT'
    else:
        rf_rotation = 'ZNE->LQT'
    # end if

    stream3c.detrend('linear')
    stream3c.taper(taper_limit, **kwargs)
    try:
        normalize = config_processing.get("normalize", True)
        if deconv_domain == 'time':
            # ZRT receiver functions must be specified
            stream3c.filter('bandpass', freqmin=filter_band_hz[0], freqmax=filter_band_hz[1],
                            corners=BANDPASS_FILTER_ORDER, zerophase=True, **kwargs).interpolate(resample_rate_hz)
            spiking = config_processing.get("spiking", DEFAULT_SPIKING)
            kwargs.update({'spiking': spiking})
            if not normalize:
                # No normalization. The "normalize" argument must be set to None.
                kwargs['normalize'] = None
            # end if
            stream3c.rf(rotate=rf_rotation, **kwargs)
        elif deconv_domain == 'freq':
            # As of https://github.com/trichter/rf/issues/15, the Gaussian parameter is directly related
            # to the cutoff freq. Requires rf version >=0.8.0
            if not normalize:
                # No normalization. The "normalize" argument must be set to None.
                kwargs['normalize'] = None
            # end if
            gauss_width = config_processing.get("gauss_width", DEFAULT_GAUSS_WIDTH)
            water_level = config_processing.get("water_level", DEFAULT_WATER_LEVEL)
            stream3c.rf(rotate=rf_rotation, deconvolve='freq', gauss=gauss_width, waterlevel=water_level, **kwargs)
            # Interpolate to requested sampling rate.
            stream3c.interpolate(resample_rate_hz)
        elif deconv_domain == 'iter':
            # For iterative deconvolution, we need to trim first before deconvolution, as the fitting to the
            # response performs much better on shorter (trimmed) data.
            stream3c.filter('bandpass', freqmin=filter_band_hz[0], freqmax=filter_band_hz[1],
                            corners=BANDPASS_FILTER_ORDER, zerophase=True, **kwargs).interpolate(resample_rate_hz)
            if not normalize:
                # No normalization. The "normalize" argument must be set to None.
                normalize = None
            else:
                normalize = 0  # Use Z-component for normalization
            # end if
            stream3c.rf(rotate=rf_rotation, trim=(trim_start_time_sec, trim_end_time_sec), deconvolve='func',
                        func=rf_iter_deconv, normalize=normalize, min_fit_threshold=75.0)
        else:
            assert False, "Not yet supported deconvolution technique '{}'".format(deconv_domain)
        # end if
    except (IndexError, ValueError) as e:
        logger.error("Failed on stream {}:\n{}\nwith error:\n{}".format(ev_id, stream3c, str(e)))
        return False
    # end try

    # Check for any empty channels after deconvolution.
    for tr in stream3c:
        if len(tr.data) == 0:
            logger.warning("No data left in channel {} of stream {} after deconv (skipping)".format(
                           tr.stats.channel, ev_id))
            return False
        # end if
    # end for

    # Perform trimming before computing max amplitude.
    if deconv_domain != 'iter':
        stream3c.trim2(trim_start_time_sec, trim_end_time_sec, reftime='onset')
    # end if

    if len(stream3c) != 3:
        logger.warning("Unexpected number of channels in stream {} after trim (skipping):\n{}"
                       .format(ev_id, stream3c))
        return False
    # end if

    assert len(stream_z) == 1, "Expected only Z channel for a single event in stream_z: {}".format(stream_z)
    for tr in stream3c:
        metadata = {
            'amp_max': np.amax(np.abs(tr.data)),  # Max absolute amplitude
            'amp_rms': np.sqrt(np.mean(np.square(tr.data))),  # RMS amplitude
            'event_id': ev_id,
            'rotation': rotation_type,
            'snr_prior': stream_z[0].stats.snr_prior,
            'z_amp_max': np.max(np.abs(stream_z[0].data)),  # Max amplitude on resampled original z-component
            'z_amp_rms': np.sqrt(np.mean(np.square(stream_z[0].data)))  # RMS thereof
        }
        tr.stats.update(metadata)
    # end for

    output_stream = stream3c
    oqueue.put(output_stream)

    return True
# end func


def event_waveforms_to_rf(input_file, output_file, config):
    """
    Main entry point for generating RFs from event traces.

    Config file consists of 3 sub-dictionaries. One named "filtering" for
    input stream filtering settings, one named "processing" for RF processing
    settings, and one named "system" for options on how the system will run the
    job. Each of these sub-dicts is described below.

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

    "system":  # job run settings
    {
      "parallel": bool # Use parallel execution
      "memmap": bool # Memmap input file for improved performance in data reading thread.
                     # Useful when data input is bottleneck, if system memory permits.
      "temp_dir": str or path # Temporary directory to use for best performance
      "aggressive_dispatch": bool # Dispatch all worker jobs as aggressively as possible to minimize
                                  # chance of worker being starved of work. Uses more memory.
    }

    :param input_file: Event waveform source file for seismograms, generated using extract_event_traces.py script
    :type input_file: str or pathlib.Path
    :param output_file: Name of hdf5 file to produce containing RFs
    :type output_file: str or pathlib.Path
    :param config: Dictionary containing job configuration parameters
    :type config: dict
    :return: None
    """

    config_filtering = config.setdefault("filtering", {})
    channel_pattern = config_filtering.get("channel_pattern")

    config_processing = config.setdefault("processing", {})

    config_system = config.setdefault("system", {})
    parallel = config_system.get("parallel", True)
    aggressive_dispatch = config_system.get("aggressive_dispatch", False)
    memmap = config_system.get("memmap", False)
    temp_dir = config_system.get("temp_dir")

    dispatch_policy = '5*n_jobs'
    if parallel:
        assert parallel_available, "Cannot run parallel as joblib import failed"
        if aggressive_dispatch:
            dispatch_policy = 'all'
        # end if
    # end if

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if channel_pattern is not None:
        channel_pattern = channel_pattern.strip()
        logger.info("Using channel matching pattern {}".format(channel_pattern))
    # end if

    # Set up asynchronous buffered writing of results to file
    mgr = Manager()
    write_queue = mgr.Queue()
    config_str = json.dumps(config)
    output_thread = Process(target=async_write, args=(write_queue, output_file, 500, config_str))
    output_thread.daemon = True
    output_thread.start()

    logger.info("Processing source file {}".format(input_file))
    if parallel:
        # Process in parallel
        logger.info("Parallel processing")
        # n_jobs is -3 to allow one dedicated processor for running main thread and one for running output thread
        status = Parallel(n_jobs=-3, verbose=5, max_nbytes='16M', temp_folder=temp_dir, pre_dispatch=dispatch_policy)\
            (delayed(transform_stream_to_rf)(write_queue, id, stream3c, config_filtering, config_processing)
             for _, id, _, stream3c in IterRfH5FileEvents(input_file, memmap, channel_pattern))
    else:
        # Process in serial
        logger.info("Serial processing")
        status = list((transform_stream_to_rf(write_queue, id, stream3c, config_filtering, config_processing)
                       for _, id, _, stream3c in IterRfH5FileEvents(input_file, memmap, channel_pattern)))
    # end if
    num_tasks = len(status)
    num_success = np.sum(status)
    num_rejected = num_tasks - num_success
    logger.info("{}/{} streams returned valid RF, {} streams rejected".format(num_success, num_tasks, num_rejected))

    # Signal completion
    logger.info("Finishing...")
    write_queue.put(None)
    write_queue.join()

    logger.info("generate_rf SUCCESS!")
# end func


@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--config-file', type=click.Path(dir_okay=False),
              help="Job configuration file in JSON format")
def _main(input_file, output_file, config_file=None):
    if config_file is not None:
        with open(config_file, 'r') as cf:
            config = json.load(cf)
        # end with
    else:
        config = {}  # all default settings
    # end if
    # Dispatch call to worker function. See worker function for documentation.
    event_waveforms_to_rf(input_file, output_file, config)
# end main


if __name__ == "__main__":
    _main()  # pylint: disable=no-value-for-parameter
# end if
