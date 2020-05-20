#!/usr/bin/env python
"""Generate Receiver Functions (RF) from collection of 3-channel seismic traces.
"""

import logging
from multiprocessing import Process, Manager

import numpy as np
import click

# from rf import IterMultipleComponents

try:
    from joblib import Parallel, delayed
    parallel_available = True
except ImportError:
    parallel_available = False

from seismic.receiver_fn.rf_process_io import async_write
from seismic.receiver_fn.rf_h5_file_event_iterator import IterRfH5FileEvents
from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.receiver_fn.rf_deconvolution import rf_iter_deconv
from seismic.stream_quality_filter import curate_stream3c


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


def transform_stream_to_rf(oqueue, ev_id, stream3c, resample_rate_hz, taper_limit, rotation_type, filter_band_hz,
                           gauss_width, water_level, trim_start_time_sec, trim_end_time_sec,
                           deconv_domain, spiking, normalize, **kwargs):
    """Generate P-phase receiver functions for a single 3-channel stream.

    :param oqueue: Output queue where filtered streams are queued
    :type oqueue: multiprocessing.Manager.Queue
    :param ev_id: The event id
    :type ev_id: int or str
    :param stream3c: Stream with 3 components of trace data
    :type stream3c: rf.RFStream
    :param resample_rate_hz: Resampling rate of generated RF
    :type resample_rate_hz: float
    :param taper_limit: Taper limit passed to obspy.core.stream.Stream.taper
    :type taper_limit: float
    :param rotation_type: Rotation type, should be one of ['ZRT', 'LQT']
    :type rotation_type: str
    :param filter_band_hz: Limits of passband for frequency filtering of raw traces.
    :type filter_band_hz: tuple(float, float)
    :param gauss_width: Gauss width for frequency-domain deconvolution
    :type gauss_width: float
    :param water_level: Water-level for frequency-domain deconvolution
    :type water_level: float
    :param trim_start_time_sec: Trim start time of RF relative to theoretical onset
    :type trim_start_time_sec: float
    :param trim_end_time_sec: Trim end time of RF relative to theoretical onset
    :type trim_end_time_sec: float
    :param deconv_domain: Domain in which to perform deconvolution
    :type deconv_domain: str
    :param spiking: Spiking value to use for time-domain deconvolution
    :type spiking: float
    :param normalize: Whether to normalize RFs
    :type normalize: bool
    :param kwargs: Keyword arguments that will be passed to filtering and deconvolution functions.
    :type kwargs: dict
    :return: Stream containing receiver function
    :rtype: rf.RFStream
    """

    logger = logging.getLogger(__name__)
    logger.info("Event #{}".format(ev_id))

    assert deconv_domain.lower() in ['time', 'freq', 'iter']

    if not curate_stream3c(ev_id, stream3c, logger):
        return False

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
        if deconv_domain == 'time':
            # ZRT receiver functions must be specified
            stream3c.filter('bandpass', freqmin=filter_band_hz[0], freqmax=filter_band_hz[1],
                            corners=BANDPASS_FILTER_ORDER, zerophase=True, **kwargs).interpolate(resample_rate_hz)
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


def event_waveforms_to_rf(input_file, output_file, resample_rate, taper_limit, filter_band,
                          trim_start_time, trim_end_time, rotation_type, deconv_domain,
                          gauss_width=None, water_level=None, spiking=None,
                          normalize=True, parallel=True, memmap=False, temp_dir=None,
                          aggressive_dispatch=False, channel_pattern=None):
    """
    Main entry point for generating RFs from event traces.

    FILL IN MISSING DOCS HERE

    :param input_file: Event waveform source file for seismograms, generated using extract_event_traces.py script
    :type input_file: str or pathlib.Path
    :param output_file: Name of hdf5 file to produce containing RFs
    :type output_file: str or pathlib.Path
    :return: None
    """
    assert resample_rate >= 2.0*filter_band[1], "Too low sample rate will alias signal"

    dispatch_policy = '5*n_jobs'
    if parallel:
        assert parallel_available, "Cannot run parallel as joblib import failed"
        if aggressive_dispatch:
            dispatch_policy = 'all'

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if channel_pattern is not None:
        channel_pattern = channel_pattern.strip()
        logger.info("Using channel matching pattern {}".format(channel_pattern))

    # Set up asynchronous buffered writing of results to file
    mgr = Manager()
    write_queue = mgr.Queue()
    output_thread = Process(target=async_write, args=(write_queue, output_file, 500))
    output_thread.daemon = True
    output_thread.start()

    logger.info("Processing source file {}".format(input_file))
    if parallel:
        # Process in parallel
        logger.info("Parallel processing")
        # n_jobs is -3 to allow one dedicated processor for running main thread and one for running output thread
        status = Parallel(n_jobs=-3, verbose=5, max_nbytes='16M', temp_folder=temp_dir, pre_dispatch=dispatch_policy)\
            (delayed(transform_stream_to_rf)(write_queue, id, stream3c, resample_rate, taper_limit, rotation_type,
                                             filter_band, gauss_width, water_level, trim_start_time, trim_end_time,
                                             deconv_domain, spiking=spiking, normalize=normalize)
             for _, id, _, stream3c in IterRfH5FileEvents(input_file, memmap, channel_pattern))
    else:
        # Process in serial
        logger.info("Serial processing")
        status = list((transform_stream_to_rf(write_queue, id, stream3c, resample_rate, taper_limit, rotation_type,
                                              filter_band, gauss_width, water_level, trim_start_time, trim_end_time,
                                              deconv_domain, spiking=spiking, normalize=normalize)
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


# --------------Main---------------------------------

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--resample-rate', type=float, default=DEFAULT_RESAMPLE_RATE_HZ, show_default=True,
              help="Resampling rate in Hz")
@click.option('--channel-pattern', type=str,
              help="Ordered list of preferred channels, e.g. 'HH*,BH*', where channel selection is ambiguous.")
@click.option('--taper-limit', type=click.FloatRange(0.0, 0.5), default=DEFAULT_TAPER_LIMIT, show_default=True,
              help="Fraction of signal to taper at end, between 0 and 0.5")
@click.option('--filter-band', type=(float, float), default=DEFAULT_FILTER_BAND_HZ, show_default=True,
              help="Filter pass band (Hz). Only required for time-domain deconvolution.")
@click.option('--gauss-width', type=float, default=DEFAULT_GAUSS_WIDTH, show_default=True,
              help="Gaussian freq domain filter width. Only required for freq-domain deconvolution")
@click.option('--water-level', type=float, default=DEFAULT_WATER_LEVEL, show_default=True,
              help="Water-level for freq domain spectrum. Only required for freq-domain deconvolution")
@click.option('--spiking', type=float, default=DEFAULT_SPIKING, show_default=True,
              help="Spiking factor (noise suppression), only required for time-domain deconvolution")
@click.option('--trim-start-time', type=float, default=DEFAULT_TRIM_START_TIME_SEC, show_default=True,
              help="Trace trim start time in sec, relative to onset")
@click.option('--trim-end-time', type=float, default=DEFAULT_TRIM_END_TIME_SEC, show_default=True,
              help="Trace trim end time in sec, relative to onset")
@click.option('--rotation-type', type=click.Choice(['zrt', 'lqt'], case_sensitive=False),
              default=DEFAULT_ROTATION_TYPE, show_default=True,
              help="Rotational coordinate system for aligning ZNE trace components with incident wave direction")
@click.option('--deconv-domain', type=click.Choice(['time', 'freq', 'iter'], case_sensitive=False),
              default=DEFAULT_DECONV_DOMAIN, show_default=True,
              help="Whether to perform deconvolution in time or freq domain, or iterative technique")
@click.option('--normalize/--no-normalize', default=True, show_default=True, help="Whether to normalize RF amplitude")
@click.option('--parallel/--no-parallel', default=True, show_default=True, help="Use parallel execution")
@click.option('--memmap/--no-memmap', default=False, show_default=True,
              help="Memmap input file for improved performance in data reading thread. Useful when data input "
                   "is bottleneck, if system memory permits.")
@click.option('--temp-dir', type=click.Path(dir_okay=True), help="Temporary directory to use for best performance")
@click.option('--aggressive-dispatch/--no-aggressive-dispatch', default=False, show_default=True,
              help="Dispatch all worker jobs as aggressively as possible to minimize chance of worker being "
                   "starved of work. Uses more memory.")
def _main(input_file, output_file, resample_rate, taper_limit, filter_band, gauss_width, water_level, spiking,
         trim_start_time, trim_end_time, rotation_type, deconv_domain, normalize=True, parallel=True, memmap=False,
         temp_dir=None, aggressive_dispatch=False, channel_pattern=None):
    # Dispatch call to worker function. See worker function for documentation.
    event_waveforms_to_rf(input_file, output_file, resample_rate, taper_limit, filter_band,
                          trim_start_time, trim_end_time,  rotation_type, deconv_domain,
                          gauss_width, water_level, spiking, normalize, parallel, memmap,
                          temp_dir, aggressive_dispatch, channel_pattern)
# end main


if __name__ == "__main__":
    _main()  # pylint: disable=no-value-for-parameter
