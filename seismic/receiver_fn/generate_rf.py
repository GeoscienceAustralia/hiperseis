#!/usr/bin/env python

import logging

import numpy as np
import click
from multiprocessing import Process, Queue, Manager

from rf import read_rf, RFStream
from rf import IterMultipleComponents

try:
    from joblib import Parallel, delayed
    parallel_available = True
except ImportError:
    parallel_available = False

from seismic.receiver_fn.rf_h5_file_iterator import IterRfH5FileEvents


logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

DEFAULT_RESAMPLE_RATE_HZ = 100
DEFAULT_FILTER_BAND_HZ = (0.03, 1.50)
DEFAULT_TAPER_LIMIT = 0.01
DEFAULT_TRIM_START_TIME_SEC = -25.0
DEFAULT_TRIM_END_TIME_SEC = 75.0
DEFAULT_DECONV_DOMAIN = 'time'
DEFAULT_GAUSS_WIDTH = 2.0
DEFAULT_WATER_LEVEL = 0.05


def transform_stream_to_rf(ioqueue, ev_id, stream3c, resample_rate_hz, taper_limit, filter_band_hz,
                           gauss_width, water_level, trim_start_time_sec, trim_end_time_sec,
                           deconv_domain=DEFAULT_DECONV_DOMAIN, **kwargs):
    """Generate receiver function for a single 3-channel stream.

    :param ev_id: The event id
    :type ev_id: int
    :param stream3c: Stream with 3 components of trace data
    :type stream3c: rf.RFStream
    :param deconv_domain: Domain in which to perform deconvolution
    :type deconv_domain: str
    :param kwargs: Keyword arguments that will be passed to filtering and deconvolution functions.
    :type kwargs: dict
    :return: Stream containing receiver function
    :rtype: rf.RFStream
    """
    logger = logging.getLogger(__name__)
    logger.info("Event #{}".format(ev_id))

    for tr in stream3c:
        if np.isnan(tr.stats.inclination):
            logger.warning("WARNING: Invalid inclination found in stream {} (skipping):\n{}".format(ev_id, stream3c))
            return

    if len(stream3c) != 3:
        logger.warning("WARNING: Unexpected number of channels in stream {} (skipping):\n{}".format(ev_id, stream3c))
        return

    if len(stream3c[0]) != len(stream3c[1]) or len(stream3c[0]) != len(stream3c[2]):
        logger.warning("WARNING: Channels in stream {} have different lengths, cannot generate RF (skipping):\n{}"
                       .format(ev_id, stream3c))
        return

    assert deconv_domain in ['time', 'freq']
    stream3c.detrend('linear').interpolate(resample_rate_hz)
    stream3c.taper(taper_limit, **kwargs)

    try:
        if deconv_domain == 'time':
            # ZRT receiver functions must be specified
            stream3c.filter('bandpass', freqmin=filter_band_hz[0], freqmax=filter_band_hz[1], corners=2, zerophase=True,
                            **kwargs)
            stream3c.rf(rotate='NE->RT', **kwargs)
        else:
            # Note the parameters of gaussian pulse and its width where
            # Value of "a" | Frequency (hz) at which G(f) = 0.1 |  Approximate Pulse Width (s)
            # 10                      4.8                                0.50
            # 5                       2.4                                0.75
            # 2.5                     1.2                                1.00
            # 1.25                    0.6                                1.50
            # 1.0                     0.5                                1.67 (5/3)
            # 0.625                   0.3                                2.10
            # 0.5                     0.24                               2.36
            # 0.4                     0.2                                2.64
            # 0.2                     0.1                                3.73
            stream3c.rf(rotate='NE->RT', deconvolve='freq', gauss=gauss_width, waterlevel=water_level, **kwargs)
        # end if
    except (IndexError, ValueError) as e:
        logger.error("ERROR: Failed on stream {}:\n{}\nwith error:\n{}".format(ev_id, stream3c, str(e)))
        return
    # end try

    event_id = {'event_id': ev_id}

    for stream_index in range(3):
        stream3c[stream_index].stats.update(event_id)
        amax = {'amax': np.amax(stream3c[stream_index].data)}
        stream3c[stream_index].stats.update(amax)
    # end for
    stream3c.trim2(trim_start_time_sec, trim_end_time_sec, reftime='onset')

    ioqueue.put(stream3c)


def async_write(rfstream_queue, outfile_name, max_buffered=100):
    buffered_streams = []
    logger = logging.getLogger(__name__)
    logger.info("Starting async write thread")
    first_write = True
    while True:
        # Passing None into the queue is taken as signal to flush buffer and terminate thread
        rfstream = rfstream_queue.get()

        terminating = (rfstream is None)
        if not terminating:
            buffered_streams.append(rfstream)
        else:
            logger.info("Flushing result buffer...")

        if len(buffered_streams) >= max_buffered or terminating:
            stream = RFStream()
            for rf in buffered_streams:
                stream.extend(rf)
            if first_write:
                mode = 'w'
                first_write = False
            else:
                mode = 'a'
            stream.write(outfile_name, 'H5', mode=mode)
            logger.info("Flushed {} streams to output file {}".format(len(buffered_streams), outfile_name))

            while buffered_streams:
                buffered_streams.pop()
                rfstream_queue.task_done()

        if terminating:
            rfstream_queue.task_done()
            break

    logger.info("Terminating async write thread")


# --------------Main---------------------------------

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--resample-rate', type=float, default=DEFAULT_RESAMPLE_RATE_HZ, show_default=True,
              help="Resampling rate in Hz")
@click.option('--taper-limit', type=click.FloatRange(0.0, 0.5), default=DEFAULT_TAPER_LIMIT, show_default=True,
              help="Fraction of signal to taper at end, between 0 and 0.5")
@click.option('--filter_band', type=(float, float), default=DEFAULT_FILTER_BAND_HZ, show_default=True,
              help="Filter pass band (Hz). Only required for time-domain deconvolution.")
@click.option('--gauss-width', type=float, default=DEFAULT_GAUSS_WIDTH, show_default=True,
              help="Gaussian freq domain filter width. Only required for freq-domain deconvolution")
@click.option('--water-level', type=float, default=DEFAULT_WATER_LEVEL, show_default=True,
              help="Water-level for freq domain spectrum. Only required for freq-domain deconvolution")
@click.option('--trim-start-time', type=float, default=DEFAULT_TRIM_START_TIME_SEC, show_default=True,
              help="Trace trim start time in sec, relative to onset")
@click.option('--trim-end-time', type=float, default=DEFAULT_TRIM_END_TIME_SEC, show_default=True,
              help="Trace trim end time in sec, relative to onset")
@click.option('--deconv-domain', type=click.Choice(['time', 'freq'], case_sensitive=False),
              default=DEFAULT_DECONV_DOMAIN, show_default=True,
              help="Whether to perform deconvolution in time or freq domain")
@click.option('--parallel/--no-parallel', default=True, show_default=True, help="Use parallel execution")
def main(input_file, output_file, resample_rate, taper_limit, filter_band, gauss_width, water_level,
         trim_start_time, trim_end_time, deconv_domain, parallel=True):
    """
    Main entry point for generating RFs from event traces. See Click documentation for details on arguments.
    """
    if parallel:
        assert parallel_available, "Cannot run parallel as joblib import failed"

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Set up asynchronous buffered writing of results to file
    mgr = Manager()
    write_queue = mgr.Queue()
    output_thread = Process(target=async_write, args=(write_queue, output_file, 1000))
    output_thread.daemon = True
    output_thread.start()

    logger.info("Processing source file {}".format(input_file))
    if parallel:
        # Process in parallel
        logger.info("Parallel processing")
        Parallel(n_jobs=-1, verbose=5, max_nbytes='8M')\
            (delayed(transform_stream_to_rf)(write_queue, id, strm3c, resample_rate, taper_limit, filter_band,
                                             gauss_width, water_level, trim_start_time, trim_end_time, deconv_domain)
             for id, strm3c in enumerate(IterRfH5FileEvents(input_file)))
    else:
        # Process in serial
        logger.info("Serial processing")
        list((transform_stream_to_rf(write_queue, id, strm3c, resample_rate, taper_limit, filter_band, gauss_width,
                                     water_level, trim_start_time, trim_end_time, deconv_domain)
              for id, strm3c in enumerate(IterRfH5FileEvents(input_file))))
    # end if

    # Signal completion
    logger.info("Finishing...")
    write_queue.put(None)
    write_queue.join()

    logger.info("generate_rf SUCCESS!")
# end func


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
