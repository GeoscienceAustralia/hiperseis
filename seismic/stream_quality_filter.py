#!/usr/bin/env python
"""
Helper functions for curating and quality controlling stream objects.
"""
import json

import numpy as np
from scipy import stats

from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.stream_processing import zrt_order, back_azimuth_filter


def curate_stream3c(ev_id, stream3c, logger=None):
    """
    Apply quality curation criteria to a stream. Modifies the stream in-place if required. Traces in stream
    must be in ZNE order. Each trace in the stream is expected to have metadata with starttime, endtime,
    channel and inclination.

    The following checks are made. If any of these checks fails, the function returns False:
    * Inclination value is not NaN
    * The stream has 3 channels
    * Each trace in the stream has the same number of samples
    * None of the traces have any NaN values in the time series data
    * None of the traces have zero variance

    The following cleanups are attempted on the stream:
    * All 3 traces have the same time range

    :param ev_id: The event id
    :type ev_id: int or str
    :param stream3c: Stream with 3 components of trace data
    :type stream3c: obspy.Stream
    :param logger: Logger in which to log messages
    :type logger: logging.Logger
    :return: True if checks pass, False otherwise
    :rtype: bool
    """

    # Apply essential sanity checks before trying to compute RFs.
    for tr in stream3c:
        if np.isnan(tr.stats.inclination):
            if logger:
                logger.warning("Invalid inclination found in stream {}:\n{}".format(ev_id, stream3c))
            return False
    # end for

    if len(stream3c) != 3:
        if logger:
            logger.warning("Unexpected number of channels in stream {}:\n{}".format(ev_id, stream3c))
        return False
    # end if

    # Strongly assert expected ordering of traces. This must be respected so that
    # RF normalization works properly.
    assert stream3c.traces[0].stats.channel[-1] == 'Z', stream3c.traces[0].stats.channel
    assert stream3c.traces[1].stats.channel[-1] == 'N', stream3c.traces[1].stats.channel
    assert stream3c.traces[2].stats.channel[-1] == 'E', stream3c.traces[2].stats.channel

    # If traces have inconsistent time ranges, clip to time range during which they overlap. Not guaranteed
    # to make time ranges consistent due to possible inconsistent sampling times across channels.
    start_times = np.array([tr.stats.starttime for tr in stream3c])
    end_times = np.array([tr.stats.endtime for tr in stream3c])
    if not (np.all(start_times == start_times[0]) and np.all(end_times == end_times[0])):
        clip_start_time = np.max(start_times)
        clip_end_time = np.min(end_times)
        stream3c.trim(clip_start_time, clip_end_time)
    # end if

    if len(stream3c[0]) != len(stream3c[1]) or len(stream3c[0]) != len(stream3c[2]):
        if logger:
            logger.warning("Channels in stream {} have different lengths:\n{}".format(ev_id, stream3c))
        return False
    # end if

    for tr in stream3c:
        # Each tr here is one component.
        # Check for any NaNs in any component. Discard such traces, as we don't want any NaNs
        # propagating through workflow.
        if np.any(np.isnan(tr.data)):
            if logger:
                logger.warning("NaN detected in trace {} of stream {}:\n{}"
                               .format(tr.stats.channel, ev_id, stream3c))
            return False
        # end if
        # Check for all zeros or all same value in any component. This is infeasible for an operational
        # station, and has been observed as a failure mode in practice.
        if np.std(tr.data) == 0:
            if logger:
                logger.warning("Invariant data detected in trace {} of stream {}:\n{}"
                               .format(tr.stats.channel, ev_id, stream3c))
            return False
        # end if
    # end for

    return True
# end func

def curate_seismograms(data_all, curation_opts, logger, rotate_to_zrt=True):
    """
    Curation function to remove bad data from streams. Note that this function
    will modify the input dataset during curation.

    :param data_all: NetworkEventDataset containing seismograms to curate. Data will be modified by this function.
    :type data_all: seismic.network_event_dataset.NetworkEventDataset
    :param curation_opts: Dict containing curation options.
    :type curation_opts: dict
    :param logger: Logger for emitting log messages
    :type logger: logging.Logger
    :param rotate_to_zrt: Whether to automatically rotate to ZRT coords.
    :type rotate_to_zrt: bool
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

    def rms_ampl_filter(_stream, rms_ampl_bounds, rotation_needed):
        if rotation_needed:
            _stream = _stream.copy()
            _stream.rotate('NE->RT')
        # end if
        _stream.traces.sort(key=zrt_order)
        tr_z = _stream[0]
        tr_r = _stream[1]
        tr_t = _stream[2]
        # Ratio of RMS amplitudes
        rms_z = np.sqrt(np.nanmean(np.square(tr_z.data)))
        rms_r = np.sqrt(np.nanmean(np.square(tr_r.data)))
        rms_t = np.sqrt(np.nanmean(np.square(tr_t.data)))
        r_on_z = rms_r/rms_z
        t_on_z = rms_t/rms_z
        t_on_r = rms_t/rms_r
        r_on_z_limit = rms_ampl_bounds.get("R/Z")
        t_on_z_limit = rms_ampl_bounds.get("T/Z")
        t_on_r_limit = rms_ampl_bounds.get("T/R")
        keep = True
        keep = (r_on_z < r_on_z_limit) if r_on_z_limit is not None else keep
        keep = keep and (t_on_z < t_on_z_limit) if t_on_z_limit is not None else keep
        keep = keep and (t_on_r < t_on_r_limit) if t_on_r_limit is not None else keep
        return keep
    # end func

    def rz_corrcoef_filter(_stream, rz_min_xcorr_coeff, rotation_needed):
        if rotation_needed:
            _stream = _stream.copy()
            _stream.rotate('NE->RT')
        # end if
        # Correlation coefficient
        tr_z = _stream.select(component='Z')[0].data
        tr_r = _stream.select(component='R')[0].data
        corr_c = np.corrcoef([tr_z, tr_r])
        z_cov_r = corr_c[0, 1]
        return z_cov_r >= rz_min_xcorr_coeff
    # end func

    logger.info("Curating input data...")
    logger.info('Curation options:\n{}'.format(json.dumps(curation_opts, indent=4)))

    # Apply curation to streams prior to rotation
    data_all.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

    if "baz_range" in curation_opts:
        # Filter by back-azimuth
        baz_range = curation_opts["baz_range"]
        assert len(baz_range) == 2
        data_all.curate(lambda _1, _2, stream: back_azimuth_filter(stream[0].stats.back_azimuth, baz_range))
    # end if

    # Keep track of whether we have already rotated
    rotated = False

    # Rotate to ZRT coordinates
    if rotate_to_zrt:
        data_all.apply(lambda stream: stream.rotate('NE->RT'))
        rotated = True
    # end if

    # Detrend the traces
    data_all.apply(lambda stream: stream.detrend('linear'))

    # Compute SNR of Z component prior to filtering to use as a quality metric
    if "min_snr" in curation_opts:
        data_all.apply(stream_snr_compute)
    # end if

    if "max_raw_amplitude" in curation_opts:
        # Filter streams with spuriously high amplitude
        max_amp = curation_opts["max_raw_amplitude"]
        data_all.curate(lambda _1, _2, stream: amplitude_nominal(stream, max_amp))
    # end if

    # Spectral filtering
    f_min = curation_opts.get("freq_min")
    f_max = curation_opts.get("freq_max")
    if f_min is not None and f_max is None:
        # Run high pass filter
        data_all.apply(lambda stream: stream.filter('highpass', freq=f_min, corners=2, zerophase=True))
    elif f_min is None and f_max is not None:
        # Run low pass filter
        data_all.apply(lambda stream: stream.filter('lowpass', freq=f_max, corners=2, zerophase=True))
    elif f_min is not None and f_max is not None:
        # Run bandpass filter
        data_all.apply(lambda stream: stream.filter('bandpass', freqmin=f_min, freqmax=f_max,
                                                    corners=2, zerophase=True))
    # end if

    # It does not make sense to filter by similarity, since these are raw waveforms, not RFs,
    # and the waveform will be dominated by the source waveform which differs for each event.

    # Filter by z-component SNR prior to filtering
    if "min_snr" in curation_opts:
        min_snr = curation_opts["min_snr"]
        data_all.curate(lambda _1, _2, stream:
                        stream.select(component='Z')[0].stats.snr_prior >= min_snr)
    # end if

    # Filter by p-arrival slope_ratio
    if "min_slope_ratio" in curation_opts:
        min_slope_ratio = curation_opts["min_slope_ratio"]
        data_all.curate(lambda _1, _2, stream:
                        stream.select(component='Z')[0].stats.slope_ratio >= min_slope_ratio)
    # end if

    # Filter by bounds on RMS channel amplitudes
    rms_ampl_bounds = curation_opts.get("rms_amplitude_bounds")
    if rms_ampl_bounds is not None:
        rotation_needed = not rotated
        data_all.curate(lambda _1, _2, stream:
                        rms_ampl_filter(stream, rms_ampl_bounds, rotation_needed))
    # end if

    # Filter by bounds on correlation coefficients between Z and R traces
    rz_min_xcorr_coeff = curation_opts.get("rz_min_corrcoef")
    if rz_min_xcorr_coeff is not None:
        rotation_needed = not rotated
        data_all.curate(lambda _1, _2, stream:
                        rz_corrcoef_filter(stream, rz_min_xcorr_coeff,  rotation_needed))
    # end if

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
# end func
