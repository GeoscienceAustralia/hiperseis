import numpy as np
from seismic.receiver_fn.rf_util import compute_vertical_snr
from seismic.receiver_fn.rf_deconvolution import rf_iter_deconv
from seismic.stream_quality_filter import curate_stream3c
from seismic.stream_processing import back_azimuth_filter
from obspy.signal.rotate import rotate_rt_ne

import logging

logging.basicConfig()

DEFAULT_RESAMPLE_RATE_HZ = 10.0
DEFAULT_FILTER_BAND_HZ = (0.02, 1.00)
DEFAULT_TAPER_LIMIT = 0.05
DEFAULT_TRIM_START_TIME_SEC = -50.0
DEFAULT_TRIM_END_TIME_SEC = 150.0
DEFAULT_ROTATION_TYPE = 'zrt'   # from ['zrt', 'lqt']
DEFAULT_DECONV_DOMAIN = 'time'  # from ['time', 'freq', 'iter']
DEFAULT_GAUSS_WIDTH = 1.0
DEFAULT_WATER_LEVEL = 0.01
DEFAULT_SPIKING = 0.5
RAW_RESAMPLE_RATE_HZ = 10.0
BANDPASS_FILTER_ORDER = 2

def transform_stream_to_rf(ev_id, stream3c, config_filtering,
                           config_processing, **kwargs):
    """Generate P-phase receiver functions for a single 3-channel stream.
    See documentation for function event_waveforms_to_rf for details of
    config dictionary contents.

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
    :return: RFstream containing receiver function if successful, None otherwise
    :rtype: rf.RFStream or NoneType
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

    assert deconv_domain.lower() in ['time', 'freq', 'iter']

    # Apply any custom preprocessing step
    custom_preproc = config_processing.get("custom_preproc")
    if custom_preproc is not None:
        # TODO: Improve security of exec and eval usage here.
        load_statement = custom_preproc.get("import")
        if load_statement is not None:
            exec(load_statement)
        # end if
        # "Func" field must be present
        func = eval(custom_preproc["func"])
        func_args = custom_preproc.get("args")
        if func_args is None:
            func_args = {}
        # Apply custom preprocessing function
        stream3c = func(ev_id, stream3c, **func_args)
        # Functor must return the preprocessed stream
        assert stream3c is not None
    # end if

    if not curate_stream3c(ev_id, stream3c, logger):
        return None

    if baz_range is not None:
        if not isinstance(baz_range[0], list):
            baz_range = [baz_range]
        baz = stream3c[0].stats.back_azimuth
        if not any([back_azimuth_filter(baz, tuple(b)) for b in baz_range]):
            return None
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
        return None
    # end try

    # Check for any empty channels after deconvolution.
    for tr in stream3c:
        if len(tr.data) == 0:
            logger.warning("No data left in channel {} of stream {} after deconv (skipping)".format(
                           tr.stats.channel, ev_id))
            return None
        # end if
    # end for

    # Perform trimming before computing max amplitude.
    if deconv_domain != 'iter':
        stream3c.trim2(trim_start_time_sec, trim_end_time_sec, reftime='onset')
    # end if

    if len(stream3c) != 3:
        logger.warning("Unexpected number of channels in stream {} after trim (skipping):\n{}"
                       .format(ev_id, stream3c))
        return None
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

    return stream3c
# end func
