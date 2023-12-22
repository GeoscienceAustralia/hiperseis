#!/usr/bin/env python
"""Export RFs to file format for external inversion code to run on.
"""

import os, sys
import numpy as np
import click
import rf
from obspy.core import Stream
from rf import RFStream
import pandas as pd
from seismic.receiver_fn import rf_util, rf_corrections
from collections import defaultdict
from scipy import stats
from seismic.stream_io import get_obspyh5_index
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from seismic.misc import setup_logger

# pylint: disable=invalid-name, logging-format-interpolation

def read_weights(file_name):
    df = None
    try:
        df = pd.read_csv(file_name, skiprows=2)
    except:
        assert 0, 'failed to read {}. Aborting..'.format(file_name)
    # end if

    result = defaultdict(float)

    for i in np.arange(len(df)):
        key = '{}.{}'.format(df.iloc[i, 0], df.iloc[i, 1]) # net.sta
        result[key] = df.iloc[i, -1]
    # end for

    if(not len(result)):
        assert 0, 'No weights found in {}'.format(file_name)
    # end if

    return result
# end func

def phase_weighted_stack(stream:RFStream, phase_weight=1.):
    if(len(stream)==0): return None

    iphase = []
    data = []
    for tr in stream:
        data.append(tr.data)
        analytic = hilbert(tr.data)
        angle = np.angle(analytic)
        iphase.append(np.exp(1j * angle))
    # end for
    data = np.array(data).T
    iphase_wt = np.abs(np.mean(np.array(iphase).T, axis=1))

    samples = np.mean(data, axis=1)
    samples *= np.power(iphase_wt, phase_weight)

    # Call RF-stack for setting appropriate headers
    dummy_stack = stream.stack()
    dummy_stack[0].data = samples

    return dummy_stack
# end func

def rf_inversion_export(input_h5_file, output_folder, network_list="*", station_list="*", station_weights_fn=None,
                        min_station_weight=-1, apply_amplitude_filter=False, apply_similarity_filter=False,
                        min_slope_ratio=-1, dereverberate=False, baz_range=(0, 360),
                        apply_phase_weighting=False, pw_exponent=1.,
                        resample_freq=None, trim_window=(-5.0, 30.0),
                        moveout=True, logger=None):
    """Export receiver function to text format for ingestion into Fortran RF inversion code.

    :param input_h5_file: Input hdf5 file containing receiver function data
    :type input_h5_file: str or Path
    :param output_folder: Folder in which to export text files, one per channel per station.
        Will be appended with network code.
    :type output_folder: str or Path
    :param network_list: List of space-separated networks to process
    :type network_list: str
    :param station_list: List of space-separated stations to process
    :type station_list: str
    :param min_station_weight: Minimum station weight
    :type float
    :param apply_amplitude_filter
    :type bool
    :param apply_similarity_filter
    :type: bool
    :param min_slope_ratio: value of minimum-slope-ratio trace attribute to use to filter bad RF traces
    :type min_slope_ratio: float
    :param dereverberate
    type: bool
    :param resample_freq: Sampling rate (Hz) of the output files. The default (None) preserves original sampling rate
    :type resample_freq: float, optional
    :param trim_window: Time window to export relative to onset, defaults to (-5.0, 30.0). If data needs
        to be resampled, the samples are anchored to the start of this time window.
    :type trim_window: tuple, optional
    :param moveout: Whether to apply moveout correction prior to exporting, defaults to True
    :type moveout: bool, optional
    :@return obspy Stream of traces exported
    """
    # Process for each station:
    # 1. Load hdf5 file containing RFs
    # 2. Filter to desired component.
    # 3. Filter stations that fail the minimum-weight criteria
    # 4. Apply amplitude filter if provided
    # 5. Apply min-slope-ratio filter if provided
    # 6. Apply dereverberation filter if specified
    # 7. Apply back-azimuth filter if specified
    # 8. Quality filter to those that meet criteria (Sippl cross-correlation similarity)
    # 9. Moveout and stack the RFs
    # 10. Resample (lanczos) and trim RF
    # 11. Export one file per station in (time, amplitude format)

    TRIM_BUFFER = 10
    if(logger is None): logger = setup_logger('__func__')

    hdfkeys = get_obspyh5_index(input_h5_file, seeds_only=True)

    # trim stations to be processed based on the user-provided network- and station-list
    hdfkeys = rf_util.trim_hdf_keys(hdfkeys, network_list, station_list)

    outputStream = Stream()
    # Iterate over all data, grouped by (net, sta, loc)
    for hdfkey in hdfkeys:
        net, sta, loc = hdfkey.split('.')

        # Load data
        data = rf_util.read_h5_rf(input_h5_file, network=net, station=sta, loc=loc)

        # Select component
        data_dict = rf_util.rf_to_dict(data)
        component = rf_util.choose_rf_source_channel(data_dict[sta])[-1]
        data = data.select(component=component)

        # Convert data to a hierarchical format, keyed by sta, cha
        data_dict = rf_util.rf_to_dict(data)
        network_code = data_dict.network
        rf_type = 'prf' if data_dict.phase == 'P' else 'srf' if data_dict.phase == 'S' else None

        if (apply_amplitude_filter and rf_type == 'prf'):
            # Label and filter quality
            rf_util.label_rf_quality_simple_amplitude('ZRT', data)
            data = rf.RFStream([tr for tr in data if tr.stats.predicted_quality == 'a'])
            if (len(data) == 0):
                logger.warn("Amplitude filter has removed all traces for {}. "
                            "Ensure rf_quality_filter was run beforehand..".format(hdfkey))
            # end if
        elif (apply_amplitude_filter and rf_type == 'srf'):
            logger.warn('Amplitude filter is only applicable for P RFs; skipping..')
        # end if

        weights_dict = None
        if(station_weights_fn): weights_dict = read_weights(station_weights_fn)

        for sta, ch_dict in data_dict:
            # Drop stations by weight
            if(weights_dict):
                wt_key = '{}.{}'.format(net, sta)
                if(weights_dict[wt_key] < min_station_weight):
                    logger.warn('Skipping station {}: weight {} falls below the minimum weight({})'.format(wt_key,
                                                                                                     weights_dict[wt_key],
                                                                                                     min_station_weight))
                    continue
                # end if
            # end if

            for cha, ch_traces in ch_dict.items():
                if len(ch_traces) < 3:
                    continue
                # end if

                # Trim data with TRIM_BUFFER seconds around trim_window relative to 'onset'.
                # This ensures slightly shorter RF traces are not culled for not being the
                # same length as the majority of the traces, while also ensuring processing
                # artefacts confined to TRIM_BUFFER are excluded in the final trim down to
                # trim_seconds around 'onset'.
                buffered_trim_window = trim_window[0] - TRIM_BUFFER, trim_window[1] + TRIM_BUFFER
                data.trim2(*buffered_trim_window, reftime='onset')

                # Drop traces that cannot be stacked
                before = len(ch_traces)
                all_trace_lens = np.array([len(tr) for tr in ch_traces])
                most_common_len, _ = stats.mode(all_trace_lens, axis=None)
                ch_traces = rf.RFStream([tr for tr in ch_traces if len(tr) == most_common_len])
                after = len(ch_traces)
                if after < before:
                    logger.info('{}.{}.{}: {}/{} traces dropped to make them stackable!'.format(network_code, sta, loc,
                                                                                       before-after, after))
                # end if

                # Apply min-slope-ratio filter
                if (min_slope_ratio > 0):
                    before = len(ch_traces)
                    ch_traces = rf.RFStream([tr for tr in ch_traces \
                                             if tr.stats.slope_ratio >= min_slope_ratio])
                    after = len(ch_traces)

                    logger.info('{}.{}.{}: {}/{} traces dropped with min-slope-ratio filter..'.format(network_code, sta,
                                                                                             loc, before - after,
                                                                                             before))
                # end if

                if (len(ch_traces) == 0):
                    logger.info('{}.{}.{}: no traces left to process..'.format(network_code, sta, loc))
                    continue
                # end if

                # Apply de-reverberation filter if specified
                if(dereverberate and rf_type == 'prf'):
                    has_reverberations = rf_corrections.has_reverberations(ch_traces)
                    if (has_reverberations):
                        logger.info('{}.{}.{}: removing reverberations..'.format(network_code, sta, loc))
                        ch_traces = rf_corrections.apply_reverberation_filter(ch_traces)
                    # end if
                elif(dereverberate and rf_type=='srf'):
                    logger.warn('Dereverberation filter is only applicable for P RFs, skipping..')
                # end if

                # Apply baz-range filter
                if(baz_range != (0, 360)):
                    before = len(ch_traces)
                    ch_traces = rf.RFStream([tr for tr in ch_traces \
                                             if((tr.stats.back_azimuth >= baz_range[0]) and
                                                (tr.stats.back_azimuth <= baz_range[1]))])
                    after = len(ch_traces)
                    logger.info('{}.{}.{}: {}/{} traces dropped with baz-range filter..'.format(network_code, sta, loc,
                                                                                       before - after,
                                                                                       before))
                # end if

                # Apply trace-similarity filter
                if(apply_similarity_filter):
                    before = len(ch_traces)
                    ch_traces = rf_util.filter_crosscorr_coeff(rf.RFStream(ch_traces), time_window=trim_window,
                                                                    apply_moveout=True)
                    after = len(ch_traces)
                    logger.info('{}.{}.{}: {}/{} traces dropped with trace-similarity filter..'.format(network_code, sta,
                                                                                              loc, before - after,
                                                                                              before))
                # end if

                # P RF amplitudes should not exceed 1.0 and should peak around onset time --
                # otherwise, such traces are deemed problematic and discarded
                if(rf_type == 'prf'):
                    before = len(ch_traces)
                    ch_traces = rf_util.filter_invalid_radial_component(ch_traces)
                    after = len(ch_traces)
                    if (before > after):
                        logger.info('{}.{}.{}: {}/{} RF traces with amplitudes > 1.0 or troughs around onset time dropped ..'.format(
                              network_code, sta, loc,
                              before - after,
                              before))
                    # end if
                # end if

                if (len(ch_traces) == 0):
                    logger.warn('{}.{}.{}: No traces left to stack. Moving along..'.format(network_code, sta, loc))
                    continue
                # end if

                if moveout:
                    logger.info('{}.{}.{}: Applying moveout..'.format(network_code, sta, loc))
                    ch_traces.moveout()
                # end if

                # report stats for traces included in stack
                logger.info('{}.{}.{}: Traces included in stack ({}): '.format(network_code, sta, loc, len(ch_traces)))
                for strc in ch_traces:
                    logger.info('\t Event id, time, lon, lat, baz: {}, {}, {:6.2f}, {:6.2f}, {:6.2f}'.format(
                                                                                        strc.stats.event_id,
                                                                                        strc.stats.event_time,
                                                                                        strc.stats.event_longitude,
                                                                                        strc.stats.event_latitude,
                                                                                        strc.stats.back_azimuth))
                # end for

                stack = None
                if(apply_phase_weighting):
                    logger.info('{}.{}.{}: Computing a phase-weighted stack..'.format(network_code, sta, loc))

                    stack = phase_weighted_stack(ch_traces, phase_weight=pw_exponent)

                    if(0): # debug plots
                        import matplotlib.pyplot as plt
                        from seismic.receiver_fn.rf_plot_utils import plot_rf_stack
                        fn1 = os.path.join(output_folder, "_".join([network_code, sta, ch_traces[0].stats.location, cha]) + "_linear.pdf")
                        fn2 = os.path.join(output_folder, "_".join([network_code, sta, ch_traces[0].stats.location, cha]) + "_pw.pdf")

                        lin = ch_traces[0].copy()
                        pw = ch_traces[0].copy()

                        lin.data = ch_traces.stack()[0].data
                        pw.data = stack[0].data
                        fig = plot_rf_stack(RFStream([lin]))
                        plt.savefig(fn1)

                        fig = plot_rf_stack(RFStream([pw]))
                        plt.savefig(fn2)
                    # end if
                else:
                    logger.info('{}.{}.{}: Computing a linear stack..'.format(network_code, sta, loc))
                    stack = ch_traces.stack()
                # end if
                trace = stack[0]

                if(resample_freq is not None):
                    exact_start_time = trace.stats.onset + trim_window[0]
                    stack.interpolate(sampling_rate=resample_freq, method='lanczos', a=10, starttime=exact_start_time)
                # end if

                stack.trim2(*trim_window, reftime='onset')

                # Apply a 2 s taper to the right side
                trace.taper(max_percentage=None, max_length=2, side='right')

                outputStream += trace

                times = trace.times() - (trace.stats.onset - trace.stats.starttime)
                column_data = np.array([times, trace.data]).T
                fname = os.path.join(output_folder, "_".join([network_code, sta, trace.stats.location, cha]) + "_rf.dat")
                np.savetxt(fname, column_data, fmt=('%5.2f', '%.8f'),
                           header='{} (lon, lat): {}, {}'.format(trace.id, trace.stats.station_longitude,
                                                                 trace.stats.station_latitude))
            # end for
        # end for
    # end for
    return outputStream
# end func

def generate_plots(stream:Stream, ofn:str):
    """
    QA/QC Plotting routine (adapted from RF.Imaging)
    @param stream:
    @param ofn:
    @return:
    """
    fig_width = 7.; trace_height = 0.5
    stack_height = 0.5; dpi = None
    scale = 0.9; trim = None

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    stats = stream[0].stats
    times = stream[0].times() - (stats.onset - stats.starttime)
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH), dpi=dpi)
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])

    def _plot(ax, t, d, i, label):
        c1, c2 = ('#000000', '#a0a0a0')
        ax.text(10, i + 0.2, label)
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k')
    # end func

    max_ = np.array([np.max(np.abs(tr.data)) for tr in stream])


    for i, tr in enumerate(stream):
        _plot(ax1, times, tr.data / max_[i] * scale, i + 1,
              "{}.{}.{}.{}".format(tr.stats.network,
                                   tr.stats.station,
                                   tr.stats.location,
                                   tr.stats.channel))
    # end for

    ymin, ymax = -0.5, N + 1.5
    #ax1.vlines(x=0, ymin=ymin, ymax=ymax, colors='grey')
    ax1.xaxis.grid()

    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_ylim(ymin, ymax)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    plt.savefig(ofn, dpi=300)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(dir_okay=True), required=True)
@click.argument('output-plot-file', type=click.Path(dir_okay=False), required=True)
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) or a text file '
                                                  'with station names in each row, w/wo location codes.', type=str,
              show_default=True)
@click.option('--station-weights', type=str, default=None, show_default=True,
              help='A comma-separated text file containing network and station in the first two columns, '
                   'respectively, with the last column being the weight')
@click.option('--min-station-weight', type=float, default=-1, show_default=True,
              help='Stations with a weight below this value are discarded. Note that this parameter '
                   'has no effect if --station-weights is not specified')
@click.option('--apply-amplitude-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF amplitude filtering to the P RFs -- not applicable for S RFs. '
                   'The default filtering logic includes: '
                   'Signal SNR >= 2.0 '
                   'RMS amplitude of signal < 0.2 '
                   'Maximum amplitude of signal < 1.0')
@click.option('--apply-similarity-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF similarity filtering to the RFs.')
@click.option('--min-slope-ratio', type=float, default=-1, show_default=True,
              help='Apply filtering to the RFs based on the "slope_ratio" metric '
                   'that indicates robustness of P-arrival. Typically, a minimum '
                   'slope-ratio of 5 is able to pick out strong arrivals. The '
                   'default value of -1 does not apply this filter')
@click.option('--baz-range', type=(float, float), default=(0, 360), show_default=True,
              help='Min and Max back-azimuth, specified as two floating point values separated by space. '
                   'RF traces with a back-azimuth outside this range are dropped')
@click.option('--dereverberate', is_flag=True, default=False,
              help='Apply de-reverberation filter. Only applicable for P RFs.')
@click.option('--apply-phase-weighting', is_flag=True, default=False, show_default=True,
              help='Compute phase-weighed RF stacks. The default is linear stacks.')
@click.option('--pw-exponent', type=float, default=1, show_default=True,
              help='Exponent used in instantaneous phase-weighting of RF amplitudes. This parameter '
                   'has no effect when --apply-phase-weighting is absent.')
@click.option('--resample-rate', default=None, type=float, show_default=True,
              help='Resampling rate (Hz) for output traces.')
def main(input_file, output_folder, output_plot_file, network_list, station_list, station_weights,
         min_station_weight, apply_amplitude_filter, apply_similarity_filter, min_slope_ratio, baz_range,
         dereverberate, apply_phase_weighting, pw_exponent, resample_rate):
    """
      INPUT_FILE : Input RFs in H5 format\n
                   (output of generate_rf.py or rf_quality_filter.py)\n
      OUTPUT_FOLDER: Path to folder to write output files to\n
      OUTPUT_FILE : Output pdf file name for plots\n
    """

    logger = setup_logger('__func__')

    assert baz_range[1] > baz_range[0], 'Invalid min/max back-azimuth; Aborting..'
    outputStream = rf_inversion_export(input_file, output_folder, network_list, station_list,
                                       station_weights_fn=station_weights,
                                       min_station_weight=min_station_weight,
                                       apply_amplitude_filter=apply_amplitude_filter,
                                       apply_similarity_filter=apply_similarity_filter,
                                       min_slope_ratio=min_slope_ratio,
                                       baz_range=baz_range,
                                       dereverberate=dereverberate,
                                       apply_phase_weighting=apply_phase_weighting,
                                       pw_exponent=pw_exponent,
                                       resample_freq=resample_rate,
                                       trim_window=(-5., 30.),
                                       moveout=True,
                                       logger=logger)
    generate_plots(outputStream, output_plot_file)
# end func

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
