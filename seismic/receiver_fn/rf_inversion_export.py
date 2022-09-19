#!/usr/bin/env python
"""Export RFs to file format for external inversion code to run on.
"""

import os
# import logging

import numpy as np
import click
import rf
import pandas as pd
from seismic.receiver_fn import rf_util, rf_corrections
from collections import defaultdict
from scipy import stats
from seismic.stream_io import get_obspyh5_index

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

def rf_inversion_export(input_h5_file, output_folder, network_list="*", station_list="*", station_weights_fn=None,
                        min_station_weight=-1, apply_amplitude_filter=False, apply_similarity_filter=False,
                        min_slope_ratio=-1, dereverberate=False, baz_range=(0, 360),
                        component='R', resample_freq=6.25, trim_window=(-5.0, 20.0), moveout=True):
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
    :param component: The channel component to export, defaults to 'R'
    :type component: str, optional
    :param resample_freq: Sampling rate (Hz) of the output files, defaults to 6.25 Hz
    :type resample_freq: float, optional
    :param trim_window: Time window to export relative to onset, defaults to (-5.0, 20.0). If data needs
        to be resampled, the samples are anchored to the start of this time window.
    :type trim_window: tuple, optional
    :param moveout: Whether to apply moveout correction prior to exporting, defaults to True
    :type moveout: bool, optional
    """
    # Process for each station:
    # 1. Load hdf5 file containing RFs
    # 2. Filter to desired component.
    # 3. Filter stations that fail the minimum-weight criteria
    # 4. Apply amplitude filter if provided
    # 4. Apply min-slope-ratio filter if provided
    # 5. Apply dereverberation filter if specified
    # 6. Quality filter to those that meet criteria (Sippl cross-correlation similarity)
    # 7. Moveout and stack the RFs
    # 8. Resample (lanczos) and trim RF
    # 9. Export one file per station in (time, amplitude format)

    hdfkeys = get_obspyh5_index(input_h5_file, seeds_only=True)

    # trim stations to be processed based on the user-provided network- and station-list
    hdfkeys = rf_util.trim_hdf_keys(hdfkeys, network_list, station_list)

    # Iterate over all data, grouped by (net, sta, loc)
    for hdfkey in hdfkeys:
        net, sta, loc = hdfkey.split('.')

        # Load data
        data = rf_util.read_h5_rf(input_h5_file, network=net, station=sta, loc=loc)

        # Select component
        data = data.select(component=component)

        #rf_util.label_rf_quality_simple_amplitude('ZRT', data, snr_cutoff=2.0, rms_amp_cutoff=0.2, max_amp_cutoff=2.0)
        #data = rf.RFStream([tr for tr in data if tr.stats.predicted_quality == 'a'])

        if apply_amplitude_filter:
            # Label and filter quality
            rf_util.label_rf_quality_simple_amplitude('ZRT', data)
            data = rf.RFStream([tr for tr in data if tr.stats.predicted_quality == 'a'])
            if (len(data) == 0):
                print ("Amplitude filter has removed all traces for {}. "
                       "Ensure rf_quality_filter was run beforehand..".format(hdfkey))
            # end if
        # end if

        # Convert data to a hierarchical format, keyed by sta, cha
        data_dict = rf_util.rf_to_dict(data)
        network_code = data_dict.network

        weights_dict = None
        if(station_weights_fn): weights_dict = read_weights(station_weights_fn)

        for sta, ch_dict in data_dict:
            # Drop stations by weight
            if(weights_dict):
                wt_key = '{}.{}'.format(net, sta)
                if(weights_dict[wt_key] < min_station_weight):
                    print('Skipping station {}: weight {} falls below the minimum weight({})'.format(wt_key,
                                                                                                     weights_dict[wt_key],
                                                                                                     min_station_weight))
                    continue
                # end if
            # end if

            for cha, ch_traces in ch_dict.items():
                if len(ch_traces) < 3:
                    continue
                # end if

                # Drop traces that cannot be stacked
                before = len(ch_traces)
                all_trace_lens = np.array([len(tr) for tr in ch_traces])
                most_common_len, _ = stats.mode(all_trace_lens, axis=None)
                ch_traces = rf.RFStream([tr for tr in ch_traces if len(tr) == most_common_len])
                after = len(ch_traces)
                if after < before:
                    print('{}.{}: {}/{} traces dropped to make them stackable!'.format(network_code, sta,
                                                                                       before-after, after))
                # end if

                # Apply min-slope-ratio filter
                if (min_slope_ratio > 0):
                    before = len(ch_traces)
                    ch_traces = rf.RFStream([tr for tr in ch_traces \
                                             if tr.stats.slope_ratio >= min_slope_ratio])
                    after = len(ch_traces)

                    print('{}.{}: {}/{} traces dropped with min-slope-ratio filter..'.format(network_code, sta,
                                                                                             before - after,
                                                                                             before))
                # end if

                if (len(ch_traces) == 0):
                    print('{}.{}: no traces left to process..'.format(network_code, sta))
                    continue
                # end if

                # Apply de-reverberation filter if specified
                if(dereverberate):
                    has_reverberations = rf_corrections.has_reverberations(ch_traces)
                    if (has_reverberations):
                        print('{}.{}: removing reverberations..'.format(network_code, sta))
                        ch_traces = rf_corrections.apply_reverberation_filter(ch_traces)
                    # end if
                # end if

                # Apply baz-range filter
                if(baz_range != (0, 360)):
                    before = len(ch_traces)
                    ch_traces = rf.RFStream([tr for tr in ch_traces \
                                             if((tr.stats.back_azimuth >= baz_range[0]) and
                                                (tr.stats.back_azimuth <= baz_range[1]))])
                    after = len(ch_traces)
                    print('{}.{}: {}/{} traces dropped with baz-range filter..'.format(network_code, sta,
                                                                                       before - after,
                                                                                       before))
                # end if

                # Apply trace-similarity filter
                if(apply_similarity_filter):
                    before = len(ch_traces)
                    ch_traces = rf_util.filter_crosscorr_coeff(rf.RFStream(ch_traces), time_window=trim_window,
                                                                    apply_moveout=True)
                    after = len(ch_traces)
                    print('{}.{}: {}/{} traces dropped with trace-similarity filter..'.format(network_code, sta,
                                                                                              before - after,
                                                                                              before))
                # end if

                if (len(ch_traces) == 0):
                    print('{}.{}: No traces left to stack. Moving along..'.format(network_code, sta))
                    continue
                # end if

                if moveout:
                    ch_traces.moveout()
                # end if

                # report stats for traces included in stack
                print('{}.{}: Traces included in stack ({}): '.format(network_code, sta, len(ch_traces)))
                for strc in ch_traces:
                    print('\t Event id, time, lon, lat, baz: {}, {}, {:6.2f}, {:6.2f}, {:6.2f}'.format(
                                                                                        strc.stats.event_id,
                                                                                        strc.stats.event_time,
                                                                                        strc.stats.event_longitude,
                                                                                        strc.stats.event_latitude,
                                                                                        strc.stats.back_azimuth))
                # end for

                stack = ch_traces.stack()
                trace = stack[0]
                exact_start_time = trace.stats.onset + trim_window[0]
                stack.interpolate(sampling_rate=resample_freq, method='lanczos', a=10, starttime=exact_start_time)
                stack.trim2(*trim_window, reftime='onset')

                times = trace.times() - (trace.stats.onset - trace.stats.starttime)
                # TODO: Remove hardwired scaling factor.
                # This scaling factor only applies to iterative deconvolution with default Gaussian width
                # factor of 2.5. Once we upgrade to rf library version >= 0.9.0, we can remove this hardwired
                # setting and instead have it determined programatically from rf processing metadata stored
                # in the trace stats structure.
                # The scaling factor originates in the amplitude attenuation effect of the filtering applied
                # in iterative deconv, see table at end of this page:
                # http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/seq01.html
                # The values in this reference table are derived as the integral of the area under the
                # Gaussian in the frequency domain. Analytically, this amounts to simply dividing by scaling
                # factor of a/sqrt(pi), where 'a' here is the Gaussian width used in iterative deconvolution.
    #            iterdeconv_scaling = 2.5/np.sqrt(np.pi)
                iterdeconv_scaling = 1
                column_data = np.array([times, trace.data/iterdeconv_scaling]).T
                fname = os.path.join(output_folder, "_".join([network_code, sta, trace.stats.location, cha]) + "_rf.dat")
                np.savetxt(fname, column_data, fmt=('%5.2f', '%.8f'))
            # end for
        # end for
    # end for
# end func

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(dir_okay=True), required=True)
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-weights', type=str, default=None, show_default=True,
              help='A comma-separated text file containing network and station in the first two columns, '
                   'respectively, with the last column being the weight')
@click.option('--min-station-weight', type=float, default=-1, show_default=True,
              help='Stations with a weight below this value are discarded. Note that this parameter '
                   'has no effect if --station-weights is not specified')
@click.option('--apply-amplitude-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF amplitude filtering to the RFs. The default filtering logic includes: '
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
@click.option('--dereverberate', is_flag=True, default=False, help='Apply de-reverberation filter')
def main(input_file, output_folder, network_list, station_list, station_weights, min_station_weight,
         apply_amplitude_filter, apply_similarity_filter, min_slope_ratio, baz_range,
         dereverberate):

    assert baz_range[1] > baz_range[0], 'Invalid min/max back-azimuth; Aborting..'

    rf_inversion_export(input_file, output_folder, network_list, station_list,
                        station_weights_fn=station_weights,
                        min_station_weight=min_station_weight,
                        apply_amplitude_filter=apply_amplitude_filter,
                        apply_similarity_filter=apply_similarity_filter,
                        min_slope_ratio=min_slope_ratio,
                        baz_range=baz_range,
                        dereverberate=dereverberate)
# end func

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
