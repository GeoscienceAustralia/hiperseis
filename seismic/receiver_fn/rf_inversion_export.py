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

# pylint: disable=invalid-name, logging-format-interpolation

def read_weights(file_name, network_name):
    df = None
    try:
        df = pd.read_csv(file_name, skiprows=2)
    except:
        assert 0, 'failed to read {}. Aborting..'.format(file_name)
    # end if

    result = defaultdict(float)

    for i in np.arange(len(df)):
        if(df.iloc[i, 0] == network_name):
            result[df.iloc[i, 1]] = df.iloc[i, -1]
        # end if
    # end for

    if(not len(result)):
        assert 0, 'No weights found in {} for network {}'.format(file_name, network_name)
    # end if

    return result
# end func

def rf_inversion_export(input_h5_file, output_folder, network_code, station_weights_fn=None, min_station_weight=-1,
                        min_slope_ratio=-1, dereverberate=False, component='R', resample_freq=6.25,
                        trim_window=(-5.0, 20.0), moveout=True):
    """Export receiver function to text format for ingestion into Fortran RF inversion code.

    :param input_h5_file: Input hdf5 file containing receiver function data
    :type input_h5_file: str or Path
    :param output_folder: Folder in which to export text files, one per channel per station.
        Will be appended with network code.
    :type output_folder: str or Path
    :param network_code: Network to which this RF data belongs, used to disambiguate and track folders.
    :type network_code: str
    :param min_slope_ratio: value of minimum-slope-ratio trace attribute to use to filter bad RF traces
    :type min_slope_ratio: float
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
    # 4. Apply min-slope-ratio filter if provided
    # 5. Apply dereverberation filter if specified
    # 6. Quality filter to those that meet criteria (Sippl cross-correlation similarity)
    # 7. Moveout and stack the RFs
    # 8. Resample (lanczos) and trim RF
    # 9. Export one file per station in (time, amplitude format)

    output_folder += "_" + network_code
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    # end if

    #data = rf_util.read_h5_rf(input_h5_file, network='OA', station='CI23', loc='0M')
    data = rf_util.read_h5_rf(input_h5_file)

    data = data.select(component=component)

    rf_util.label_rf_quality_simple_amplitude('ZRT', data, snr_cutoff=2.0, rms_amp_cutoff=0.2, max_amp_cutoff=2.0)
    data = rf.RFStream([tr for tr in data if tr.stats.predicted_quality == 'a'])

    data_dict = rf_util.rf_to_dict(data)
    weights_dict = None
    if(station_weights_fn): weights_dict = read_weights(station_weights_fn, network_code)

    for sta, ch_dict in data_dict:
        if(weights_dict):
            if(weights_dict[sta] < min_station_weight):
                print('Skipping station {}: weight {} falls below the minimum weight({})'.format(sta,
                                                                                                 weights_dict[sta],
                                                                                                 min_station_weight))
                continue
            # end if
        # end if

        for cha, ch_traces in ch_dict.items():
            if len(ch_traces) < 3:
                continue
            # end if

            # Apply min-slope-ratio filter
            if (min_slope_ratio > 0):
                before = len(ch_traces)
                ch_traces= rf.RFStream([tr for tr in ch_traces \
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

            # Apply trace-similarity filter
            before = len(ch_traces)
            similar_traces = rf_util.filter_crosscorr_coeff(rf.RFStream(ch_traces), time_window=trim_window,
                                                            apply_moveout=True)
            after = len(similar_traces)
            print('{}.{}: {}/{} traces dropped with trace-similarity filter..'.format(network_code, sta,
                                                                                      before - after,
                                                                                      before))

            if not similar_traces:
                continue
            if moveout:
                similar_traces.moveout()
            # end if

            # report stats for traces included in stack
            print('{}.{}: Traces included in stack ({}): '.format(network_code, sta, len(similar_traces)))
            for strc in similar_traces:
                print('\t Event id, time, lon, lat: {}, {}, {}, {}'. format(strc.stats.event_id,
                                                                         strc.stats.event_time,
                                                                         strc.stats.event_longitude,
                                                                         strc.stats.event_latitude))
            # end for

            stack = similar_traces.stack()
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

# end func

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(dir_okay=True), required=True)
@click.option('--network-code', type=str, required=True,
              help='Network code that this input file is for. Appended to output folder name.')
@click.option('--station-weights', type=str, default=None, show_default=True,
              help='A comma-separated text file containing network and station in the first two columns, '
                   'respectively, with the last column being the weight')
@click.option('--min-station-weight', type=float, default=-1, show_default=True,
              help='Stations with a weight below this value are discarded. Note that this parameter '
                   'has no effect if --station-weights is not specified')
@click.option('--min-slope-ratio', type=float, default=-1, show_default=True,
              help='Apply filtering to the RFs based on the "slope_ratio" metric '
                   'that indicates robustness of P-arrival. Typically, a minimum '
                   'slope-ratio of 5 is able to pick out strong arrivals. The '
                   'default value of -1 does not apply this filter')
@click.option('--dereverberate', is_flag=True, default=False, help='Apply de-reverberation filter')
def main(input_file, output_folder, network_code, station_weights, min_station_weight,
         min_slope_ratio, dereverberate):  # pylint: disable=missing-docstring
    rf_inversion_export(input_file, output_folder, network_code, station_weights_fn=station_weights,
                        min_station_weight=min_station_weight, min_slope_ratio=min_slope_ratio,
                        dereverberate=dereverberate)
# end func


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
