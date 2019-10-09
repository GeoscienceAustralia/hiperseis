#!/usr/bin/env python
"""Export RFs to file format for external inversion code to run on.
"""

import os
# import logging

import numpy as np
import click
import rf

from seismic.receiver_fn import rf_util

# pylint: disable=invalid-name, logging-format-interpolation


def rf_inversion_export(input_h5_file, output_folder, network_code, component='R',
                        resample_freq=6.25, trim_window=(-5.0, 20.0), moveout=True):
    # Process for each station:
    # 1. Load hdf5 file containing RFs
    # 2. Filter to desired component.
    # 3. Quality filter to those that meet criteria (Sippl cross-correlation similarity)
    # 4. Moveout and stack the RFs
    # 5. Resample (lanczos) and trim RF
    # 6. Export one file per station in (time, amplitude format)

    output_folder += "_" + network_code
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    # end if

    data = rf_util.read_h5_rf(input_h5_file)

    data = data.select(component=component)

    rf_util.label_rf_quality_simple_amplitude('ZRT', data, snr_cutoff=2.0, rms_amp_cutoff=0.2, max_amp_cutoff=2.0)
    data = rf.RFStream([tr for tr in data if tr.stats.predicted_quality == 'a'])

    data_dict = rf_util.rf_to_dict(data)

    for sta, ch_dict in data_dict.items():
        for cha, ch_traces in ch_dict.items():
            similar_traces = rf_util.filter_crosscorr_coeff(rf.RFStream(ch_traces))
            if not similar_traces:
                continue
            if moveout:
                similar_traces.moveout()
            # end if
            stack = similar_traces.stack()
            stack.interpolate(sampling_rate=resample_freq, method='lanczos', a=10)
            stack.trim2(*trim_window, reftime='onset')

            trace = stack[0]
            times = trace.times() - (trace.stats.onset - trace.stats.starttime)
            column_data = np.array([times, trace.data]).T
            fname = os.path.join(output_folder, "_".join([network_code, sta, cha]) + "_rf.dat")
            np.savetxt(fname, column_data, fmt=('%5.2f', '%.8f'))
        # end for
    # end for

# end func


@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(dir_okay=True), required=True)
# @click.option('network-code', type=str, required=True)
def main(input_file, output_folder, network_code='7X'):
    rf_inversion_export(input_file, output_folder, network_code)
# end func


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
