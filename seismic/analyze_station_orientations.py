#!/usr/bin/env python
"""
Analyze a data set of seismic arrival events on a per-station basis and try to
detect and estimate any station orientation error.

In future, consider moving this script to the `inventory` module and applying
corrections to the station inventory xml (to the azimuth tag).

Reference:
Wilde-Piórko, M., Grycuk, M., Polkowski, M. et al.
On the rotation of teleseismic seismograms based on the receiver function technique.
J Seismol 21, 857-868 (2017). https://doi.org/10.1007/s10950-017-9640-x
"""

import json
import click
import logging
import copy
from joblib import Parallel, delayed

import numpy as np
from scipy import optimize, interpolate
from rf import RFStream

import matplotlib.pyplot as plt

from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.runners import curate_seismograms
from seismic.receiver_fn.generate_rf import transform_stream_to_rf

# Take care not to use any curation options that would vary if there were a station orientation error.
DEFAULT_CURATION_OPTS = {
    "min_snr": 2.0,
    "max_raw_amplitude": 20000.0,
    "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
}

DEFAULT_CONFIG_FILTERING = {
    "resample_rate": 10.0,
    "taper_limit": 0.05,
    "filter_band": [0.02, 1.0],
}

DEFAULT_CONFIG_PROCESSING = {
    "rotation_type": "ZRT",
    "deconv_domain": "time",
    "normalize": True,
    "trim_start_time": -10,
    "trim_end_time": 30,
    "spiking": 1.0
}


def _run_single_station(db_evid, angles, config_filtering, config_processing):
    """
    Internal processing function for running sequence of candidate angles
    over a single station.

    :param db_evid: Dictionary of event streams (3-channel ZNE) keyed by event ID.
        Best obtained using class NetworkEventDataset
    :param angles: Sequence of candidate correction angles to try (degrees)
    :param config_filtering: Waveform filtering options for RF processing
    :param config_processing: RF processing options
    :return:
    """
    ampls = []
    for correction in angles:
        rf_stream_all = RFStream()
        for evid, stream in db_evid.items():
            stream_rot = copy.deepcopy(stream)
            for tr in stream_rot:
                tr.stats.back_azimuth += correction
                while tr.stats.back_azimuth < 0:
                    tr.stats.back_azimuth += 360
                while tr.stats.back_azimuth >= 360:
                    tr.stats.back_azimuth -= 360
            # end for

            rf_3ch = transform_stream_to_rf(evid, RFStream(stream_rot),
                                            config_filtering, config_processing)
            if rf_3ch is None:
                continue

            rf_stream_all += rf_3ch
        # end for
        if len(rf_stream_all) > 0:
            rf_stream_R = rf_stream_all.select(component='R')
            rf_stream_R.trim2(-5, 5, reftime='onset')
            rf_stream_R.detrend('linear')
            rf_stream_R.taper(0.1)
            R_stack = rf_stream_R.stack().trim2(-1, 1, reftime='onset')[0].data
            ampl_mean = np.mean(R_stack)
        else:
            ampl_mean = np.nan
        # endif
        ampls.append(ampl_mean)
    # end for
    return ampls
# end func


def analyze_station_orientations(ned, curation_opts, config_filtering,
                                 config_processing, save_plot=False):
    """
    Main processing function for analyzing station orientation using 3-channel
    event waveforms. Uses method of Wilde-Piorko https://doi.org/10.1007/s10950-017-9640-x

    One should not worry about estimates that come back with error of less than
    about 20 degrees from zero, since this analysis provides only an estimate.

    :param ned: NetworkEventDataset containing waveforms to analyze
    :param curation_opts: Seismogram curation options.
        Safe default to use is `DEFAULT_CURATION_OPTS`.
    :param config_filtering: Seismogram filtering options for RF computation.
        Safe default to use is `DEFAULT_CONFIG_FILTERING`.
    :param config_processing: Seismogram RF processing options.
        Safe default to use is `DEFAULT_CONFIG_PROCESSING`.
    :param save_plot: Flag to save plot of analysis results on each station
    :return: Dict of estimated orientation error with net.sta code as the key.
    """
    assert isinstance(ned, NetworkEventDataset), 'Pass NetworkEventDataset as input'
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    evids_orig = set([evid for _, evid, _ in ned])

    # Trim streams to time window
    logger.info('Trimming dataset')
    ned.apply(lambda stream: stream.trim(stream[0].stats.onset - 10, stream[0].stats.onset + 30))

    # Downsample.
    logger.info('Downsampling dataset')
    fs = 20.0
    ned.apply(lambda stream: stream.filter('lowpass', freq=fs/2.0, corners=2, zerophase=True) \
              .interpolate(fs, method='lanczos', a=10))

    logger.info('Curating dataset')
    curate_seismograms(ned, curation_opts, logger, rotate_to_zrt=False)
    evids_to_keep = set([evid for _, evid, _ in ned])
    evids_discarded = evids_orig - evids_to_keep
    logger.info('Discarded {}/{} events'.format(len(evids_discarded), len(evids_orig)))

    logger.info('Analysing arrivals')
    angles = np.linspace(-180, 180, num=20, endpoint=False)
    job_runner = Parallel(n_jobs=-2, verbose=5, max_nbytes='16M')
    jobs = []
    stations = []
    for sta, db_evid in ned.by_station():
        stations.append((sta, len(db_evid)))
        job = delayed(_run_single_station)(db_evid, angles, config_filtering, config_processing)
        jobs.append(job)
    # end for
    sta_ampls = job_runner(jobs)
    sta_ori_metrics = [(sta, N, ampls) for (sta, N), ampls in zip(stations, sta_ampls)]

    x = np.hstack((angles, angles[0] + 360))
    angles_fine = np.linspace(-180, 180, num=361)
    results = {}
    for sta, N, ampls in sta_ori_metrics:
        y = np.array(ampls + [ampls[0]])
        mask = np.isfinite(y)
        n_valid = np.sum(mask)
        if n_valid < 5:
            logger.info('{}: Insufficient data ({} valid points)'.format(sta, n_valid))
            continue
        # end if
        x_valid = x[mask]
        y_valid = y[mask]
        #  Fit cubic spline through all points, used for visualization
        if np.isfinite(y[0]):
            bc_type = 'periodic'
        else:
            bc_type = 'not-a-knot'
        # end if
        interp_cubsp = interpolate.CubicSpline(x_valid, y_valid, bc_type=bc_type)
        yint_cubsp = interp_cubsp(angles_fine)
        # Fit cosine curve to all data points, used to estimate correction angle
        fitted, cov = optimize.curve_fit(lambda _x, _amp, _ph: _amp*np.cos(np.deg2rad(_x + _ph)),
                                         x_valid, y_valid, p0=[0.2, 0.0], bounds=([-1, -180], [1, 180]))
        A_fit, ph_fit = fitted
        y_fitted = A_fit*np.cos(np.deg2rad(angles_fine + ph_fit))
        ph_uncertainty = np.sqrt(cov[1, 1])
        # Estimate correction angle
        angle_max = angles_fine[np.argmax(y_fitted)]
        results['.'.join([ned.network, sta])] = angle_max
        logger.info('{}: {:2.3f}°, stddev {:2.3f}° (N = {:3d})'.format(
            sta, angle_max, ph_uncertainty, N))
        if save_plot:
            _f = plt.figure(figsize=(16,9))
            plt.plot(x_valid, y_valid, 'x', label='Computed P arrival strength')
            plt.plot(angles_fine, yint_cubsp, '--', alpha=0.7, label='Cubic spline fit')
            plt.plot(angles_fine, y_fitted, '--', alpha=0.7, label='Cosine fit')
            plt.grid(linestyle=':', color="#80808080")
            plt.xlabel('Orientation correction (deg)', fontsize=12)
            plt.ylabel('P phase ampl. mean (0-1 sec)', fontsize=12)
            plt.title('{}.{}'.format(ned.network, sta), fontsize=14)
            plt.text(0.9, 0.9, 'N = {}'.format(N), ha='right', va='top',
                     transform=plt.gca().transAxes)
            plt.legend(framealpha=0.7)
            plt.savefig('_'.join([sta, config_processing["deconv_domain"], 'ori.png']), dpi=300)
            plt.close()
        # end if
    # end for

    return results
# end func


def process_event_file(src_h5_event_file, dest_file=None, save_plot=False):
    """
    Use event dataset from an HDF5 file to analyze station for orientation errors.

    :param src_h5_event_file: HDF5 file to load. Typically one created by `extract_event_traces.py` script
    :param dest_file: File in which to save results in JSON format
    :param save_plot: Flag to save plot of analysis results on each station
    :return: None
    """

    ned = NetworkEventDataset(src_h5_event_file)

    results = analyze_station_orientations(ned, curation_opts=DEFAULT_CURATION_OPTS,
                                           config_filtering=DEFAULT_CONFIG_FILTERING,
                                           config_processing=DEFAULT_CONFIG_PROCESSING,
                                           save_plot=save_plot)

    if dest_file is not None:
        with open(dest_file, 'w') as f:
            json.dump(results, f, indent=4)
    # end if
# end func


@click.command()
@click.option('--dest-file', type=click.Path(dir_okay=False),
              help='Output file in which to store results in JSON text format')
@click.option('--save-plot', is_flag=True, help='Flag to save plot per station')
@click.argument('src-h5-event-file', type=click.Path(exists=True, dir_okay=False),
                required=True)
def main(src_h5_event_file, dest_file=None, save_plot=False):
    """
    Run station orientation checks.

    Example usage:
    python seismic/analyze_station_orientations.py --dest-file 7X_ori_estimates.json \
      /g/data/ha3/am7399/shared/7X_RF_analysis/7X_event_waveforms_for_rf_20090616T034200-20110401T231849_rev2.h5

    :param src_h5_event_file: Event waveform file whose waveforms are used to perform checks
    :param dest_file: Output file in which to store results in JSON text format
    :param save_plot: Set option to save plot per station of mean arrival RF amplitude as function of correction angle
    """
    process_event_file(src_h5_event_file, dest_file, save_plot=save_plot)
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
