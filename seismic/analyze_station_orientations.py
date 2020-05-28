#!/usr/bin/env python
"""
Analyze a data set of seismic arrival events on a per-station basis and try to
detect and estimate any station orientation error.

In future, consider moving this script to the `inventory` module and applying
corrections to the station inventory xml (to the azimuth tag).

Methods explored:

1. Wang, Xin & Chen, Qi-fu & Li, Juan & Wei, Shengji. (2016).
Seismic Sensor Misorientation Measurement Using P -Wave Particle Motion: An Application to the NECsaids Array.
Seismological Research Letters. 87. 901-911. doi: 10.1785/0220160005

2. Wilde-Piórko, M., Grycuk, M., Polkowski, M. et al.
On the rotation of teleseismic seismograms based on the receiver function technique.
J Seismol 21, 857-868 (2017). https://doi.org/10.1007/s10950-017-9640-x
"""

import click
import logging
import copy
from joblib import Parallel, delayed

import numpy as np
from scipy import stats, optimize, interpolate
from sklearn.decomposition import PCA
from rf import RFStream

import matplotlib.pyplot as plt

from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.runners import curate_seismograms
from seismic.receiver_fn.generate_rf import transform_stream_to_rf


def pdf_prune_outliers(data, cull_n_stddev=2.5):
    n_outside = len(data)
    while n_outside > 0:
        mean, stddev = stats.norm.fit(data)
        mask = ((data < mean - cull_n_stddev*stddev) | (data > mean + cull_n_stddev*stddev))
        n_outside = np.sum(mask)
        data = data[~mask]
    # end while
    return data
# end func


def resids_stats(resids):
    N = len(resids)
    mean = np.mean(resids)
    stddev = np.std(resids)
    stderr = stddev / np.sqrt(N)  # standard error of the mean
    return N, mean, stderr, stddev
# end func


def method_wang(src_h5_event_file, dest_file=None):
    # EXPERIMENT: Wang PCA method.

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.info('Loading dataset')
    ned = NetworkEventDataset(src_h5_event_file)
    evids_orig = set([evid for _, evid, _ in ned])

    # Trim streams to time window
    logger.info('Trimming dataset')
    ned.apply(lambda stream: stream.trim(stream[0].stats.onset - 10, stream[0].stats.onset + 30))

    # Downsample.
    logger.info('Downsampling dataset')
    fs = 20.0
    ned.apply(lambda stream: stream.filter('lowpass', freq=fs/2.0, corners=2, zerophase=True) \
              .interpolate(fs, method='lanczos', a=10))

    curation_opts = {
        "min_snr": 2.0,
        "max_raw_amplitude": 20000.0,
        "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
    }
    logger.info('Curating dataset')
    # This also rotates to ZRT coordinates
    curate_seismograms(ned, curation_opts, logger)
    evids_to_keep = set([evid for _, evid, _ in ned])
    evids_discarded = evids_orig - evids_to_keep
    logger.info('Discarded {}/{} events'.format(len(evids_discarded), len(evids_orig)))

    logger.info('Analysing arrivals')
    for sta, db_evid in ned.by_station():
        resids = []
        for stream in db_evid.values():
            s = stream.copy()
            s.detrend(type='linear')
            s.taper(0.05)
            s.resample(10.0, no_filter=False)
            tr_r = s.select(component='R')[0]
            tr_r.trim(tr_r.stats.onset - 1.0, tr_r.stats.onset + 3.0)
            if tr_r.stats.event_magnitude < 5.5:
                continue
            tr_t = s.select(component='T')[0]
            tr_t.trim(tr_t.stats.onset - 1.0, tr_t.stats.onset + 3.0)
            if len(tr_r) != len(tr_t):
                continue

            data = np.array([tr_r.data, tr_t.data])

            # sklearn method
            sk_pca = PCA()
            sk_pca.fit(data.T)
            # sk_evals = sk_pca.singular_values_
            sk_evecs = sk_pca.components_
            if sk_pca.explained_variance_ratio_[0] < 0.80:
                continue
            sk_baz_error = np.rad2deg(np.arctan2(sk_evecs[0, 1], sk_evecs[0, 0]))
            while sk_baz_error < -90:
                sk_baz_error += 180
            while sk_baz_error > 90:
                sk_baz_error -= 180
            resids.append(sk_baz_error)

        # end for
        resids = np.array(sorted(resids))
        N, mean, stderr, stddev = resids_stats(resids)
        # Detect outliers to Gaussian fit and remove from set before computing stats
        resids_culled = pdf_prune_outliers(resids)
        N_c, mean_c, stderr_c, stddev_c = resids_stats(resids_culled)
        if N >= 5:
            logger.info('{}: {:2.3f}° ± {:2.3f}°, stddev {:2.3f}° (N = {:3d})  '
                        '[CULLED: {:2.3f}° ± {:2.3f}°, stddev {:2.3f}° (N = {:3d})]'
                        .format(sta, mean, stderr, stddev, N,
                                mean_c, stderr_c, stddev_c, N_c))
        else:
            logger.info('{}:  Insufficient data (N = {})'.format(sta, N))
        # end if
        # print(resids)
    # end for

    # if dest_file is not None:
    #     ned.write(dest_file)
    # # end if
# end func


def _run_single_station(db_evid, angles, config_filtering, config_processing):
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
        rf_stream_R = rf_stream_all.select(component='R')
        rf_stream_R.trim2(-5, 5, reftime='onset')
        rf_stream_R.detrend('linear')
        rf_stream_R.taper(0.1)
        R_stack = rf_stream_R.stack().trim2(-1, 1, reftime='onset')[0].data
        ampl_mean = np.mean(R_stack)
        ampls.append(ampl_mean)
    # end for
    return ampls
# end func


def method_wilde_piorko(src_h5_event_file, dest_file=None, save_plot=False):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.info('Loading dataset')
    ned = NetworkEventDataset(src_h5_event_file)
    evids_orig = set([evid for _, evid, _ in ned])

    # Trim streams to time window
    logger.info('Trimming dataset')
    ned.apply(lambda stream: stream.trim(stream[0].stats.onset - 10, stream[0].stats.onset + 30))

    # Downsample.
    logger.info('Downsampling dataset')
    fs = 20.0
    ned.apply(lambda stream: stream.filter('lowpass', freq=fs/2.0, corners=2, zerophase=True) \
              .interpolate(fs, method='lanczos', a=10))

    curation_opts = {
        "min_snr": 2.0,
        "max_raw_amplitude": 20000.0,
        "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
    }
    logger.info('Curating dataset')
    curate_seismograms(ned, curation_opts, logger, rotate_to_zrt=False)
    evids_to_keep = set([evid for _, evid, _ in ned])
    evids_discarded = evids_orig - evids_to_keep
    logger.info('Discarded {}/{} events'.format(len(evids_discarded), len(evids_orig)))

    logger.info('Analysing arrivals')
    config_filtering = {
        "resample_rate": 10.0,
        "taper_limit": 0.05,
        "filter_band": [0.02, 1.0],
    }
    config_processing = {
        "rotation_type": "ZRT",
        "deconv_domain": "iter",
        "normalize": True,
        "trim_start_time": -10,
        "trim_end_time": 30
    }
    angles = np.linspace(-180, 180, num=30, endpoint=False)
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
    for sta, N, ampls in sta_ori_metrics:
        y = np.array(ampls + [ampls[0]])
        #  Fit cubic spline through all points, used for visualization
        interp_cubsp = interpolate.CubicSpline(x, y, bc_type='periodic')
        yint_cubsp = interp_cubsp(angles_fine)
        # Fit cosine curve to all data points, used to estimate correction angle
        fitted, cov = optimize.curve_fit(lambda _x, _amp, _ph: _amp*np.cos(np.deg2rad(_x + _ph)),
                                         x, y, p0=[0.2, 0.0], bounds=([-1, -180], [1, 180]))
        A_fit, ph_fit = fitted
        y_fitted = A_fit*np.cos(np.deg2rad(angles_fine + ph_fit))
        ph_uncertainty = np.sqrt(cov[1, 1])
        # Estimate correction angle
        angle_max = angles_fine[np.argmax(y_fitted)]
        logger.info('{}: {:2.3f}°, stddev {:2.3f}° (N = {:3d})'.format(
            sta, angle_max, ph_uncertainty, N))
        if save_plot:
            f = plt.figure(figsize=(16,9))
            plt.plot(x, y, 'x', label='Computed P arrival strength')
            plt.plot(angles_fine, yint_cubsp, '--', alpha=0.7, label='Cubic spline fit')
            plt.plot(angles_fine, y_fitted, '--', alpha=0.7, label='Cosine fit')
            plt.grid(linestyle=':', color="#80808080")
            plt.xlabel('Orientation correction (deg)', fontsize=12)
            plt.ylabel('P phase ampl. mean (0-1 sec)', fontsize=12)
            plt.title('{}.{}'.format(ned.network, sta), fontsize=14)
            plt.text(0.9, 0.9, 'N = {}'.format(N), ha='right', va='top',
                     transform=plt.gca().transAxes)
            plt.legend(framealpha=0.7)
            plt.savefig(sta + '_ori.png', dpi=300)
            plt.close()
        # end if
    # end for
# end if


@click.command()
@click.option('--dest-file', type=click.Path(dir_okay=False), help='NYI')
@click.argument('src-h5-event-file', type=click.Path(exists=True, dir_okay=False),
                required=True)
def main(src_h5_event_file, dest_file=None):
    method_wang(src_h5_event_file, dest_file)
    method_wilde_piorko(src_h5_event_file, dest_file, save_plot=True)
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
