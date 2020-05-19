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

2. Abt, D. L., Fischer, K. M., French, S. W., Ford, H. A., Yuan, H., and Romanowicz, B. ( 2010),
North American lithospheric discontinuity structure imaged by Ps and Sp receiver functions,
J. Geophys. Res., 115, B09301, doi:10.1029/2009JB006914.

3. Wilde-Piórko, M., Grycuk, M., Polkowski, M. et al.
On the rotation of teleseismic seismograms based on the receiver function technique.
J Seismol 21, 857-868 (2017). https://doi.org/10.1007/s10950-017-9640-x
"""

import click
import logging

import numpy as np
from scipy import stats
from sklearn.decomposition import PCA

from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.runners import curate_seismograms


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

    if dest_file is not None:
        ned.write(dest_file)
    # end if
# end func


@click.command()
@click.option('--dest-file', type=click.Path(dir_okay=False))
@click.argument('src-h5-event-file', type=click.Path(exists=True, dir_okay=False),
                required=True)
def main(src_h5_event_file, dest_file=None):
    method_wang(src_h5_event_file, dest_file)
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
