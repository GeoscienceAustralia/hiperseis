#!/usr/bin/env python
"""
Analyze a data set of seismic arrival events on a per-station basis and try to
detect and estimate any station orientation error.

In future, consider moving this script to the `inventory` module and applying
corrections to the station inventory xml (to the azimuth tag).

Reference:
Wilde-Piórko, M., Grycuk, M., Polkowski, M. et al. On the rotation of teleseismic seismograms based on the receiver function technique.
J Seismol 21, 857-868 (2017). https://doi.org/10.1007/s10950-017-9640-x
"""

import logging

import numpy as np
from numpy.linalg import eig
from sklearn.decomposition import PCA

from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.runners import curate_seismograms


def main(src_h5_event_file, network=None, station=None, dest_file=None):

    ned = NetworkEventDataset(src_h5_event_file, network, station)
    evids_orig = set([evid for _, evid, _ in ned])

    logger = logging.getLogger(__name__)
    curation_opts = {
        "min_snr": 2.0,
        "max_raw_amplitude": 20000.0,
        "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
    }
    curate_seismograms(ned, curation_opts, logger)
    evids_to_keep = set([evid for _, evid, _ in ned])

    evids_discarded = evids_orig - evids_to_keep
    print('Discarded {}/{} events'.format(len(evids_discarded), len(evids_orig)))

    # Reload and apply event filter
    ned = NetworkEventDataset(src_h5_event_file, network, station)
    ned.curate(lambda _0, evid, _2: evid in evids_to_keep)

    # Methodology from M. Wilde-Piórko et al. paper (ZRT method):
    # 1.

    resids = []
    for sta, evid, stream in ned:
        s = stream.copy()
        s.detrend(type='linear')
        s.taper(0.05)
        s.resample(10.0, no_filter=False)
        tr_n = s.select(component='N')[0]
        tr_n.trim(tr_n.stats.onset - 1.0, tr_n.stats.onset + 3.0)
        if tr_n.stats.event_magnitude < 5.5:
            continue
        # print('Mag:', tr_n.stats.event_magnitude)
        tr_e = s.select(component='E')[0]
        tr_e.trim(tr_e.stats.onset - 1.0, tr_e.stats.onset + 3.0)
        if len(tr_n) != len(tr_e):
            continue
        event_baz = tr_n.stats.back_azimuth
        # print('baz:', event_baz)
        data = np.array([tr_n.data, tr_e.data])

        # -----------
        # EXPERIMENT: NOT Wilde-Piórko method.
        # sklearn method
        sk_pca = PCA()
        sk_pca.fit(data.T)
        sk_evals = sk_pca.singular_values_
        sk_evecs = sk_pca.components_
        if sk_pca.explained_variance_ratio_[0] < 0.80:
            continue
        sk_baz = np.rad2deg(np.arctan2(sk_evecs[0,1], sk_evecs[0,0]))
        sk_residual = sk_baz - event_baz
        while sk_residual < -90:
            sk_baz += 180
            sk_residual = sk_baz - event_baz
        while sk_residual > 90:
            sk_baz -= 180
            sk_residual = sk_baz - event_baz
        resids.append(sk_residual)
        # -----------

        pass
    # end for
    resids = np.array(sorted(resids))
    mean = np.mean(resids)
    stddev = np.std(resids)
    stderr = stddev/np.sqrt(len(resids))  # standard error of the mean
    print('{}:  {:.4f}° ± {:.4f}°, stddev {:.4f}° (N = {})'.format(station, mean, stderr, stddev, len(resids)))
    print(resids)

    if dest_file is not None:
        ned.write(dest_file)
    # end if
# end func


if __name__ == '__main__':
    main('/g/data1a/ha3/am7399/shared/7X_RF_analysis/7X_event_waveforms_for_rf_20090616T034200-20110401T231849_rev2.h5',
         '7X', 'MA01')
# end if
