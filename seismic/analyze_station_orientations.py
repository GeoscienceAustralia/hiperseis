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

import numpy as np
from numpy.linalg import eig
from sklearn.decomposition import PCA

from seismic.network_event_dataset import NetworkEventDataset


def main(src_h5_event_file, network=None, station=None, dest_file=None):
    ned = NetworkEventDataset(src_h5_event_file, network, station)

    # Methodology from M. Wilde-Piórko et al. paper (ZRT method):
    # 1.

    for sta, evid, stream in ned:
        s = stream.copy()
        s.detrend()
        s.taper(0.05)
        s.resample(20.0, no_filter=False)
        tr_n = s.select(component='N')[0]
        tr_n.trim(tr_n.stats.onset - 1.0, tr_n.stats.onset + 2.0)
        tr_e = s.select(component='E')[0]
        tr_e.trim(tr_e.stats.onset - 1.0, tr_e.stats.onset + 2.0)
        if len(tr_n) != len(tr_e):
            continue
        event_baz = tr_n.stats.back_azimuth
        data = np.array([tr_n.data, tr_e.data])

        #-----------
        # EXPERIMENT: NOT Wilde-Piórko method.
        # Might still be applicable to R-T receiver functions to determine
        # numpy direct method
        cov = np.cov(data)
        evals, evecs = eig(cov)
        # np_residual = np_baz - event_baz

        # sklearn method
        sk_pca = PCA()
        sk_pca.fit(data.T)
        sk_evals = sk_pca.singular_values_
        sk_evecs = sk_pca.components_
        print(sk_pca.explained_variance_ratio_)
        # sk_residual = sk_baz - event_baz
        #-----------

        pass
    # end for

    if dest_file is not None:
        ned.write(dest_file)
    # end if
# end func


if __name__ == '__main__':
    main('/g/data1a/ha3/am7399/shared/7X_RF_analysis/7X_event_waveforms_for_rf_20090616T034200-20110401T231849_rev2.h5',
         '7X', 'MA01')
# end if
