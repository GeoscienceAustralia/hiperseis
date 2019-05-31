#!/usr/bin/env python
from __future__ import absolute_import

import numpy as np
# import pandas as pd

from seismic.gps_corrections import relative_tt_residuals_plotter as rttr


def test_global_filtering():
    pass


def test_get_iris_station_codes(iris_stntxt_file):
    df_au = rttr.get_iris_station_codes(iris_stntxt_file, 'AU')
    assert list(df_au.index) == ['ARMA', 'EIDS', 'KMBL', 'MILA', 'NRFK', 'RABL', 'XMI']
    assert np.allclose(df_au['lat'], [-30.4198, -25.369101, -31.366899, -37.054699, -29.040001, -4.127967, -10.4498],
                       rtol=1.0e-6)
    assert np.allclose(df_au['lon'], [151.628006, 151.081696, 121.882103, 149.154999, 167.962997, 152.108765,
                                      105.688950], rtol=1.0e-6)
    df_ge = rttr.get_iris_station_codes(iris_stntxt_file, 'GE')
    assert list(df_ge.index) == ['BKNI', 'DAG', 'GSI', 'KAAM', 'KWP', 'MARJ', 'MORC', 'TRTE']
    assert np.allclose(df_ge['lat'], [0.3262, 76.771301, 1.3039, 0.49264, 49.630501, 32.522598, 49.83105, 58.378601],
                       rtol=1.0e-6)
    assert np.allclose(df_ge['lon'], [101.039597, -18.655001, 97.5755, 72.994858, 22.7078, 20.8776, 17.577573,
                                      26.720501], rtol=1.0e-6)
