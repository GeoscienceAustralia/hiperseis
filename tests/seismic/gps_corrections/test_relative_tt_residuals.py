#!/usr/bin/env python
from __future__ import absolute_import

# pylint: disable=invalid-name, missing-docstring

import numpy as np
# import pandas as pd

from seismic.gps_corrections import relative_tt_residuals_plotter as rttr


def test_global_filtering(df_picks):
    ## Test filter_limit_channels
    channels_broad = rttr.filter_limit_channels(df_picks, rttr.CHANNEL_PREF_NO_SHZ)
    assert len(channels_broad) == 40
    channels_broad_shz = rttr.filter_limit_channels(df_picks, rttr.CHANNEL_PREF_BALANCED)
    assert len(channels_broad_shz) == 40 + 6
    channels_all = rttr.filter_limit_channels(df_picks, rttr.CHANNEL_PREF_GREEDY)
    assert len(channels_all) == 40 + 6 + 5

    ## Test filter_to_teleseismic
    df_ts = rttr.filter_to_teleseismic(df_picks, 0, 10)
    assert df_ts.empty
    df_ts = rttr.filter_to_teleseismic(df_picks, 0, 180)
    assert len(df_ts) == 51
    assert df_ts['distance'].min() >= 0
    assert df_ts['distance'].max() <= 180
    df_ts = rttr.filter_to_teleseismic(df_picks, 0, 20)
    assert len(df_ts) == 10
    assert df_ts['distance'].min() >= 0
    assert df_ts['distance'].max() <= 20
    df_ts = rttr.filter_to_teleseismic(df_picks, 20, 35)
    assert len(df_ts) == 34
    assert df_ts['distance'].min() >= 20
    assert df_ts['distance'].max() <= 35
    df_ts = rttr.filter_to_teleseismic(df_picks, 35, 90)
    assert len(df_ts) == 7
    assert df_ts['distance'].min() >= 35
    assert df_ts['distance'].max() <= 90
    df_ts = rttr.filter_to_teleseismic(df_picks, 50, 90)
    assert df_ts.empty


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


# Time-boxing the work to just the following remaining functions:
#   apply_event_quality_filtering
#   broadcast_ref_residual_per_event
# Maybe:
#   analyze_target_relative_to_ref
