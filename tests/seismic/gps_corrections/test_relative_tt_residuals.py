#!/usr/bin/env python
from __future__ import absolute_import

# pylint: disable=invalid-name, missing-docstring

import numpy as np
import pytest
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


def test_apply_event_quality_filtering(df_picks):
    '''
        self.strict_filtering = DEFAULT_STRICT_FILTERING
        self.min_event_snr = DEFAULT_MIN_EVENT_SNR
        self.cwt_cutoff = DEFAULT_CWT_CUTOFF
        self.slope_cutoff = DEFAULT_SLOPE_CUTOFF
        self.nsigma_cutoff = DEFAULT_NSIGMA_CUTOFF
        self.min_event_mag = DEFAULT_MIN_EVENT_MAG
        self.channel_preference = CHANNEL_PREF
    '''
    filter_options = rttr.FilterOptions()
    # Set filter to loosest settings first, then test adjusting one filter variable at a time.
    filter_options.strict_filtering = False
    filter_options.min_event_mag = 5.0
    filter_options.min_event_snr = 0.0
    filter_options.cwt_cutoff = 0.0
    filter_options.slope_cutoff = 0.0
    filter_options.nsigma_cutoff = 0

    net_sta = {'net': ['AU'], 'sta': ['MEEK']}

    # Test magnitude filter
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 51
    filter_options.min_event_mag = 6.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 49
    filter_options.min_event_mag = 7.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 48
    filter_options.min_event_mag = 8.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 46
    filter_options.min_event_mag = 9.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 44

    # Test SNR filter
    filter_options.min_event_mag = 5.0
    filter_options.min_event_snr = 5.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 47
    filter_options.min_event_snr = 10.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 34
    filter_options.min_event_snr = 15.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 23
    # This encompasses in AU.MEEK, which gets ignored since filtering is non-strict.
    filter_options.min_event_snr = 20.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 12
    # Turn on strict filtering and check filter criteria is now applied to AU.MEEK.
    filter_options.strict_filtering = True
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 11
    filter_options.min_event_snr = 25.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 7
    filter_options.strict_filtering = False
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 8

    # Test CWT cutoff
    filter_options.min_event_snr = 0.0
    filter_options.cwt_cutoff = 5
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 39
    # This encompasses in AU.MEEK, which gets ignored since filtering is non-strict.
    filter_options.cwt_cutoff = 10
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 30
    # Turn on strict filtering and check filter criteria is now applied to AU.MEEK.
    filter_options.strict_filtering = True
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 29
    filter_options.cwt_cutoff = 20
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 22
    filter_options.cwt_cutoff = 50
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 12
    filter_options.strict_filtering = False
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 13

    # Test slope cutoff
    filter_options.cwt_cutoff = 0.0
    filter_options.slope_cutoff = 1.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 43
    filter_options.slope_cutoff = 2.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 28
    # This encompasses in AU.MEEK, which gets ignored since filtering is non-strict.
    filter_options.slope_cutoff = 5.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 16
    # Turn on strict filtering and check filter criteria is now applied to AU.MEEK.
    filter_options.strict_filtering = True
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 15

    # Test nSigma cutoff
    filter_options.strict_filtering = False
    filter_options.slope_cutoff = 0.0
    filter_options.nsigma_cutoff = 5
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 45
    filter_options.nsigma_cutoff = 6
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 43
    filter_options.nsigma_cutoff = 7
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 41

    ## Test combined filtering
    filter_options.strict_filtering = True
    filter_options.min_event_mag = 6.0
    filter_options.min_event_snr = 8.0
    filter_options.cwt_cutoff = 2.0
    filter_options.slope_cutoff = 1.2
    filter_options.nsigma_cutoff = 5
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 5 + 24  # 5 from the zero quality metrics group + 24 from the rest
    mask_zero_metrics = (df_filt[['snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']] == 0).all(axis=1)
    df_zero_mets = df_filt[mask_zero_metrics]
    df_other = df_filt[~mask_zero_metrics]
    assert np.all(df_zero_mets['mag'] >= filter_options.min_event_mag)
    assert np.all((df_other['snr'] >= filter_options.min_event_snr) &
                  (df_other['qualityMeasureCWT'] >= filter_options.cwt_cutoff) &
                  (df_other['qualityMeasureSlope'] >= filter_options.slope_cutoff) &
                  (df_other['nSigma'] >= filter_options.nsigma_cutoff))

    filter_options.slope_cutoff = 3.0
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 18  # 5 from the zero quality metrics group + 13 from the rest

    # Not applying filtering to AU.MEEK
    filter_options.strict_filtering = False
    df_filt = rttr.apply_event_quality_filtering(df_picks, net_sta, filter_options)
    assert len(df_filt) == 19

    # Test passing more than one target station
    net_sta['net'].append('GE')
    net_sta['sta'].append('WSI')
    pytest.raises(ValueError, rttr.apply_event_quality_filtering, df_picks, net_sta, filter_options)


def test_broadcast_ref_residual_per_event(df_picks):
    filter_options = rttr.FilterOptions()
    df_ref = rttr.broadcast_ref_residual_per_event(df_picks, 'AU', 'MEEK', filter_options)
    # AU.MEEK is only there in event4, so only that event will have a non-Nan reference residual.
    expected_residual = 1.029852
    assert np.all(df_ref.loc[df_ref['#eventID'] == 'event4', 'ttResidualRef'] == expected_residual)
    assert np.all(df_ref.loc[df_ref['#eventID'] != 'event4', 'ttResidualRef'].isnull())

    # Repeat with a station that appears in multiple events
    df_ref = rttr.broadcast_ref_residual_per_event(df_picks, 'AU', 'KNA', filter_options)
    expected_event0_residual = -0.134955
    expected_event4_residual = 1.192448
    assert np.allclose(df_ref.loc[df_ref['#eventID'] == 'event0', 'ttResidualRef'], expected_event0_residual)
    assert np.allclose(df_ref.loc[df_ref['#eventID'] == 'event4', 'ttResidualRef'], expected_event4_residual)
    assert np.all(df_ref.loc[(df_ref['#eventID'] != 'event0') & (df_ref['#eventID'] != 'event4'),
                             'ttResidualRef'].isnull())

    # Repeat with NWAO which is common to event1 and event4
    df_ref = rttr.broadcast_ref_residual_per_event(df_picks, 'IR', 'NWAO', filter_options)
    expected_event1_residual = -1.786222
    expected_event4_residual = -0.692545
    assert np.allclose(df_ref.loc[df_ref['#eventID'] == 'event1', 'ttResidualRef'], expected_event1_residual)
    assert np.allclose(df_ref.loc[df_ref['#eventID'] == 'event4', 'ttResidualRef'], expected_event4_residual)
    assert np.all(df_ref.loc[(df_ref['#eventID'] != 'event1') & (df_ref['#eventID'] != 'event4'),
                             'ttResidualRef'].isnull())


def test_analyze_target_relative_to_ref(df_picks):
    filter_options = rttr.FilterOptions()


if __name__ == "__main__":
    # Select explicit test to run.
    import conftest
    picks = conftest._read_picks()
    test_broadcast_ref_residual_per_event(picks)
    pass
