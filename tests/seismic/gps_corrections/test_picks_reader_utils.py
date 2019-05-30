#!/usr/bin/env python
from __future__ import absolute_import

# import os
import numpy as np
# import pandas as pd

from seismic.gps_corrections import picks_reader_utils as pru


def test_read_picks_ensemble(df_picks):
    assert df_picks is not None
    assert not df_picks.empty
    assert sorted(list(pru.PICKS_TABLE_SCHEMA.keys())) == sorted(list(df_picks.columns))
    assert len(df_picks['#eventID'].unique()) == 5


def test_get_network_stations(df_picks):
    net_stations = pru.get_network_stations(df_picks, 'AU')
    assert net_stations == ['AS01', 'FITZ', 'KNA', 'QIS', 'RMQ']
    net_stations = pru.get_network_stations(df_picks, 'G')
    assert net_stations == ['CAN']
    net_stations = pru.get_network_stations(df_picks, 'GE')
    assert net_stations == ['ARMA', 'CHN', 'CNB', 'FUQ', 'KAPI', 'KMBL', 'MEEK', 'NAH1',
                            'NGO1', 'NJ2', 'NZJ', 'PMG', 'SSE', 'TIA', 'WSI']


def test_get_network_location_mean(df_picks):
    # Test getting mean of network with a single station record, should just match exact lat-long of that station.
    g_mean = pru.get_network_location_mean(df_picks, 'G')
    assert g_mean == (-35.3203, 148.999)
    # Test that changing case doesn't affect result.
    g_mean_2 = pru.get_network_location_mean(df_picks, 'g')
    assert g_mean_2 == g_mean
    # Test mean of network for which a station appears multiple times. Should only count once towards the mean location.
    au_mean = pru.get_network_location_mean(df_picks, 'AU')
    assert np.allclose(au_mean, (-20.91066, 135.3442))


if __name__ == "__main__":
    # Select explicit test to run
    # import conftest
    # test_read_picks_ensemble(conftest._read_picks())
    pass
