#!/usr/bin/env python
from __future__ import absolute_import

# pylint: disable=invalid-name, missing-docstring

import numpy as np
import pytest

from obspy import UTCDateTime

from seismic.gps_corrections import picks_reader_utils as pru


def test_read_picks_ensemble(df_picks):
    assert df_picks is not None
    assert not df_picks.empty
    assert sorted(list(pru.PICKS_TABLE_SCHEMA.keys())) == sorted(list(df_picks.columns))
    assert len(df_picks['#eventID'].unique()) == 5


def test_get_network_stations(df_picks):
    # Test result for non-existent codes
    net_stations = pru.get_network_stations(df_picks, 'NONEXIST')
    assert not net_stations
    net_stations = pru.get_network_stations(df_picks, 'AU')
    assert net_stations == ['AS01', 'FITZ', 'KNA', 'MEEK', 'QIS', 'RMQ']
    net_stations = pru.get_network_stations(df_picks, 'G')
    assert net_stations == ['CAN']
    net_stations = pru.get_network_stations(df_picks, 'GE')
    assert net_stations == ['ARMA', 'CHN', 'CNB', 'FUQ', 'KAPI', 'KMBL', 'MEEK', 'NAH1',
                            'NGO1', 'NJ2', 'NZJ', 'SSE', 'TIA', 'WSI']


def test_get_network_location_mean(df_picks):
    # Test result for non-existent codes
    g_mean = pru.get_network_location_mean(df_picks, 'NONEXIST')
    assert np.all(np.isnan(g_mean))

    # Test getting mean of network with a single station record, should just match exact lat-long of that station.
    g_mean = pru.get_network_location_mean(df_picks, 'G')
    assert g_mean == (-35.3203, 148.999)

    # Test that changing case doesn't affect result.
    g_mean_2 = pru.get_network_location_mean(df_picks, 'g')
    assert g_mean_2 == g_mean

    # Test mean of network for which a station appears multiple times. Should only count once towards the mean location.
    au_mean = pru.get_network_location_mean(df_picks, 'AU')
    assert np.allclose(au_mean, (-21.86503, 132.556))


def test_get_network_date_range(df_picks):
    # Test result for non-existent codes
    date_range = pru.get_network_date_range(df_picks, 'NONEXIST')
    assert date_range == (None, None)

    # Get dates of picks for network with a single record.
    date_range = pru.get_network_date_range(df_picks, 'G')
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(777031773.378))

    # Get dates for other networks.
    date_range = pru.get_network_date_range(df_picks, 'GE')
    assert date_range == (UTCDateTime(742222992.735), UTCDateTime(1197051882.89))
    date_range = pru.get_network_date_range(df_picks, 'AU')
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(1175528404.0))

    # Repeat with lowercase netcode, results should be the same
    date_range_2 = pru.get_network_date_range(df_picks, 'au')
    assert date_range_2 == date_range


def test_get_station_date_range(df_picks):
    # Test result for non-existent codes
    date_range = pru.get_station_date_range(df_picks, 'NONEXIST', 'WHAT')
    assert date_range == (None, None)
    date_range = pru.get_station_date_range(df_picks, 'AU', 'NONEXIST')
    assert date_range == (None, None)

    # Get dates of picks for network.station with a single record
    date_range = pru.get_station_date_range(df_picks, 'BL', 'ATDB')
    assert date_range == (UTCDateTime(829873172.636), UTCDateTime(829873172.636))

    # Get dates of picks spanning more than one event
    date_range = pru.get_station_date_range(df_picks, 'AU', 'KNA')
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(1175528404.0))
    date_range_2 = pru.get_station_date_range(df_picks, 'au', 'KNA')
    date_range_3 = pru.get_station_date_range(df_picks, 'AU', 'kna')
    assert date_range == date_range_2 == date_range_3


def test_get_overlapping_date_range(df_picks):
    # The expected results of queries here are determined by inspection of the test data file.

    # Pass in empty query and make sure we get empty result
    empty = {}
    no_query = {'net': [], 'sta': []}
    date_range = pru.get_overlapping_date_range(df_picks, empty, empty)
    assert date_range == (None, None)
    date_range = pru.get_overlapping_date_range(df_picks, no_query, empty)
    assert date_range == (None, None)
    date_range = pru.get_overlapping_date_range(df_picks, empty, no_query)
    assert date_range == (None, None)
    date_range = pru.get_overlapping_date_range(df_picks, no_query, no_query)
    assert date_range == (None, None)

    # Query with invalid dict keys
    invalid = {'ne': ['AU'], 'sta': ['KNA']}
    pytest.raises(KeyError, pru.get_overlapping_date_range, df_picks, invalid, invalid)
    invalid = {'net': ['AU'], 'st': ['KNA']}
    pytest.raises(KeyError, pru.get_overlapping_date_range, df_picks, invalid, invalid)

    # Invalid query due to mismatching net.sta pairs
    invalid = {'net': ['AU', 'GE'], 'sta': ['KNA']}
    pytest.raises(ValueError, pru.get_overlapping_date_range, df_picks, invalid, invalid)
    invalid = {'net': ['AU'], 'sta': ['KNA', 'RMQ']}
    pytest.raises(ValueError, pru.get_overlapping_date_range, df_picks, invalid, invalid)

    # Query with invalid codes should return empty result
    non_existent = {'net': ['NONEXIST'], 'sta': ['WHAT']}
    set1 = {'net': ['AU'], 'sta': ['RMQ']}
    date_range = pru.get_overlapping_date_range(df_picks, non_existent, non_existent)
    assert date_range == (None, None)
    date_range = pru.get_overlapping_date_range(df_picks, non_existent, set1)
    assert date_range == (None, None)
    date_range = pru.get_overlapping_date_range(df_picks, set1, non_existent)
    assert date_range == (None, None)

    # Query degenerate cases querying what events are common to same station.
    date_range = pru.get_overlapping_date_range(df_picks, set1, set1)
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(777031773.378))
    set2 = {'net': ['AU'], 'sta': ['KNA']}
    date_range = pru.get_overlapping_date_range(df_picks, set2, set2)
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(1175528404.0))

    # Check that when the network code exists and the station code exists, but not in the
    # same record, then it is not generating a matched record.
    # First case, crossed within the same event.
    set_crossed = {'net': ['AU'], 'sta': ['FORT']}
    date_range = pru.get_overlapping_date_range(df_picks, set_crossed, set_crossed)
    assert date_range == (None, None)
    # Second case, crossed between different events.
    set_crossed = {'net': ['AU'], 'sta': ['HIA']}
    date_range = pru.get_overlapping_date_range(df_picks, set_crossed, set_crossed)
    assert date_range == (None, None)
    # Third case, where the second station matches the first network and the first station matches the second network,
    # within the same event. This should NOT match, since the sequence of station and network codes are supposed to
    # line up.
    set_crossed = {'net': ['AU', 'IR'], 'sta': ['FORT', 'KNA']}
    date_range = pru.get_overlapping_date_range(df_picks, set_crossed, set_crossed)
    assert date_range == (None, None)
    # Same as third case, but across two different events.
    set_crossed = {'net': ['AU', 'CD'], 'sta': ['HIA', 'KNA']}
    date_range = pru.get_overlapping_date_range(df_picks, set_crossed, set_crossed)
    assert date_range == (None, None)

    # Query set of stations with no events in common
    set_A = {'net': ['IR'], 'sta': ['MBWA']}
    set_B = {'net': ['CD'], 'sta': ['HIA']}
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (None, None)

    ## Start expanding the query in increments and make sure expected results come back.
    ## Do spot checks on the symmetry of the results (i.e. same results if order of arguments is swapped)
    # Add station GE.WSI to set_B, which makes event0 the only event in common.
    set_B['net'].append('GE')   # --> set_B = {'net': ['CD', 'GE'], 'sta': ['HIA']}
    set_B['sta'].append('WSI')  # --> set_B = {'net': ['CD', 'GE'], 'sta': ['HIA', 'WSI']}
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(777031773.378))
    # Add station IR.SKR to set_A, which brings in event3 as another common event.
    set_A['net'].append('IR')   # --> set_A = {'net': ['IR', 'IR'], 'sta': ['MBWA']}
    set_A['sta'].append('SKR')  # --> set_A = {'net': ['IR', 'IR'], 'sta': ['MBWA', 'SKR']}
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(1197051882.89))
    date_range_2 = pru.get_overlapping_date_range(df_picks, set_B, set_A)
    assert date_range_2 == date_range
    # Add station XC.QRZ to set_A. This does not change the events in common with set_A and set_B.
    set_A['net'].append('XC')   # --> set_A = {'net': ['IR', 'IR', 'XC'], 'sta': ['MBWA', 'SKR']}
    set_A['sta'].append('QRZ')  # --> set_A = {'net': ['IR', 'IR', 'XC'], 'sta': ['MBWA', 'SKR', 'QRZ']}
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (UTCDateTime(777031773.378), UTCDateTime(1197051882.89))
    # Add station GE.CNB to set_B, which brings in event1 as another common event.
    set_B['net'].append('GE')   # --> set_B = {'net': ['CD', 'GE', 'GE'], 'sta': ['HIA', 'WSI']}
    set_B['sta'].append('CNB')  # --> set_B = {'net': ['CD', 'GE', 'GE'], 'sta': ['HIA', 'WSI', 'CNB']}
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (UTCDateTime(742222992.735), UTCDateTime(1197051882.89))
    date_range_2 = pru.get_overlapping_date_range(df_picks, set_B, set_A)
    assert date_range_2 == date_range
    # Adding new common event with date between current range results in no change in query result.
    # Add IR.BAO to both sets to bring in event2 with timestamp 829873172.636.
    set_A['net'].append('IR')
    set_A['sta'].append('BAO')
    set_B['net'].append('IR')
    set_B['sta'].append('BAO')
    date_range = pru.get_overlapping_date_range(df_picks, set_A, set_B)
    assert date_range == (UTCDateTime(742222992.735), UTCDateTime(1197051882.89))


def test_generate_large_events_catalog(df_picks):
    event_cat = pru.generate_large_events_catalog(df_picks, min_magnitude=8.0, label_historical_events=False)
    # No result, since min record count is too big for test dataset.
    assert event_cat.empty
    # Repeat with min record count set to 1.
    event_cat = pru.generate_large_events_catalog(df_picks, min_magnitude=8.0, min_record_count=1,
                                                  label_historical_events=False)
    assert len(event_cat) == 1
    event_cat = pru.generate_large_events_catalog(df_picks, min_magnitude=7.0, min_record_count=1,
                                                  label_historical_events=False)
    assert len(event_cat) == 3
    event_cat = pru.generate_large_events_catalog(df_picks, min_magnitude=7.0, min_record_count=1)
    assert len(event_cat) > 3


def test_compute_matching_network_mask(df_picks):
    set_crossed = {'net': ['AU'], 'sta': ['FORT']}
    mask = pru.compute_matching_network_mask(df_picks, set_crossed)
    assert len(mask) == len(df_picks)
    assert not np.any(mask)

    set_crossed = {'net': ['AU', 'IR'], 'sta': ['FORT', 'KNA']}
    mask = pru.compute_matching_network_mask(df_picks, set_crossed)
    assert not np.any(mask)

    set_aligned = {'net': ['AU'], 'sta': ['MEEK']}
    mask = pru.compute_matching_network_mask(df_picks, set_aligned)
    assert np.sum(mask) == 1

    set_mismatched = {'net': ['AU', 'GE'], 'sta': ['MEEK']}
    pytest.raises(ValueError, pru.compute_matching_network_mask, df_picks, set_mismatched)

    set_aligned = {'net': ['AU', 'GE'], 'sta': ['MEEK', 'MEEK']}
    mask = pru.compute_matching_network_mask(df_picks, set_aligned)
    assert np.sum(mask) == 3


if __name__ == "__main__":
    # Select explicit test to run.
    # import conftest
    # picks = conftest._read_picks()
    # test_get_overlapping_date_range(picks)
    pass
