#!/usr/bin/env python

import io

import numpy as np
from obspy import read_inventory

from seismic.inventory.pdconvert import inventory_to_dataframe
from seismic.inventory.table_format import TABLE_COLUMNS


def test_inventory_to_dataframe(iris_mocker):
    """Test conversion of an obspy Inventory object to dataframe representation with
       schema custom to seismic.inventory module.
    """
    # Test with minimal input
    with io.BytesIO(iris_mocker.get_minimal_response().encode('utf-8')) as buffer_file:
        obspy_inv = read_inventory(buffer_file)
        test_df = inventory_to_dataframe(obspy_inv)
        num_channel_records = iris_mocker.get_minimal_response().count(u"<Channel code=")
        assert len(test_df) == num_channel_records
        assert np.all(test_df.columns.values == list(TABLE_COLUMNS))
        assert test_df.iloc[0]['NetworkCode'] == 'GE'
        assert test_df.iloc[0]['StationCode'] == 'MAHO'
        assert np.isclose(test_df.iloc[0]['Latitude'], 39.895901)
        assert np.isclose(test_df.iloc[0]['Longitude'], 4.2665)
        assert np.isclose(test_df.iloc[0]['Elevation'], 15.0)
        assert test_df.iloc[0]['StationStart'] == np.datetime64("1999-04-10T00:00:00")
        assert test_df.iloc[0]['StationEnd'] == np.datetime64("2001-02-13T00:00:00")
        assert test_df.iloc[0]['ChannelCode'] == 'BHZ'
        assert test_df.iloc[0]['ChannelStart'] == np.datetime64("1999-04-10T00:00:00")
        assert test_df.iloc[0]['ChannelEnd'] == np.datetime64("2001-02-13T00:00:00")

    # Test with larger scale input and perform bulk checks on the results
    with io.BytesIO(iris_mocker.get_full_response().encode('utf-8')) as buffer_file:
        obspy_inv = read_inventory(buffer_file)
        test_df = inventory_to_dataframe(obspy_inv)
        num_channel_records = iris_mocker.get_full_response().count(u"<Channel code=")
        assert len(test_df) == num_channel_records
        assert np.all(test_df.columns.values == list(TABLE_COLUMNS))
