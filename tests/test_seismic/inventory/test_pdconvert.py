#!/usr/bin/env python

import os
import io

import numpy as np
import requests
from obspy import read_inventory

from seismic.inventory.pdconvert import inventory_to_dataframe, dataframe_to_fdsn_station_xml
from seismic.inventory.table_format import TABLE_COLUMNS
from seismic.inventory.inventory_util import extract_unique_sensors_responses


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


def test_dataframe_to_inventory(iris_mocker, tmp_path):
    """Test conversion of a Pandas dataframe representation conforming to TABLE_SCHEMA back to station XML file
    """
    outfile = str(tmp_path / 'df_to_inv.xml')
    assert not os.path.exists(outfile)
    # Generate test dataframe by using function inventory_to_dataframe().
    with io.BytesIO(iris_mocker.get_minimal_response().encode('utf-8')) as buffer_file:
        obspy_inv = read_inventory(buffer_file)
        test_df = inventory_to_dataframe(obspy_inv)
        expected_query = "https://service.iris.edu/fdsnws/station/1/query?net=GE&sta=*&cha=*HZ&level=response" \
                         "&format=xml&includerestricted=false&includecomments=false&nodata=404"
        iris_mocker.get(expected_query, text=iris_mocker.get_minimal_response())
        instruments = extract_unique_sensors_responses(obspy_inv, requests, show_progress=True, test_mode=True)
        dataframe_to_fdsn_station_xml(test_df, instruments, outfile, show_progress=True)
    assert os.path.exists(outfile)
    # Test that the expected network/station/channel codes are there and that each channel
    # has an instrument response.
    with open(outfile, 'r') as f:
        contents = f.read()
        assert contents.count('<Network code="GE"') == 1
        assert contents.count('<Station code="MAHO"') == 1
        assert contents.count('<Channel code="BHZ"') == 1
        assert contents.count('<Response>') == 1
    os.remove(outfile)
