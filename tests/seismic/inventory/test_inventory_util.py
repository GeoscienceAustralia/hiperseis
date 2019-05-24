#!/usr/bin/env python

import io
from builtins import dict

import numpy as np
import requests
import requests_mock

from obspy import read_inventory

from seismic.inventory.inventory_util import (load_station_xml,
                                              obtain_nominal_instrument_response,
                                              extract_unique_sensors_responses)
from seismic.inventory.iris_query import form_response_request_url


def test_load_station_xml(iris_mocker):
    """Basic test of loading a station xml file
    """
    obspy_inv = load_station_xml(iris_mocker.iris_data_file)
    assert [n.code for n in obspy_inv.networks] == ['GE']
    print(len(obspy_inv.networks))
    print(len(obspy_inv.networks[0].stations))
    assert [s.code for n in obspy_inv.networks for s in n.stations] == [
        'MAHO', 'MAHO', 'MALT', 'MARJ', 'MATE', 'MAUI', 'MAUI', 'MELI', 'MELI', 'MHV', 'MHV', 'MLR', 'MLR', 'MMRI',
        'MNAI', 'MORC', 'MORC', 'MRNI', 'MRNI', 'MTE', 'MTE']
    assert np.all([c.code in ('BHE', 'BHN', 'BHZ') for n in obspy_inv.networks for s in n.stations for c in s.channels])


def test_obtain_nominal_instr_resp(iris_mocker):
    """Test finding a dummy instrument response
    """
    # Exercise the retry logic in case connection is difficult
    query_fields = ('GE', 'M*', '*')
    query = form_response_request_url(*query_fields)
    # First two attempts fail, then third attempt succeeds
    iris_mocker.get(query, [{'exc': requests.exceptions.RequestException},
                            {'exc': requests.exceptions.RequestException},
                            {'text': iris_mocker.get_full_response()}])
    responses = obtain_nominal_instrument_response(*query_fields, req=requests)
    # Check that all the expected responses are there
    assert sorted(responses.keys()) == ['BHE', 'BHN', 'BHZ']


def test_extract_unique_sensor_responses(iris_mocker):
    """Exercise function that finds an instrument response for every channel in an obspy Inventory.
    """
    # Test with larger scale input and perform bulk checks on the results
    iris_mocker.setup_internal_inv(iris_mocker.get_full_response)
    with io.BytesIO(iris_mocker.get_full_response().encode('utf-8')) as buffer_file:
        obspy_inv = read_inventory(buffer_file)
        # Setup response. All queries here handled by the dynamic response generator
        iris_mocker.get(requests_mock.ANY, text=iris_mocker.dynamic_query)
        # Call extract_unique_sensors_responses
        responses = extract_unique_sensors_responses(obspy_inv, requests, test_mode=True)
        expected_response_keys = ['BHE', 'BHN', 'BHZ', 'LAST_RESORT']
        assert list(sorted(responses.keys())) == expected_response_keys
        for instrument in responses.values():
            assert instrument
            assert instrument.response
