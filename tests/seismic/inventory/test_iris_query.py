#!/usr/bin/env python

import os
import requests

from seismic.inventory import iris_query

# TODO: Factor this out into MockIrisResponse class that derives from requests_mock.Mocker,
# and provides access to this 'full' IRIS response and minimal response from test_update_iris_inventory.
iris_data_file = os.path.join(os.path.dirname(__file__), "iris_db_test_chan_GE.xml")


def test_url_formatter():
    """Test IRIS web services query URL generation functions
    """
    url = iris_query.form_channel_request_url()
    assert url == "https://service.iris.edu/fdsnws/station/1/query?net=*&sta=*&cha=*&level=channel&"\
                    "format=xml&includerestricted=false&includecomments=false&nodata=404"
    url = iris_query.form_channel_request_url(netmask='AU', statmask='C*', chanmask='BH*')
    assert url == "https://service.iris.edu/fdsnws/station/1/query?net=AU&sta=C*&cha=BH*&level=channel&"\
                    "format=xml&includerestricted=false&includecomments=false&nodata=404"
    url = iris_query.form_response_request_url(netmask='AU', statmask='C*', chanmask='BH*')
    assert url == "https://service.iris.edu/fdsnws/station/1/query?net=AU&sta=C*&cha=BH*&level=response&"\
                    "format=xml&includerestricted=false&includecomments=false&nodata=404"


def test_set_text_encoding(requests_mock):
    """Test detection of encoding from XML and setting it in the response object.

    :param requests_mock: Automatic requests_mock object passed by pytest
    :type requests_mock: requests_mock
    """
    # Syntax of this function customized to pytest framework.
    # See https://requests-mock.readthedocs.io/en/latest/pytest.html
    default_channel_query = iris_query.form_channel_request_url()
    with open(iris_data_file, 'rb') as f:
        original_response = f.read().decode('utf-8')
    requests_mock.get(default_channel_query, text=original_response)
    response = requests.get(default_channel_query)
    assert response.encoding != "ISO-8859-1"
    iris_query.set_text_encoding(response)
    assert response.encoding == "ISO-8859-1"
