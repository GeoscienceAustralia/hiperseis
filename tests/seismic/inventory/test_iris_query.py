#!/usr/bin/env python

import requests

from seismic.inventory import iris_query


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


def test_set_text_encoding(iris_mocker):
    """Test detection of encoding from XML and setting it in the response object.

    :param iris_mocker: Automatic fixture object passed by pytest which duck types requests
    :type iris_mocker: tests.mocks.inventory.mock_fdsn_xml.MockIrisResponse
    """
    default_channel_query = iris_query.form_channel_request_url()
    iris_mocker.get(default_channel_query, text=iris_mocker.get_full_response())
    response = requests.get(default_channel_query)
    assert response.encoding != "ISO-8859-1"
    iris_query.set_text_encoding(response)
    assert response.encoding == "ISO-8859-1"
