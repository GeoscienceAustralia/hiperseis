#!/usr/bin/env python
"""
Helper functions for making and managing web queries to IRIS web service.
"""

import re
import sys
if sys.version_info[0] < 3:
    import cStringIO as sio  # pylint: disable=import-error
else:
    import io as sio


def form_channel_request_url(netmask="*", statmask="*", chanmask="*"):
    """
    Form request URL to download station inventory in stationxml format, down to channel level,
    with the given filters applied to network codes, station codes and channel codes.

    :param netmask: Pattern of network codes to match, comma separated with wildcards, defaults to "*"
    :param netmask: str, optional
    :param statmask: Pattern of station codes to match, comma separated with wildcards, defaults to "*"
    :param statmask: str, optional
    :param chanmask: Pattern of channel codes to match, comma separated with wildcards, defaults to "*"
    :param chanmask: str, optional
    :return: Fully formed URL to perform IRIS query and get back FDSN station XML result.
    :rtype: str
    """
    # Hardwired to exclude restricted channels and exclude comments to reduce file size.
    return "https://service.iris.edu/fdsnws/station/1/query?net=" + netmask + \
           "&sta=" + statmask + \
           "&cha=" + chanmask + \
           "&level=channel&format=xml&includerestricted=false&includecomments=false&nodata=404"


def form_response_request_url(netmask, statmask, chanmask):
    """
    Form request URL to download station inventory in stationxml format, down to response level,
    for the given network, station and channel codes.

    :param netmask: Pattern of network codes to match, comma separated with wildcards
    :param netmask: str, optional
    :param statmask: Pattern of station codes to match, comma separated with wildcards
    :param statmask: str, optional
    :param chanmask: Pattern of channel codes to match, comma separated with wildcards
    :param chanmask: str, optional
    :return: Fully formed URL to perform IRIS query and get back FDSN station XML result.
    :rtype: str
    """
    # Hardwired to exclude restricted channels and exclude comments to reduce file size.
    return "https://service.iris.edu/fdsnws/station/1/query?net=" + netmask + \
           "&sta=" + statmask + \
           "&cha=" + chanmask + \
           "&level=response&format=xml&includerestricted=false&includecomments=false&nodata=404"


def set_text_encoding(resp, quiet=False):
    """
    For the given response object, set its encoding from the contents of the text returned from server.

    :param resp: Query response object returned by response.get()
    :type resp: requests.Response
    """
    encoding_pattern = r"^<\?xml .* encoding=[\"'](.+)[\"']"
    matcher = re.compile(encoding_pattern)
    first_line = sio.StringIO(resp.text).readline().rstrip()
    match = matcher.search(first_line)
    assert match, "Encoding match missing in response\n{}".format(resp.text[0:1000])
    encoding = match.group(1)
    if not quiet:
        print("Detected text encoding {}".format(encoding))
    resp.encoding = encoding
    assert resp.encoding is not None
