#!/usr/bin/env python
"""Utility functions and constants shared amongst inventory management modules.
"""

from __future__ import print_function

import sys
from collections import namedtuple, defaultdict
import time

from obspy import read_inventory

from seismic.inventory.iris_query import set_text_encoding, form_response_request_url

PY2 = sys.version_info[0] < 3  # pylint: disable=invalid-name

if PY2:
    import cStringIO as sio  # pylint: disable=import-error
    import io as bio
else:
    import io as sio  # pylint: disable=ungrouped-imports
    bio = sio
    basestring = str  # pylint: disable=invalid-name


NOMINAL_EARTH_RADIUS_KM = 6378.1370  # pylint: disable=invalid-name

SORT_ORDERING = ['NetworkCode', 'StationCode', 'StationStart', 'StationEnd',  # pylint: disable=invalid-name
                 'ChannelCode', 'ChannelStart', 'ChannelEnd']

# Bundled container for related sensor and response.
Instrument = namedtuple("Instrument", ['sensor', 'response'])


def load_station_xml(inventory_file):
    """Load a stationxml file

    :param inventory_file: [description]
    :type inventory_file: str or path
    """
    if PY2:
        import io
        with io.open(inventory_file, mode='r', encoding='utf-8') as f:
            obspy_inv = read_inventory(f)
    else:
        with open(inventory_file, mode='r', encoding='utf-8') as f:
            obspy_inv = read_inventory(f)

    return obspy_inv


def obtain_nominal_instrument_response(netcode, statcode, chcode, req):
    """
    For given network, station and channel code, find suitable response(s) in IRIS database and
    return as dict of obspy instrument responses.

    :param netcode: Network code (may include wildcards)
    :type netcode: str
    :param statcode: Station code (may include wildcards)
    :type statcode: str
    :param chcode: Channel code (may include wildcards)
    :type chcode: str
    :param req: Request object to use for URI query
    :type req: Object conforming to interface of 'requests' library
    :return: Dictionary of instrument responses from IRIS for given network(s), station(s) and channel(s).
    :rtype: dict of {str, Instrument(obspy.core.inventory.util.Equipment, obspy.core.inventory.response.Response)}
    """
    from obspy.core.util.obspy_types import FloatWithUncertaintiesAndUnit

    query_url = form_response_request_url(netcode, statcode, chcode)
    tries = 10
    while tries > 0:
        try:
            tries -= 1
            response_xml = req.get(query_url)
            first_line = sio.StringIO(response_xml.text).readline().rstrip()
            assert 'Error 404' not in first_line
            break
        except req.exceptions.RequestException as e:  # pylint: disable=unused-variable
            time.sleep(1)
    assert tries > 0
    set_text_encoding(response_xml, quiet=True)
    # This line decodes when .text attribute is extracted, then encodes to utf-8
    obspy_input = bio.BytesIO(response_xml.text.encode('utf-8'))
    try:
        channel_data = read_inventory(obspy_input)
        responses = {cha.code: Instrument(cha.sensor, cha.response) for net in channel_data.networks
                     for sta in net.stations for cha in sta.channels
                     if cha.code is not None and cha.response is not None}
        # Make responses valid for Seiscomp3
        for inst in responses.values():
            assert inst.response
            for rs in inst.response.response_stages:
                if rs.decimation_delay is None:
                    rs.decimation_delay = FloatWithUncertaintiesAndUnit(0)
                if rs.decimation_correction is None:
                    rs.decimation_correction = FloatWithUncertaintiesAndUnit(0)
    except ValueError:
        responses = {}
    return responses


def extract_unique_sensors_responses(inv, req, show_progress=True, blacklisted_networks=None, test_mode=False):
    """
    For the channel codes in the given inventory, determine a nominal instrument response suitable
    for that code. Note that no attempt is made here to determine an ACTUAL correct response for
    a given network and station. The only requirement here is to populate a plausible, non-empty
    response for a given channel code, to placate Seiscomp3 which requires that an instrument
    response always be present.

    :param inv: Seismic station inventory
    :type inv: obspy.Inventory
    :param req: Request object to use for URI query
    :type req: Object conforming to interface of 'requests' library
    :return: Python dict of (obspy.core.inventory.util.Equipment, obspy.core.inventory.response.Response)
        indexed by str representing channel code
    :rtype: {str: Instrument(obspy.core.inventory.util.Equipment, obspy.core.inventory.response.Response) }
        where Instrument is a namedtuple("Instrument", ['sensor', 'response'])
    """
    if blacklisted_networks is None:
        blacklisted_networks = []

    # Create like this so if later indexed with invalid key, returns None instead of exception.
    nominal_instruments = defaultdict(lambda: None)
    if test_mode:
        reference_networks = (('GE', '*', '*HZ'),)  # trailing comma required to make it a tuple
    else:
        reference_networks = (('GE', '*', '*'), ('IU', '*', '*'), ('BK', '*', '*'))
    print("Preparing common instrument response database from networks {} "
          "(this may take a while)...".format(reference_networks))
    for query in reference_networks:
        print("  querying {} as {}.{}.{}".format(query[0], *query))
        nominal_instruments.update(obtain_nominal_instrument_response(*query, req=req))

    if show_progress:
        import tqdm
        num_entries = sum(len(sta.channels) for net in inv.networks for sta in net.stations)
        pbar = tqdm.tqdm(total=num_entries, ascii=True, desc="Finding additional instrument responses")
        std_print = tqdm.tqdm.write
    else:
        std_print = print

    failed_codes = set()
    for net in inv.networks:
        if net.code in blacklisted_networks:
            if show_progress:
                for sta in net.stations:
                    pbar.update(len(sta.channels))
            continue
        for sta in net.stations:
            if show_progress:
                pbar.update(len(sta.channels))
            for cha in sta.channels:
                if cha.code is None or cha.code in nominal_instruments:
                    continue
                assert isinstance(cha.code, basestring), type(cha.code)
                # For each channel code, obtain a nominal instrument response by IRIS query.
                if cha.code not in nominal_instruments:
                    response = obtain_nominal_instrument_response(net.code, sta.code, cha.code, req)
                    nominal_instruments.update(response)
                    if cha.code in response:
                        std_print("Found nominal instrument response for channel code {} in "
                                  "{}.{}".format(cha.code, net.code, sta.code))
                    else:
                        std_print("Failed to acquire instrument response for channel code {} "
                                  "in {}.{}".format(cha.code, net.code, sta.code))
                        failed_codes.add(cha.code)
    if show_progress:
        pbar.close()

    # Flag a fallback response of last resort for channel codes for which we find no valid response from IRIS.
    last_resort_responses = ['BHE', 'BHN', 'BHZ', 'HHE', 'HHN', 'HHZ', 'LHE', 'LHN', 'LHZ', 'SHE', 'SHN', 'SHZ', None]
    for last_resort_response in last_resort_responses:
        assert last_resort_response is not None
        if last_resort_response in nominal_instruments:
            nominal_instruments['LAST_RESORT'] = nominal_instruments[last_resort_response]
            print("Last resort response code: {}".format(last_resort_response))
            assert nominal_instruments['LAST_RESORT'] is not None
            break

    # Report on channel codes for which no response could be found
    failed_codes = sorted(list(failed_codes - set(nominal_instruments.keys())))
    if failed_codes:
        print("WARNING: No instrument response could be determined for these channel codes:\n{}".format(failed_codes))
        print("         {} selected as response of last resort.".format(last_resort_response))

    return nominal_instruments
