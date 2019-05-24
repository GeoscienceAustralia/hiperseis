#!/usr/bin/env python

import os
import sys
import re

# import requests
import requests_mock

PY2 = sys.version_info[0] < 3  # pylint: disable=invalid-name

if PY2:
    import cStringIO as sio  # pylint: disable=import-error
    import io as bio
else:
    import io as sio  # pylint: disable=ungrouped-imports
    bio = sio

from obspy import read_inventory


class MockIrisResponse(requests_mock.Mocker):
    """Class to conveniently generate consistent mock request responses for FDSN station xml
       data across the test suites.
    """

    def __init__(self):
        super(MockIrisResponse, self).__init__()

        self.iris_data_file = os.path.join(os.path.dirname(__file__), "iris_db_test_chan_GE.xml")
        with open(self.iris_data_file, 'rb') as f:
            self.full_test_response = f.read().decode('utf-8')

        invalid_response_file = os.path.join(os.path.dirname(__file__), "iris_db_invalid_response.xml")
        with open(invalid_response_file, 'rb') as f:
            self.invalid_response_template = f.read().decode('utf-8')

        self.obspy_inv = None

        iris_record_pattern = r'net=(.+?)&sta=(.+?)&cha=(.+?)&'
        self.url_query_matcher = re.compile(iris_record_pattern)


    def get_full_response(self):
        """Getter for full station XML response

        Returns:
            str -- Detailed test station XMl with multiple stations and channels
        """
        return self.full_test_response


    def get_invalid_response(self, netcode, statcode='*', chcode='*'):
        """Get the response string for a channel level query for records that don't exist.

        :param netcode: The network code to query
        :type netcode: str
        :param statcode: The station code to query, defaults to '*'
        :type statcode: str, optional
        :param chcode: The channel code to query, defaults to '*'
        :type chcode: str, optional
        :return: Simulated IRIS response when the query cannot be fulfilled (404 error)
        :rtype: str
        """
        return self.invalid_response_template.format(netcode, statcode, chcode)


    def setup_internal_inv(self, response_functor):
        """Using given callable response generator response_function, obtain the
           response and convert it into an internal obspy inventory for later
           dynamic query.

        :param response_functor: Callable object to generate query response string (emulating IRIS web query)
            for a station inventory.
        :type response_functor: Any callable generating readable inventory string for ingestion into obspy.
        """
        resp = response_functor().encode('utf-8')
        self.obspy_inv = read_inventory(bio.BytesIO(resp))


    def dynamic_query(self, req, context):
        """Dynamically generate a simulated request response.

        :param req: The dynamic request.
        :type req: requests.Request
        :param context: Container for response data
        :type context: requests_mock.response._Context mimmicking fields of requests.Response
        :return: Emulated query response string
        :rtype: str
        """
        assert self.obspy_inv is not None, "Initialize self.obspy_inv by calling setup_internal_inv"
        request_str = req.url
        # Get network, station and channel query fields out of request_str and
        # pass to select(), then use that result to generate FDSN stationxml string response.
        match_result = self.url_query_matcher.search(request_str)
        assert match_result
        net, sta, cha = match_result[1], match_result[2], match_result[3]
        query_result = self.obspy_inv.select(network=net, station=sta, channel=cha)
        if query_result.networks:
            context.status_code = 200
            byte_buffer = bio.BytesIO()
            query_result.write(byte_buffer, "stationxml")
            return byte_buffer.getvalue().decode('utf-8')
        # No resource found matching the query
        context.status_code = 404
        return self.get_invalid_response(net, sta, cha)


    def get_minimal_response(self):
        """Get a very small station xml snippet for testing purposes

        Returns:
            unicode str -- FDSN station XML snippet
        """
        # Keep this as instance method for now, this allows flexibility to encapsulate mutation of the response
        # to generate error conditions in future.
        return u'''<?xml version="1.0" encoding="ISO-8859-1"?>
            <FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" schemaVersion="1.0" xsi:schemaLocation="http://www.fdsn.org/xml/station/1 http://www.fdsn.org/xml/station/fdsn-station-1.0.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:iris="http://www.fdsn.org/xml/station/1/iris">
            <Source>IRIS-DMC</Source>
            <Sender>IRIS-DMC</Sender>
            <Module>IRIS WEB SERVICE: fdsnws-station | version: 1.1.36</Module>
            <ModuleURI>http://service.iris.edu/fdsnws/station/1/query?net=GE&amp;sta=M*&amp;cha=BH*&amp;level=channel&amp;format=xml&amp;includerestricted=false&amp;includecomments=false&amp;nodata=404</ModuleURI>
            <Created>2019-05-22T03:17:10</Created>
            <Network code="GE" startDate="1991-01-01T00:00:00" endDate="2599-12-31T23:59:59" restrictedStatus="open">
            <Description>GEOFON</Description>
            <Station code="MAHO" startDate="1999-04-10T00:00:00" endDate="2001-02-13T00:00:00" restrictedStatus="open" iris:alternateNetworkCodes="_FDSN-ALL,.UNRESTRICTED">
            <Latitude>39.895901</Latitude>
            <Longitude>4.2665</Longitude>
            <Elevation>15.0</Elevation>
            <Site>
                <Name>ROA/UCM/GEOFON Station Mahon, Menorca, Spain</Name>
            </Site>
            <CreationDate>1999-04-10T00:00:00</CreationDate>
            <Channel code="BHZ" endDate="2001-02-13T00:00:00" locationCode="" restrictedStatus="open" startDate="1999-04-10T00:00:00">
                <Latitude>39.895901</Latitude>
                <Longitude>4.2665</Longitude>
                <Elevation>15</Elevation>
                <Depth>0</Depth>
                <Azimuth>0</Azimuth>
                <Dip>-90</Dip>
                <Type>TRIGGERED</Type>
                <Type>GEOPHYSICAL</Type>
                <SampleRate>2E01</SampleRate>
                <ClockDrift>2E-02</ClockDrift>
                <CalibrationUnits>
                <Name>A</Name>
                <Description>Amperes</Description>
                </CalibrationUnits>
                <Sensor>
                <Description>Streckeisen STS-2/N seismometer</Description>
                </Sensor>
                <Response>
                <InstrumentSensitivity>
                <Value>6.5006E8</Value>
                <Frequency>2E-2</Frequency>
                <InputUnits>
                <Name>M/S</Name>
                <Description>Velocity in Meters per Second</Description>
                </InputUnits>
                <OutputUnits>
                <Name>COUNTS</Name>
                <Description>Digital Counts</Description>
                </OutputUnits>
                </InstrumentSensitivity>
                </Response>
            </Channel>
            </Station>
            </Network>
            </FDSNStationXML>'''
