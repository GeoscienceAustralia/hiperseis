"""
Description:
    Implements a Class for reading/storing responses from a number of sources.
    The ResponseFactory class is used for attaching 'bogus' responses to
    station inventories that lack them.

References:

CreationDate:   14/02/19

Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     14/02/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from collections import defaultdict
from io import StringIO

from obspy.core import UTCDateTime
from obspy import read_inventory
import sqlite3
import atexit
from io import BytesIO
from types import SimpleNamespace

class ResponseFactory:
    """
    The ResponseFactory class encapsulates the loading and retrieval of Instrument Response Objects from a
    variety of sources:
    i) from a stationXML file or obspy.core.Inventory containing a single Response Object
    ii) from Poles and Zeros
    iii) from a database, in which each row identified by net, sta, cha contains a response-level inventory
    """

    def __init__(self):
        self.response_cache = defaultdict(list)
        self.db_source = None
    # end func

    class ResponseFromDB(object):
        def __init__(self, db_fn):
            try:
                # Connect to SQLite database
                self.conn = sqlite3.connect(db_fn)
                self.cursor = self.conn.cursor()

            except sqlite3.Error as e:
                print("An error occurred: {}".format(e))
            # end try

            atexit.register(lambda: self.conn.close())
        # end func

        def getResponse(self, net, sta, loc, cha):
            try:
                q = "select sta_xml from responses where net='{}' and sta='{}' and cha='{}';".format(net, sta, cha)
                self.cursor.execute(q)

                row = self.cursor.fetchall()
                sta_xml = row[0][0]
                inv = read_inventory(BytesIO(sta_xml))

                # select desired location
                inv = inv.select(network=net, station=sta, location=loc, channel=cha)
                if(len(inv)):
                    try:
                        return inv.networks[0].stations[0].channels[0].response
                    except:
                        return None
                    # end try
                # end func
            except sqlite3.Error as e:
                print("An error occurred: {}".format(e))
            # end try
            return None
        # end func
    # end class

    class ResponseFromInventory(object):
        """Helper class to get Obspy Response object from an Inventory

        :raises RuntimeError: Raises error if response not found
        """
        def __init__(self, source_inventory):
            """Constructor

            :param source_inventory: Inventory from which to extract response
            :type source_inventory: obspy.core.inventory.inventory.Inventory
            """
            self.inventory = source_inventory
            self.response = None
            self._get_response_from_inventory()
        #end func

        def _get_response_from_inventory(self):
            n = None
            s = None
            c = None
            found = 0
            # Extract network, station and channel codes
            if self.inventory.networks:
                n = self.inventory.networks[0]
                found += 1
                if n.stations:
                    s = n.stations[0]
                    found += 1
                    if n.stations[0].channels:
                        c = n.stations[0].channels[0]
                        found += 1
                    # end if
                # end if
            # end if

            if (found < 1):
                msg = 'Network, station or channel information missing in RESP file.'
                raise RuntimeError(msg)
            else:
                seedid = self.inventory.get_contents()['channels'][0]
                self.response = self.inventory.get_response(seedid, c.start_date)
            # end if
        #end func
    # end class

    class ResponseFromStationXML(ResponseFromInventory):
        """Helper class to get Obspy Response object from a station xml file
        """
        def __init__(self, respFileName):
            """Constructor

            :param respFileName: XML file to load
            :type respFileName: str
            """
            xml_inventory = read_inventory(respFileName)
            super(ResponseFactory.ResponseFromStationXML, self).__init__(xml_inventory)
    # end class

    class ResponseFromResp(ResponseFromInventory):
        """Helper class to get Obspy Response object from a station xml file
        """
        def __init__(self, respFileName):
            """Constructor

            :param respFileName: XML file to load
            :type respFileName: str
            """
            resp_inventory = read_inventory(respFileName)
            super(ResponseFactory.ResponseFromResp, self).__init__(resp_inventory)
    # end class

    class ResponseFromPAZ:
        def __init__(self, pzTransferFunctionType='LAPLACE (RADIANS/SECOND)',
                     normFactor=8e4,
                     normFreq=1e-2,
                     stageGain=2e3,
                     stageGainFreq=1e-2,
                     poles=[0 + 0j],
                     zeros=[0 + 0j]):
            self.base = '''<?xml version="1.0" standalone="yes"?>
            <FDSNStationXML
                xmlns="http://www.fdsn.org/xml/station/1" schemaVersion="1">
                <Source>-i</Source>
                <Module>fdsn-stationxml-converter/1.0.9</Module>
                <ModuleURI>http://www.iris.edu/fdsnstationconverter</ModuleURI>
                <Created>2017-09-13T16:06:47.256+10:00</Created>
                <Network code="IU">
                    <Description></Description>
                    <Station code="ANMO" startDate="2002-11-19T21:07:00.000" endDate="2008-06-30T00:00:00.000">
                        <Latitude>0.0</Latitude>
                        <Longitude>0.0</Longitude>
                        <Elevation>0.0</Elevation>
                        <Site>
                            <Name>0</Name>
                        </Site>
                        <CreationDate>2002-11-19T21:07:00.000</CreationDate>
                        <Channel locationCode="00" code="BHZ" startDate="2002-11-19T21:07:00.000" endDate="2008-06-30T00:00:00.000">
                            <Latitude>0.0</Latitude>
                            <Longitude>0.0</Longitude>
                            <Elevation>0.0</Elevation>
                            <Depth>0.0</Depth>
                            <Azimuth>0.0</Azimuth>
                            <Dip>0.0</Dip>
                            <SampleRate>0.0</SampleRate>
                            <ClockDrift>0.0</ClockDrift>
                            <Response>
                                <InstrumentSensitivity>
                                    <Value>8.11597E8</Value>
                                    <Frequency>0.02</Frequency>
                                    <InputUnits>
                                        <Name>M/S</Name>
                                        <Description>Velocity in Meters Per Second</Description>
                                    </InputUnits>
                                    <OutputUnits>
                                        <Name>COUNTS</Name>
                                        <Description>Digital Counts</Description>
                                    </OutputUnits>
                                </InstrumentSensitivity>
                                <Stage number="1">
                                    <PolesZeros>
                                        <InputUnits>
                                            <Name>M/S</Name>
                                            <Description>Velocity in Meters Per Second</Description>
                                        </InputUnits>
                                        <OutputUnits>
                                            <Name>V</Name>
                                            <Description>Volts</Description>
                                        </OutputUnits>
                                        <PzTransferFunctionType>LAPLACE (RADIANS/SECOND)</PzTransferFunctionType>
                                        <NormalizationFactor>86083.0</NormalizationFactor>
                                        <NormalizationFrequency>0.02</NormalizationFrequency>
                                        <Zero number="0">
                                            <Real plusError="0.0" minusError="0.0">0.0</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">0.0</Imaginary>
                                        </Zero>
                                        <Zero number="1">
                                            <Real plusError="0.0" minusError="0.0">0.0</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">0.0</Imaginary>
                                        </Zero>
                                        <Pole number="0">
                                            <Real plusError="0.0" minusError="0.0">-59.4313</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">0.0</Imaginary>
                                        </Pole>
                                        <Pole number="1">
                                            <Real plusError="0.0" minusError="0.0">-22.7121</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">27.1065</Imaginary>
                                        </Pole>
                                        <Pole number="2">
                                            <Real plusError="0.0" minusError="0.0">-22.7121</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">-27.1065</Imaginary>
                                        </Pole>
                                        <Pole number="3">
                                            <Real plusError="0.0" minusError="0.0">-0.0048004</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">0.0</Imaginary>
                                        </Pole>
                                        <Pole number="4">
                                            <Real plusError="0.0" minusError="0.0">-0.073199</Real>
                                            <Imaginary plusError="0.0" minusError="0.0">0.0</Imaginary>
                                        </Pole>
                                    </PolesZeros>
                                    <StageGain>
                                        <Value>1935.0</Value>
                                        <Frequency>0.02</Frequency>
                                    </StageGain>
                                </Stage>
                            </Response>
                        </Channel>
                    </Station>
                </Network>
            </FDSNStationXML>
            '''
            self.response = None
            self.pzTransferFunctionType = pzTransferFunctionType
            self.normFactor = normFactor
            self.normFreq = normFreq
            self.stageGain = stageGain
            self.stageGainFreq = stageGainFreq
            self.poles = poles
            self.zeros = zeros

            # Generate a valid response inventory based on the fdsnstationxml file, which
            # was downloaded from IRIS as an example.
            inv = read_inventory(StringIO(self.base))
            datetime = UTCDateTime("2002-11-19T21:07:00.000")

            # Fetch the response object and adapt its parameters based on user-input.
            self.response = inv.get_response('IU.ANMO.00.BHZ', datetime)

            self.response.response_stages[0].pz_transfer_function_type = self.pzTransferFunctionType
            self.response.response_stages[0].normalization_factor = self.normFactor
            self.response.response_stages[0].normalization_frequency = self.normFreq
            self.response.response_stages[0].stage_gain = self.stageGain
            self.response.response_stages[0].stage_gain_frequency = self.stageGainFreq
            self.response.response_stages[0].poles = poles
            self.response.response_stages[0].zeros = zeros
        # end func
    # end class

    def createFromInventory(self, name, obspy_inventory):
        """Create response from an Inventory

        :param name: Name of the response for later retrieval
        :type name: str
        :param obspy_inventory: Inventory from which to extract response
        :type obspy_inventory: obspy.core.inventory.inventory.Inventory
        """
        self.response_cache[name] = ResponseFactory.ResponseFromInventory(obspy_inventory)
    # end func

    def createFromStationXML(self, name, staionXMLFileName):
        """Create response from an XML file

        :param name: Name of the response for later retrieval
        :type name: str
        :param staionXMLFileName: XML file to load
        :type staionXMLFileName: str
        """
        self.response_cache[name] = ResponseFactory.ResponseFromStationXML(staionXMLFileName)
    # end func

    def createFromRespFile(self, name, respFileName):
        """Create response from an Resp file

        :param name: Name of the response for later retrieval
        :type name: str
        :param respFileName: Resp file to load
        :type respFileName: str
        """
        self.response_cache[name] = ResponseFactory.ResponseFromStationXML(respFileName)
    # end func

    def createFromPAZ(self, name, pzTransferFunctionType,
                      normFactor,
                      normFreq,
                      stageGain,
                      stageGainFreq,
                      poles,
                      zeros):

        self.response_cache[name] = ResponseFactory.ResponseFromPAZ(pzTransferFunctionType,
                                                                    normFactor,
                                                                    normFreq,
                                                                    stageGain,
                                                                    stageGainFreq,
                                                                    poles,
                                                                    zeros)
    #end func

    def createFromDB(self, db_fn):
        self.db_source = self.ResponseFromDB(db_fn)
    # end func

    def getResponse(self, name):
        """Retrieve response by name

        :param name: Name given to response at creation time, or a string containing
                     'network.station.location.channel' to query a database source, if it exists
        :type name: str
        :raises RuntimeError: Raises error if name is not recognized
        :return: The requested response
        :rtype: obspy.core.inventory.response.Response
        """
        if (name in self.response_cache.keys()):
            return self.response_cache[name].response
        elif(self.db_source):
            # try breaking down name into net, sta, loc, cha
            net, sta, loc, cha = name.split('.')
            resp = self.db_source.getResponse(net, sta, loc, cha)

            if(resp): self.response_cache[name] = SimpleNamespace(**{'response': resp})

            return resp
        else:
            msg = "Response with name: %s not found.." % (name)
            raise RuntimeError(msg)
        # end if
    # end func
# end class
