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


class ResponseFactory:
    """
    The ResponseFactory class encapsulates the generation of a collection of named Instrument Response Objects from a
    variety of sources. Currently it provides the facility to create Response objects from two sources, namely,
    Poles and Zeroes supplied by the user and from StationXML files generated from RESP files using the PDCC tool
    (link below). The conversion of RESP files into a corresponding StationXML file, at this stage, must take place
    externally, because ObsPy lacks that functionality.
    The intended usage of this class during the creation of an ASDF dataset is as follows::

      1. User creates a number of uniquely named Response objects (see associated tests as well) pertaining to different
         channels in a given survey.
      2. User fetches these Response objects from an instance of ResponseFactory as needed, while creating ObsPy Channel
         objects, during which Response objects can be passed in as an argument.
      3. User builds a hierarchy of channel->station->network inventories, with the appropriate instrument response
         information embedded
      4. The master FDSN StaionXML file output after step 3 can then be converted into an SC3ML file (which can be
         ingested by SeisComp3) using the fdsnxml2inv tool.

    PDCC tool: https://ds.iris.edu/ds/nodes/dmc/software/downloads/pdcc/
    """

    def __init__(self):
        self.m_responseInventory = defaultdict(list)
    # end func

    class ResponseFromInventory(object):
        """Helper class to get Obspy Response object from an Inventory

        :raises RuntimeError: Raises error if response not found
        """
        def __init__(self, source_inventory):
            """Constructor

            :param source_inventory: Inventory from which to extract response
            :type source_inventory: obspy.core.inventory.inventory.Inventory
            """
            self.m_inventory = source_inventory
            self.m_response = None
            self._get_response_from_inventory()
        #end func

        def _get_response_from_inventory(self):
            n = None
            s = None
            c = None
            found = 0
            # Extract network, station and channel codes
            if self.m_inventory.networks:
                n = self.m_inventory.networks[0]
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
                seedid = self.m_inventory.get_contents()['channels'][0]
                self.m_response = self.m_inventory.get_response(seedid, c.start_date)
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

    class ResponseFromPAZ:
        def __init__(self, pzTransferFunctionType='LAPLACE (RADIANS/SECOND)',
                     normFactor=8e4,
                     normFreq=1e-2,
                     stageGain=2e3,
                     stageGainFreq=1e-2,
                     poles=[0 + 0j],
                     zeros=[0 + 0j]):
            self.m_base = u'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
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

            self.m_response = None
            self.m_pzTransferFunctionType = pzTransferFunctionType
            self.m_normFactor = normFactor
            self.m_normFreq = normFreq
            self.m_stageGain = stageGain
            self.m_stageGainFreq = stageGainFreq
            self.m_poles = poles
            self.m_zeros = zeros

            # Generate a valid response inventory based on the fdsnstationxml file, which
            # was downloaded from IRIS as an example.
            inv = read_inventory(StringIO(self.m_base))
            datetime = UTCDateTime("2002-11-19T21:07:00.000")

            # Fetch the response object and adapt its parameters based on user-input.
            self.m_response = inv.get_response('IU.ANMO.00.BHZ', datetime)

            self.m_response.instrument_sensitivity = None  # Get rid of redundant information
            self.m_response.response_stages[0].pz_transfer_function_type = self.m_pzTransferFunctionType
            self.m_response.response_stages[0].normalization_factor = self.m_normFactor
            self.m_response.response_stages[0].normalization_frequency = self.m_normFreq
            self.m_response.response_stages[0].stage_gain = self.m_stageGain
            self.m_response.response_stages[0].stage_gain_frequency = self.m_stageGainFreq
            self.m_response.response_stages[0].poles = poles
            self.m_response.response_stages[0].zeros = zeros
        # end func
    # end class

    def CreateFromInventory(self, name, obspy_inventory):
        """Create response from an Inventory

        :param name: Name of the response for later retrieval
        :type name: str
        :param obspy_inventory: Inventory from which to extract response
        :type obspy_inventory: obspy.core.inventory.inventory.Inventory
        """
        self.m_responseInventory[name] = ResponseFactory.ResponseFromInventory(obspy_inventory)
    # end func

    def CreateFromStationXML(self, name, respFileName):
        """Create response from an XML file

        :param name: Name of the response for later retrieval
        :type name: str
        :param respFileName: XML file to load
        :type respFileName: str
        """
        self.m_responseInventory[name] = ResponseFactory.ResponseFromStationXML(respFileName)
    # end func

    def CreateFromPAZ(self, name, pzTransferFunctionType,
                      normFactor,
                      normFreq,
                      stageGain,
                      stageGainFreq,
                      poles,
                      zeros):

        self.m_responseInventory[name] = ResponseFactory.ResponseFromPAZ(pzTransferFunctionType,
                                                                         normFactor,
                                                                         normFreq,
                                                                         stageGain,
                                                                         stageGainFreq,
                                                                         poles,
                                                                         zeros)

    #end func

    def getResponse(self, name):
        """Retrieve response by name

        :param name: Name given to response at creation time
        :type name: str
        :raises RuntimeError: Raises error if name is not recognized
        :return: The requested response
        :rtype: obspy.core.inventory.response.Response
        """
        if (name in self.m_responseInventory.keys()):
            return self.m_responseInventory[name].m_response
        else:
            msg = "Response with name: %s not found.." % (name)
            raise RuntimeError(msg)
        # end if
    # end func
# end class
