from lxml import etree as ET
import re
from collections import defaultdict
from obspy import UTCDateTime


class Origin:
    def __init__(self, utctime, lat, lon, depthkm):
        self.utctime = utctime
        self.lat = lat
        self.lon = lon
        self.depthkm = depthkm
        self.magnitude_list = []
# end func


# end class

class Event:
    def __init__(self):
        self.preferred_origin = None
        self.preferred_magnitude = None
        self.origin_list = []
    # end func
# end class

class Magnitude:
    def __init__(self, mag, mag_type):
        self.magnitude_value = mag
        self.magnitude_type = mag_type
    # end func
# end class

class EventParser:
    def __init__(self, xml_filename):
        self.xml_filename = xml_filename
        self.events  = []
        self.origins = []
        self.magnitudes = []

        for xmlevent, elem in ET.iterparse(xml_filename):
            tag = ET.QName(elem.tag).localname
            if (tag.upper() == 'EVENT'):
                self.events.append(elem)
            elif (tag.upper() == 'ORIGIN'):
                self.origins.append(elem)
            elif (tag.upper() == 'MAGNITUDE'):
                self.magnitudes.append(elem)
            # end if
        # end for

        self.origin_dict = defaultdict(list)
        self.magnitude_dict = defaultdict(list)

        # parse magnitudes
        for m in self.magnitudes:
            pid = None
            try:
                pid = int(re.split("[/=]", m.get('publicID'))[-1])
            except:
                pass
            # end try
            if(pid): self.magnitude_dict[pid] = self.parseMagnitude(m)
        # end for

        # parse origins
        for o in self.origins:
            pid = int(re.split("[/=]", o.get('publicID'))[-1])
            self.origin_dict[pid] = self.parseOrigin(o)
            # end for
        assert len(self.origins) == len(self.origin_dict), 'Invalid originIDs found..'

        self.event_dict = defaultdict(list)
        # parse events
        for e in self.events:
            pid = int(re.split("[/=]", e.get('publicID'))[-1])
            self.event_dict[pid] = self.parseEvent(e)
        # end for
        assert len(self.events) == len(self.event_dict), 'Invalid eventIDs found..'

    # end func

    def parseEvent(self, e):
        resultEvent = Event()
        origList = []
        for element in e:
            if ('preferredOriginID' in element.tag):
                oid = int(re.split("[/=]", element.text)[-1])
                resultEvent.preferred_origin = self.origin_dict[oid]
            elif('preferredMagnitudeID' in element.tag):
                mid = int(re.split("[/=]", element.text)[-1])
                resultEvent.preferred_magnitude = self.magnitude_dict[mid]
            elif ('origin' in element.tag and \
                  'Reference' not in element.tag):
                oid = int(re.split("[/=]", element.get('publicID'))[-1])
                origList.append(self.origin_dict[oid])
            # end if
        # end for
        resultEvent.origin_list = origList
        return resultEvent
    # end func

    def parseOrigin(self, o):
        utctime = None
        lat = None
        lon = None
        depthkm = None
        pref_mag = None
        magList = []
        for element in o:
            if ('time' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        utctime = UTCDateTime(subelement.text)
                        # end if
                        # end for
            elif ('latitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        lat = float(subelement.text)
                    # end if
                # end for
            # end if
            elif ('longitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        lon = float(subelement.text)
                    # end if
                # end for
            # end if
            elif ('depth' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        depthkm = float(subelement.text)
                    # end if
                # end for
            elif ('magnitude' in element.tag and \
                          'station' not in element.tag):
                mid = None
                try:
                    mid = int(re.split("[/=]", element.get('publicID'))[-1])
                except:
                    pass
                # end try
                if(mid): magList.append(self.magnitude_dict[mid])
            # end if
        # end for

        assert None not in [utctime, lat, lon], 'Failed to find required values for Origin'

        o = Origin(utctime, lat, lon, depthkm)
        o.magnitude_list = magList
        return o
    # end func

    def parseMagnitude(self, m):
        mag_value = None
        mag_type = None
        for element in m:
            if ('magnitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        mag_value = float(subelement.text)
                    # end if
                # end for
            elif ('type' in element.tag):
                mag_type = element.text
            # end if
        # end for

        return Magnitude(mag_value, mag_type)
    # end func

    def getEvents(self):
        return self.event_dict
    # end func
# end class
