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
    # end func
# end class

class Event:
    def __init__(self):
        self.preferred_origin = None
        self.origins = []
    # end func
# end class

class EventParser:
    def __init__(self, xml_filename):
        self.events  = []
        self.origins = []

        for xmlevent, elem in ET.iterparse(xml_filename):
            tag = ET.QName(elem.tag).localname
            if (tag.upper() == 'EVENT'):
                self.events.append(elem)
            # end if
            if (tag.upper() == 'ORIGIN'):
                self.origins.append(elem)
            # end if
        # end for

        self.origin_dict = defaultdict(list)
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
        for element in e:
            if ('preferredOriginID' in element.tag):
                oid = int(re.split("[/=]", element.text)[-1])
                resultEvent.preferred_origin = self.origin_dict[oid]
                # end if
        # end for

        return resultEvent

    # end func

    def parseOrigin(self, o):
        utctime = None
        lat = None
        lon = None
        depthkm = None
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
                        # end if
        # end for

        assert None not in [utctime, lat, lon], 'Failed to find required values for Origin'

        return Origin(utctime, lat, lon, depthkm)

    # end func

    def getEvents(self):
        return self.event_dict
    # end func
# end class
