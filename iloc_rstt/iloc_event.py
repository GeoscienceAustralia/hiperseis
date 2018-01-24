import pandas as pd
from collections import defaultdict
from subprocess import check_call
from obspy import UTCDateTime
from obspy.core.event import (Event, Catalog, Origin, Comment, Pick, Arrival,
                              CreationInfo, ResourceIdentifier,
                              WaveformStreamID)
from obspy import read_events
from obspy.geodetics import locations2degrees
from inventory.parse_inventory import read_stations, sc3_inventory
from iloc_rstt.config import event_file, max_station_dist

stations = pd.read_csv(sc3_inventory)
stations_dict = read_stations(sc3_inventory)
DELTA = 'delta'  # distance in degrees between stations and source
STATION_CODE = 'station_code'
LATITUDE = 'latitude'
LONGITUDE = 'longitude'
ELEVATION = 'elevation'
NETWORK_CODE = 'network_code'
STATION_COLS = [STATION_CODE, LATITUDE , LONGITUDE,
                ELEVATION, NETWORK_CODE]
PHASEHINT = 'phase_hint'
ARRIVAL = 'arrival'
DBFLAG = "mysql://sysop:sysop@localhost/seiscomp3"


class ILocEvent(object):
    """
    Just a convenience class for event enhancement using iLoc
    """
    def __init__(self, event):
        self._pref_origin = event.preferred_origin() or \
                            event.origins[0]  # best guess origin

        # preferred_origin_id to point to old_origin
        self.preferred_origin_id = self._pref_origin.resource_id
        self._event = event
        self._iloc_origin = None
        self.picks = event.picks
        self.magnitudes = event.magnitudes
        self.amplitudes = event.amplitudes
        self.station_magnitudes = event.station_magnitudes
        self.resource_id = event.resource_id
        self.event_descriptions = event.event_descriptions
        self.comments = event.comments
        self.creation_info = event.creation_info
        self.origins = event.origins
        self.focal_mechanisms = event.focal_mechanisms
        self._max_stations_dist = max_station_dist
        self.__all_stations = None

    def __call__(self, *args, **kwargs):
        _ = self.all_stations
        return self

    def add_picks(self):
        pass
        # if self.all_stations is None:
        #     self.all_stations = self.add_stations()
        # return self.all_stations
        # mseed = self._get_miniseed()

    def add_magnitudes(self):
        return

    def _get_miniseed(self):
        # can make system call to seiscomp3 based on stations and event time
        # return mseed
        pass

    def add_origin(self):
        """
        Return the origin found by iLoc
        """
        if self._iloc_origin is None:
            # run iloc and return iloc determined origin
            return self.origins.append(Origin())
        else:
            return self.origins

    @property
    def all_stations(self):
        """
        :param max_station_dist: percentage distance to scan for extra stations

        Distance is calculated as percentage of the farthest station from the
        original preferred origin. A value of 100 implies an euclidean
        distance same as the of the farthest arriving station will be used.

        :return: pd.DataFrame of stations that satisfy range criteria
        """
        if self.__all_stations is None:
            self.__all_stations = pd.concat([self._add_primary_stations(),
                                             self._add_temporary_stations()])
        return self.__all_stations

    def _swap_preferred_origin_with_iloc(self):
        self._event.preferred_origin = self.add_origin()

    @staticmethod
    def _delta(x, origin):
        return locations2degrees(x['latitude'], x['longitude'],
                                 origin.latitude, origin.longitude)

    @staticmethod
    def _match_sta(x, sta_phs_dict):
        return sta_phs_dict[x['station_code']] \
            if x['station_code'] in sta_phs_dict else None

    def _add_primary_stations(self):
        stations[DELTA] = stations.apply(self._delta, axis=1,
                                         origin=self._pref_origin)
        max_dist, sta_phs_dict = self._farthest_station_dist()
        stations[ARRIVAL] = stations.apply(self._match_sta, axis=1,
                                           sta_phs_dict=sta_phs_dict)
        sel_stations = stations[stations[DELTA] <
                                self._max_stations_dist / 100 * max_dist]
        return sel_stations

    def _farthest_station_dist(self):
        max_dist = -1.0
        sta_phs_dict = defaultdict(dict)

        for arr in self._pref_origin.arrivals:
            pick = arr.pick_id.get_referred_object()
            sta = stations_dict[pick.waveform_id.station_code]
            cha = pick.waveform_id.channel_code
            dist = locations2degrees(self._pref_origin.latitude,
                                     self._pref_origin.longitude,
                                     sta.latitude, sta.longitude)
            if dist > max_dist:
                max_dist = dist

            # better to use phase as key than channel as same channel can be
            # picked for multiple phases
            sta_phs_dict[sta.station_code].update(
                    {arr.phase: {'channel': cha, 'time': pick.time}}
            )

        return max_dist, sta_phs_dict

    def _add_temporary_stations(self):
        return pd.DataFrame()


class ILocCatalog(Catalog):
    """
    Class that supports event/catalog enhancement using iLoc.
    iLoc ref: http://www.seismology.hu/index.php/en/home/iloc
    """
    def __init__(self, event_xml, **kwargs):
        self.event_xml = event_xml
        self.cat = read_events(event_xml)
        self.orig_events = self.cat.events
        self.events = None

        new_comments = kwargs.get("comments", [])
        self.comments = new_comments + self.cat.comments
        self._set_resource_id(kwargs.get("resource_id", None))
        new_description = kwargs.get("description", "PST ILocCatalog Modified")
        self.description = (('orig_description: '
                             + self.cat.description + '; ')
                            if self.cat.description is not None else '') + \
                           ('description: ' + new_description)
        old_ci = self.cat.creation_info

        self.creation_info = CreationInfo(
            agency_id=old_ci.agency_id,
            agency_uri=old_ci.agency_uri,
            author=(('orig_author: ' + old_ci.author + '; ')
                    if old_ci.author is not None else '') + 'PST ILocCatalog',
            creation_info=UTCDateTime()
        )

        super(ILocCatalog, self).__init__(
            events=self.events,
            comments=self.comments,
            creation_info=self.creation_info,
            description=self.description,
            resource_id=self.resource_id
        )

    def insert_into_sc3db(self):
        """
        needs seiscomp3 installed.
        """
        cmd = 'scdb -i'.split()
        cmd.append(self.event_xml)
        cmd += ['-d', DBFLAG]

        try:
            check_call(cmd)
        except OSError:
            raise OSError('seiscomp3 is not available in your system.')

    def update(self):
        for e in self.orig_events:
            iloc_ev = ILocEvent(e)()
            self.events.append(iloc_ev)


if __name__ == "__main__":
    cat = ILocCatalog(event_file)
    cat.update()
