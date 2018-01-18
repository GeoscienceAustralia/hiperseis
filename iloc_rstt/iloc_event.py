import pandas as pd
from obspy.core.event import Event, Catalog, Origin, Comment, Pick
from obspy import read_events
from obspy.geodetics import locations2degrees
from inventory.parse_inventory import read_stations, sc3_inventory
from iloc_rstt.config import event_file, max_station_dist

stations = pd.read_csv(sc3_inventory)
stations_dict = read_stations(sc3_inventory)
DELTA = 'delta'  # distance in degrees between stations and source
STATION_COLS = ['station_code', 'latitude', 'longitude',
                'elevation', 'network_code']


class ILocEvent:
    """
    Just a convenience class for event enhancement using iLoc
    """
    def __init__(self, event):
        self.event = event
        self.iloc_origin = None
        self.old_origin = self.event.preferred_origin() or \
            self.event.origins[0]  # best guess origin
        self.new_stations = None

    def update(self):
        self.add_picks()
        # self.event.picks += self.add_picks()
        # self.event.origins += [self.add_origin()]

    def add_picks(self):
        if self.new_stations is None:
            self.new_stations = self.add_stations()
        return self.new_stations

    def add_origin(self):
        """
        Return the origin found by iLoc
        """
        if self.iloc_origin is None:
            # run iloc and return iloc determined origin
            return Origin()
        else:
            return self.iloc_origin

    def add_stations(self, max_station_dist=max_station_dist):
        """
        :param max_station_dist: percentage distance to scan for extra stations

        Distance is calculated as percentage of the farthest station from the
        original preferred origin. A value of 100 implies an euclidean
        distance same as the of the farthest arriving station will be used.

        :return: pd.DataFrame of stations that satisfy range criteria
        """
        return pd.concat([self._add_primary_stations(max_station_dist),
                          self._add_temporary_stations()])

    def _swap_preferred_origin_with_iloc(self):
        self.event.preferred_origin = self.add_origin()

    @staticmethod
    def _delta(x, origin):
        return locations2degrees(x['latitude'], x['longitude'],
                                 origin.latitude, origin.longitude)

    def _add_primary_stations(self, max_range):
        stations[DELTA] = stations.apply(self._delta, axis=1,
                                         origin=self.old_origin)
        max_dist, orig_stas = self._farthest_station_dist()
        new_stations = stations.loc[~stations['station_code'].isin(orig_stas)]
        return new_stations[new_stations[DELTA] < max_range / 100 * max_dist]

    def _farthest_station_dist(self):
        max_dist = -1.0
        orig_stas_list = []
        for arr in self.old_origin.arrivals:
            sta = stations_dict[
                arr.pick_id.get_referred_object().waveform_id.station_code
            ]
            dist = locations2degrees(self.old_origin.latitude,
                                     self.old_origin.longitude,
                                     sta.latitude, sta.longitude)
            if dist > max_dist:
                max_dist = dist
            orig_stas_list.append(sta.station_code)

        return max_dist, orig_stas_list

    def _add_temporary_stations(self):
        return pd.DataFrame()


class ILocCatalog:
    """
    Class that supports event/catalog enhancement using iLoc.
    iLoc ref: http://www.seismology.hu/index.php/en/home/iloc
    """
    def __init__(self, event_xml):
        self.event_xml = event_xml
        self.cat = read_events(event_xml)
        self.events = self.cat.events
        self.iloc_events = []

    def update(self):
        for e in self.events:
            iloc_ev = ILocEvent(e)
            iloc_ev.update()
            self.iloc_events.append(iloc_ev)

    def write(self, new_sc3ml):
        self.cat.write(new_sc3ml, format='SC3ML')


if __name__ == "__main__":
    cat = ILocCatalog(event_file)
    cat.update()
