import pandas as pd
from obspy.core.event import Event, Catalog, Origin, Comment, Pick
from obspy import read_events
from obspy.geodetics import locations2degrees
from inventory.parse_inventory import read_stations, sc3_inventory
from iloc_rstt.config import event_file

stations = pd.read_csv(sc3_inventory)
stations_dict = read_stations(sc3_inventory)
DELTA = 'delta'  # distance in degrees between stations and source


class ILocEvent:
    """
    Just a convenience class for event enhancement using iLoc
    """
    def __init__(self, event):
        self.event = event
        self.iloc_origin = None
        self.old_origin = self.event.preferred_origin() or \
            self.event.origins[0]  # best guess origin

    def update(self):
        self.event.picks += self.add_picks()
        self.event.origins += [self.add_origin()]

    def add_picks(self):
        stations = self.add_stations()
        return [Pick()]

    def add_origin(self):
        """
        Return the origin found by iLoc
        """
        if self.iloc_origin is None:
            # run iloc and return iloc determined origin
            return Origin()
        else:
            return self.iloc_origin

    def add_stations(self, range=120.0):
        """
        :param range: percentage distance to scan for extra stations.

        Distance is calculated as percentage of the farthest station from the
        original preferred origin. A value of 100 implies an euclidean
        distance same as the of the farthest arriving station will be used.

        :return: list of stations that satisfy range criteria
        """
        return pd.concat(self._add_primary_stations(range),
                         self._add_temporary_stations())

    def _swap_preferred_origin_with_iloc(self):
        self.event.preferred_origin = self.add_origin()

    @staticmethod
    def _delta(x, origin):
        return locations2degrees(x['latitude'], x['longitude'],
                                 origin.latitude, origin.longitude)

    def _add_primary_stations(self, range):
        stations[DELTA] = stations.apply(self._delta, axis=1,
                                         origin=self.old_origin)

        return stations[stations[DELTA] < range/100 *
                        self._farthest_station_dist()]

    def _farthest_station_dist(self):
        max_dist = -1.0
        for arr in self.old_origin.arrivals:
            sta = stations_dict[
                arr.pick_id.get_referred_object().waveform_id.station_code
            ]
            dist = locations2degrees(self.old_origin.latitude,
                                     self.old_origin.longitude,
                                     sta.latitude, sta.longitude)
            if dist > max_dist:
                max_dist = dist


        return max_dist

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
            self.iloc_events.append(ILocEvent(e).update())

    def write(self, new_sc3ml):
        self.cat.write(new_sc3ml, format='SC3ML')


if __name__ == "__main__":
    cat = ILocCatalog(event_file)
    cat.update()
