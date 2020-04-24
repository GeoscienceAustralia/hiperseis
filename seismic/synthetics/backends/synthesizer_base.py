#!/usr/bin/env python
"""
TBD
"""

import abc

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client as ClientF
import geohash as gh

from seismic.units_utils import KM_PER_DEG


class Synthesizer(object):
    """
    Base class for seismogram synthesizers.
    """

    def __init__(self, station_latlon):
        """

        :param station_latlon: Either a tuple of (lat, lon) coordinates, or a station
            code in the format 'NET.STA' string.
        """
        if isinstance(station_latlon, str):
            net, sta = station_latlon.split('.')
            client_real = ClientF()
            station_metadata = client_real.get_stations(network=net, station=sta)
            station_metadata = station_metadata.select(channel='*Z')
            receiver_lat = station_metadata.networks[0].stations[0].latitude
            receiver_lon = station_metadata.networks[0].stations[0].longitude
            station_latlon = (receiver_lat, receiver_lon)
        # end if
        self.station_latlon = station_latlon
    # end func

    @property
    def station_latlon(self):
        return self._station_latlon
    # end func

    @station_latlon.setter
    def station_latlon(self, new_latlon):
        assert isinstance(new_latlon, (tuple, list))
        assert len(new_latlon) == 2
        assert abs(new_latlon[0] <= 90)
        assert abs(new_latlon[1] <= 180)
        self._station_latlon = new_latlon
    # end func

    @abc.abstractmethod
    def synthesize(self, src_latlon, fs, time_window):
        """

        :param src_latlon: Iterable of source (lat, lon) locations
        :param fs: Sampling rate in Hz
        :param time_window: Pair of time values relative to onset
        :return: obspy.Stream containing ZNE velocity seismogram
        """
        pass  # Implement this interface in the derived class
    # end func

    def compute_event_stats(self, src_lat, src_lon, eventid_base, src_depth_m=0,
                            earth_model='iasp91', phase='P', origin_time=None):
        """
        TBD

        :param src_lat:
        :param src_lon:
        :param eventid_base:
        :param src_depth_m:
        :param earth_model:
        :param phase:
        :param origin_time:
        :return:
        """

        if origin_time is None:
            origin_time = UTCDateTime.now()
        # end if
        receiver_lat, receiver_lon = self.station_latlon
        dist_m, baz, _ = gps2dist_azimuth(receiver_lat, receiver_lon, src_lat, src_lon)
        dist_deg = dist_m / 1000 / KM_PER_DEG
        tt_model = TauPyModel(model=earth_model)
        event_depth_km = src_depth_m / 1000
        arrivals = tt_model.get_travel_times(event_depth_km, dist_deg, (phase,))
        arrival = arrivals[0]
        ray_param = arrival.ray_param_sec_degree
        onset = origin_time + arrival.time
        inc = arrival.incident_angle
        event_id = eventid_base + '_'.join([
            gh.encode(receiver_lat, receiver_lon),
            gh.encode(src_lat, src_lon),
            origin_time.format_fissures()
        ])
        stats = {'distance': dist_deg, 'back_azimuth': baz, 'inclination': inc,
                 'onset': onset, 'slowness': ray_param, 'phase': phase,
                 'origin_time': origin_time, 'tt_model': earth_model,
                 'event_id': event_id}
        return stats
    # end func

# end class
