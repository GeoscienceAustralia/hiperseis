#!/usr/bin/env python
"""
Base class for seismogram synthesis class
"""

import abc

import numpy as np
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
        Initialization

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
        Function signature for function to compute synthetic dataset of obspy streams.

        :param src_latlon: Iterable of source (lat, lon) locations
        :type src_latlon: iterable of pairs
        :param fs: Sampling rate in Hz
        :type fs: float
        :param time_window: Pair of time values relative to onset
        :type time_window: tuple(float, float)
        :return: obspy.Stream containing ZNE velocity seismogram
        :rtype: obspy.Stream
        """
        pass  # Implement this interface in the derived class
    # end func

    def compute_event_stats(self, src_lat, src_lon, eventid_base, src_depth_m=0,
                            earth_model='iasp91', phase='P', origin_time=None):
        """
        Compute trace stats fields for a source single event.

        :param src_lat: Source latitude
        :type src_lat: float
        :param src_lon: Source longitude
        :type src_lon: float
        :param eventid_base: Base string for event id
        :type eventid_base: str
        :param src_depth_m: Source depth in metres
        :type src_depth_m: float
        :param earth_model: String name of earth model to use for ray tracing
        :type earth_model: str
        :param phase: Which phase is being modelled
        :type phase: str
        :param origin_time: Timestamp of the source event. If empty, will be a random offset from now.
        :type origin_time: obspy.UTCDateTime
        :return: Stats dictionary
        :rtype: dict
        """

        if origin_time is None:
            origin_time = UTCDateTime.now() - 60*np.random.rand()
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
        inc = arrival.incident_angle  # degrees
        event_id = eventid_base + '_'.join([
            gh.encode(receiver_lat, receiver_lon),
            gh.encode(src_lat, src_lon),
            origin_time.format_fissures()
        ])
        stats = {'distance': dist_deg, 'back_azimuth': baz, 'inclination': inc,
                 'onset': onset, 'slowness': ray_param, 'phase': phase,
                 'event_time': origin_time, 'tt_model': earth_model,
                 'event_id': event_id}
        return stats
    # end func

# end class
