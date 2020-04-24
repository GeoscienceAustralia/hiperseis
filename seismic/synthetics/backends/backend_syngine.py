#!/usr/bin/env python
"""
TBD
"""

import obspy

from seismic.synthetics.backends.synthesizer_base import Synthesizer
from seismic.units_utils import KM_PER_DEG


def synthesizer():
    return SynthesizerSyngine
# end func


class SynthesizerSyngine(Synthesizer):
    """
    TBD
    """

    def __init__(self, station_latlon, earthmodel):
        """

        :param station_latlon:
        :param earthmodel: String naming which standard earth model to use.
        """
        super().__init__(station_latlon)
        pass
    # end func

    def synthesize(self, src_latlon, fs, time_window):
        pass
    # end func

# end class


def synthesize_ideal_seismogram(network, station, units, sourcelatitude, sourcelongitude, sourcedepthmetres=0,
                                timewindow=(-20, 60), components='ZRT', origintime=None, f_s=None):
    """
    Given a receiving station and basic seismic source parameters, generate apure synthetic seismic
    waveform at the surface at the receiver location, without any instrument response.

    Uses IASP91 model.

    To write resultant stream to HDF5 format, add 'ignore' option:
        `synth_stream.write('test_synth.h5', 'h5', ignore=('mseed',))`

    :return: Stream containing synthetic waveform.
    """
    from obspy.clients.syngine import Client as ClientS
    from obspy.clients.fdsn import Client as ClientF
    from obspy.taup import TauPyModel
    from obspy.geodetics import gps2dist_azimuth

    earth_model = 'iasp91'
    default_model = earth_model + '_2s'

    if origintime is None:
        origintime = obspy.UTCDateTime.now()
    # end if

    if f_s is not None:
        dt = 1.0/f_s
    else:
        dt = None
    # end if

    client_synth = ClientS()
    assert timewindow[0] <= 0
    assert timewindow[1] >= 0
    starttime = 'P' + '-{}'.format(abs(timewindow[0]))
    endtime = 'P' + '+{}'.format(timewindow[1])
    synth_stream = client_synth.get_waveforms(model=default_model, network=network, station=station,
                                              starttime=starttime, endtime=endtime, components=components,
                                              units=units, sourcelatitude=sourcelatitude, dt=dt,
                                              sourcelongitude=sourcelongitude, sourcedepthinmeters=sourcedepthmetres,
                                              origintime=origintime)

    origintime = obspy.UTCDateTime(origintime)

    client_real = ClientF()
    station_metadata = client_real.get_stations(network=network, station=station)
    station_metadata = station_metadata.select(channel='*Z')
    receiver_lat = station_metadata.networks[0].stations[0].latitude
    receiver_lon = station_metadata.networks[0].stations[0].longitude
    dist, baz, _ = gps2dist_azimuth(receiver_lat, receiver_lon, sourcelatitude, sourcelongitude)
    dist_deg = dist / 1000 / KM_PER_DEG
    event_depth_km = sourcedepthmetres/1000
    tt_model = TauPyModel(model=earth_model)
    phase = 'P'
    arrivals = tt_model.get_travel_times(event_depth_km, dist_deg, (phase,))
    arrival = arrivals[0]
    onset = origintime + arrival.time
    inc = arrival.incident_angle
    slowness = arrival.ray_param_sec_degree
    stats = {'distance': dist_deg, 'back_azimuth': baz, 'inclination': inc,
             'onset': onset, 'slowness': slowness, 'phase': phase, 'tt_model': earth_model}
    for tr in synth_stream:
        tr.stats.update(stats)
    # end for

    return synth_stream
# end func


if __name__ == "__main__":
    # s = synthesize_ideal_seismogram('AU', 'QIS', 'velocity', 40, 140)
    pass
# end if
