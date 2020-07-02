#!/usr/bin/env python
"""
Backend for making synthetic seismograms using Syngine.
"""

import obspy
from obspy.clients.syngine import Client as ClientS

from seismic.synthetics.backends.synthesizer_base import Synthesizer
from seismic.stream_processing import zne_order

# pylint: disable=invalid-name


def synthesizer():
    """Getter for backend Synthesizer class

    :return: Class name
    :rtype: SynthesizerSyngine
    """
    return SynthesizerSyngine
# end func


class SynthesizerSyngine(Synthesizer):
    """
    Class to synthesize seismograms using online Syngine service.

    To write resultant stream to HDF5 format, add 'ignore' option::

        synth_stream.write('test_synth.h5', 'h5', ignore=('mseed',))

    """

    def __init__(self, station_latlon, earth_model='iasp91'):
        """
        Initialization

        :param station_latlon: See documentation for :func:`~seismic.synthetics.backends.synthesizer_base.Synthesizer.synthesize`
        :param earth_model: String naming which standard earth model to use.
        :type earth_model: str
        """
        super().__init__(station_latlon)
        self.earth_model = earth_model
    # end func

    def synthesize(self, src_latlon, fs, time_window):
        """
        See documentation for :func:`~seismic.synthetics.backends.synthesizer_base.Synthesizer.synthesize`
        """
        dt = 1.0/fs
        stream_all = obspy.Stream()
        for src_lat, src_lon in src_latlon:
            stream = self._synthone(src_lat, src_lon, time_window, dt)
            stream.traces = sorted(stream.traces, key=zne_order)
            stream_all += stream
        # end for
        return stream_all
    # end func

    def _synthone(self, src_lat, src_lon, time_window, dt, src_depth_m=0, origintime=None):
        """
        Synthesize single event
        """
        components = 'ZNE'
        units = 'velocity'
        default_model = self.earth_model + '_2s'

        if origintime is None:
            origintime = obspy.UTCDateTime.now()
        # end if

        client_synth = ClientS()
        assert time_window[0] <= 0
        assert time_window[1] >= 0
        receiver_lat, receiver_lon = self.station_latlon
        starttime = 'P' + '-{}'.format(abs(time_window[0]))
        endtime = 'P' + '+{}'.format(time_window[1])
        synth_stream = client_synth.get_waveforms(model=default_model,
                                                  receiverlatitude=receiver_lat,
                                                  receiverlongitude=receiver_lon,
                                                  starttime=starttime,
                                                  endtime=endtime,
                                                  dt=dt,
                                                  components=components,
                                                  units=units,
                                                  scale=1e9,  # nanometres
                                                  sourcelatitude=src_lat,
                                                  sourcelongitude=src_lon,
                                                  sourcedepthinmeters=src_depth_m,
                                                  origintime=origintime)

        event_id_base = 'hiperseis:SYN/evid='
        stats = self.compute_event_stats(src_lat, src_lon, event_id_base,
                                         earth_model=self.earth_model,
                                         origin_time=origintime)
        for tr in synth_stream:
            tr.stats.update(stats)
        # end for

        return synth_stream
    # end func

# end class


if __name__ == '__main__':
    example = SynthesizerSyngine('AU.ARMA')
    test_strm = example.synthesize(((30, 150),), 100.0, (-10, 30))
    print(test_strm)
# end if
