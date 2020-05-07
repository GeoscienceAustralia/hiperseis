#!/usr/bin/env python
"""
Backend for making synthetic seismograms using Telewavesim.
"""

import numpy as np
from telewavesim.utils import Model, run_plane
import obspy

from seismic.synthetics.backends.synthesizer_base import Synthesizer
from seismic.model_properties import LayerProps
from seismic.units_utils import KM_PER_DEG
from seismic.stream_processing import zne_order


def synthesizer():
    return SynthesizerMatrixPropagator
# end func


class SynthesizerMatrixPropagator(Synthesizer):
    """
    Class to synthesize seismograms from a 1D model description using Kennett's
    matrix propagator method.
    """

    def __init__(self, station_latlon, layerprops):
        """
        Initialization

        :param station_latlon: See documentation for :func:`~seismic.synthetics.backends.synthesizer_base.Synthesizer.synthesize`
        :param layerprops: List of LayerProps. Last layer should be mantle properties.
        """
        super().__init__(station_latlon)
        Vp = [layer.Vp for layer in layerprops]
        Vs = [layer.Vs for layer in layerprops]
        self.model = Model([layer.H for layer in layerprops],
                           [layer.rho for layer in layerprops],
                           Vp, Vs, ['iso']*len(layerprops))
        self._kappa = np.array(Vp)/np.array(Vs)
    # end func

    @property
    def kappa(self):
        # Return ratio of Vp/Vs for each layer
        return self._kappa
    # end func

    def synthesize(self, src_latlon, fs, time_window):
        """
        See documentation for :func:`~seismic.synthetics.backends.synthesizer_base.Synthesizer.synthesize`
        """
        duration = time_window[1] - time_window[0]
        npts = int(np.ceil(duration*fs))
        dt = 1.0/fs
        stream_all = obspy.Stream()
        for src_lat, src_lon in src_latlon:
            stream = self._synthone(src_lat, src_lon, dt, npts)
            stream.traces = sorted(stream.traces, key=zne_order)
            for tr in stream:
                tr.stats.starttime = tr.stats.onset + time_window[0]
            # end for
            stream_all += stream
        # end for
        return stream_all
    # end func

    def _synthone(self, src_lat, src_lon, dt, npts):
        """
        Synthesize single event
        """
        event_id_base = 'hiperseis:PRM/evid='
        stats = self.compute_event_stats(src_lat, src_lon, event_id_base)
        stats.pop('tt_model')
        ray_param_sec_per_km = stats['slowness']/KM_PER_DEG
        traces_zne = run_plane(self.model, ray_param_sec_per_km, npts, dt, baz=stats['back_azimuth'])
        traces_zne.differentiate()
        for tr in traces_zne:
            stats['channel'] = '**' + tr.stats['channel'][-1]
            tr.stats.update(stats)
        # end for
        return traces_zne
    # end func

# end class


if __name__ == '__main__':
    example = SynthesizerMatrixPropagator((-20, 160),
                                          [LayerProps(6.4, 4.2, 2.7, 35.0),
                                           LayerProps(8.2, 6.8, 3.3, np.nan)])
    test_strm = example.synthesize(((30, 150),), 100.0, (-10, 30))
    print(test_strm)
# end if
