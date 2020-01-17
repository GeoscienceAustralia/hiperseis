#!/usr/bin/env python
"""Helper functions for producing synthetic pseudo-Receiver function traces
"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d

import obspy
import rf

import seismic.receiver_fn.rf_util as rf_util

# pylint: disable=invalid-name

def generate_synth_rf(arrival_times, arrival_amplitudes, fs_hz=100.0, window_sec=(-10, 30), f_cutoff_hz=2.0):
    """Simple generator of synthetic R component receiver function with pulses at given arrival times.

    :param arrival_times: Iterable of arrival times as numerical values in seconds
    :type arrival_times: iterable of float
    :param arrival_amplitudes: Iterable of arrival amplitudes
    :type arrival_amplitudes: iterable of float
    :param fs_hz: Sampling rate (Hz) of output signal, defaults to 100.0
    :type fs_hz: float, optional
    :param window_sec: Time window over which to create signal (sec), defaults to (-10, 30)
    :type window_sec: tuple, optional
    :param f_cutoff_hz: Cutoff frequency (Hz) for low-pass filtering to generate realistic result, defaults to 2.0
    :type f_cutoff_hz: float, optional
    :return: Array of times and corresponding signal amplitudes
    :rtype: numpy.array, numpy.array
    """
    # Compute array of time values and indexes of arrivals
    duration = window_sec[1] - window_sec[0]
    N = int(fs_hz*duration)
    times = np.linspace(window_sec[0], window_sec[1], N)
    arrivals_index = np.round((np.array(arrival_times) - times[0])*fs_hz).astype(int)

    # Generate kernel of delta functions at specified arrival times
    kernel = np.zeros_like(times)
    kernel[arrivals_index] = np.array(arrival_amplitudes)  # pylint: disable=unsupported-assignment-operation

    # Filter to pass low frequencies
    if f_cutoff_hz:
        waveform = signal.butter(4, f_cutoff_hz/fs_hz)
        signal_filt = signal.filtfilt(waveform[0], waveform[1], kernel)
    else:
        signal_filt = kernel
    # end if

    return times, signal_filt
# end func


def synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, ds, log=None, include_t3=False,
                          amplitudes=None, baz=0.0):
    """Synthesize RF R-component data set over range of inclinations and distances
    and get result as a rf.RFStream instance.

    :param H: Moho depth (km)
    :type H: float
    :param V_p: P body wave velocity in uppermost layer
    :type V_p: float
    :param V_s: S body wave velocity in uppermost layer
    :type V_s: float
    :param inclinations: Array of inclinations for which to create RFs
    :type inclinations: numpy.array(float)
    :param distances: Array of teleseismic distances corresponding to inclinations
    :type distances: numpy.array(float)
    :param ds: Final sampling rate (Hz) for the downsampled output signal
    :type ds: float
    :param log: Logger to send output to, defaults to None
    :type log: logger, optional
    :param include_t3: If True, include the third expected multiple PpSs+PsPs
    :type include_t3: bool, optional
    :param amplitudes: Custom amplitudes to apply to the multiples
    :type amplitudes: list(float), optional
    :param baz: Back azimuth for metadata
    :type baz: float, optional
    :return: Stream containing synthetic RFs
    :rtype: rf.RFStream
    """
    assert len(inclinations) == len(distances), "Must provide 1:1 inclination and distance pairs"

    k = V_p/V_s
    traces = []
    for i, inc_deg in enumerate(inclinations):
        theta_p = np.deg2rad(inc_deg)
        p = np.sin(theta_p)/V_p

        t1 = H*(np.sqrt((k*k/V_p/V_p) - p*p) - np.sqrt(1.0/V_p/V_p - p*p))
        t2 = H*(np.sqrt((k*k/V_p/V_p) - p*p) + np.sqrt(1.0/V_p/V_p - p*p))
        arrivals = [t1, t2]
        if include_t3:
            t3 = t1 + t2
            arrivals.append(t3)
        if log is not None:
            log.info("Inclination {:3g} arrival times: {}".format(inc_deg, arrivals))

        arrivals = [0] + arrivals
        if amplitudes is None:
            amplitudes = [1, 0.5, 0.4]
            if include_t3:
                amplitudes.append(-0.3)
            # end if
        else:
            assert len(amplitudes) == 3 + int(include_t3)
            # t3 amplitude should be negative
            assert (not include_t3) or (amplitudes[3] <= 0)
        # end if
        window = (-5.0, 50.0)  # sec
        fs = 100.0  # Hz
        _, synth_signal = generate_synth_rf(arrivals, amplitudes, fs_hz=fs, window_sec=window)

        now = obspy.UTCDateTime.now()
        # Make sure time difference of events is at least 1 second, since onset time is used as part of
        # logic for identifying related channels in rf.RFStream.
        now += float(i)
        dt = float(window[1] - window[0])
        end = now + dt
        onset = now - window[0]
        header = {'network': 'SY', 'station': 'TST', 'location': 'GA', 'channel': 'HHR', 'sampling_rate': fs,
                  'starttime': now, 'endtime': end, 'onset': onset,
                  'station_latitude': -19.0, 'station_longitude': 137.0,  # arbitrary (approx location of OA deployment)
                  'slowness': p*rf_util.KM_PER_DEG, 'inclination': inc_deg,
                  'back_azimuth': baz, 'distance': float(distances[i])}
        tr = rf.rfstream.RFTrace(data=synth_signal.copy(), header=header)
        tr = tr.decimate(int(np.round(fs/ds)), no_filter=True)
        traces.append(tr)
    # end for

    stream = rf.RFStream(traces)

    return stream, arrivals
# end func


def convert_inclination_to_distance(inclinations, model="iasp91", nominal_source_depth_km=10.0):
    """Helper function to convert range of inclinations to teleseismic distance in degrees.

    :param inclinations: Array of inclination angles in degrees
    :type inclinations: numpy.array(float)
    :param model: Name of model to use for ray tracing, defaults to "iasp91"
    :type model: str, optional
    :param nominal_source_depth_km: Assumed depth of source events, defaults to 10.0
    :type nominal_source_depth_km: float, optional
    :return: Array of teleseismic distances in degrees corresponding to input inclination angles.
    :rtype: numpy.array(float)
    """
    # Generate function mapping ray parameter to teleseismic distance.
    # The distances are not strictly required for H-k stacking, but rf behaves better when they are there.
    ts_distance = np.linspace(25, 95, 71)
    inc = np.zeros_like(ts_distance)
    model = obspy.taup.TauPyModel(model=model)
    source_depth_km = nominal_source_depth_km
    for i, d in enumerate(ts_distance):
        ray = model.get_ray_paths(source_depth_km, d, phase_list=['P'])
        inc[i] = ray[0].incident_angle  # pylint: disable=unsupported-assignment-operation
    # end for

    # Create interpolator to convert inclination to teleseismic distance
    interp_dist = interp1d(inc, ts_distance, bounds_error=False,
                           fill_value=(np.max(ts_distance), np.min(ts_distance)))

    # Loop over valid range of inclinations and generate synthetic RFs
    distances = interp_dist(inclinations)

    return distances
# end func


def synthesize_ideal_seismogram(network, station, units, sourcelatitude, sourcelongitude, sourcedepthmetres=80000,
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
    dist_deg = dist / 1000 / rf_util.KM_PER_DEG
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
