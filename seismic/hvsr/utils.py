import math
import numpy as np
from obspy.core import UTCDateTime, Stream
from scipy.signal.signaltools import detrend
import stockwell as st
import pywt
from collections import defaultdict
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

def getIntervalsInAreas(areas, min_interval=None):
    """
    Searches areas and gets all intervals in it.

    Areas is a list with tuples. The first item is the beginning of the area
    and the second the end. Both in samples.
    """
    # The smallest area will be the interval if none is given.
    if min_interval is None:
        min_interval = min([i[1] - i[0] for i in areas])
        # Choose 99 Percent.
        min_interval *= 0.99
    intervals = []
    # Loop over each area.
    for area in areas:
        # The span of one area.
        span = area[1] - area[0]
        # How many intervals fit in the area?
        count = int(span / min_interval)
        # Determine the rest to be able to distribute the intervals evenly in
        # each area. Remove one just for savety reasons.
        rest = span - count * min_interval - 1
        # Distribute evenly.
        if count > 1:
            space_inbetween = rest // (count - 1)
            for _i in range(count):
                start = area[0] + _i * (min_interval + space_inbetween)
                end = start + min_interval
                intervals.append((start, end))
        # Center.
        elif count == 1:
            start = area[0] + rest / 2
            end = start + min_interval
            intervals.append((start, end))
    return intervals


def getAreasWithinThreshold(c_funct, threshold, min_width, feather=0):
    """
    Parses any characteristic function and returns all areas in samples since
    start and length where the values of the function are below (or above) a
    certain threshold.

    :type c_funct: Tuple with two numpy.ndarrays.
    :param c_funct: The first array are the x-values and the second array are
                    the y-values. If it just is one array, then it will
                    be assumed to be the y-values and the x-values will be
                    created using numpy.arange(len(y-values)).
    :type threshold: Integer, Float
    :param threshold: Threshold value.
    :type min_width: Integer
    :param min_width: Minimum width of the returned areas in samples. Any
                      smaller areas will be discarded.
    """
    if type(c_funct) == np.ndarray:
        y_values = c_funct
        x_values = np.arange(len(y_values))
    else:
        x_values = c_funct[0]
        y_values = c_funct[1]
    if len(x_values) != len(y_values):
        raise
    areas = []
    # Init values for loop.
    start = 0
    last = False
    for _i, _j in zip(x_values, y_values):
        if _j < threshold:
            last = True
            continue
        # Already larger than threshold.
        if last:
            if _i - start < min_width:
                start = _i
                last = False
                continue
            areas.append((start + feather, _i - feather))
        start = _i
        last = False
    if last and x_values[-1] - start >= min_width:
        areas.append((start + feather, x_values[-1] - feather))
    return np.array(areas)


def cwt_TFA(data, delta, nf, f_min=0, f_max=0, w0=8, freq_spacing='log'):
    """
    :param data: time dependent signal.
    :param delta: time step between two samples in data (in seconds)
    :param nf: number of logarithmically spaced frequencies between fmin and
               fmax
    :param f_min: minimum frequency (in Hz)
    :param f_max: maximum frequency (in Hz)
    :param wf: wavelet to use
    :param w0: parameter w0 for morlet wavelets
    :param useMlpy: use the continuous wavelet transform from MLPY by default, otherwise
           use alternative implementation from ObsPy
    :return: 1. time frequency representation of data, type numpy.ndarray of complex
                values, shape = (nf, len(data)).
             2. frequency bins
             3. wavelet scales
    """

    wl = data.shape[0]
    if f_min == 0:
        f_min = 1.0 / wl
    if f_max == 0:
        f_max = 0.4 / delta

    freqs = np.logspace(np.log10(f_min), np.log10(f_max), nf)
    # spacing option only applies to mlpy
    if (freq_spacing == 'linear'):
        oldnf = nf
        sr = 1.0/delta
        spacing = sr/wl
        f_min_idx = math.ceil(f_min / spacing)
        f_max_idx = math.ceil(f_max / spacing)
        nf = f_max_idx - f_min_idx + 1

        freqs = np.linspace(f_min_idx, f_max_idx, nf) * spacing
    # end if

    scales = pywt.scale2frequency('cmor2-1.25', freqs * delta)
    cfs, r_freqs = pywt.cwt(data, scales, 'cmor2-1.25', delta)

    assert np.allclose(freqs, r_freqs), 'Unknown error in pywt. Aborting..'

    return cfs, freqs, scales*delta
# end func

def st_TFA(data, delta, nf, f_min=0, f_max=0):
    cfs = st.st(data)
    npts = int(data.shape[0])

    # function that maps frequencies to indices
    def st_freq(f):
        return int(np.floor(f * npts / (1. / delta) + .5))
    # end func

    freqsOut = np.logspace(np.log10(f_min), np.log10(f_max), nf)
    indices = list(map(st_freq, freqsOut))

    return cfs[indices, :], freqsOut
#end func

def single_taper_spectrum(data, delta, taper_name=None):
    """
    Returns the spectrum and the corresponding frequencies for data with the
    given taper.
    """
    length = len(data)
    good_length = length // 2 + 1
    # Create the frequencies.
    # XXX: This might be some kind of hack
    freq = abs(np.fft.fftfreq(length, delta)[:good_length])
    # Create the tapers.
    if taper_name == 'bartlett':
        taper = np.bartlett(length)
    elif taper_name == 'blackman':
        taper = np.blackman(length)
    elif taper_name == 'boxcar':
        taper = np.ones(length)
    elif taper_name == 'hamming':
        taper = np.hamming(length)
    elif taper_name == 'hanning':
        taper = np.hanning(length)
    elif 'kaiser' in taper_name:
        taper = np.kaiser(length, beta=14)
    # Should never happen.
    else:
        msg = 'Something went wrong.'
        raise Exception(msg)
    # Detrend the data.
    data = detrend(data)
    # Apply the taper.
    data *= taper
    spec = abs(np.fft.rfft(data)) ** 2
    return spec, freq
#end func

def waveform_iterator3c(fds:FederatedASDFDataSet, net:str, sta:str,
                        st:UTCDateTime, et:UTCDateTime, step=86400):
    # get stations
    rows = fds.get_stations(st, et, network=net, station=sta)
    uloc_dict = defaultdict(set)  # dict of channels keyed by location code
    for row in rows: uloc_dict[row[2]].add(row[3])

    wst = st
    wet = st
    while (wet < et):
        if (wst + step > et):
            wet = et
        else:
            wet = wst + step
        # end if

        wd = Stream()
        for loc in uloc_dict.keys():
            channels = uloc_dict[loc]
            if(len(channels)>3):
                print('Expecting 3 but found {} channels for {}.{}.{}. '
                      'Moving along..'.format(len(channels), net, sta, loc))
                continue
            # end if
            for cha in channels:
                wd += fds.get_waveforms(net, sta, loc, cha, wst, wet)
            # end for

            # merge filling gaps
            try:
                wd.merge(method=1, fill_value=0)
            except:
                print(('Merge failed for station %s.%s.%s..' % (net, sta, loc)))
                wd = Stream()
            # end try
            yield wd
        # end for

        wst += step
    # wend
# end func