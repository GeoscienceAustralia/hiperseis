from obspy.core import Trace, Stats
import pytest
import numpy as np
from scipy import signal
from obspy.core import Trace, UTCDateTime, Stats
from obspy.core.event import Event, Origin, Magnitude
from rf import RFStream, rfstats
from rf.util import DEG2KM
from obspy.taup import TauPyModel
from seismic.receiver_fn.rf_stacking import compute_hk_stack

def synthetic_radial_component(trc_npts, offset, sr, Vp, h, k, p):
    Vp_inv = 1./Vp
    Vs_inv = k * Vp_inv

    term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
    term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

    t1 = h * (term1 - term2)
    t2 = h * (term1 + term2)
    t3 = h * 2 * term1

    t = np.linspace(-1, 1, 5 * sr)
    _, _, phase = signal.gausspulse(t, fc=sr/2, retquad=True, retenv=True)

    trc_data = np.random.random(trc_npts)*0.01
    amp_factor = 1.
    for t in [t1, t2, t3]:
        idx = int((t+offset) * sr)
        trc_data[idx - len(phase)//2:idx + len(phase)//2] += phase*amp_factor
        amp_factor /= 2.
    # end for
    return trc_data
# end func

def get_rf_stream(nevents, Vp, h, k):
    """
    Generate RF-stream based on test parameters
    """
    dist_range = [31, 89]
    trc_npts = 700
    trc_sr = 10
    sta_lon = 0
    sta_lat = 0
    sta_elev = 0
    event_depth = 100e3
    offset = 10
    tt_model = TauPyModel(model='iasp91')

    trc_list = []
    for dist in np.linspace(*dist_range, nevents):
        e = Event(origins=[Origin(
                  time=UTCDateTime(2000, 1, 1), latitude=0,
                  longitude=dist, depth=event_depth)],
                  magnitudes=[Magnitude(mag=7.4)])

        arrivals = tt_model.get_travel_times(event_depth/1e3,
                                             dist, 'P',)
        arrival = arrivals[0]

        trc = Trace(data = np.zeros(trc_npts),
                    header=Stats({'sampling_rate':trc_sr,
                                  'npts':trc_npts,
                                  'network':'SYN',
                                  'station':'A',
                                  'channel':'BHR',
                                  'starttime':e.origins[0].time+
                                  arrival.time-offset}))
        trc.stats = rfstats(trc.stats, event=e,
                            station={'latitude':sta_lat,
                                     'longitude':sta_lon,
                                     'elevation':sta_elev})

        p = trc.stats.slowness / DEG2KM
        trc.data = synthetic_radial_component(trc_npts,
                                              offset,
                                              trc_sr, Vp,
                                              h, k, p)

        trc_list.append(trc)
    # end for
    rfstream = RFStream(trc_list)

    return rfstream
# end func

@pytest.fixture(params=[5, 10])
def nevents(request):
    return request.param

@pytest.fixture(params=np.linspace(6, 7, 3))
def Vp(request):
    return request.param

@pytest.fixture(params=np.linspace(30, 60, 3))
def h(request):
    return request.param

@pytest.fixture(params=np.linspace(1.7, 1.9, 3))
def k(request):
    return request.param

@pytest.fixture(params=[True, False])
def semblance_weighting(request):
    return request.param

def test_stacking(nevents, Vp, h, k, semblance_weighting):
    rfstream = get_rf_stream(nevents, Vp, h, k)

    kg, hg, stack = compute_hk_stack(rfstream, semblance_weighted=semblance_weighting)
    hc, kc = hg.ravel()[np.argmax(stack)], kg.ravel()[np.argmax(stack)]

    htol = 10
    ktol = 0.2
    assert np.fabs(h - hc) < htol
    assert np.fabs(k - kc) < ktol
# end func