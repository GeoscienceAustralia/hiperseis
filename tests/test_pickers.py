import os
import pytest
from obspy.core import read as obspy_read, Stream
# from obspy.core.event import Event, read_events
import sys

is_windows = sys.platform.startswith('win')
algos = ['aicdpicker', 'fbpicker', 'ktpicker']


@pytest.fixture(params=algos)
def algorithm(request):
    return request.param


# def test_pickermaps(algorithm):
#     """
#     basic operation test for now
#     """
#     picker = pickermaps[algorithm]()
#     st = obspy_read(mseed)
#     for s in st[:1]:
#         picker.picks(s)

@pytest.mark.skipif(is_windows, reason='Availability of compilers cannot be guaranteed')
def test_pick_amplitude_assocs(miniseed_conf, algorithm, mseed):
    from legacy.pickers_integration.pickers import pickermaps
    picker = pickermaps[algorithm]()
    st = obspy_read(mseed)
    st2 = Stream(st[0:1])
    event = picker.event(st2, config=miniseed_conf)
    assert len(event.picks) == len(event.amplitudes)
    for p, a in zip(event.picks, event.amplitudes):
        assert a.pick_id == p.resource_id
