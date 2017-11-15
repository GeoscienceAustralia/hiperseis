import os
import csv
import glob
import numpy
import pytest
from obspy import read_events
from seismic.cluster.cluster import (process_event,
                                     read_stations,
                                     process_many_events)

TESTS = os.path.dirname(__file__)
EVENTS = os.path.join(TESTS, 'mocks', 'events')
xmls = glob.glob(os.path.join(EVENTS, '*.xml'))
engdhal_xmls = glob.glob(os.path.join(EVENTS, 'engdahl_sample', '*.xml'))
stations_file = os.path.join(TESTS, 'mocks', 'inventory', 'stations.csv')
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')


@pytest.fixture(params=xmls + engdhal_xmls)
def event_xml(request):
    return request.param


@pytest.fixture(params=['P S'])
def arr_type(request):
    return request.param


def test_single_event_output(xml, random_filename):
    p_file = random_filename(ext='_p.csv')
    s_file = random_filename(ext='_s.csv')
    event = read_events(xml).events[0]
    origin = event.preferred_origin()
    with open(stations_file, 'r') as sta_f:
        with open(p_file, 'w') as p_writer:
            with open(s_file, 'w') as s_writer:
                p_writer = csv.writer(p_writer)
                process_event(read_events(xml)[0],
                              stations=read_stations(sta_f),
                              p_writer=p_writer,
                              s_writer=s_writer,
                              nx=1440, ny=720, dz=25.0,
                              wave_type='P S')
    inputs = numpy.genfromtxt(saved_out, delimiter=',')
    outputs = numpy.genfromtxt(p_file, delimiter=',')
    numpy.testing.assert_array_almost_equal(inputs, outputs)

    # s_file is created
    assert os.path.exists(s_file)

    # make sure number of arrivals match that of output lines
    assert len(origin.arrivals) == len(outputs)  # no s arrivals for this event


def test_single_event_arrivals(event_xml, random_filename, arr_type):
    p_file = random_filename(ext='_p.csv')
    s_file = random_filename(ext='_s.csv')
    event = read_events(event_xml).events[0]
    origin = event.preferred_origin()
    with open(stations_file, 'r') as sta_f:
        with open(p_file, 'w') as p_wrt:
            with open(s_file, 'w') as s_wrt:
                p_writer = csv.writer(p_wrt)
                s_writer = csv.writer(s_wrt)
                process_event(read_events(event_xml)[0],
                              stations=read_stations(sta_f),
                              p_writer=p_writer,
                              s_writer=s_writer,
                              nx=1440, ny=720, dz=25.0,
                              wave_type=arr_type)
    outputs_p = numpy.genfromtxt(p_file, delimiter=',')
    outputs_s = numpy.genfromtxt(s_file, delimiter=',')

    stations = read_stations(open(stations_file, 'r'))

    p_arrivals = []
    s_arrivals = []
    p, s = arr_type.split()
    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
            ).waveform_id.station_code
        if sta_code in stations:
            if arr.phase == p:
                p_arrivals.append(arr)
            if arr.phase == s:
                s_arrivals.append(arr)

    # make sure number of arrivals match that of output lines
    assert len(p_arrivals) == len(outputs_p)
    assert len(s_arrivals) == len(outputs_s)


def test_sorted():
    pass


def test_filtered():
    pass


def test_multiple_event_output(random_filename):

    events = read_events(os.path.join(EVENTS, '*.xml')).events
    outfile = random_filename()

    with open(stations_file, 'r') as sta_f:
        process_many_events(events,
                            stations=read_stations(sta_f),
                            nx=1440, ny=720, dz=25.0,
                            wave_type='P S',
                            output_file=outfile)

    # check files created
    assert os.path.exists(outfile + '_' + 'P' + '.csv')
    assert os.path.exists(outfile + '_' + 'S' + '.csv')

    # check all arrivals are present
    arrivals = []
    for e in events:
        origin = e.preferred_origin()
        arrivals += origin.arrivals

    p_arr = numpy.genfromtxt(outfile + '_' + 'P' + '.csv', delimiter=',')
    s_arr = numpy.genfromtxt(outfile + '_' + 'S' + '.csv', delimiter=',')

    assert len(arrivals) == len(p_arr) + len(s_arr)


def test_matched_files():
    pass
